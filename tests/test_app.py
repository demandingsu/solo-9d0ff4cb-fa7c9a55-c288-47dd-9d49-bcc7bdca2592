"""关键情形自测:拆分正确性、纠错规则、损坏处理、原子发布、恢复、不可改写。"""
import os
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.demux import (  # noqa: E402
    FastqFormatError,
    assign_index,
    build_matcher,
    iter_fastq,
    read_id,
)
from app.samples import SampleSheetError, load_samples  # noqa: E402

SAMPLES_CSV = "sample_id,index\nS1,AAAACCCC\nS2,AAAACCGG\nS3,GGGGTTTT\n"

# S1=AAAACCCC 与 S2=AAAACCGG 相距 2,便于构造歧义错配;S3 远离二者
READS = [
    ("r001", "AAAACCCC"),  # S1 精确
    ("r002", "AAAACCGG"),  # S2 精确
    ("r003", "AAAACCCA"),  # 距 S1 唯一单碱基错配 -> 纠错到 S1
    ("r004", "AAAACCGC"),  # 距 S1、S2 各 1 -> 歧义 -> 未分配
    ("r005", "NNNNNNNN"),  # 无候选 -> 未分配
    ("r006", "GGGGTTTA"),  # 距 S3 唯一单碱基错配 -> 纠错到 S3
    ("r007", "AAAACCCC"),  # S1 精确(检验同组内顺序)
]

SEQ1, SEQ2 = "ACGT" * 10, "TGCA" * 10
QUAL = "I" * 40


def fq(records):
    return "".join(f"@{h}\n{s}\n+\n{q}\n" for h, s, q in records)


def make_inputs(reads=READS):
    r1, r2, i1 = [], [], []
    for n, (rid, idx) in enumerate(reads):
        head = f"INST:1:FC:1:1101:{1000 + n}:{2000 + n}"
        r1.append((f"{head} 1:N:0:0", SEQ1, QUAL))
        r2.append((f"{head} 2:N:0:0", SEQ2, QUAL))
        i1.append((f"{head} 1:N:0:0", idx, "I" * len(idx)))
    return fq(r1).encode(), fq(r2).encode(), fq(i1).encode(), SAMPLES_CSV.encode()


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATA_DIR", str(tmp_path / "data"))
    from fastapi.testclient import TestClient
    from app.main import app

    with TestClient(app) as c:
        c.data_dir = tmp_path / "data"
        yield c


def submit(client, run_number, r1, r2, i1, samples):
    return client.post(
        "/runs",
        data={"run_number": run_number},
        files={
            "r1": ("r1.fastq", r1),
            "r2": ("r2.fastq", r2),
            "i1": ("i1.fastq", i1),
            "samples": ("samples.csv", samples),
        },
        follow_redirects=False,
    )


def submit_good(client, run_number="RUN1"):
    r1, r2, i1, csv = make_inputs()
    return submit(client, run_number, r1, r2, i1, csv)


def run_id_of(resp):
    assert resp.status_code == 303, resp.text
    return int(resp.headers["location"].rsplit("/", 1)[1])


def wait_done(client, run_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = client.get(f"/api/runs/{run_id}").json()
        if st["status"] != "processing":
            return st
        time.sleep(0.05)
    raise TimeoutError("处理超时")


def parse_fastq_ids(text):
    lines = text.splitlines()
    return [lines[i][1:].split()[0] for i in range(0, len(lines), 4)]


def head_of(n):
    """第 n 条 READS 记录对应的读段号(见 make_inputs)。"""
    return f"INST:1:FC:1:1101:{1000 + n}:{2000 + n}"


def counts_of(status):
    return {c["sample_id"]: (c["pairs"], c["corrected"]) for c in status["counts"]}


# ------------------------------------------------------------ 单元级

def test_assign_index_rules():
    samples = [("S1", "AAAACCCC"), ("S2", "AAAACCGG"), ("S3", "GGGGTTTT")]
    exact, neighbors = build_matcher(samples)
    assert assign_index("AAAACCCC", exact, neighbors) == ("S1", False)   # 精确优先
    assert assign_index("AAAACCCA", exact, neighbors) == ("S1", True)    # 唯一错配
    assert assign_index("AAAACCGC", exact, neighbors) == (None, False)   # 歧义
    assert assign_index("NNNNNNNN", exact, neighbors) == (None, False)   # 无候选
    assert assign_index("AAAA", exact, neighbors) == (None, False)       # 长度不符


def test_read_id():
    assert read_id("INST:1:FC:1:1101:100:100 1:N:0:0") == "INST:1:FC:1:1101:100:100"
    assert read_id("PLAINID") == "PLAINID"
    assert read_id("") == ""


def test_iter_fastq_truncated(tmp_path):
    p = tmp_path / "bad.fastq"
    p.write_text("@a\nACGT\n+\n")  # 缺质量行
    with pytest.raises(FastqFormatError, match="截断"):
        list(iter_fastq(str(p), "R1"))


def test_iter_fastq_len_mismatch(tmp_path):
    p = tmp_path / "bad.fastq"
    p.write_text("@a\nACGT\n+\nIII\n")  # 质量短 1
    with pytest.raises(FastqFormatError, match="长度"):
        list(iter_fastq(str(p), "R1"))


def test_samples_validation(tmp_path):
    def load(text):
        p = tmp_path / "s.csv"
        p.write_text(text)
        return load_samples(str(p))

    assert load("sid,idx\nA,ACGT\n") == [("A", "ACGT")]
    assert load("A,acgt\n") == [("A", "ACGT")]  # 小写归一、无表头
    with pytest.raises(SampleSheetError, match="重复"):
        load("A,ACGT\nB,ACGT\n")
    with pytest.raises(SampleSheetError, match="长度不一致"):
        load("A,ACGT\nB,ACGTA\n")
    with pytest.raises(SampleSheetError, match="非法字符"):
        load("A,ACGT\nB,ACGB\n")
    with pytest.raises(SampleSheetError, match="保留字"):
        load("undetermined,ACGT\n")


# ------------------------------------------------------------ 端到端

def test_happy_path(client):
    rid = run_id_of(submit_good(client))
    st = wait_done(client, rid)
    assert st["status"] == "published", st["error"]
    assert st["total_pairs"] == 7
    assert counts_of(st) == {
        "S1": (3, 1), "S2": (1, 0), "S3": (1, 1), "undetermined": (2, 0),
    }

    s1_r1 = client.get(f"/runs/{rid}/download/S1_R1.fastq")
    assert s1_r1.status_code == 200
    assert parse_fastq_ids(s1_r1.text) == [head_of(0), head_of(2), head_of(6)]  # 原顺序
    s1_r2 = client.get(f"/runs/{rid}/download/S1_R2.fastq")
    assert parse_fastq_ids(s1_r2.text) == [head_of(0), head_of(2), head_of(6)]
    assert SEQ2 in s1_r2.text  # R2 内容确实是 R2
    und = client.get(f"/runs/{rid}/download/undetermined_R1.fastq")
    assert parse_fastq_ids(und.text) == [head_of(3), head_of(4)]

    log = client.get(f"/runs/{rid}/download/demux.log")
    assert log.status_code == 200 and "总读段对: 7" in log.text

    page = client.get(f"/runs/{rid}")
    assert page.status_code == 200 and "AAAACCCC" in page.text and "已发布" in page.text


def test_truncated_input_fails_then_retry_succeeds(client):
    r1, r2, i1, csv = make_inputs()
    bad_r2 = r2[: r2.rfind(b"\n", 0, r2.rfind(b"\n"))]  # 末条记录缺质量行
    rid = run_id_of(submit(client, "RUN2", r1, bad_r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed"
    assert "截断" in st["error"]
    # 不发布半套文件
    assert not (client.data_dir / "published" / "RUN2").exists()
    assert client.get(f"/runs/{rid}/download/S1_R1.fastq").status_code == 409

    # 同一运行号重试可完成
    rid2 = run_id_of(submit_good(client, "RUN2"))
    st2 = wait_done(client, rid2)
    assert st2["status"] == "published"
    assert st2["total_pairs"] == 7
    assert client.get(f"/runs/{rid2}/download/S1_R1.fastq").status_code == 200


def test_read_id_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    bad_r2 = r2.replace(b"1101:1001:2001", b"1101:9999:2001", 1)
    rid = run_id_of(submit(client, "RUN3", r1, bad_r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "读段号不一致" in st["error"]
    assert not (client.data_dir / "published" / "RUN3").exists()


def test_seq_qual_len_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    bad_i1 = i1.replace(b"IIIIIIII\n", b"IIIIIII\n", 1)  # 首条 I1 质量少 1
    rid = run_id_of(submit(client, "RUN4", r1, r2, bad_i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "长度" in st["error"]


def test_record_count_mismatch_fails(client):
    r1, r2, i1, csv = make_inputs()
    extra = f"@INST:1:FC:1:1101:9:9 1:N:0:0\n{SEQ1}\n+\n{QUAL}\n".encode()
    rid = run_id_of(submit(client, "RUN5", r1 + extra, r2, i1, csv))
    st = wait_done(client, rid)
    assert st["status"] == "failed" and "记录数不一致" in st["error"]


def test_published_run_immutable(client):
    rid = run_id_of(submit_good(client, "RUN6"))
    assert wait_done(client, rid)["status"] == "published"
    before = client.get(f"/runs/{rid}/download/S1_R1.fastq").text

    resp = submit_good(client, "RUN6")  # 同号再提交
    assert resp.status_code == 409
    assert "不可改写" in resp.json()["detail"]
    after = client.get(f"/runs/{rid}/download/S1_R1.fastq").text
    assert before == after  # 文件未被改写


def test_bad_sample_sheet_rejected(client):
    r1, r2, i1, _ = make_inputs()
    dup = b"A,ACGT\nB,ACGT\n"
    resp = submit(client, "RUN7", r1, r2, i1, dup)
    assert resp.status_code == 400
    assert "重复" in resp.json()["detail"]


def test_interrupted_run_recovered_on_startup(client):
    """模拟进程中断:遗留 processing 记录 + 半成品目录,重启后标失败、清理干净。"""
    from app.main import recover_interrupted
    from app.db import db_conn
    from app.config import db_path, staging_dir, published_dir

    staging = staging_dir() / "deadbeef"
    (staging / "out").mkdir(parents=True)
    (staging / "out" / "S1_R1.fastq").write_text("半成品")
    stray_pub = published_dir() / "RUN8"
    stray_pub.mkdir(parents=True)
    (stray_pub / "S1_R1.fastq").write_text("改名后来不及落库的文件")
    with db_conn(db_path()) as conn:
        conn.execute(
            "INSERT INTO runs(run_number, status, staging_dir, created_at, updated_at) "
            "VALUES ('RUN8', 'processing', ?, 't', 't')",
            (str(staging),),
        )

    recover_interrupted()  # 模拟重启

    with db_conn(db_path()) as conn:
        row = conn.execute("SELECT * FROM runs WHERE run_number='RUN8'").fetchone()
    assert row["status"] == "failed" and "中断" in row["error"]
    assert not staging.exists()
    assert not stray_pub.exists()
    # 恢复后同一运行号可重试并完成
    rid = run_id_of(submit_good(client, "RUN8"))
    assert wait_done(client, rid)["status"] == "published"


def test_index_page(client):
    resp = client.get("/")
    assert resp.status_code == 200
    assert "新建拆样运行" in resp.text
