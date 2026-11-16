"""泳道拆样台:上传 R1/R2/I1 与样本索引 CSV,按样本拆出成对读段。

发布模型:
- 处理中:上传件与半成品输出都在 staging/<token>/,不对外可见。
- 成功:输出目录整体 os.rename 到 published/<运行号>/,同事务写入计数,
  之后状态置为 published。任何一步失败都不会留下半套已发布文件。
- 失败/中断:状态置 failed,暂存与可能的半成品目录被清理;同一运行号可重试。
- 已发布:发布时登记每组 R1/R2 的字节数与 SHA-256;启动及下载前流式重算,
  发现缺失/截断/等长内容不符即把整次运行标 failed 并移出供数区;
  同一运行号拒绝改写。
"""
import os
import re
import shutil
import sqlite3
import threading
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, RedirectResponse
from fastapi.templating import Jinja2Templates

from .config import db_path, published_dir, quarantine_dir, staging_dir
from .db import db_conn, init_db
from .demux import demultiplex
from .integrity import IntegrityCheckError, hash_file, verify_published
from .samples import load_samples

RUN_NUMBER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
STATUS_LABELS = {"processing": "处理中", "published": "已发布", "failed": "失败"}
UNDETERMINED = "undetermined"

BASE_DIR = Path(__file__).resolve().parent
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _remove(path) -> None:
    p = Path(path)
    if p.is_dir() and not p.is_symlink():
        shutil.rmtree(p, ignore_errors=True)
    else:
        try:
            p.unlink()
        except (FileNotFoundError, NotADirectoryError):
            pass


def quarantine_run(run_number: str) -> None:
    """把核验失败的已发布目录移出供数区(保留现场),移不动则删除。"""
    src = published_dir() / run_number
    if not src.exists() and not src.is_symlink():
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    dest = quarantine_dir() / f"{run_number}-{stamp}-{uuid.uuid4().hex[:8]}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.rename(src, dest)
    except OSError:
        _remove(src)


def fail_integrity(run_id: int, run_number: str, detail: str) -> str:
    """把一次运行标为失败并停止供数;返回写入页面的原因。"""
    reason = f"已发布数据完整性核验失败,已停止整次运行供数: {detail}"
    with db_conn(db_path()) as conn:
        conn.execute(
            "UPDATE runs SET status='failed', error=?, integrity_mode=NULL, "
            "verified_at=NULL, updated_at=? WHERE id=? AND status='published'",
            (reason, utcnow(), run_id),
        )
    quarantine_run(run_number)
    return reason


def _published_snapshot():
    """一次性取出待核验的已发布运行、分组计数和内容摘要依据。"""
    with db_conn(db_path()) as conn:
        runs = conn.execute(
            "SELECT * FROM runs WHERE status='published' ORDER BY id"
        ).fetchall()
        snapshot = []
        for run in runs:
            counts = [
                (r["sample_id"], r["pairs"])
                for r in conn.execute(
                    "SELECT sample_id, pairs FROM sample_counts WHERE run_id=?",
                    (run["id"],),
                )
            ]
            evidence = {
                r["filename"]: {
                    "size_bytes": r["size_bytes"],
                    "sha256": r["sha256"],
                    "records": r["records"],
                }
                for r in conn.execute(
                    "SELECT filename, size_bytes, sha256, records "
                    "FROM published_files WHERE run_id=?",
                    (run["id"],),
                )
            }
            snapshot.append((dict(run), counts, evidence))
    return snapshot


def recover_interrupted() -> None:
    """启动恢复:凡未走到 published 的运行一律标失败并清掉半成品;
    已发布运行逐个重新核验,任何一组文件缺失、截断或内容不符只停该次运行,
    不影响其他已发布运行。"""
    with db_conn(db_path()) as conn:
        stale = conn.execute(
            "SELECT id, run_number, staging_dir FROM runs WHERE status = 'processing'"
        ).fetchall()
        for row in stale:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                ("进程中断,未完成即停止;可用同一运行号重新提交", utcnow(), row["id"]),
            )
        published = {
            r["run_number"]
            for r in conn.execute("SELECT run_number FROM runs WHERE status='published'")
        }
    for row in stale:
        _remove(row["staging_dir"])
        _remove(published_dir() / row["run_number"])
    if staging_dir().exists():
        for child in staging_dir().iterdir():  # 无主的暂存目录一并清理
            _remove(child)
    if published_dir().exists():
        for child in published_dir().iterdir():
            if child.name not in published:
                _remove(child)

    for run, counts, evidence in _published_snapshot():
        try:
            mode = verify_published(
                published_dir() / run["run_number"], counts,
                run["total_pairs"], evidence, run.get("integrity_mode"),
            )
            with db_conn(db_path()) as conn:
                conn.execute(
                    "UPDATE runs SET integrity_mode=?, verified_at=? WHERE id=?",
                    (mode, utcnow(), run["id"]),
                )
        except Exception as exc:  # 单次运行损坏不得影响其他已发布运行
            try:
                fail_integrity(run["id"], run["run_number"], str(exc))
            except Exception:
                pass


def setup() -> None:
    staging_dir().mkdir(parents=True, exist_ok=True)
    published_dir().mkdir(parents=True, exist_ok=True)
    init_db(db_path())
    recover_interrupted()


@asynccontextmanager
async def lifespan(_app: FastAPI):
    setup()
    yield


app = FastAPI(title="泳道拆样台", lifespan=lifespan)


# ---------------------------------------------------------------- 后台处理

def write_log(path: Path, run, sample_list, stats) -> None:
    lines = [
        f"运行号: {run['run_number']}",
        f"提交时间: {run['created_at']}",
        f"完成时间: {utcnow()}",
        f"输入文件: R1={run['r1_name']} R2={run['r2_name']} "
        f"I1={run['i1_name']} 样本表={run['samples_name']}",
        "样本索引:",
    ]
    lines += [f"  {sid}\t{idx}" for sid, idx in sample_list]
    lines.append(f"总读段对: {stats['total_pairs']}")
    lines.append("分组对数(对数\t其中单碱基纠错):")
    for sid in stats["order"]:
        pairs, corrected = stats["counts"][sid]
        lines.append(f"  {sid}\t{pairs}\t{corrected}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_run(run_id: int) -> None:
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        return
    run_number = run["run_number"]
    staging = Path(run["staging_dir"])
    uploads = staging / "uploads"
    out = staging / "out"
    final = published_dir() / run_number
    try:
        sample_list = load_samples(uploads / "samples.csv")
        out.mkdir(parents=True)
        stats = demultiplex(
            str(uploads / "r1.fastq"), str(uploads / "r2.fastq"),
            str(uploads / "i1.fastq"), sample_list, str(out),
        )
        write_log(out / "demux.log", run, sample_list, stats)

        # 发布前为每组 R1/R2 登记字节数、记录数和流式 SHA-256;这是重启/下载
        # 核验的内容依据,能识别截断以及无法靠结构发现的等长内容改写。
        evidence = {}
        for sid in stats["order"]:
            pairs = stats["counts"][sid][0]
            for suffix in ("R1", "R2"):
                filename = f"{sid}_{suffix}.fastq"
                size_bytes, sha256 = hash_file(out / filename)
                evidence[filename] = {
                    "size_bytes": size_bytes,
                    "sha256": sha256,
                    "records": pairs,
                }
        log_size, log_sha = hash_file(out / "demux.log")
        evidence["demux.log"] = {
            "size_bytes": log_size, "sha256": log_sha, "records": 0,
        }

        if final.exists():
            raise RuntimeError(f"发布目录已存在: {final}")
        os.rename(out, final)  # 同一数据卷内原子移动:先完整落盘,再发布
        with db_conn(db_path()) as conn:
            updated = conn.execute(
                "UPDATE runs SET status='published', total_pairs=?, published_at=?, "
                "updated_at=?, integrity_mode='digest', verified_at=? "
                "WHERE id=? AND status='processing'",
                (stats["total_pairs"], utcnow(), utcnow(), utcnow(), run_id),
            ).rowcount
            if updated != 1:
                raise RuntimeError("运行状态已变化,放弃发布")
            index_of = dict(sample_list)
            for sid in stats["order"]:
                pairs, corrected = stats["counts"][sid]
                conn.execute(
                    "INSERT INTO sample_counts(run_id, sample_id, index_seq, pairs, "
                    "corrected) VALUES (?,?,?,?,?)",
                    (run_id, sid, index_of.get(sid, ""), pairs, corrected),
                )
            for filename, item in evidence.items():
                conn.execute(
                    "INSERT INTO published_files(run_id, filename, size_bytes, sha256, "
                    "records) VALUES (?,?,?,?,?)",
                    (run_id, filename, item["size_bytes"], item["sha256"],
                     item["records"]),
                )
        _remove(staging)
    except Exception as exc:  # 任何失败:标失败 + 清半成品,不发布
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (str(exc), utcnow(), run_id),
            )
        _remove(staging)
        _remove(final)


# ---------------------------------------------------------------- 页面与接口

@app.exception_handler(HTTPException)
async def http_exception_handler(request: Request, exc: HTTPException):
    if "text/html" in request.headers.get("accept", ""):
        return templates.TemplateResponse(
            request, "error.html",
            {"status": exc.status_code, "detail": exc.detail},
            status_code=exc.status_code,
        )
    return JSONResponse({"detail": exc.detail}, status_code=exc.status_code)


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    with db_conn(db_path()) as conn:
        runs = conn.execute("SELECT * FROM runs ORDER BY id DESC").fetchall()
    return templates.TemplateResponse(
        request, "index.html", {"runs": runs, "labels": STATUS_LABELS}
    )


async def save_upload(upload: UploadFile, dest: Path) -> None:
    with open(dest, "wb") as fh:
        while chunk := await upload.read(1 << 20):
            fh.write(chunk)


@app.post("/runs")
async def create_run(
    run_number: str = Form(...),
    r1: UploadFile = File(...),
    r2: UploadFile = File(...),
    i1: UploadFile = File(...),
    samples: UploadFile = File(...),
):
    run_number = run_number.strip()
    if not RUN_NUMBER_RE.fullmatch(run_number):
        raise HTTPException(400, "运行号格式不正确:字母或数字开头,可含 - _ .,最长 64 字符")
    for label, up in (("R1", r1), ("R2", r2), ("I1", i1), ("样本表", samples)):
        if not up.filename:
            raise HTTPException(400, f"请上传{label}文件")

    token = uuid.uuid4().hex
    staging = staging_dir() / token
    uploads = staging / "uploads"
    names = {k: Path(u.filename or "").name for k, u in
             (("r1", r1), ("r2", r2), ("i1", i1), ("samples", samples))}
    try:
        with db_conn(db_path()) as conn:
            row = conn.execute(
                "SELECT id, status FROM runs WHERE run_number=?", (run_number,)
            ).fetchone()
            if row and row["status"] == "published":
                raise HTTPException(409, f"运行 {run_number} 已发布,不可改写")
            if row and row["status"] == "processing":
                raise HTTPException(409, f"运行 {run_number} 正在处理中,请稍后查看")
            if row:  # failed:允许同号重试
                conn.execute("DELETE FROM runs WHERE id=?", (row["id"],))
                retry_stale_dir = published_dir() / run_number
            else:
                retry_stale_dir = None
            cur = conn.execute(
                "INSERT INTO runs(run_number, status, staging_dir, r1_name, r2_name, "
                "i1_name, samples_name, created_at, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                (run_number, "processing", str(staging), names["r1"], names["r2"],
                 names["i1"], names["samples"], utcnow(), utcnow()),
            )
            run_id = cur.lastrowid
    except sqlite3.IntegrityError:
        raise HTTPException(409, f"运行 {run_number} 已存在,请刷新后查看")

    if retry_stale_dir is not None:
        _remove(retry_stale_dir)  # 失败重试:清掉上次可能残留的发布目录

    try:
        uploads.mkdir(parents=True)
        await save_upload(r1, uploads / "r1.fastq")
        await save_upload(r2, uploads / "r2.fastq")
        await save_upload(i1, uploads / "i1.fastq")
        await save_upload(samples, uploads / "samples.csv")
        load_samples(uploads / "samples.csv")  # 先校验样本表,尽早报错
    except HTTPException:
        raise
    except Exception as exc:
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET status='failed', error=?, updated_at=? WHERE id=?",
                (f"上传或样本表错误: {exc}", utcnow(), run_id),
            )
        _remove(staging)
        _remove(published_dir() / run_number)
        raise HTTPException(400, f"上传或样本表错误: {exc}")

    threading.Thread(target=process_run, args=(run_id,), daemon=True).start()
    return RedirectResponse(f"/runs/{run_id}", status_code=303)


def get_run_or_404(run_id: int):
    with db_conn(db_path()) as conn:
        run = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    if run is None:
        raise HTTPException(404, "运行不存在")
    return run


def get_counts(run_id: int):
    with db_conn(db_path()) as conn:
        return conn.execute(
            "SELECT * FROM sample_counts WHERE run_id=? "
            "ORDER BY sample_id = 'undetermined', sample_id",
            (run_id,),
        ).fetchall()


def get_evidence(run_id: int):
    with db_conn(db_path()) as conn:
        return {
            r["filename"]: {
                "size_bytes": r["size_bytes"],
                "sha256": r["sha256"],
                "records": r["records"],
            }
            for r in conn.execute(
                "SELECT filename, size_bytes, sha256, records "
                "FROM published_files WHERE run_id=?",
                (run_id,),
            )
        }


_verify_locks = [threading.Lock() for _ in range(32)]


def _verify_lock(run_id: int) -> threading.Lock:
    """按运行号分片的锁:同一运行只核验一次,内存不随历史运行数增长。"""
    return _verify_locks[run_id % len(_verify_locks)]


def verify_run_or_fail(run):
    """供数前核验整次运行;任何一组文件有问题就停掉整次运行,抛 409。

    同一运行并发下载时只让一个线程做核验,其余线程等待结果。
    """
    with _verify_lock(run["id"]):
        fresh = get_run_or_404(run["id"])
        if fresh["status"] != "published":
            raise HTTPException(409, fresh["error"] or "运行未发布,没有完整文件可下载")
        counts = [(c["sample_id"], c["pairs"]) for c in get_counts(fresh["id"])]
        try:
            mode = verify_published(
                published_dir() / fresh["run_number"], counts,
                fresh["total_pairs"], get_evidence(fresh["id"]),
                fresh["integrity_mode"],
            )
        except IntegrityCheckError as exc:
            reason = fail_integrity(fresh["id"], fresh["run_number"], str(exc))
            raise HTTPException(409, reason) from exc
        except Exception as exc:  # 核验本身出错同样停止供数,不冒险提供文件
            reason = fail_integrity(
                fresh["id"], fresh["run_number"], f"核验无法完成: {exc}"
            )
            raise HTTPException(409, reason) from exc
        with db_conn(db_path()) as conn:
            conn.execute(
                "UPDATE runs SET integrity_mode=?, verified_at=? WHERE id=?",
                (mode, utcnow(), fresh["id"]),
            )
        return fresh


@app.get("/runs/{run_id}", response_class=HTMLResponse)
def run_detail(request: Request, run_id: int):
    run = get_run_or_404(run_id)
    groups = [
        {
            "sample_id": c["sample_id"],
            "index_seq": c["index_seq"],
            "pairs": c["pairs"],
            "corrected": c["corrected"],
            "r1": f"{c['sample_id']}_R1.fastq",
            "r2": f"{c['sample_id']}_R2.fastq",
        }
        for c in get_counts(run_id)
    ]
    return templates.TemplateResponse(
        request, "run.html",
        {"run": run, "groups": groups, "labels": STATUS_LABELS,
         # 未登记摘要的已发布运行按"仅结构核验"展示,绝不冒充有内容依据。
         "integrity_mode": (run["integrity_mode"] or "structure")
         if run["status"] == "published" else None},
    )


@app.get("/api/runs/{run_id}")
def run_status(run_id: int):
    run = get_run_or_404(run_id)
    published = run["status"] == "published"
    counts = get_counts(run_id) if published else []
    return {
        "run_number": run["run_number"],
        "status": run["status"],
        "error": run["error"],
        # 未发布(含核验失败)时不提供任何可被当作可信计数的数字。
        "total_pairs": run["total_pairs"] if published else None,
        "integrity_mode": run["integrity_mode"] if published else None,
        "verified_at": run["verified_at"] if published else None,
        "counts": [
            {"sample_id": c["sample_id"], "index_seq": c["index_seq"],
             "pairs": c["pairs"], "corrected": c["corrected"]}
            for c in counts
        ],
    }


@app.get("/runs/{run_id}/download/{filename}")
def download(run_id: int, filename: str):
    run = get_run_or_404(run_id)
    if run["status"] != "published":
        raise HTTPException(409, run["error"] or "运行未发布,没有完整文件可下载")
    valid = {"demux.log"}
    for c in get_counts(run_id):
        valid.add(f"{c['sample_id']}_R1.fastq")
        valid.add(f"{c['sample_id']}_R2.fastq")
    if filename not in valid:
        raise HTTPException(404, "文件不存在")
    run = verify_run_or_fail(run)  # 核验整次运行后,才提供其中任何文件
    path = published_dir() / run["run_number"] / filename
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, filename=filename, media_type="application/octet-stream")
