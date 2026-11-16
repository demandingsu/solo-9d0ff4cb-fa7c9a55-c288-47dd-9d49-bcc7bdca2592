"""已发布数据的完整性核验。

新运行发布时会把每组 R1/R2 的字节数与 SHA-256 登记到数据库,之后在启动
和下载前流式重算摘要,可发现缺失、截断和**等长内容改写**。旧库运行没有
摘要依据,只能按 FASTQ 结构、双端配对和入库对数核验;这种情况不会被说成
能识别无法追溯的等长改动。

所有文件读取都按固定大小的块进行,不会把整份 FASTQ 读入内存。
"""
from __future__ import annotations

import hashlib
from itertools import zip_longest
from pathlib import Path

from .demux import FastqFormatError, iter_fastq, read_id

CHUNK_SIZE = 1 << 20
_MISSING = object()


class IntegrityCheckError(Exception):
    """已发布文件缺失、截断、结构损坏或内容与发布登记不符。"""


def hash_file(path) -> tuple[int, str]:
    """流式返回 (字节数, SHA-256 十六进制摘要),内存占用与文件大小无关。"""
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK_SIZE):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _check_against_evidence(run_dir: Path, filename: str, size_bytes: int,
                            expected_sha256: str) -> None:
    path = run_dir / filename
    if not path.is_file():
        raise IntegrityCheckError(f"已发布文件缺失: {filename}")
    try:
        actual_size = path.stat().st_size
    except OSError as exc:
        raise IntegrityCheckError(f"无法读取 {filename}: {exc}") from exc
    if actual_size < size_bytes:
        raise IntegrityCheckError(
            f"已发布文件被截断: {filename}(实际 {actual_size} 字节,"
            f"发布时 {size_bytes} 字节)"
        )
    if actual_size > size_bytes:
        raise IntegrityCheckError(
            f"已发布文件大小异常: {filename}(实际 {actual_size} 字节,"
            f"发布时 {size_bytes} 字节)"
        )
    try:
        hashed_size, actual_sha256 = hash_file(path)
    except OSError as exc:
        raise IntegrityCheckError(f"无法读取 {filename}: {exc}") from exc
    if hashed_size != size_bytes or actual_sha256 != expected_sha256:
        raise IntegrityCheckError(
            f"已发布文件等长内容与发布时登记不符: {filename}"
            "(大小未变但 SHA-256 不一致,可能被改写)"
        )


def _verify_structure(run_dir: Path, counts, total_pairs: int) -> None:
    """旧库回退核验:FASTQ 结构 + R1/R2 配对 + 入库对数。"""
    seen_total = 0
    for sample_id, pairs in counts:
        r1_name = f"{sample_id}_R1.fastq"
        r2_name = f"{sample_id}_R2.fastq"
        for filename in (r1_name, r2_name):
            if not (run_dir / filename).is_file():
                raise IntegrityCheckError(f"已发布文件缺失: {filename}")

        seen = 0
        try:
            streams = zip_longest(
                iter_fastq(str(run_dir / r1_name), r1_name),
                iter_fastq(str(run_dir / r2_name), r2_name),
                fillvalue=_MISSING,
            )
            for rec1, rec2 in streams:
                if rec1 is _MISSING or rec2 is _MISSING:
                    raise IntegrityCheckError(
                        f"{sample_id}: R1/R2 记录数不一致(第 {seen + 1} 对起)"
                    )
                seen += 1
                if read_id(rec1[0]) != read_id(rec2[0]):
                    raise IntegrityCheckError(
                        f"{sample_id}: 第 {seen} 对读段号不一致: "
                        f"R1={rec1[0]!r} R2={rec2[0]!r}"
                    )
        except FastqFormatError as exc:
            raise IntegrityCheckError(str(exc)) from exc
        except OSError as exc:
            raise IntegrityCheckError(f"无法读取 {sample_id} 的输出文件: {exc}") from exc

        if seen != pairs:
            raise IntegrityCheckError(
                f"{sample_id}: 输出文件记录数 {seen} 与入库对数 {pairs} 不一致"
            )
        seen_total += seen

    if seen_total != total_pairs:
        raise IntegrityCheckError(
            f"各组记录数合计 {seen_total} 与运行入库总对数 {total_pairs} 不一致"
        )


def verify_published(run_dir, counts, total_pairs: int, evidence,
                     declared_mode=None) -> str:
    """核验一次已发布运行;失败抛 IntegrityCheckError,成功返回核验模式。

    counts: [(样本号, 入库对数), ...](含 undetermined)
    evidence: {文件名: {"size_bytes", "sha256", "records"}, ...}
    declared_mode: 发布时登记的模式;'digest' 的运行若摘要丢失则直接判失败,
                   不能悄悄降级成只能做结构核验。
    返回 "digest"(有内容摘要依据)或 "structure"(旧库无摘要,仅结构核验)。
    """
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise IntegrityCheckError(f"已发布运行目录缺失: {run_dir.name}")
    if not counts:
        raise IntegrityCheckError("运行缺少入库分组计数,无法核对发布文件")
    if not (run_dir / "demux.log").is_file():
        raise IntegrityCheckError("已发布处理日志缺失: demux.log")

    count_sum = sum(pairs for _sid, pairs in counts)
    if count_sum != total_pairs:
        raise IntegrityCheckError(
            f"入库分组对数合计 {count_sum} 与运行总对数 {total_pairs} 不一致"
        )

    expected = {}
    for sample_id, pairs in counts:
        expected[f"{sample_id}_R1.fastq"] = pairs
        expected[f"{sample_id}_R2.fastq"] = pairs

    if evidence:
        missing_evidence = sorted(set(expected) - set(evidence))
        if "demux.log" not in evidence:
            missing_evidence.append("demux.log")
        if missing_evidence:
            raise IntegrityCheckError(
                "发布内容摘要不完整,缺少: " + ", ".join(missing_evidence)
            )
        for filename, pairs in expected.items():
            item = evidence[filename]
            if item["records"] != pairs:
                raise IntegrityCheckError(
                    f"{filename}: 发布登记的记录数 {item['records']} "
                    f"与入库对数 {pairs} 不一致"
                )
            _check_against_evidence(
                run_dir, filename, item["size_bytes"], item["sha256"]
            )
        log_item = evidence["demux.log"]  # 日志同样按发布摘要核验
        _check_against_evidence(
            run_dir, "demux.log", log_item["size_bytes"], log_item["sha256"]
        )
        return "digest"

    # 旧库运行:没有内容摘要,只能做结构与计数核验。
    if declared_mode == "digest":
        raise IntegrityCheckError(
            "该运行发布时登记了内容摘要,但摘要记录缺失,不能降级为仅结构核验"
        )
    _verify_structure(run_dir, counts, total_pairs)
    return "structure"
