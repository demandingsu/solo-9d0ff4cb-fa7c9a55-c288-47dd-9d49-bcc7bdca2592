"""逐对流式拆分核心:一次只处理一条记录,不整份读入内存。

输入为三份未压缩的 4 行 FASTQ(R1/R2/I1)。任何结构损坏、读段号不一致、
序列与质量长度不符都会抛 FastqFormatError,由上层把运行标记为失败,
保证不发布半套文件。
"""
from __future__ import annotations

import os
from itertools import zip_longest

BASES = "ACGT"
_MISSING = object()  # zip_longest 哨兵:某份文件先读完


class FastqFormatError(Exception):
    """FASTQ 结构损坏或三份文件读对不一致。"""


def read_id(header: str) -> str:
    """读段号 = 头部第一个空白之前的部分,如 '@INST:1:FC:1:1101:100:100 1:N:0:0'。"""
    parts = header.split(None, 1)
    return parts[0] if parts else ""


def iter_fastq(path, label):
    """按 4 行一条产出 (header, seq, plus, qual);结构损坏即抛错。"""
    rec = 0
    with open(path, "r", encoding="ascii", newline="") as fh:
        while True:
            try:
                block = [fh.readline() for _ in range(4)]
            except UnicodeDecodeError as exc:
                raise FastqFormatError(
                    f"{label}: 第 {rec + 1} 条附近不是纯文本(需要未压缩 FASTQ): {exc}"
                ) from exc
            if block[0] == "":
                return  # 恰好结束在记录边界
            rec += 1
            if block[1] == "" or block[2] == "" or block[3] == "":
                raise FastqFormatError(f"{label}: 第 {rec} 条记录在文件末尾被截断")
            header, seq, plus, qual = (line.rstrip("\r\n") for line in block)
            if not header.startswith("@"):
                raise FastqFormatError(f"{label}: 第 {rec} 条记录头缺少 '@'")
            if not plus.startswith("+"):
                raise FastqFormatError(f"{label}: 第 {rec} 条记录分隔行缺少 '+'")
            if len(seq) != len(qual):
                raise FastqFormatError(
                    f"{label}: 第 {rec} 条({header})序列长度 {len(seq)} 与质量长度 {len(qual)} 不一致"
                )
            yield header[1:], seq, plus, qual


def build_matcher(samples):
    """预计算精确匹配表与单碱基错配候选表(错配串 -> 候选样本集合)。"""
    exact = {}
    neighbors: dict[str, set[str]] = {}
    for sid, index in samples:
        exact[index] = sid
    for sid, index in samples:
        for pos in range(len(index)):
            for base in BASES:
                if base != index[pos]:
                    mut = index[:pos] + base + index[pos + 1 :]
                    neighbors.setdefault(mut, set()).add(sid)
    return exact, neighbors


def assign_index(index, exact, neighbors):
    """精确匹配优先;否则仅当单碱基错配候选唯一时纠错;其余返回未分配。"""
    sid = exact.get(index)
    if sid is not None:
        return sid, False
    cands = neighbors.get(index)
    if cands is not None and len(cands) == 1:
        return next(iter(cands)), True
    return None, False


def demultiplex(r1_path, r2_path, i1_path, samples, out_dir):
    """逐对拆分。samples 为 [(样本号, 索引序列)],输出写入 out_dir。

    返回 {'total_pairs', 'counts': {组: [对数, 纠错对数]}, 'order': [组...]}。
    输出文件按输入顺序追加,保证各组内读段顺序与原始一致。
    """
    exact, neighbors = build_matcher(samples)
    order = [sid for sid, _ in samples] + ["undetermined"]
    writers = {}
    counts = {sid: [0, 0] for sid in order}
    try:
        for sid in order:
            writers[sid] = (
                open(os.path.join(out_dir, f"{sid}_R1.fastq"), "w", encoding="ascii", buffering=1 << 20),
                open(os.path.join(out_dir, f"{sid}_R2.fastq"), "w", encoding="ascii", buffering=1 << 20),
            )
        total = 0
        streams = zip_longest(
            iter_fastq(r1_path, "R1"),
            iter_fastq(r2_path, "R2"),
            iter_fastq(i1_path, "I1"),
            fillvalue=_MISSING,
        )
        for rec1, rec2, reci in streams:
            if rec1 is _MISSING or rec2 is _MISSING or reci is _MISSING:
                raise FastqFormatError(f"第 {total + 1} 条起,三份文件的记录数不一致")
            total += 1
            h1, s1, p1, q1 = rec1
            h2, s2, p2, q2 = rec2
            hi, si, _pi, _qi = reci
            id1, id2, idi = read_id(h1), read_id(h2), read_id(hi)
            if not id1 or id1 != id2 or id1 != idi:
                raise FastqFormatError(
                    f"第 {total} 条读段号不一致: R1={h1!r} R2={h2!r} I1={hi!r}"
                )
            sid, corrected = assign_index(si.upper(), exact, neighbors)
            key = sid if sid is not None else "undetermined"
            w1, w2 = writers[key]
            w1.write(f"@{h1}\n{s1}\n{p1}\n{q1}\n")
            w2.write(f"@{h2}\n{s2}\n{p2}\n{q2}\n")
            counts[key][0] += 1
            if corrected:
                counts[key][1] += 1
    finally:
        for pair in writers.values():
            for fh in pair:
                fh.close()
    return {"total_pairs": total, "counts": counts, "order": order}
