"""样本索引 CSV 解析与校验:两列(样本号,索引序列),索引等长且唯一。"""
import csv
import re


class SampleSheetError(Exception):
    """样本表内容不合法。"""


SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
INDEX_RE = re.compile(r"^[ACGTNacgtn]+$")
RESERVED = {"undetermined"}  # 未分配组的固定名,样本号不得占用


def load_samples(path):
    """返回 [(样本号, 大写索引序列)];首行若不像数据行则视为表头跳过。"""
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for raw in csv.reader(fh):
            cells = [c.strip() for c in raw]
            while cells and cells[-1] == "":
                cells.pop()  # 容忍行尾多余逗号
            if cells:
                rows.append(cells)
    if rows:
        first = rows[0]
        looks_like_data = len(first) == 2 and bool(INDEX_RE.fullmatch(first[1]))
        if not looks_like_data:
            rows = rows[1:]
    if not rows:
        raise SampleSheetError("样本表为空或只有表头")

    samples, seen_ids, seen_idx = [], set(), set()
    idx_len = None
    for lineno, cells in enumerate(rows, start=1):
        if len(cells) != 2:
            raise SampleSheetError(f"第 {lineno} 行应为两列:样本号,索引序列")
        sid, index = cells
        if not SAMPLE_ID_RE.fullmatch(sid):
            raise SampleSheetError(f"第 {lineno} 行样本号无效: {sid!r}")
        if sid.lower() in RESERVED:
            raise SampleSheetError(f"样本号 {sid!r} 为保留字,请改用其他名字")
        if not INDEX_RE.fullmatch(index):
            raise SampleSheetError(f"第 {lineno} 行索引序列含非法字符(仅允许 ACGTN): {index!r}")
        index = index.upper()
        if idx_len is None:
            idx_len = len(index)
        elif len(index) != idx_len:
            raise SampleSheetError(f"索引序列长度不一致: {sid} 为 {len(index)},应为 {idx_len}")
        if sid in seen_ids:
            raise SampleSheetError(f"样本号重复: {sid}")
        if index in seen_idx:
            raise SampleSheetError(f"索引序列重复: {index}")
        seen_ids.add(sid)
        seen_idx.add(index)
        samples.append((sid, index))
    return samples
