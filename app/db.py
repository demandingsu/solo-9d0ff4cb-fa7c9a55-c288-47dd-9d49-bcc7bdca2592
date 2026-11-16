"""SQLite 访问:每次操作独立连接,提交即持久化;WAL 模式支持读写并发。"""
import sqlite3
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  run_number    TEXT NOT NULL UNIQUE,
  status        TEXT NOT NULL CHECK (status IN ('processing', 'published', 'failed')),
  error         TEXT,
  staging_dir   TEXT NOT NULL,
  r1_name       TEXT NOT NULL DEFAULT '',
  r2_name       TEXT NOT NULL DEFAULT '',
  i1_name       TEXT NOT NULL DEFAULT '',
  samples_name  TEXT NOT NULL DEFAULT '',
  total_pairs   INTEGER NOT NULL DEFAULT 0,
  created_at    TEXT NOT NULL,
  updated_at    TEXT NOT NULL,
  published_at  TEXT,
  integrity_mode TEXT,                     -- digest=有内容摘要;structure=旧库仅结构核验
  verified_at    TEXT
);
CREATE TABLE IF NOT EXISTS sample_counts (
  run_id     INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
  sample_id  TEXT NOT NULL,                -- 'undetermined' 表示未分配组
  index_seq  TEXT NOT NULL DEFAULT '',
  pairs      INTEGER NOT NULL,
  corrected  INTEGER NOT NULL DEFAULT 0,   -- 其中经单碱基纠错的对数
  PRIMARY KEY (run_id, sample_id)
);
CREATE TABLE IF NOT EXISTS published_files (
  run_id     INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
  filename   TEXT NOT NULL,                -- 如 S1_R1.fastq
  size_bytes INTEGER NOT NULL,             -- 发布瞬间的字节数
  sha256     TEXT NOT NULL,                -- 发布瞬间的流式内容摘要
  records    INTEGER NOT NULL,             -- 该文件应有的 FASTQ 记录数
  PRIMARY KEY (run_id, filename)
);
"""

# 老库升级:CREATE TABLE IF NOT EXISTS 不会补列,这里按需 ALTER。
_ADDED_COLUMNS = (
    ("runs", "integrity_mode", "TEXT"),
    ("runs", "verified_at", "TEXT"),
)


@contextmanager
def db_conn(path):
    """成功提交、异常回滚、用完关闭的连接上下文。"""
    conn = sqlite3.connect(str(path), timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def init_db(path):
    with db_conn(path) as conn:
        conn.executescript(SCHEMA)
        for table, column, decl in _ADDED_COLUMNS:
            existing = {r["name"] for r in conn.execute(f"PRAGMA table_info({table})")}
            if column not in existing:
                conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {decl}")
        conn.execute("PRAGMA journal_mode = WAL")
