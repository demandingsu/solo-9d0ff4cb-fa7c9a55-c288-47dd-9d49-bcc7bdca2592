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
  published_at  TEXT
);
CREATE TABLE IF NOT EXISTS sample_counts (
  run_id     INTEGER NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
  sample_id  TEXT NOT NULL,                -- 'undetermined' 表示未分配组
  index_seq  TEXT NOT NULL DEFAULT '',
  pairs      INTEGER NOT NULL,
  corrected  INTEGER NOT NULL DEFAULT 0,   -- 其中经单碱基纠错的对数
  PRIMARY KEY (run_id, sample_id)
);
"""


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
        conn.execute("PRAGMA journal_mode = WAL")
