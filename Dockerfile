FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

RUN useradd --system --uid 10001 demux \
    && mkdir -p /app/data \
    && chown -R demux:demux /app/data
USER demux

# 所有状态(SQLite、暂存区、发布区)都在这一个数据目录里
ENV DATA_DIR=/app/data

EXPOSE 8000
# 单 worker:后台线程在本进程内处理,多 worker 会与中断恢复逻辑冲突
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
