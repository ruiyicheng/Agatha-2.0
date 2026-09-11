FROM python:3.13-slim@sha256:6771159cd4fa5d9bba1258caf0b82e6b73458c694d178ad97c5e925c2d0e1a91
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN apt-get update && apt-get install -y --no-install-recommends poppler-utils tmux procps && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir PyMuPDF==1.28.2 jsonschema==4.26.0
CMD ["sleep", "infinity"]
