FROM python:3.11-slim

RUN set -e; \
    apt-get update -y && \
    apt-get install -y curl && \
    curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash -

# npm install -g @openai/codex@v0.1.2504172351 \
# v0.1.2504172351
# v0.1.2504251709

RUN bash -c "\
      . \"/root/.nvm/nvm.sh\" && \
      nvm install 22 && \
      node -v && nvm -v && \
      npm install -g @openai/codex@v0.1.2504251709 \
    "
RUN pip install --no-cache-dir scanpy==1.11.1

WORKDIR /workspace

# ENTRYPOINT ["/bin/bash"]
