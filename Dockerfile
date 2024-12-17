# Stage 1: General Python environment
FROM python:3.12-slim AS base

RUN <<EOT
apt-get update -qy
apt-get install -qyy \
    --no-install-recommends \
    build-essential \
    ca-certificates
EOT

ENV UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=python3.12 \
    UV_PROJECT_ENVIRONMENT=/venv

ENV PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"

# Stage 2: Python dependencies
FROM base AS build

# Copy uv binary from the official repository
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Copy lockfiles for dependency installation
COPY pyproject.toml /_lock/
COPY uv.lock /_lock/

RUN --mount=type=cache,target=/root/.cache <<EOT
cd /_lock
uv sync \
    --locked \
    --no-dev \
    --no-install-project
EOT

# Install the project
COPY . /src
RUN --mount=type=cache,target=/root/.cache \
    uv pip install \
    --python=$UV_PROJECT_ENVIRONMENT \
    --no-deps \
    /src

# Stage 3: Runtime Environment
FROM python:3.12-slim AS runtime

ENV UV_PROJECT_ENVIRONMENT="/venv"
ENV PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"

# Copy the virtual environment from the build stage
COPY --from=build $UV_PROJECT_ENVIRONMENT $UV_PROJECT_ENVIRONMENT

# Set entrypoint for CLI
CMD ["alvoc"]
