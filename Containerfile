FROM rust:1.76 as builder

LABEL org.opencontainers.image.source=https://github.com/berquist/rchem
LABEL org.opencontainers.image.description="quantum chemistry in Rust"
LABEL org.opencontainers.image.licenses=BSD-3-Clause

# hadolint ignore=DL3008
RUN \
    --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update \
    && apt-get install -y --no-install-recommends \
      cmake \
      gfortran \
      libclang-dev \
      libgsl-dev \
      liblapacke-dev \
    && rm -rf /var/lib/apt/lists/*

ARG PYENV_ROOT=/pyenv
ENV PYENV_ROOT "${PYENV_ROOT}"
ENV PATH "${PYENV_ROOT}/shims:${PYENV_ROOT}/bin:${PATH}"
ARG PYTHON_VERSION=3.11.8
RUN git clone https://github.com/pyenv/pyenv.git "${PYENV_ROOT}" \
    && pyenv install "${PYTHON_VERSION}" \
    && pyenv rehash \
    && pyenv global "${PYTHON_VERSION}"

WORKDIR /code/rchem
COPY requirements.txt requirements.txt
# hadolint ignore=DL3013
RUN python -m venv ${HOME}/venv \
    && . "${HOME}"/venv/bin/activate \
    && python -m pip install --no-cache-dir -U pip setuptools \
    && python -m pip install --no-cache-dir -r requirements.txt

RUN cargo install cargo-tarpaulin

COPY . .

ENV LD_LIBRARY_PATH "${PYENV_ROOT}/versions/${PYTHON_VERSION}/lib:${LD_LIBRARY_PATH}"

RUN \
    --mount=type=cache,target=/code/rchem/target,sharing=locked \
    cargo tarpaulin --workspace --all-features --out Xml Html \
    && cargo build --release
