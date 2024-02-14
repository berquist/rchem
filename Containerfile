ARG RUST_VERSION=1.76
FROM rust:${RUST_VERSION} as builder

# hadolint ignore=DL3008
RUN \
    --mount=type=cache,target=/var/cache/apt \
    apt-get update \
    && apt-get install -y --no-install-recommends libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN cargo install cargo-tarpaulin


WORKDIR /code/rchem
COPY . .
RUN \
    --mount=type=cache,target=/code/rchem/target \
    cargo tarpaulin --workspace --all-features --out Xml Html \
    && cargo build --release
