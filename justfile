# flexalign task runner — `just <recipe>` (run `just` with no args to list).
# This crate requires the Rust nightly toolchain (unstable features + the
# `profile-rustflags` cargo feature), so every cargo call goes through nightly.

cargo := "cargo +nightly"

# Show available recipes.
default:
    @just --list

# Build the optimized release binary (target/release/flexalign).
build:
    {{cargo}} build --release

# Build the debug binary (faster compile, slower binary).
build-debug:
    {{cargo}} build

# Run unit tests (the #[test] modules in src/).
unit:
    {{cargo}} test --lib

# Integration/regression check: align a fixed read subset, diff vs committed baseline.
integration: build
    ./tests/regression/run.sh

# Run everything: unit tests then the integration check.
test: unit integration

# Refresh the baseline (archives the old one + logs it). Add a note: `just regression-update "why"`.
regression-update NOTE="": build
    ./tests/regression/run.sh --update "{{NOTE}}"

# Run the aligner, passing args through, e.g. `just run -r ref.fna -1 r1.fq -2 r2.fq`.
run *ARGS:
    {{cargo}} run --release -- {{ARGS}}

# Lint with clippy.
lint:
    {{cargo}} clippy --release

# Install the flexalign binary to ~/.cargo/bin.
install:
    {{cargo}} install --path .

# Install flexalign on a cluster, from a clean checkout, into PREFIX (default ~/bin).
#
# Builds against the PINNED GitHub dependencies -- it moves .cargo/config.toml aside first, because
# that file redirects bioreader/kmerrs/flexmap to sibling working copies that exist only on a
# developer machine. Without that step a cluster build either fails or silently builds different
# code than the pinned revisions.
#
# Needs: rustup (for nightly) and a host C toolchain (cc, make, cmake) for libwfa2. Load those with
# `module load` first if your cluster provides them that way.
#
#   just install-cluster                 # -> ~/bin
#   just install-cluster ~/opt/flexalign # -> that prefix
install-cluster PREFIX="~/bin":
    #!/usr/bin/env bash
    set -euo pipefail
    command -v rustup >/dev/null || { echo "rustup not found -- install from https://rustup.rs"; exit 1; }
    for t in cc make cmake; do command -v $t >/dev/null || { echo "missing $t (needed by libwfa2)"; exit 1; }; done
    rustup toolchain install nightly --profile minimal
    # A local path override would defeat the point of a reproducible cluster build.
    if [ -f .cargo/config.toml ]; then
        echo ">> setting aside .cargo/config.toml (local path override) for this build"
        mv .cargo/config.toml .cargo/config.toml.clusterbuild
        trap 'mv .cargo/config.toml.clusterbuild .cargo/config.toml 2>/dev/null || true' EXIT
    fi
    {{cargo}} build --release
    PREFIX=$(eval echo {{PREFIX}}); mkdir -p "$PREFIX"
    install -m755 target/release/flexalign "$PREFIX/flexalign"
    echo ">> installed to $PREFIX:"; "$PREFIX/flexalign" --version || true
    echo ">> add to PATH:  export PATH=\"$PREFIX:\$PATH\""

# Verify the pinned dependencies resolve from GitHub (no local path override). Fast, no compile.
check-deps:
    #!/usr/bin/env bash
    set -euo pipefail
    [ -f .cargo/config.toml ] && { mv .cargo/config.toml .cargo/config.toml.depcheck; trap 'mv .cargo/config.toml.depcheck .cargo/config.toml' EXIT; }
    CARGO_TARGET_DIR=$(mktemp -d) {{cargo}} metadata --format-version 1 >/dev/null
    echo "pinned dependencies resolve from GitHub"

# Remove build artifacts.
clean:
    {{cargo}} clean

# NOTE: the alignment benchmark (flexalign vs strobealign, correctness + speed) lives in its own
# repo — ../flexalign_benchmark — managed with pixi. Run `pixi run setup` then `pixi run bench`
# there. This crate only provides the aligner binary it drives (flexalign; sharded runs are the
# same binary, which detects a sliced index from the reference path).
