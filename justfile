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

# Remove build artifacts.
clean:
    {{cargo}} clean

# NOTE: the alignment benchmark (flexalign vs strobealign, correctness + speed) lives in its own
# repo — ../flexalign_benchmark — managed with pixi. Run `pixi run setup` then `pixi run bench`
# there. This crate only provides the aligner binaries it drives (flexalign, examples/shard_align).
