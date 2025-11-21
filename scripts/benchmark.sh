#!/bin/bash
set -e

# Ensure cargo build --release has been run
echo "Building Rust binary..."
cargo build --release

# Check for hyperfine
if ! command -v hyperfine &> /dev/null; then
    echo "Error: hyperfine is not installed. Please install it to run benchmarks."
    echo "  brew install hyperfine  # macOS"
    echo "  cargo install hyperfine # Rust"
    exit 1
fi

# Define paths
RUST_BIN="./target/release/pairfq"
PERL_SCRIPT="./scripts/pairfq_lite.pl"
DATA_DIR="tests/data"
FWD="$DATA_DIR/bench_1.fq.gz"
REV="$DATA_DIR/bench_2.fq.gz"

# Create temp directory for outputs to avoid clutter
OUT_DIR=$(mktemp -d)
trap 'rm -rf "$OUT_DIR"' EXIT

echo "Running benchmarks with:"
echo "  Forward: $FWD"
echo "  Reverse: $REV"
echo "  Output:  $OUT_DIR"
echo ""

hyperfine --warmup 3 --export-markdown benchmark_results.md \
    -n "Rust (pairfq)" \
    "$RUST_BIN makepairs -f $FWD -r $REV -p $OUT_DIR/fp.fq -P $OUT_DIR/rp.fq -s $OUT_DIR/fs.fq -S $OUT_DIR/rs.fq" \
    -n "Perl (pairfq_lite.pl)" \
    "perl $PERL_SCRIPT makepairs -f $FWD -r $REV -fp $OUT_DIR/fp_pl.fq -rp $OUT_DIR/rp_pl.fq -fs $OUT_DIR/fs_pl.fq -rs $OUT_DIR/rs_pl.fq"

echo ""
echo "Benchmark complete. Results saved to benchmark_results.md"
cat benchmark_results.md
