#!/bin/bash
set -e

# Check for seqtk
if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk is not installed."
    exit 1
fi

DATA_DIR="tests/data"
SRC_1="$DATA_DIR/test_1.fq.gz"
SRC_2="$DATA_DIR/test_2.fq.gz"
DEST_1="$DATA_DIR/bench_1.fq.gz"
DEST_2="$DATA_DIR/bench_2.fq.gz"

SEED=100
N_READS=1000000
# Keep 90% of the reads in the second file to simulate missing pairs
# 1000000 reads * 4 lines/read * 0.9 = 3600000 lines
N_LINES_2=3600000 

echo "Generating benchmark data..."

echo "Sampling $N_READS reads from $SRC_1 to $DEST_1..."
seqtk sample -s$SEED $SRC_1 $N_READS | gzip > $DEST_1

echo "Sampling $N_READS reads from $SRC_2 and truncating to 90% to $DEST_2..."
seqtk sample -s$SEED $SRC_2 $N_READS | head -n $N_LINES_2 | gzip > $DEST_2

echo "Done. Created:"
ls -lh $DEST_1 $DEST_2
