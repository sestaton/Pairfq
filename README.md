# Pairfq

Sync paired-end FASTA/FASTQ files and keep singleton reads.

## Build Status

[![CI](https://github.com/sestaton/Pairfq/actions/workflows/main.yml/badge.svg)](https://github.com/sestaton/Pairfq/actions/workflows/main.yml) | [![GitHub version](https://badge.fury.io/gh/sestaton%2FPairfq.svg)](https://badge.fury.io/gh/sestaton%2FPairfq)

## Installation

### From Source (Rust)

To build and install `pairfq` from source, you need to have Rust installed. You can install Rust using [rustup](https://rustup.rs/).

```bash
git clone https://github.com/sestaton/Pairfq.git
cd Pairfq
cargo build --release
```

The binary will be located in `target/release/pairfq`. You can copy it to a directory in your PATH, e.g.:

```bash
cp target/release/pairfq /usr/local/bin/
```

## Usage

`pairfq` provides several subcommands to manipulate paired-end data.

### `makepairs`

Syncs paired-end reads from two separate files or an interleaved file.

```bash
pairfq makepairs -f forward.fastq -r reverse.fastq -fp forward_paired.fastq -rp reverse_paired.fastq -fs forward_unpaired.fastq -rs reverse_unpaired.fastq
```

**Options:**
*   `-f, --forward`: Input forward reads.
*   `-r, --reverse`: Input reverse reads.
*   `-i, --infile`: Input interleaved file (alternative to -f/-r).
*   `-p, --forw_paired`: Output paired forward reads.
*   `-P, --rev_paired`: Output paired reverse reads.
*   `-s, --forw_unpaired`: Output unpaired forward reads.
*   `-S, --rev_unpaired`: Output unpaired reverse reads.
*   `--index`: Use on-disk indexing for large files. This uses `sled`, a high-performance embedded database, to store reads on disk instead of in memory. This allows processing of massive datasets (tens of millions of reads) with constant low memory usage, similar to the original Perl version's SQLite implementation but faster and without external dependencies.
*   `--stats`: Print statistics.

### `joinpairs`

Interleaves two paired-end files into a single file.

```bash
pairfq joinpairs -f forward.fastq -r reverse.fastq -o interleaved.fastq
```

### `splitpairs`

Splits an interleaved file into two separate files.

```bash
pairfq splitpairs -i interleaved.fastq -f forward.fastq -r reverse.fastq
```

### `addinfo`

Adds pairing information (e.g., `/1`, `/2`) to read headers.

```bash
pairfq addinfo -i input.fastq -o output.fastq -p 1
```

## For Developers

### Building

```bash
cargo build
```

### Testing

Run the test suite:

```bash
cargo test
```

The test suite includes integration tests that verify the functionality of all subcommands against the expected behavior defined in the original Perl test suite.

## Legacy Lite Script

There is a standalone script in the 'scripts' directory that has no dependencies and will work with Perl version 5.6 or newer. This script has fewer features (mainly, it lacks the indexing function for working with large data) than the main application but it may be useful in an environment where installing libraries is not convenient. Obtaining this version can be done with curl:

    curl -sL git.io/pairfq_lite > pairfq_lite

You can then make the script executable and check the usage:

    chmod +x pairfq_lite
    ./pairfq_lite -h

Alternatively, you can use this version without storing it locally.

    curl -sL git.io/pairfq_lite | perl -

The above command will show the options. To see a specific subcommand menu, for example the `makepairs` command, just type that subcommand with no options.

    curl -sL git.io/pairfq_lite | perl - makepairs

## License

The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013-2025 S. Evan Staton
