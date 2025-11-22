# 🧬 Pairfq

> **Sync paired-end FASTA/FASTQ files and keep singleton reads.** 🚀

[![CI](https://github.com/sestaton/Pairfq/actions/workflows/main.yml/badge.svg)](https://github.com/sestaton/Pairfq/actions/workflows/main.yml)
[![GitHub version](https://badge.fury.io/gh/sestaton%2FPairfq.svg)](https://badge.fury.io/gh/sestaton%2FPairfq)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Pairfq** is a high-performance tool designed to handle paired-end sequencing data. It provides blazing fast speed, memory safety, and efficient handling of massive datasets.

---

## ✨ Features

*   **⚡️ Blazing Fast**: Uses the Rust library `needletail` for rapid FASTX parsing.
*   **🗄️ Low Memory Footprint**: Optional on-disk indexing with `sled` allows processing of huge datasets (tens of millions of reads) with constant low memory usage.
*   **📦 Zero Dependencies**: The main binary is self-contained (no external DB drivers needed).
*   **🔧 Versatile**: Handles FASTA and FASTQ formats, gzip/bzip2 compression, and interleaved or separate files.

---

## 📦 Installation

### From Source

You need to have [Rust installed](https://rustup.rs/).

```bash
# Clone the repository
git clone https://github.com/sestaton/Pairfq.git
cd Pairfq

# Build for release
cargo build --release

# Install (optional)
cp target/release/pairfq /usr/local/bin/
```

---

## 🚀 Usage

`pairfq` provides a suite of subcommands to deal with paired-end FASTA/FASTQ files.

### `makepairs`
**Sync paired-end reads.**
Matches forward and reverse reads, keeping them in sync and separating singletons.

```bash
pairfq makepairs \
  -f forward.fastq -r reverse.fastq \
  -fp forward_paired.fastq -rp reverse_paired.fastq \
  -fs forward_unpaired.fastq -rs reverse_unpaired.fastq
```

**Key Options:**
*   `--index`: **Recommended for large files!** Uses `sled` (embedded DB) to index reads on disk, keeping memory usage low. 📉
*   `--stats`: Print detailed statistics after processing. 📊

### `joinpairs`
**Interleave paired files.**
Combines separate forward and reverse files into a single interleaved file.

```bash
pairfq joinpairs -f forward.fastq -r reverse.fastq -o interleaved.fastq
```

### `splitpairs`
**De-interleave files.**
Splits a single interleaved file back into separate forward and reverse files.

```bash
pairfq splitpairs -i interleaved.fastq -f forward.fastq -r reverse.fastq
```

### `checkpairs`

Check the integrity and pairing of forward and reverse files.

```bash
pairfq checkpairs -f <forward_reads> -r <reverse_reads>
```

**Output:**
A tab-delimited table showing the status of each file:
- **integrity**: Checks if the file can be parsed (validates gzip/bzip2 compression if applicable).
- **paired**: Checks if the file has the same number of records as its pair.
- **paired_reads**: Count of reads that are paired.
- **unpaired_reads**: Count of reads that are unpaired (difference in counts).

Example:
```
file	integrity	paired	paired_reads	unpaired_reads
file1.fq	✅	✅	100	0
file2.fq	✅	✅	100	0
```

### `addinfo`
**Fix headers.**
Adds standard pairing information (e.g., `/1`, `/2`) to read headers.

```bash
pairfq addinfo -i input.fastq -o output.fastq -p 1
```

---

## 🛠️ For Developers

Want to contribute? Great!

```bash
# Run the test suite (includes ported Perl tests)
cargo test
```

---

## 📜 Legacy Lite Script

For environments where you cannot install the Rust binary, we preserve the legacy Perl script. It has **no dependencies** and works with Perl 5.6+.

> ⚠️ **Note:** This version lacks the high-performance indexing of the Rust version.

**Quick Install:**
```bash
curl -sL git.io/pairfq_lite > pairfq_lite
chmod +x pairfq_lite
./pairfq_lite -h
```

Alternatively, you can use this version without storing it locally.

```bash
    curl -sL git.io/pairfq_lite | perl -
```

The above command will show the options. To see a specific subcommand menu, for example the `pairfq makepairs` command, just type that subcommand with no options.

```bash
    curl -sL git.io/pairfq_lite | perl - makepairs
```

However, repeatedly running the above command with `curl` is not efficient. You can save it to a file and make it executable as shown above for the quick install.

---

## 📊 Benchmark Results

Pairfq is designed to be efficient and portable. Here are some benchmark results with hyperfine comparing Pairfq to the legacy Perl implementation with a test set of 1M reads in the R1 file and 900k reads in the R2 file. For transparency, I have included the scripts to create the test data and run the benchmarks in the `scripts/` directory.

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `Rust (pairfq)` | 2.702 ± 0.084 | 2.610 | 2.891 | 1.00 |
| `Perl (pairfq_lite.pl)` | 9.806 ± 0.326 | 9.461 | 10.446 | 3.63 ± 0.17 |

---


## 📄 License

This project is licensed under the **MIT License**.

Copyright (C) 2013-2025 S. Evan Staton
