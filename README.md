# 🧬 Pairfq

> **Sync paired-end FASTA/FASTQ files and keep singleton reads.** 🚀

[![CI](https://github.com/sestaton/Pairfq/actions/workflows/main.yml/badge.svg)](https://github.com/sestaton/Pairfq/actions/workflows/main.yml)
[![GitHub version](https://badge.fury.io/gh/sestaton%2FPairfq.svg)](https://badge.fury.io/gh/sestaton%2FPairfq)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Pairfq** is a high-performance tool designed to handle paired-end sequencing data. Originally written in Perl, it has been **ported to Rust** 🦀 to provide blazing fast speed, memory safety, and efficient handling of massive datasets.

---

## ✨ Features

*   **🦀 Rust Power**: Built with Rust for maximum performance and reliability.
*   **⚡️ Blazing Fast**: Uses `needletail` for rapid FASTX parsing.
*   **🗄️ Low Memory Footprint**: Optional on-disk indexing with `sled` allows processing of huge datasets (tens of millions of reads) with constant low memory usage.
*   **📦 Zero Dependencies**: The main binary is self-contained (no external DB drivers needed).
*   **🔧 Versatile**: Handles FASTA and FASTQ formats, gzip/bzip2 compression, and interleaved or separate files.

---

## 📦 Installation

### From Source (Rust)

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

`pairfq` provides a suite of subcommands to manipulate your data.

### 1️⃣ `makepairs`
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

### 2️⃣ `joinpairs`
**Interleave paired files.**
Combines separate forward and reverse files into a single interleaved file.

```bash
pairfq joinpairs -f forward.fastq -r reverse.fastq -o interleaved.fastq
```

### 3️⃣ `splitpairs`
**De-interleave files.**
Splits a single interleaved file back into separate forward and reverse files.

```bash
pairfq splitpairs -i interleaved.fastq -f forward.fastq -r reverse.fastq
```

### 4️⃣ `addinfo`
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

---

## 📄 License

This project is licensed under the **MIT License**.

Copyright (C) 2013-2025 S. Evan Staton
