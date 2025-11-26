# Pairfq Benchmarks and Performance Optimization

This document details the performance optimization journey of porting `Pairfq` from Perl to Rust.

## Summary

The final Rust implementation achieves a **~3.8x speedup** over the legacy Perl script on a dataset of 1 million reads.

| Implementation | Time (Mean) | Speedup |
| :--- | :--- | :--- |
| Perl (`pairfq_lite.pl`) | 9.29s | 1.0x |
| Rust (`pairfq` v1.2.0) | 2.42s | 3.84x |

## Optimization Journey

### 1. Baseline Implementation (SipHash)
The initial Rust port used the standard library's `HashMap`, which uses **SipHash** by default. SipHash is a cryptographically secure hash function designed to prevent HashDoS attacks. While safe, it is relatively slow for non-adversarial use cases like hashing read IDs.

-   **Rust Time**: ~3.38s
-   **Perl Time**: ~9.62s
-   **Speedup**: ~2.85x

### 2. Switching to AHashMap
We switched to [`AHashMap`](https://github.com/tkaitchuck/aHash), which uses AES instruction sets (where available) or a specialized fallback. Since read IDs are not malicious user input, the slight reduction in DoS protection was an acceptable trade-off for significant performance gains.

### 3. Zero-Copy Writes
Profiling revealed that formatting FASTQ records using `format!` was allocating a new `String` on the heap for every single record. We implemented a "zero-copy" writing strategy:
-   Instead of creating intermediate strings, we write byte slices (`&[u8]`) directly to the buffered writer.
-   This eliminates millions of small heap allocations and copies.

### 4. I/O Buffer Tuning
The default buffer size for `BufReader` and `BufWriter` in Rust is typically 8KB. For high-throughput bioinformatics tools processing large files, this results in frequent system calls (`read` and `write`).
-   We increased the buffer size to **64KB**.
-   This reduces the overhead of context switching between user space and kernel space.

### 5. Flexible Pairing Impact
In v1.2.0, we introduced flexible pairing logic (splitting IDs at whitespace) and header retention (storing the full reverse header).
-   **Concern**: Storing full headers increases memory usage and copying.
-   **Result**: Despite the added logic, the performance remained excellent (**2.42s**). The efficiency of `needletail`'s parsing and our zero-copy output pipeline outweighed the cost of handling larger header strings.

## Methodology

Benchmarks were run using [`hyperfine`](https://github.com/sharkdp/hyperfine) on a MacBook Pro (M3 Max).
-   **Data**: 1 million reads (forward) and 900k reads (reverse) to simulate unpaired data.
-   **Warmup**: 3 runs.
-   **Runs**: 10 runs.
