# Distributing Pairfq

This document outlines the steps to make `pairfq` available on various package managers.

## Bioconda

To add `pairfq` to Bioconda:

1.  **Fork the Bioconda recipes repository**: https://github.com/bioconda/bioconda-recipes
2.  **Create a new recipe**:
    *   Create a directory `recipes/pairfq`.
    *   Add a `meta.yaml` file.
    *   Example `meta.yaml` for a Rust project:

```yaml
{% set name = "pairfq" %}
{% set version = "0.1.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://github.com/sestaton/Pairfq/archive/v{{ version }}.tar.gz
  sha256: <SHA256_HASH_OF_RELEASE_TARBALL>

build:
  number: 0

requirements:
  build:
    - {{ compiler('rust') }}
  host:
    - openssl
    - pkg-config
  run:
    - openssl

test:
  commands:
    - pairfq --help

about:
  home: https://github.com/sestaton/Pairfq
  license: MIT
  license_file: LICENSE
  summary: Sync paired-end FASTA/FASTQ files and keep singleton reads
```

3.  **Submit a Pull Request**: Push your changes and submit a PR to the bioconda-recipes repo.

## Cargo

To publish to crates.io:

1.  Ensure `Cargo.toml` has the correct metadata (authors, description, license, repository).
2.  Login to crates.io: `cargo login <API_TOKEN>`
3.  Publish: `cargo publish`

## Homebrew (Optional)

To add to Homebrew:

1.  Create a formula file `pairfq.rb`.
2.  Submit a PR to `homebrew-core` or create a custom tap.
