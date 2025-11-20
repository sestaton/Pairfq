# TODO

- [ ] Add 'checkpairs' command to check integrity and pairing of input files. Good validation to have before running large-scale projects.
- [ ] ~~Refactor code to use Moo and App::Cmd for easier extension and development.~~ (Not relevant for Rust port)
- [x] Use low-level, faster interface for handling of FASTQ data. (Done - ported to Rust in v1.0.0)
- [ ] Add docker build.
- [ ] Add nextflow recipe so users can easily run utilities at scale across many files or machines.