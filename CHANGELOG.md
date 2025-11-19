# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-11-19

### Added
- Ported the entire application from Perl to Rust for improved performance and maintainability.
- Added `sled` dependency for high-performance, embedded key-value storage (replaces SQLite).
- Added `needletail` for fast FASTX parsing.
- Added `clap` for robust CLI argument parsing.
- Added `flate2` and `bzip2` for native compression support.
- Added `env_logger` for configurable logging.
- Added GitHub Actions CI workflow for Rust.

### Changed
- Replaced `DBD::SQLite` with `sled` for the `--index` option, removing the need for external database drivers.
- Updated `makepairs` to use streaming processing by default, with optional on-disk indexing.
- Updated build system to use `cargo` instead of `Makefile.PL`.

### Removed
- Removed Perl dependencies (`DBI`, `DBD::SQLite`, `Tie::Hash::DBD`, `Getopt::Long`, etc.).
- Removed `Makefile.PL` and `MANIFEST` related to Perl distribution.

## [0.18.0] - 2024-08-13

### Changed
- Allow pairing without having pair info in the name for the 'makepairs' command.
- Expand FASTQ header formats for working with SRA data that may contain pair info at the end of the comment.

## [0.17.0] - 2017-03-07

### Fixed
- Bug fix for handling opened processes (named pipes or piped processes).
- Bug fix for handling options to main program (help/man/version options).

## [0.16.1] - 2016-12-09

### Changed
- Modify how options are processed in main application.

## [0.16.0] - 2016-07-24

### Fixed
- Clean up Unix paths correctly (fixes #10).

## [0.15.0] - 2015-12-22

### Changed
- Major change in how compressing output is performed (piped directly to compression program).
- Improved file extension handling.
- Reduced Perl dependencies for compression.

## [0.14.7] - 2015-11-18

### Fixed
- Fix bug in compression method.

## [0.14.6] - 2015-10-03

### Changed
- Use dispatch tables in main application.

## [0.14.5] - 2015-08-08

### Changed
- Allow build system to find Perl instead of using environment variables.

## [0.14.4] - 2015-06-10

### Added
- Add option to write to STDOUT for 'addinfo' and 'joinpairs' commands.
- Allow to read from stdin with a '-' or case-insensitive 'stdin'.

## [0.14.3] - 2015-05-11

### Added
- Add conditional to check if comment is missing after pair id with Illumina 1.8+ identifiers.

## [0.14.2] - 2015-04-18

### Fixed
- Fix logging script name when reading from stdin.

## [0.14.1] - 2014-11-18

### Removed
- Remove dependencies used only for testing (List::MoreUtils and IPC::System::Simple).

## [0.14] - 2014-10-27

### Changed
- Replace BerkeleyDB indexing method with SQLite.

## [0.13.2] - 2014-10-13

### Fixed
- Fix bug in 'pairfq_lite.pl' script (temp file creation).

## [0.13.1] - 2014-08-15

### Added
- Add DOI to make repo citable.

## [0.13] - 2014-08-11

### Added
- Add method for finding pairs from an interleaved file.

## [0.12] - 2014-06-03

### Changed
- Add method to not load modules until they are required.
- Remove List::Util dependency from main application.

## [0.11] - 2014-04-14

### Fixed
- Bug fix for reading files containing a dash.

## [0.10] - 2014-03-20

### Added
- Add standalone script with no dependencies in `/scripts`.

### Changed
- Improve doc readability.
- Modify standalone script to work with v5.6.
- Fix help message option listings.

## [0.09.1] - 2014-03-18

### Added
- Add check for ExtUtils::MakeMaker version.

### Changed
- Remove tying objects for dbm files.

## [0.09] - 2014-02-20

### Changed
- Rework mk_key and mk_vec functions for storing data (no UTF-8).
- Indexing method uses less memory and is faster.
- Lowered version requirement to 5.10.
- Make printing stats for 'makepairs' method an option.

## [0.08] - 2014-02-10

### Fixed
- Bugfix for encoding read names with a comment line.

### Changed
- Update README to point to wiki.

## [0.07] - 2014-01-27

### Removed
- Remove use of DB_File.

## [0.06] - 2014-01-25

### Added
- Added support for using BerkeleyDB if available.

## [0.05] - 2014-01-22

### Changed
- Make in-memory processing default.
- Make --index an option.
- Add explicit exit or return from subroutines.
- More sig handlers for indexing method.

## [0.04] - 2013-11-18

### Fixed
- Bug fix for taking positional arguments.

### Added
- Add option to get version.

## [0.03] - 2013-11-18

### Fixed
- Bug fixes for splitting pairs, storing files in memory, and naming DBM files.

## [0.02] - 2013-11-11

### Changed
- Removed all utility scripts and rolled methods into a single executable.
- Make individual methods positional arguments.

## [0.01] - 2013-11-01

### Added
- Initial release.
