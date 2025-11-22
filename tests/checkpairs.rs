use assert_cmd::Command;
use predicates::prelude::*;
use std::io::Write;
use tempfile::NamedTempFile;

mod common;

#[test]
fn test_checkpairs_valid() {
    let content1 = "@seq1/1\nACGT\n+\nIIII\n@seq2/1\nGGGG\n+\nIIII\n";
    let content2 = "@seq1/2\nTGCA\n+\nIIII\n@seq2/2\nCCCC\n+\nIIII\n";

    let fq1 = common::create_fastq_file(content1);
    let fq2 = common::create_fastq_file(content2);

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("checkpairs")
        .arg("-f")
        .arg(fq1.path())
        .arg("-r")
        .arg(fq2.path())
        .assert();

    assert
        .success()
        .stdout(predicate::str::contains("✅").count(4)); // 2 files * (integrity + paired) = 4 checks
}

#[test]
fn test_checkpairs_unpaired() {
    let content1 = "@seq1/1\nACGT\n+\nIIII\n@seq2/1\nGGGG\n+\nIIII\n";
    let content2 = "@seq1/2\nTGCA\n+\nIIII\n"; // Missing second read

    let fq1 = common::create_fastq_file(content1);
    let fq2 = common::create_fastq_file(content2);

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("checkpairs")
        .arg("-f")
        .arg(fq1.path())
        .arg("-r")
        .arg(fq2.path())
        .assert();

    assert
        .success()
        .stdout(predicate::str::contains("✅").count(2)) // Integrity OK for both
        .stdout(predicate::str::contains("❌").count(2)); // Paired Failed for both
}
