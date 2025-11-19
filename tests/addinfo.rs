use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::NamedTempFile;
use std::io::{Write, Read};

mod common;

#[test]
fn test_addinfo() {
    let content = "\
@seq1
ACGT
+
IIII
";
    let mut infile = NamedTempFile::new().unwrap();
    write!(infile, "{}", content).unwrap();
    let outfile = NamedTempFile::new().unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("addinfo")
        .arg("-i").arg(infile.path())
        .arg("-o").arg(outfile.path())
        .arg("-p").arg("1")
        .assert();

    assert.success();

    let mut out_content = String::new();
    std::fs::File::open(outfile.path()).unwrap().read_to_string(&mut out_content).unwrap();
    assert!(out_content.contains("@seq1/1"));
}

#[test]
fn test_addinfo_stdout() {
    let content = "\
@seq1
ACGT
+
IIII
";
    let mut infile = NamedTempFile::new().unwrap();
    write!(infile, "{}", content).unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("addinfo")
        .arg("-i").arg(infile.path())
        .arg("-o").arg("-")
        .arg("-p").arg("2")
        .assert();

    assert
        .success()
        .stdout(predicate::str::contains("@seq1/2"));
}
