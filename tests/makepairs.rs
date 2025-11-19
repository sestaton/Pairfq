use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::NamedTempFile;

mod common;

#[test]
fn test_makepairs_inmemory() {
    let (fq1, fq2) = common::build_fq_data();
    let fp = NamedTempFile::new().unwrap();
    let rp = NamedTempFile::new().unwrap();
    let fs = NamedTempFile::new().unwrap();
    let rs = NamedTempFile::new().unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("makepairs")
        .arg("-f").arg(fq1.path())
        .arg("-r").arg(fq2.path())
        .arg("-p").arg(fp.path())
        .arg("-P").arg(rp.path())
        .arg("-s").arg(fs.path())
        .arg("-S").arg(rs.path())
        .arg("--stats")
        .assert();

    assert
        .success()
        .stdout(predicate::str::contains("Total forward reads").and(predicate::str::contains("8")))
        .stdout(predicate::str::contains("Total reverse reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total forward paired reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total reverse paired reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total forward unpaired reads").and(predicate::str::contains("2")))
        .stdout(predicate::str::contains("Total reverse unpaired reads").and(predicate::str::contains("0")))
        .stdout(predicate::str::contains("Total paired reads").and(predicate::str::contains("12")))
        .stdout(predicate::str::contains("Total unpaired reads").and(predicate::str::contains("2")));
}

#[test]
fn test_makepairs_ondisk() {
    let (fq1, fq2) = common::build_fq_data();
    let fp = NamedTempFile::new().unwrap();
    let rp = NamedTempFile::new().unwrap();
    let fs = NamedTempFile::new().unwrap();
    let rs = NamedTempFile::new().unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("makepairs")
        .arg("-f").arg(fq1.path())
        .arg("-r").arg(fq2.path())
        .arg("-p").arg(fp.path())
        .arg("-P").arg(rp.path())
        .arg("-s").arg(fs.path())
        .arg("-S").arg(rs.path())
        .arg("--index")
        .arg("--stats")
        .assert();

    assert
        .success()
        .stdout(predicate::str::contains("Total forward reads").and(predicate::str::contains("8")))
        .stdout(predicate::str::contains("Total reverse reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total forward paired reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total reverse paired reads").and(predicate::str::contains("6")))
        .stdout(predicate::str::contains("Total forward unpaired reads").and(predicate::str::contains("2")))
        .stdout(predicate::str::contains("Total reverse unpaired reads").and(predicate::str::contains("0")))
        .stdout(predicate::str::contains("Total paired reads").and(predicate::str::contains("12")))
        .stdout(predicate::str::contains("Total unpaired reads").and(predicate::str::contains("2")));
}
