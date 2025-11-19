use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::NamedTempFile;
use std::io::Read;

mod common;

#[test]
fn test_joinpairs() {
    // Create paired data (2 pairs)
    let fq1_content = "\
@seq1/1
ACGT
+
IIII
@seq2/1
GGGG
+
IIII
";
    let fq2_content = "\
@seq1/2
TGCA
+
IIII
@seq2/2
CCCC
+
IIII
";
    let fq1 = common::create_fastq_file(fq1_content);
    let fq2 = common::create_fastq_file(fq2_content);
    let outfile = NamedTempFile::new().unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("joinpairs")
        .arg("-f").arg(fq1.path())
        .arg("-r").arg(fq2.path())
        .arg("-o").arg(outfile.path())
        .assert();

    assert.success();

    // 2 pairs * 2 reads/pair * 4 lines/read = 16 lines
    let mut file = std::fs::File::open(outfile.path()).unwrap();
    let mut content = String::new();
    file.read_to_string(&mut content).unwrap();
    let line_count = content.lines().count();
    
    assert_eq!(line_count, 16);
}
