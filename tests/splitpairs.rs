use assert_cmd::Command;

use tempfile::NamedTempFile;
use std::io::{Write, Read};

mod common;

#[test]
fn test_splitpairs() {
    // Create interleaved data
    // 2 pairs
    let content = "\
@seq1/1
ACGT
+
IIII
@seq1/2
TGCA
+
IIII
@seq2/1
GGGG
+
IIII
@seq2/2
CCCC
+
IIII
";
    let mut infile = NamedTempFile::new().unwrap();
    write!(infile, "{}", content).unwrap();
    
    let fwd = NamedTempFile::new().unwrap();
    let rev = NamedTempFile::new().unwrap();

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("splitpairs")
        .arg("-i").arg(infile.path())
        .arg("-f").arg(fwd.path())
        .arg("-r").arg(rev.path())
        .assert();

    assert.success();

    let mut f_content = String::new();
    std::fs::File::open(fwd.path()).unwrap().read_to_string(&mut f_content).unwrap();
    assert!(f_content.contains("@seq1/1"));
    assert!(f_content.contains("@seq2/1"));
    assert!(!f_content.contains("@seq1/2"));

    let mut r_content = String::new();
    std::fs::File::open(rev.path()).unwrap().read_to_string(&mut r_content).unwrap();
    assert!(r_content.contains("@seq1/2"));
    assert!(r_content.contains("@seq2/2"));
    assert!(!r_content.contains("@seq1/1"));
}
