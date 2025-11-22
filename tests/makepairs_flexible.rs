use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::NamedTempFile;
use std::io::Write;

mod common;

#[test]
fn test_makepairs_sra_style() {
    // SRA style: Identical IDs, no /1 /2
    let f_content = "@SRR12345 length=10\nACGTACGTAC\n+\nIIIIIIIIII\n";
    let r_content = "@SRR12345 length=10\nTGCATGCATG\n+\nIIIIIIIIII\n";
    
    let fq1 = common::create_fastq_file(f_content);
    let fq2 = common::create_fastq_file(r_content);
    
    let out_dir = tempfile::tempdir().unwrap();
    let fp = out_dir.path().join("fp.fq");
    let rp = out_dir.path().join("rp.fq");
    let fs = out_dir.path().join("fs.fq");
    let rs = out_dir.path().join("rs.fq");

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    let assert = cmd
        .arg("makepairs")
        .arg("-f").arg(fq1.path())
        .arg("-r").arg(fq2.path())
        .arg("-p").arg(&fp)
        .arg("-P").arg(&rp)
        .arg("-s").arg(&fs)
        .arg("-S").arg(&rs)
        .assert();

    assert.success();

    // Check if paired
    let fp_content = std::fs::read_to_string(&fp).unwrap();
    let rp_content = std::fs::read_to_string(&rp).unwrap();
    
    assert!(fp_content.contains("@SRR12345 length=10"), "Forward header preserved");
    // Currently, this might fail if it tries to add /2 or if it fails to match
    // If it matches (because IDs are identical), it might write @SRR12345 (without /2) which is correct for this input.
    // But let's check if it matches at all.
    assert!(fp_content.contains("ACGTACGTAC"));
    assert!(rp_content.contains("TGCATGCATG"));
    
    // Check that we didn't add /2 if it wasn't there (or did we?)
    // The current logic: if F ends with /1, add /2. Else use base_id.
    // Here F is SRR12345. base_id is SRR12345. r_id is SRR12345.
    // So this test might actually PASS with current logic IF get_base_id works on full header.
    // But get_base_id on "SRR12345 length=10" -> "SRR12345 length=10".
    // Matches "SRR12345 length=10".
    // So this specific case passes.
}

#[test]
fn test_makepairs_comments_preservation() {
    // F: @seq1/1 comment1
    // R: @seq1/2 comment2
    // We want R output to have "comment2", not "comment1" or reconstructed "/2".
    
    let f_content = "@seq1/1 comment1\nACGT\n+\nIIII\n";
    let r_content = "@seq1/2 comment2\nTGCA\n+\nIIII\n";
    
    let fq1 = common::create_fastq_file(f_content);
    let fq2 = common::create_fastq_file(r_content);
    
    let out_dir = tempfile::tempdir().unwrap();
    let fp = out_dir.path().join("fp.fq");
    let rp = out_dir.path().join("rp.fq");
    let fs = out_dir.path().join("fs.fq");
    let rs = out_dir.path().join("rs.fq");

    let mut cmd = Command::cargo_bin("pairfq").unwrap();
    cmd.arg("makepairs")
        .arg("-f").arg(fq1.path())
        .arg("-r").arg(fq2.path())
        .arg("-p").arg(&fp)
        .arg("-P").arg(&rp)
        .arg("-s").arg(&fs)
        .arg("-S").arg(&rs)
        .assert()
        .success();

    let rp_content = std::fs::read_to_string(&rp).unwrap();
    
    // Current logic:
    // F id: "seq1/1 comment1". base_id: "seq1/1 comment1" (fails to strip /1 because of comment).
    // R id: "seq1/2 comment2". base_id: "seq1/2 comment2".
    // No match.
    // So they will end up in singles.
    
    // We want them paired.
    // And we want rp to contain "comment2".
    
    assert!(rp_content.contains("comment2"), "Reverse comment preserved");
    assert!(!rp_content.contains("comment1"), "Reverse comment is not forward comment");
}
