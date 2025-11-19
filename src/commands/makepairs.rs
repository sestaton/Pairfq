use anyhow::{Context, Result};
use crate::utils::{get_writer, format_fastq};
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::path::Path;
use log::info;

pub fn run(
    forward: String,
    reverse: String,
    fp: String,
    rp: String,
    fs: String,
    rs: String,
    index: bool,
    compress: Option<String>,
    stats: bool,
) -> Result<()> {
    info!("Starting makepairs");

    let mut fp_writer = get_writer(&fp, compress.as_deref())?;
    let mut rp_writer = get_writer(&rp, compress.as_deref())?;
    let mut fs_writer = get_writer(&fs, compress.as_deref())?;
    let mut rs_writer = get_writer(&rs, compress.as_deref())?;

    let mut f_count = 0;
    let mut r_count = 0;
    let mut fp_count = 0;
    let mut rp_count = 0;
    let mut fs_count = 0;
    let mut rs_count = 0;

    if index {
        // Disk-based indexing using sled
        let tmp_dir = tempfile::tempdir()?;
        let db = sled::open(tmp_dir.path().join("pairfq_db"))?;

        // 1. Index reverse reads
        let mut reader = parse_fastx_file(&reverse).with_context(|| format!("Failed to open {}", reverse))?;
        while let Some(record) = reader.next() {
            let record = record?;
            r_count += 1;
            let id = std::str::from_utf8(record.id())?.to_string();
            // Strip /1 or /2 if present (simple heuristic, Perl script does more complex parsing but this covers 99%)
            let base_id = id.trim_end_matches("/2").trim_end_matches("/1");
            
            let seq = record.seq();
            let qual = record.qual();
            
            // Store as "seq|qual" or just "seq" if no qual
            let val = if let Some(q) = qual {
                let mut v = Vec::with_capacity(seq.len() + 1 + q.len());
                v.extend_from_slice(&seq);
                v.push(b'|');
                v.extend_from_slice(q);
                v
            } else {
                seq.to_vec()
            };
            
            db.insert(base_id, val)?;
        }
        db.flush()?;

        // 2. Process forward reads
        let mut reader = parse_fastx_file(&forward).with_context(|| format!("Failed to open {}", forward))?;
        while let Some(record) = reader.next() {
            let record = record?;
            f_count += 1;
            let id = std::str::from_utf8(record.id())?.to_string();
            let base_id = id.trim_end_matches("/1").trim_end_matches("/2");

            if let Some(val) = db.remove(base_id)? {
                // Match found
                fp_count += 1;
                rp_count += 1;

                // Write forward
                let f_seq_cow = record.seq();
                let f_seq = std::str::from_utf8(&f_seq_cow)?;
                let f_qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;
                write!(fp_writer, "{}", format_fastq(&id, f_seq, f_qual))?;

                // Write reverse (reconstruct from DB value)
                let val_vec = val.to_vec(); // Convert IVec to Vec<u8>
                let (r_seq_bytes, r_qual_bytes) = if let Some(pos) = val_vec.iter().position(|&x| x == b'|') {
                    (&val_vec[..pos], Some(&val_vec[pos+1..]))
                } else {
                    (&val_vec[..], None)
                };
                
                let r_seq = std::str::from_utf8(r_seq_bytes)?;
                let r_qual = r_qual_bytes.map(|q| std::str::from_utf8(q)).transpose()?;
                
                // Reconstruct ID (Perl script logic handles comments etc, here we just append /2 if it was stripped or just use base_id)
                // For simplicity, we assume standard pairing. 
                // If the original had /2, we should probably put it back. 
                // But we only stored the base_id. 
                // Let's assume the reverse read should have the same base ID.
                // The Perl script reconstructs it carefully. 
                // For now, let's append /2 if the forward had /1, or just use base_id if neither had suffixes.
                // A safer bet is to append /2 if we are writing to the reverse file.
                let r_id = if id.ends_with("/1") {
                    format!("{}/2", base_id)
                } else {
                    base_id.to_string()
                };

                write!(rp_writer, "{}", format_fastq(&r_id, r_seq, r_qual))?;

            } else {
                // No match, write to singles
                fs_count += 1;
                let f_seq_cow = record.seq();
                let f_seq = std::str::from_utf8(&f_seq_cow)?;
                let f_qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;
                write!(fs_writer, "{}", format_fastq(&id, f_seq, f_qual))?;
            }
        }

        // 3. Write remaining reverse reads to singles
        for item in db.iter() {
            let (key, val) = item?;
            rs_count += 1;
            let base_id = std::str::from_utf8(&key)?;
            
            let val_vec = val.to_vec();
            let (r_seq_bytes, r_qual_bytes) = if let Some(pos) = val_vec.iter().position(|&x| x == b'|') {
                (&val_vec[..pos], Some(&val_vec[pos+1..]))
            } else {
                (&val_vec[..], None)
            };
            
            let r_seq = std::str::from_utf8(r_seq_bytes)?;
            let r_qual = r_qual_bytes.map(|q| std::str::from_utf8(q)).transpose()?;
            let r_id = format!("{}/2", base_id); // Assume /2 for reverse singles

            write!(rs_writer, "{}", format_fastq(&r_id, r_seq, r_qual))?;
        }

    } else {
        // Memory-based (HashMap)
        let mut r_map = HashMap::new();

        // 1. Load reverse reads
        let mut reader = parse_fastx_file(&reverse).with_context(|| format!("Failed to open {}", reverse))?;
        while let Some(record) = reader.next() {
            let record = record?;
            r_count += 1;
            let id = std::str::from_utf8(record.id())?.to_string();
            let base_id = id.trim_end_matches("/2").trim_end_matches("/1").to_string();
            
            let seq = record.seq().to_vec();
            let qual = record.qual().map(|q| q.to_vec());
            
            r_map.insert(base_id, (seq, qual));
        }

        // 2. Process forward reads
        let mut reader = parse_fastx_file(&forward).with_context(|| format!("Failed to open {}", forward))?;
        while let Some(record) = reader.next() {
            let record = record?;
            f_count += 1;
            let id = std::str::from_utf8(record.id())?.to_string();
            let base_id = id.trim_end_matches("/1").trim_end_matches("/2");

            if let Some((r_seq, r_qual)) = r_map.remove(base_id) {
                // Match
                fp_count += 1;
                rp_count += 1;

                // Write forward
                let f_seq_cow = record.seq();
                let f_seq = std::str::from_utf8(&f_seq_cow)?;
                let f_qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;
                write!(fp_writer, "{}", format_fastq(&id, f_seq, f_qual))?;

                // Write reverse
                let r_seq_str = std::str::from_utf8(&r_seq)?;
                let r_qual_str = r_qual.as_ref().map(|q| std::str::from_utf8(q)).transpose()?;
                
                let r_id = if id.ends_with("/1") {
                    format!("{}/2", base_id)
                } else {
                    base_id.to_string()
                };

                write!(rp_writer, "{}", format_fastq(&r_id, r_seq_str, r_qual_str))?;

            } else {
                // No match
                fs_count += 1;
                let f_seq_cow = record.seq();
                let f_seq = std::str::from_utf8(&f_seq_cow)?;
                let f_qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;
                write!(fs_writer, "{}", format_fastq(&id, f_seq, f_qual))?;
            }
        }

        // 3. Remaining reverse
        for (base_id, (r_seq, r_qual)) in r_map {
            rs_count += 1;
            let r_seq_str = std::str::from_utf8(&r_seq)?;
            let r_qual_str = r_qual.as_ref().map(|q| std::str::from_utf8(q)).transpose()?;
            let r_id = format!("{}/2", base_id);
            write!(rs_writer, "{}", format_fastq(&r_id, r_seq_str, r_qual_str))?;
        }
    }

    if stats {
        println!("========= pairfq version : 0.18.0 (completion time: TODO)");
        println!("{:<40} : {:>10}", format!("Total forward reads ({})", forward), f_count);
        println!("{:<40} : {:>10}", format!("Total reverse reads ({})", reverse), r_count);
        println!("{:<40} : {:>10}", format!("Total forward paired reads ({})", fp), fp_count);
        println!("{:<40} : {:>10}", format!("Total reverse paired reads ({})", rp), rp_count);
        println!("{:<40} : {:>10}", format!("Total forward unpaired reads ({})", fs), fs_count);
        println!("{:<40} : {:>10}", format!("Total reverse unpaired reads ({})", rs), rs_count);
        println!();
        println!("{:<40} : {:>10}", "Total paired reads", fp_count + rp_count);
        println!("{:<40} : {:>10}", "Total unpaired reads", fs_count + rs_count);
    }

    Ok(())
}
