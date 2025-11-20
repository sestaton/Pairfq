use anyhow::{Context, Result};
use crate::utils::{get_writer, format_fastq, get_reader};
use needletail::parse_fastx_file;
use std::collections::HashMap;
use log::info;
use std::io::Write;

struct Stats {
    forward_reads: usize,
    reverse_reads: usize,
    forward_paired: usize,
    reverse_paired: usize,
    forward_unpaired: usize,
    reverse_unpaired: usize,
    total_paired: usize,
    total_unpaired: usize,
}

pub fn run(
    forward: Option<String>,
    reverse: Option<String>,
    infile: Option<String>,
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

    let mut stats_counts = Stats {
        forward_reads: 0,
        reverse_reads: 0,
        forward_paired: 0,
        reverse_paired: 0,
        forward_unpaired: 0,
        reverse_unpaired: 0,
        total_paired: 0,
        total_unpaired: 0,
    };

    if let Some(infile_path) = infile {
        // Interleaved input mode
        let reader = get_reader(&infile_path)?;
        let mut parser = needletail::parse_fastx_reader(reader)?;
        
        struct BufferedRecord {
            id: Vec<u8>,
            seq: Vec<u8>,
            qual: Option<Vec<u8>>,
        }

        let mut record_buffer: Option<BufferedRecord> = None;

        while let Some(record) = parser.next() {
            let record = record?;
            if let Some(first_record) = record_buffer.take() {
                // We have a pair
                stats_counts.forward_reads += 1;
                stats_counts.reverse_reads += 1;
                stats_counts.forward_paired += 1;
                stats_counts.reverse_paired += 1;
                stats_counts.total_paired += 2;

                // Write first record (forward)
                let f_id = std::str::from_utf8(&first_record.id)?;
                let f_seq = std::str::from_utf8(&first_record.seq)?;
                let f_qual = first_record.qual.as_deref().map(std::str::from_utf8).transpose()?;
                write!(fp_writer, "{}", format_fastq(f_id, f_seq, f_qual))?;

                // Write second record (reverse)
                write_record(&mut rp_writer, &record)?;
            } else {
                // Buffer this record
                record_buffer = Some(BufferedRecord {
                    id: record.id().to_vec(),
                    seq: record.seq().to_vec(),
                    qual: record.qual().map(|q| q.to_vec()),
                });
            }
        }
        
        // If there's a record left in the buffer, it's a singleton (unpaired forward)
        if let Some(record) = record_buffer {
             stats_counts.forward_reads += 1;
             stats_counts.forward_unpaired += 1;
             stats_counts.total_unpaired += 1;
             
             let f_id = std::str::from_utf8(&record.id)?;
             let f_seq = std::str::from_utf8(&record.seq)?;
             let f_qual = record.qual.as_deref().map(std::str::from_utf8).transpose()?;
             write!(fs_writer, "{}", format_fastq(f_id, f_seq, f_qual))?;
        }

    } else if let (Some(f_path), Some(r_path)) = (forward, reverse) {
        if index {
            run_ondisk(&f_path, &r_path, &mut fp_writer, &mut rp_writer, &mut fs_writer, &mut rs_writer, &mut stats_counts)?;
        } else {
            run_inmemory(&f_path, &r_path, &mut fp_writer, &mut rp_writer, &mut fs_writer, &mut rs_writer, &mut stats_counts)?;
        }
    } else {
        anyhow::bail!("Must provide either --infile or both --forward and --reverse");
    }

    if stats {
        print_stats(&stats_counts, "TODO"); // TODO: Add timing
    }

    Ok(())
}

fn run_ondisk<W: Write>(
    forward: &str,
    reverse: &str,
    fp_writer: &mut W,
    rp_writer: &mut W,
    fs_writer: &mut W,
    rs_writer: &mut W,
    stats: &mut Stats,
) -> Result<()> {
    // Disk-based indexing using sled
    let tmp_dir = tempfile::tempdir()?;
    let db = sled::open(tmp_dir.path().join("pairfq_db"))?;

    // 1. Index reverse reads
    let mut reader = parse_fastx_file(reverse).with_context(|| format!("Failed to open {}", reverse))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.reverse_reads += 1;
        let id = std::str::from_utf8(record.id())?.to_string();
        let base_id = id.trim_end_matches("/2").trim_end_matches("/1");
        
        let seq = record.seq();
        let qual = record.qual();
        
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
    let mut reader = parse_fastx_file(forward).with_context(|| format!("Failed to open {}", forward))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.forward_reads += 1;
        let id = std::str::from_utf8(record.id())?.to_string();
        let base_id = id.trim_end_matches("/1").trim_end_matches("/2");

        if let Some(val) = db.remove(base_id)? {
            // Match found
            stats.forward_paired += 1;
            stats.reverse_paired += 1;
            stats.total_paired += 2;

            // Write forward
            write_record(fp_writer, &record)?;

            // Write reverse (reconstruct from DB value)
            let val_vec = val.to_vec();
            let (r_seq_bytes, r_qual_bytes) = if let Some(pos) = val_vec.iter().position(|&x| x == b'|') {
                (&val_vec[..pos], Some(&val_vec[pos+1..]))
            } else {
                (&val_vec[..], None)
            };
            
            let r_seq = std::str::from_utf8(r_seq_bytes)?;
            let r_qual = r_qual_bytes.map(|q| std::str::from_utf8(q)).transpose()?;
            
            let r_id = if id.ends_with("/1") {
                format!("{}/2", base_id)
            } else {
                base_id.to_string()
            };

            write!(rp_writer, "{}", format_fastq(&r_id, r_seq, r_qual))?;

        } else {
            // No match, write to singles
            stats.forward_unpaired += 1;
            stats.total_unpaired += 1;
            write_record(fs_writer, &record)?;
        }
    }

    // 3. Write remaining reverse reads to singles
    for item in db.iter() {
        let (key, val) = item?;
        stats.reverse_unpaired += 1;
        stats.total_unpaired += 1;
        let base_id = std::str::from_utf8(&key)?;
        
        let val_vec = val.to_vec();
        let (r_seq_bytes, r_qual_bytes) = if let Some(pos) = val_vec.iter().position(|&x| x == b'|') {
            (&val_vec[..pos], Some(&val_vec[pos+1..]))
        } else {
            (&val_vec[..], None)
        };
        
        let r_seq = std::str::from_utf8(r_seq_bytes)?;
        let r_qual = r_qual_bytes.map(|q| std::str::from_utf8(q)).transpose()?;
        let r_id = format!("{}/2", base_id);

        write!(rs_writer, "{}", format_fastq(&r_id, r_seq, r_qual))?;
    }

    Ok(())
}

fn run_inmemory<W: Write>(
    forward: &str,
    reverse: &str,
    fp_writer: &mut W,
    rp_writer: &mut W,
    fs_writer: &mut W,
    rs_writer: &mut W,
    stats: &mut Stats,
) -> Result<()> {
    let mut r_map = HashMap::new();

    // 1. Load reverse reads
    let mut reader = parse_fastx_file(reverse).with_context(|| format!("Failed to open {}", reverse))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.reverse_reads += 1;
        let id = std::str::from_utf8(record.id())?.to_string();
        let base_id = id.trim_end_matches("/2").trim_end_matches("/1").to_string();
        
        let seq = record.seq().to_vec();
        let qual = record.qual().map(|q| q.to_vec());
        
        r_map.insert(base_id, (seq, qual));
    }

    // 2. Process forward reads
    let mut reader = parse_fastx_file(forward).with_context(|| format!("Failed to open {}", forward))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.forward_reads += 1;
        let id = std::str::from_utf8(record.id())?.to_string();
        let base_id = id.trim_end_matches("/1").trim_end_matches("/2");

        if let Some((r_seq, r_qual)) = r_map.remove(base_id) {
            // Match
            stats.forward_paired += 1;
            stats.reverse_paired += 1;
            stats.total_paired += 2;

            // Write forward
            write_record(fp_writer, &record)?;

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
            stats.forward_unpaired += 1;
            stats.total_unpaired += 1;
            write_record(fs_writer, &record)?;
        }
    }

    // 3. Remaining reverse
    for (base_id, (r_seq, r_qual)) in r_map {
        stats.reverse_unpaired += 1;
        stats.total_unpaired += 1;
        let r_seq_str = std::str::from_utf8(&r_seq)?;
        let r_qual_str = r_qual.as_ref().map(|q| std::str::from_utf8(q)).transpose()?;
        let r_id = format!("{}/2", base_id);
        write!(rs_writer, "{}", format_fastq(&r_id, r_seq_str, r_qual_str))?;
    }

    Ok(())
}

fn write_record<W: Write>(writer: &mut W, record: &needletail::parser::SequenceRecord) -> Result<()> {
    record.write(writer, None)?;
    Ok(())
}

fn print_stats(stats: &Stats, _time: &str) {
    println!("========= pairfq version : 0.18.0 (completion time: TODO)");
    println!("{:<40} : {:>10}", "Total forward reads", stats.forward_reads);
    println!("{:<40} : {:>10}", "Total reverse reads", stats.reverse_reads);
    println!("{:<40} : {:>10}", "Total forward paired reads", stats.forward_paired);
    println!("{:<40} : {:>10}", "Total reverse paired reads", stats.reverse_paired);
    println!("{:<40} : {:>10}", "Total forward unpaired reads", stats.forward_unpaired);
    println!("{:<40} : {:>10}", "Total reverse unpaired reads", stats.reverse_unpaired);
    println!();
    println!("{:<40} : {:>10}", "Total paired reads", stats.total_paired);
    println!("{:<40} : {:>10}", "Total unpaired reads", stats.total_unpaired);
}
