use anyhow::{Context, Result};
use crate::utils::{get_writer, get_reader, write_fastq};
use needletail::parse_fastx_file;
use ahash::AHashMap;
use log::info;
use std::io::Write;
use std::time::Instant;

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
    let start_time = Instant::now();
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
                write_fastq(&mut fp_writer, &first_record.id, &first_record.seq, first_record.qual.as_deref())?;

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
             
             write_fastq(&mut fs_writer, &record.id, &record.seq, record.qual.as_deref())?;
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
        let duration = start_time.elapsed();
        print_stats(&stats_counts, duration);
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
        let id = record.id();
        let base_id = get_base_id(id);
        
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
        let id = record.id();
        let base_id = get_base_id(id);

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
            
            let r_id = if id.ends_with(b"/1") {
                let mut rid = base_id.to_vec();
                rid.extend_from_slice(b"/2");
                rid
            } else {
                base_id.to_vec()
            };

            write_fastq(rp_writer, &r_id, r_seq_bytes, r_qual_bytes)?;

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
        let base_id = key;
        
        let val_vec = val.to_vec();
        let (r_seq_bytes, r_qual_bytes) = if let Some(pos) = val_vec.iter().position(|&x| x == b'|') {
            (&val_vec[..pos], Some(&val_vec[pos+1..]))
        } else {
            (&val_vec[..], None)
        };
        
        let mut r_id = base_id.to_vec();
        r_id.extend_from_slice(b"/2");

        write_fastq(rs_writer, &r_id, r_seq_bytes, r_qual_bytes)?;
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
    let mut r_map = AHashMap::new();

    // 1. Load reverse reads
    let mut reader = parse_fastx_file(reverse).with_context(|| format!("Failed to open {}", reverse))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.reverse_reads += 1;
        let id = record.id();
        let base_id = get_base_id(id).to_vec();
        
        let seq = record.seq().to_vec();
        let qual = record.qual().map(|q| q.to_vec());
        
        r_map.insert(base_id, (seq, qual));
    }

    // 2. Process forward reads
    let mut reader = parse_fastx_file(forward).with_context(|| format!("Failed to open {}", forward))?;
    while let Some(record) = reader.next() {
        let record = record?;
        stats.forward_reads += 1;
        let id = record.id();
        let base_id = get_base_id(id);

        if let Some((r_seq, r_qual)) = r_map.remove(base_id) {
            // Match
            stats.forward_paired += 1;
            stats.reverse_paired += 1;
            stats.total_paired += 2;

            // Write forward
            write_record(fp_writer, &record)?;

            // Write reverse
            let r_id = if id.ends_with(b"/1") {
                let mut rid = base_id.to_vec();
                rid.extend_from_slice(b"/2");
                rid
            } else {
                base_id.to_vec()
            };

            write_fastq(rp_writer, &r_id, &r_seq, r_qual.as_deref())?;

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
        
        let mut r_id = base_id;
        r_id.extend_from_slice(b"/2");
        
        write_fastq(rs_writer, &r_id, &r_seq, r_qual.as_deref())?;
    }

    Ok(())
}

fn write_record<W: Write>(writer: &mut W, record: &needletail::parser::SequenceRecord) -> Result<()> {
    write_fastq(writer, record.id(), &record.seq(), record.qual())
}

fn get_base_id(id: &[u8]) -> &[u8] {
    if id.ends_with(b"/1") || id.ends_with(b"/2") {
        &id[..id.len()-2]
    } else {
        id
    }
}

fn print_stats(stats: &Stats, duration: std::time::Duration) {
    println!("========= pairfq version : 1.1.0 (completion time: {:.2?})", duration);
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
