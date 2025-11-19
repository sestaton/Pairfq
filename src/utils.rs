use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use bzip2::read::BzDecoder;
use bzip2::write::BzEncoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn get_reader(path: &str) -> Result<Box<dyn BufRead + Send>> {
    let reader: Box<dyn BufRead + Send> = if path == "-" {
        Box::new(BufReader::new(io::stdin()))
    } else {
        let file = File::open(path).with_context(|| format!("Failed to open file: {}", path))?;
        if path.ends_with(".gz") {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else if path.ends_with(".bz2") {
            Box::new(BufReader::new(BzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    };
    Ok(reader)
}

pub fn get_writer(path: &str, compress: Option<&str>) -> Result<Box<dyn Write + Send>> {
    let compression_type = if let Some(c) = compress {
        c
    } else if path.ends_with(".gz") {
        "gzip"
    } else if path.ends_with(".bz2") {
        "bzip2"
    } else {
        "none"
    };

    let writer: Box<dyn Write + Send> = if path == "-" {
        match compression_type {
            "gzip" => Box::new(GzEncoder::new(io::stdout(), Compression::default())),
            "bzip2" => Box::new(BzEncoder::new(io::stdout(), bzip2::Compression::default())),
            _ => Box::new(BufWriter::new(io::stdout())),
        }
    } else {
        let file = File::create(path).with_context(|| format!("Failed to create file: {}", path))?;
        match compression_type {
            "gzip" => Box::new(GzEncoder::new(file, Compression::default())),
            "bzip2" => Box::new(BzEncoder::new(file, bzip2::Compression::default())),
            _ => Box::new(BufWriter::new(file)),
        }
    };
    Ok(writer)
}

// Helper to format FASTQ record
pub fn format_fastq(id: &str, seq: &str, qual: Option<&str>) -> String {
    if let Some(q) = qual {
        format!("@{}\n{}\n+\n{}\n", id, seq, q)
    } else {
        format!(">{}\n{}\n", id, seq)
    }
}
