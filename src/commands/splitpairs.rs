use crate::utils::{format_fastq, get_writer};
use anyhow::{Context, Result};
use log::info;
use needletail::parse_fastx_file;

pub fn run(
    infile: String,
    forward: String,
    reverse: String,
    compress: Option<String>,
) -> Result<()> {
    info!("Starting splitpairs");

    let mut f_writer = get_writer(&forward, compress.as_deref())?;
    let mut r_writer = get_writer(&reverse, compress.as_deref())?;

    let mut reader =
        parse_fastx_file(&infile).with_context(|| format!("Failed to open {}", infile))?;

    let mut count = 0;
    while let Some(record) = reader.next() {
        let record = record?;
        let id = std::str::from_utf8(record.id())?;
        let seq_cow = record.seq();
        let seq = std::str::from_utf8(&seq_cow)?;
        let qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;

        // Check suffix or just alternate
        // Perl script checks for /1 or /2 or comments.
        // If explicit suffix found, use it. Else alternate.

        let is_forward = if id.ends_with("/1") {
            true
        } else if id.ends_with("/2") {
            false
        } else {
            count % 2 == 0
        };

        if is_forward {
            write!(f_writer, "{}", format_fastq(id, seq, qual))?;
        } else {
            write!(r_writer, "{}", format_fastq(id, seq, qual))?;
        }

        count += 1;
    }

    Ok(())
}
