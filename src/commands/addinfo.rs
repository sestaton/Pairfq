use anyhow::{Context, Result};
use crate::utils::{get_writer, format_fastq};
use needletail::parse_fastx_file;
use log::info;

pub fn run(
    infile: String,
    outfile: String,
    pairnum: u8,
    compress: Option<String>,
    uppercase: bool,
) -> Result<()> {
    info!("Starting addinfo");

    if pairnum != 1 && pairnum != 2 {
        return Err(anyhow::anyhow!("pairnum must be 1 or 2"));
    }

    let mut writer = get_writer(&outfile, compress.as_deref())?;
    let mut reader = parse_fastx_file(&infile).with_context(|| format!("Failed to open {}", infile))?;
    
    let suffix = format!("/{}", pairnum);

    while let Some(record) = reader.next() {
        let record = record?;
        let id = std::str::from_utf8(record.id())?;
        
        let new_id = if id.ends_with(&suffix) {
            id.to_string()
        } else {
            format!("{}{}", id, suffix)
        };

        let seq_cow = record.seq();
        let seq = std::str::from_utf8(&seq_cow)?;
        let seq = if uppercase {
            seq.to_uppercase()
        } else {
            seq.to_string()
        };

        let qual = record.qual().map(|q| std::str::from_utf8(q)).transpose()?;
        
        write!(writer, "{}", format_fastq(&new_id, &seq, qual))?;
    }

    Ok(())
}
