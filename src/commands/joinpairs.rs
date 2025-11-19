use anyhow::{Context, Result};
use crate::utils::{get_writer, format_fastq};
use needletail::parse_fastx_file;
use log::info;

pub fn run(
    forward: String,
    reverse: String,
    outfile: String,
    _index: bool, // Not used in this implementation as we stream both
    compress: Option<String>,
) -> Result<()> {
    info!("Starting joinpairs");

    let mut writer = get_writer(&outfile, compress.as_deref())?;
    
    let mut f_reader = parse_fastx_file(&forward).with_context(|| format!("Failed to open {}", forward))?;
    let mut r_reader = parse_fastx_file(&reverse).with_context(|| format!("Failed to open {}", reverse))?;

    loop {
        let f_rec = f_reader.next();
        let r_rec = r_reader.next();

        match (f_rec, r_rec) {
            (Some(f), Some(r)) => {
                let f = f?;
                let r = r?;

                let f_id = std::str::from_utf8(f.id())?;
                let r_id = std::str::from_utf8(r.id())?;
                
                // Basic check for ID match (ignoring suffixes)
                let f_base = f_id.trim_end_matches("/1");
                let r_base = r_id.trim_end_matches("/2");
                
                if f_base != r_base {
                    log::warn!("IDs do not match: {} vs {}", f_id, r_id);
                }

                let f_seq_cow = f.seq();
                let f_seq = std::str::from_utf8(&f_seq_cow)?;
                let f_qual = f.qual().map(|q| std::str::from_utf8(q)).transpose()?;
                
                let r_seq_cow = r.seq();
                let r_seq = std::str::from_utf8(&r_seq_cow)?;
                let r_qual = r.qual().map(|q| std::str::from_utf8(q)).transpose()?;

                write!(writer, "{}", format_fastq(f_id, f_seq, f_qual))?;
                write!(writer, "{}", format_fastq(r_id, r_seq, r_qual))?;
            }
            (None, None) => break,
            _ => {
                return Err(anyhow::anyhow!("Files have different number of records"));
            }
        }
    }

    Ok(())
}
