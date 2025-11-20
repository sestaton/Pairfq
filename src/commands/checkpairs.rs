use anyhow::{Context, Result};
use crate::utils::get_reader;
use needletail::parse_fastx_reader;
use log::info;

struct FileCheckResult {
    path: String,
    integrity_ok: bool,
    count: usize,
}

pub fn run(forward: String, reverse: String) -> Result<()> {
    info!("Starting checkpairs");

    let f_res = check_file(&forward)?;
    let r_res = check_file(&reverse)?;

    let paired_ok = f_res.integrity_ok && r_res.integrity_ok && f_res.count == r_res.count;
    
    // Header
    println!("file\tintegrity\tpaired\tpaired_reads\tunpaired_reads");

    print_row(&f_res, paired_ok, &r_res);
    print_row(&r_res, paired_ok, &f_res);

    Ok(())
}

fn check_file(path: &str) -> Result<FileCheckResult> {
    let reader = get_reader(path).with_context(|| format!("Failed to open {}", path))?;
    let mut parser = parse_fastx_reader(reader)?;
    
    let mut count = 0;
    let mut integrity_ok = true;

    // Iterate through the file to check integrity and count
    while let Some(record) = parser.next() {
        match record {
            Ok(_) => count += 1,
            Err(_) => {
                integrity_ok = false;
                break;
            }
        }
    }

    Ok(FileCheckResult {
        path: path.to_string(),
        integrity_ok,
        count,
    })
}

fn print_row(res: &FileCheckResult, paired_ok: bool, other: &FileCheckResult) {
    let integrity_symbol = if res.integrity_ok { "\u{2705}" } else { "\u{274C}" }; // Check mark or Cross
    let paired_symbol = if paired_ok { "\u{2705}" } else { "\u{274C}" };

    let paired_reads = if paired_ok {
        res.count
    } else {
        std::cmp::min(res.count, other.count)
    };

    let unpaired_reads = if paired_ok {
        0
    } else {
        if res.count > other.count {
            res.count - other.count
        } else {
            0 // If this file has fewer reads, it doesn't have "unpaired" reads in the sense of extras, 
              // but the other file does. Or should we report the difference? 
              // User said "unpaired_reads". If F=100, R=90. F has 10 unpaired. R has 0 unpaired (all 90 are paired with first 90 of F).
              // So: max(0, self.count - other.count)
        }
    };
    
    // Calculate unpaired logic:
    // If F > R: F has (F-R) unpaired. R has 0.
    // If R > F: R has (R-F) unpaired. F has 0.
    let unpaired_count = if res.count > other.count {
        res.count - other.count
    } else {
        0
    };

    println!("{}\t{}\t{}\t{}\t{}", 
        res.path, 
        integrity_symbol, 
        paired_symbol, 
        paired_reads, 
        unpaired_count
    );
}
