use clap::{Parser, Subcommand};
use anyhow::Result;


mod commands;
mod utils;

#[derive(Parser)]
#[command(name = "pairfq")]
#[command(version = "1.0.0")]
#[command(about = "Sync paired-end sequences from separate FASTA/Q files", long_about = "Re-pair paired-end sequences that may have been separated by quality trimming.\nThis script also writes the unpaired forward and reverse sequences to separate\nfiles so that they may be used for assembly or mapping. The input may be FastA\nor FastQ format in either Illumina 1.3+ or Illumina 1.8 format. The input files\nmay be compressed with gzip or bzip2. Optionally, the script can interleave paired\nfiles, separate interleaved files into separate forward and reverse files, and\nfix paired-end files which have lost the pair information.")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Pair the forward and reverse reads and write singletons for both forward and reverse reads to separate files.
    Makepairs {
        /// File of interleaved forward and reverse reads that has been trimmed. As below, the forward and reverse reads must be labeled in the name or comment.
        #[arg(short = 'i', long = "infile", conflicts_with_all = ["forward", "reverse"])]
        infile: Option<String>,

        /// File of foward reads (usually with "/1" or " 1" in the header).
        #[arg(short = 'f', long = "forward", required_unless_present = "infile", requires = "reverse")]
        forward: Option<String>,

        /// File of reverse reads (usually with "/2" or " 2" in the header).
        #[arg(short = 'r', long = "reverse", required_unless_present = "infile", requires = "forward")]
        reverse: Option<String>,

        /// Name for the file of paired forward reads.
        #[arg(long = "forw_paired", short = 'p', alias = "fp")]
        fp: String,

        /// Name for the file of paired reverse reads.
        #[arg(long = "rev_paired", short = 'P', alias = "rp")]
        rp: String,

        /// Name for the file of singleton forward reads.
        #[arg(long = "forw_unpaired", short = 's', alias = "fs")]
        fs: String,

        /// Name for the file of singleton reverse reads.
        #[arg(long = "rev_unpaired", short = 'S', alias = "rs")]
        rs: String,

        /// Use disk-based index (slower but less memory)
        #[arg(long, short = 'x', alias = "idx")]
        index: bool,

        /// Compress output (gzip or bzip2)
        #[arg(long, short = 'c')]
        compress: Option<String>,

        /// Print statistics on the pairing results to STDOUT (Default: No).
        #[arg(long, short = 't', alias = "stats")]
        stats: bool,
    },
    /// Interleave the paired forward and reverse files.
    Joinpairs {
        /// File of foward reads (usually with "/1" or " 1" in the header).
        #[arg(short = 'f', long = "forward")]
        forward: String,

        /// File of reverse reads (usually with "/2" or " 2" in the header).
        #[arg(short = 'r', long = "reverse")]
        reverse: String,

        /// File of interleaved reads.
        #[arg(short = 'o', long = "outfile")]
        outfile: String,

        /// Use disk-based index
        #[arg(long, short = 'x', alias = "idx")]
        index: bool,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,
    },
    /// Split the interleaved file into separate files for the forward and reverse reads.
    Splitpairs {
        /// File of interleaved forward and reverse reads.
        #[arg(short = 'i', long = "infile")]
        infile: String,

        /// File to place the foward reads.
        #[arg(short = 'f', long = "forward")]
        forward: String,

        /// File to place the reverse reads.
        #[arg(short = 'r', long = "reverse")]
        reverse: String,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,
    },
    /// Add the pair info back to the FASTA/Q header.
    Addinfo {
        /// The file of sequences without the pair information in the sequence name.
        #[arg(short = 'i', long = "infile")]
        infile: String,

        /// The file of sequences that will contain the sequence names with the pair information.
        #[arg(short = 'o', long = "outfile")]
        outfile: String,

        /// The number to append to the sequence name. Integer (Must be 1 or 2).
        #[arg(short = 'p', long = "pairnum")]
        pairnum: u8,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,

        /// Convert the sequence to uppercase.
        #[arg(long, short = 'u', alias = "uc")]
        uppercase: bool,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Makepairs { forward, reverse, infile, fp, rp, fs, rs, index, compress, stats } => {
            commands::makepairs::run(forward, reverse, infile, fp, rp, fs, rs, index, compress, stats)
        }
        Commands::Joinpairs { forward, reverse, outfile, index, compress } => {
            commands::joinpairs::run(forward, reverse, outfile, index, compress)
        }
        Commands::Splitpairs { infile, forward, reverse, compress } => {
            commands::splitpairs::run(infile, forward, reverse, compress)
        }
        Commands::Addinfo { infile, outfile, pairnum, compress, uppercase } => {
            commands::addinfo::run(infile, outfile, pairnum, compress, uppercase)
        }
    }
}
