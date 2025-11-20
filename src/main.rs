use clap::{Parser, Subcommand};
use anyhow::Result;


mod commands;
mod utils;

#[derive(Parser)]
#[command(name = "pairfq")]
#[command(version = "0.18.0")]
#[command(about = "Sync paired-end FASTX files", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Pair the forward and reverse reads and write singletons
    Makepairs {
        /// Input file (interleaved)
        #[arg(short = 'i', long = "infile", conflicts_with_all = ["forward", "reverse"])]
        infile: Option<String>,

        /// Forward reads input
        #[arg(short = 'f', long = "forward", required_unless_present = "infile", requires = "reverse")]
        forward: Option<String>,

        /// Reverse reads input
        #[arg(short = 'r', long = "reverse", required_unless_present = "infile", requires = "forward")]
        reverse: Option<String>,

        /// Forward paired output
        #[arg(long = "forw_paired", short = 'p', alias = "fp")] // alias for -fp
        fp: String,

        /// Reverse paired output
        #[arg(long = "rev_paired", short = 'P', alias = "rp")] // alias for -rp
        rp: String,

        /// Forward unpaired output
        #[arg(long = "forw_unpaired", short = 's', alias = "fs")] // alias for -fs
        fs: String,

        /// Reverse unpaired output
        #[arg(long = "rev_unpaired", short = 'S', alias = "rs")] // alias for -rs
        rs: String,

        /// Use disk-based index (slower but less memory)
        #[arg(long, short = 'x', alias = "idx")] // Changed short flag to 'x' to avoid conflict with 'i' (infile)
        index: bool,

        /// Compress output (gzip or bzip2)
        #[arg(long, short = 'c')]
        compress: Option<String>,

        /// Print stats
        #[arg(long, short = 't', alias = "stats")] // alias for -s
        stats: bool,
    },
    /// Interleave the paired forward and reverse files
    Joinpairs {
        /// Forward reads input
        #[arg(short = 'f', long = "forward")]
        forward: String,

        /// Reverse reads input
        #[arg(short = 'r', long = "reverse")]
        reverse: String,

        /// Output interleaved file
        #[arg(short = 'o', long = "outfile")]
        outfile: String,

        /// Use disk-based index
        #[arg(long, short = 'x', alias = "idx")]
        index: bool,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,
    },
    /// Split the interleaved file into separate files
    Splitpairs {
        /// Input interleaved file
        #[arg(short = 'i', long = "infile")]
        infile: String,

        /// Forward reads output
        #[arg(short = 'f', long = "forward")]
        forward: String,

        /// Reverse reads output
        #[arg(short = 'r', long = "reverse")]
        reverse: String,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,
    },
    /// Add the pair info back to the FASTA/Q header
    Addinfo {
        /// Input file
        #[arg(short = 'i', long = "infile")]
        infile: String,

        /// Output file
        #[arg(short = 'o', long = "outfile")]
        outfile: String,

        /// Pair number (1 or 2)
        #[arg(short = 'p', long = "pairnum")]
        pairnum: u8,

        /// Compress output
        #[arg(long, short = 'c')]
        compress: Option<String>,

        /// Uppercase sequence
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
