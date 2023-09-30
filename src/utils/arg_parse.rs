use std::ops::RangeInclusive;
use std::path::PathBuf;
use clap::Parser;

/*
    Requirements:
        Size of trigger: optional, default 36
        Input genomes: Path to file(s) with at least one genome in fasta format
        Exclude [^ATGC]: optional, default yes
        Number of threads: optional, default max
        path to library: optional, no default
        save intermediate: optional, override to on, if path to library is not given
        save final: optional, override to on, if path to library is given
        file_suffix: optional, default _out
*/
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Size of trigger, default 36 nucleotides. Range 10 to 100.
    #[arg(short='l', long="length_trigger", default_value_t=36, value_parser=trg_validate)]
    pub trg_len: usize,

    /// Path to file(s) with at least one genome in fasta format
    #[arg(short='i', long="infile", required = true, action=clap::ArgAction::Append)]
    pub input_sequence: Vec<PathBuf>,

    /// Switch to include non standard IUPAC nucleic acids. Default false.
    #[arg(long="include-non-standard", default_value_t=false)]
    pub include_non_standard: bool,

    /// Number of threads to be used for the program. Optional.
    /// If not specified maximum threads available will be used.
    #[arg(long="threads")]
    pub nproc: Option<usize>,

    /// Path to library to search against. Optional
    #[arg(long="library", id="library")]
    pub library_path: Option<PathBuf>,

    /// save intermediate results, list of all triggers.
    #[arg(long="save-trg-list", default_value_t=false)]
    pub save_trg_list: bool,

    /// save final results, catalogue of genomes covered per trigger.
    #[arg(short='o', long="output", default_value_t=false, requires="library")]
    pub save_trg_genome_tab: bool,

    /// file suffix to be added to output file. Default "_output"
    #[arg(short='s', long="output-suffix", default_value_t=String::from("_output"))]
    pub output_suffix: String,
}

fn trg_validate(trg_len: &str) -> Result<usize, String> {
    let trg_len: usize = trg_len.parse().unwrap_or(36);
    const TRG_LEN_RANGE: RangeInclusive<usize> = 10..=100;
    if TRG_LEN_RANGE.contains(&trg_len) {
        Ok(trg_len)
    } else {
        Err(format!("Trigger length not in range considered: {} - {}",
            TRG_LEN_RANGE.start(), TRG_LEN_RANGE.end()))
    }
}