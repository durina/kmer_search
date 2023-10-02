/*
    Generate n-nucleotide windows
        Include non-ATGC if flag
    Combine all triggers and output all possible n-nucleotide triggers
*/

/*
    Relevant flags
        trg_len
        input_Sequences
        include_non_standard
        nproc
*/

use std::collections::HashSet;
use std::io::{BufRead, Seek, SeekFrom};
use std::sync::{Arc, Mutex};
use check_fasta::lib_utils::struct_helper::FileBufferHelper;
use log::{debug, info, trace};
use regex::Regex;
use threadpool::ThreadPool;
use crate::utils::arg_parse::Cli;

const STANDARD_DNA: &str = "ATGC";
pub const IUPAC_DNA: &str = "ATGCUWSMKRYBDHVN";
pub fn collect_windows( infile: FileBufferHelper, cli: &Cli) -> Vec<String> {
    let dna_letters: Vec<char>;
    let dna_filter: Regex;
    match cli.include_non_standard {
        true => {
            dna_letters = IUPAC_DNA.chars().collect();
            dna_filter = Regex::new(format!("[^{}]", IUPAC_DNA).as_str()).unwrap();
        },
        false => {
            dna_letters = STANDARD_DNA.chars().collect();
            dna_filter = Regex::new(format!("[^{}]", STANDARD_DNA).as_str()).unwrap();
        }
    }
    debug!("Using {:?} characters", dna_letters);
    let trg_len = cli.trg_len;
    debug!("Trigger length set to: {}", trg_len);
    let nproc = match cli.nproc {
        Some(x) => x,
        None => num_cpus::get()
    };
    rayon::ThreadPoolBuilder::new().num_threads(nproc).build_global().unwrap();
    analyse_input( &trg_len, infile, dna_filter, nproc)
}

fn analyse_input(
                    trg_len: &usize,
                    mut infile: FileBufferHelper,
                    dna_filter: Regex,
                    nproc: usize
                ) -> Vec<String> {
    debug!("Analysing input");
    let trg_vec: Vec<String> = Vec::new();
    let arc_trg_vec: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(trg_vec));
    let pool = ThreadPool::new(nproc);
    // reset file buffer to start
    infile.buffer_reader.seek(SeekFrom::Start(0)).unwrap();
    let mut genome = String::new();
    infile.line.clear();
    while infile.buffer_reader.read_line(&mut infile.line).unwrap_or(0) >= 1 {
        trace!("Line starts with {:?}", infile.line.chars().nth(0).unwrap());
        if &infile.line[0..1] == ">" && !genome.is_empty() {
            info!("Entering Header block");
            let temp_genome = genome.drain(..).collect();
            let arc_clone = Arc::clone(&arc_trg_vec);
            let copy_trg_len = trg_len.clone();
            pool.execute(move ||
                { make_windows(temp_genome, arc_clone, copy_trg_len) }
            );
        } else if IUPAC_DNA.contains(&infile.line[0..1]) {
            genome = genome + infile.line.trim();
        }
        infile.line.clear();
    }
    make_windows(genome, Arc::clone(&arc_trg_vec), *trg_len);
    pool.join();

    let trg_vec = Arc::try_unwrap(arc_trg_vec)
        .unwrap()
        .into_inner()
        .unwrap();
    let mut trg_vec = trg_vec.into_iter()
        .collect::<HashSet<_>>()
        .into_iter()
        .collect::<Vec<String>>();
    trg_vec.retain(|trg| !dna_filter.is_match(trg));
    trg_vec
}

fn make_windows(genome: String, arc_trg_vec: Arc<Mutex<Vec<String>>>,
                trg_len: usize) {
    trace!("Making windows");
    if let Ok(mut trg_vec) = arc_trg_vec.lock() {
        for trg_idx in 0..=(genome.len() - trg_len) {
            let trg = &genome[trg_idx..=trg_idx+trg_len-1];
            trace!("{} - {}", trg, trg.len());
            trg_vec.push(trg.to_string().to_uppercase());
        }
    }
}