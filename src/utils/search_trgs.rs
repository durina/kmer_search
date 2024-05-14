/*
    Search the triggers collected against a library of genomes
    Tabulate number of genomes for each trigger
*/

/*
    Relevant parameters
        library path
        nproc
        save_trg_genome_tab
        output_suffix
*/

use std::io::{BufRead, Seek, SeekFrom};
use std::sync::{Arc, Mutex};
use std::thread::sleep;
use std::time::Duration;
use ahash::{AHashMap};
use check_fasta::lib_utils::struct_helper::FileBufferHelper;
use log::{debug, info, trace};
use threadpool::ThreadPool;
use crate::utils::arg_parse::Cli;
use super::window_generate::IUPAC_DNA;
use rayon::prelude::*;
pub fn check_coverage(trg_windows: Vec<String>, cli: &Cli, library: FileBufferHelper)
                -> (usize, Vec<(String, usize)>) {
    debug!("Checking coverage");
    let mut trg_map: AHashMap<String, usize> = AHashMap::new();
    trg_windows.into_iter().for_each(|trg| { trg_map.insert(trg, 0); });
    let nproc = match cli.nproc {
        Some(x) => x,
        None => num_cpus::get()
    };
    let pool = ThreadPool::new(nproc);
    count_genomes(trg_map, library, pool)
}

fn count_genomes(trg_map: AHashMap<String, usize>,
                 mut library: FileBufferHelper, pool: ThreadPool,
                 ) -> (usize, Vec<(String, usize)>) {
    debug!("Counting number of genomes per trigger");
    let arc_trg_map: Arc<Mutex<AHashMap<String, usize>>> = Arc::new(Mutex::new(trg_map));
    library.buffer_reader.seek(SeekFrom::Start(0)).unwrap();
    library.line.clear();
    let mut total_genome_count : usize = 0;
    let mut genome = String::new();
    while library.buffer_reader.read_line(&mut library.line).unwrap_or(0) >= 1 {
        if &library.line[0..1] == ">" && !genome.is_empty() {
            trace!("Entering header block");
            let temp_genome: String = genome.drain(..).collect();
            let arc_clone = Arc::clone(&arc_trg_map);
            total_genome_count += 1;
            pool.execute(move || tabulate(temp_genome, arc_clone));
        } else if IUPAC_DNA.contains(&library.line[0..1]) {
            genome += library.line.trim();
        }
        library.line.clear();
    }
    total_genome_count += 1;
    let arc_clone = Arc::clone(&arc_trg_map);
    pool.execute(move || tabulate(genome, arc_clone));
    while pool.queued_count() > 0 {
        sleep(Duration::from_secs(5));
        info!("Queue {}, max {}", pool.queued_count(), pool.max_count());
    }
    pool.join();
    let trg_map = Arc::try_unwrap(arc_trg_map)
        .expect("Arc failed to unwrap")
        .into_inner()
        .expect("Mutex failed to unwrap");
    let trg_vec: Vec<(String, usize)> = trg_map.into_iter().collect();
    (total_genome_count, trg_vec)
}
fn tabulate(genome: String, trg_map: Arc<Mutex<AHashMap<String,usize>>>) {
    trace!("Tabulating number of genomes");
    if let Ok(mut trg_map) = trg_map.lock() {
        trg_map.par_iter_mut().for_each(
            |(trg, trg_count)| {
                if genome.to_uppercase().find(trg).is_some() {
                    *trg_count += 1;
                }
            }
        )
    }
}