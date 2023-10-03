/*
    Rationale:
        Input list of genomes from a virus of choice
        Extract all possible triggers from all the genomes
        Compute the number of genomes the triggers are found
*/

/*
    Requirements:
        Size of trigger: optional, default 36
        Input genomes: required, at least one genome required
        Exclude [^ATGC]: optional, default yes
        Number of threads: optional, default max
        path to library: optional, no default
        save intermediate: optional, override to on, if path to library is not given
        save final: optional, override to on, if path to library is given
        file_suffix: optional, default _out
*/

mod utils;

use std::fs::File;
use std::io::{BufWriter, Write};
use log::{debug, info, trace};
use env_logger;
use utils::arg_parse::Cli;
use check_fasta::check_fasta;
use clap::Parser;
use utils::window_generate::collect_windows;
use utils::search_trgs::check_coverage;

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    debug!("Parsing commandline arguments");
    cli.input_sequence.iter().for_each(|file|
      {
        debug!("Processing input alignment files");
        match check_fasta(file, false) {
            Ok(input_file) => {
                info!("File is in fasta format {:?}", file);
                let windows = collect_windows(input_file, &cli);
                debug!("# of unique triggers: {}", windows.len());
                if cli.save_trg_list {
                    // create file to save trigger list output
                    let trg_list_out_file_name = format!("{}_{}_trglist",
                                                         file.to_str().unwrap(), cli.output_suffix);
                    let trg_list_save_file =File::create(&trg_list_out_file_name)
                        .expect("Could not create trigger list file.");
                    // BufWriter to save the list of triggers
                    let mut save_buffer = BufWriter::new(trg_list_save_file);
                    debug!("Saving list of triggers to: {}", trg_list_out_file_name);
                    windows.iter().for_each(|trg| {
                        writeln!(save_buffer, "{}", trg)
                            .expect(&format!("Unable to write to: {} ",
                                            trg_list_out_file_name));
                    });
                    debug!("Output created successfully: {}", trg_list_out_file_name);
                } else {
                    debug!("Printing list of triggers to trace. Enable trace to see triggers\\
                    on terminal.");
                    windows.iter().for_each(|trg| trace!("{}", trg))
                }
                if cli.library_path.is_some() {
                    debug!("Checking the given library file");
                    match check_fasta(&cli.library_path.clone().unwrap(), false) {
                        Ok(library) => {
                            debug!("Given library is in expected fasta format.\\
                            Proceeding to check the coverage of triggers identified.");
                            let (total_genomes, mut window_coverage) =
                                check_coverage(windows, &cli, library);
                            debug!("Analysed {} genomes in total", total_genomes);
                            if cli.save_trg_genome_tab {
                                let trg_count_save_file_name = format!("{}_{}_trgcount",
                                file.to_str().unwrap(), cli.output_suffix);
                                let trg_count_save_file = File::create(&trg_count_save_file_name)
                                    .expect("Unable to create trigger count output file");
                                debug!("Saving results to {}", trg_count_save_file_name);
                                let mut save_buffer = BufWriter::new(trg_count_save_file);
                                let total_genome_float = f64::from(total_genomes as u32);
                                writeln!(save_buffer, "Trigger,Count,%Count")
                                                .expect("Unable to file");
                                window_coverage.iter().for_each(|(trg, count)|{
                                    writeln!(save_buffer, "{}_{}={}",
                                             trg,
                                             count,
                                             f64::from(*count as u32)/total_genome_float*100.0)
                                        .expect("Unable to write to file");
                                });
                            } else {
                                debug!("Printing trigger coverage across library to trace.\\
                                Enable trace to see output on terminal.");
                                let total_genome_float = f64::from(total_genomes as u32);
                                trace!("Trigger,Count,%Count");
                                window_coverage.iter().for_each(|(trg, count)|{
                                trace!("{},{},{}",
                                             trg,
                                             count,
                                             f64::from(*count as u32)/total_genome_float*100.0)
                                });
                            }
                        },
                        Err(e) => eprintln!("{}", e)
                    }
                } else {
                    debug!("End: no library specified.");
                }
            }
            Err(e) => {
                eprintln!("{}", e)
            }
        }
      }
    );
}