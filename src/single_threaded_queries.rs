use std::{ops::Range, path::Path};

use jseqio::reader::DynamicFastXReader;

use crate::color_storage::SimpleColorStorage;
use crate::single_colored_kmers::SingleColoredKmers;
use crate::traits::*;

fn print_run(seq_id: usize, run_color: Option<usize>, range: Range<usize>, root_id: usize) {
    // This code is almost duplicated in parallel_queries.rs
    if !range.is_empty() {
        match run_color {
            Some(c) if c == root_id => println!("{seq_id}\t{}\t{}\t*", range.start, range.end),
            Some(c) => println!("{seq_id}\t{}\t{}\t{}", range.start, range.end, c),
            None => (), // Do not print runs of misses
        }
    }
}

pub fn lookup_single_threaded<L,C>(query_path: &Path, index: &SingleColoredKmers<L,C>, k: usize)
where L: sbwt::ContractLeft + Clone + MySerialize + From<sbwt::LcsArray> + LcsAccess,
      C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>
{

    let root_id = index.color_hierarchy().root();
    let mut reader = DynamicFastXReader::from_file(&query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
    let mut seq_id = 0_usize;

    println!("seq_rank\tfrom_kmer\tto_kmer\tcolor");
    while let Some(rec) = reader.read_next().unwrap() {
        let mut run_start: usize = 0;
        let mut run_color: Option<usize> = None;
        let mut n_kmers = 0;
        for (i, color) in index.lookup_kmers(rec.seq, k).enumerate() {
            if i == 0 { // Start a new run
                run_start = i;
                run_color = color;
            }
            else if color != run_color {
                // Run ends
                print_run(seq_id, run_color, run_start..i, root_id);
                run_start = i;
                run_color = color;
            } else {
                // Run continues
            }
            n_kmers += 1;
        }

        // Done processing the sequence. Close the last run
        if n_kmers > 0 {
            print_run(seq_id, run_color, run_start..n_kmers, root_id);
        }

        seq_id += 1;
    }

}