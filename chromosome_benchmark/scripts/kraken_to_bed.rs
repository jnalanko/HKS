#!/usr/bin/env rust-script
//! Convert Kraken2 output to BED format with integer taxid labels.
//!
//! Each classified kmer run becomes a BED record:
//!   col 1  = seq_id (query sequence identifier)
//!   col 2  = start position
//!   col 3  = end position
//!   col 4  = taxid (integer) or 'A' for ambiguous kmers
//!
//! Unclassified sequences (status 'U') are skipped.
//!
//! Usage:
//!   kraken_to_bed <kraken_output> [output.bed]
//!
//! ```cargo
//! [dependencies]
//! ```

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: kraken_to_bed <kraken_output> [output.bed]");
        std::process::exit(1);
    }

    let kraken_path = &args[1];
    let out_path    = args.get(2);

    let input = BufReader::new(
        File::open(kraken_path).unwrap_or_else(|e| panic!("Cannot open {kraken_path}: {e}")),
    );

    let stdout = io::stdout();
    let mut out: Box<dyn Write> = match out_path {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).unwrap_or_else(|e| panic!("Cannot create {p}: {e}")),
        )),
        None => Box::new(BufWriter::new(stdout.lock())),
    };

    for line in input.lines() {
        let line = line.expect("read error");
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.splitn(5, '\t').collect();
        if cols.len() < 5 {
            continue;
        }
        if cols[0] == "U" {
            continue;
        }
        let seq_id = cols[1];
        let hits   = cols[4];

        let mut pos: u64 = 0;
        for pair in hits.split_whitespace() {
            let (taxid_str, count_str) = pair.split_once(':').expect("malformed hit pair");
            let count: u64 = count_str.parse().expect("non-integer count");
            writeln!(out, "{seq_id}\t{pos}\t{}\t{taxid_str}", pos + count).expect("write error");
            pos += count;
        }
    }
}
