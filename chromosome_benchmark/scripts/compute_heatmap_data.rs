#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! flate2 = "1.1.9"
//! ```

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use flate2::read::GzDecoder;

// ---- file loading ----

type Counts = HashMap<String, HashMap<String, i64>>;

fn load_counts(path: &str) -> Counts {
    let mut reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(
            File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}")),
        )))
    } else {
        Box::new(BufReader::new(
            File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}")),
        ))
    };

    let mut counts: Counts = HashMap::new();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line).expect("read error") == 0 { break; }
        let s = line.trim_end_matches(|c| c == '\n' || c == '\r');
        if s.is_empty() || !s.contains('\t') { continue; }

        let mut parts = s.split('\t');
        let row_str = match parts.next() { Some(s) => s, None => continue };
        let start: i64 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let end:   i64 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let col_str = s.rsplit('\t').next().map(str::trim).unwrap_or("");
        if col_str.is_empty() { continue; }

        let delta = end - start;

        // Only allocate a new String when the key is new.
        if !counts.contains_key(row_str) {
            counts.insert(row_str.to_string(), HashMap::new());
        }
        let row_map = counts.get_mut(row_str).unwrap();

        if let Some(c) = row_map.get_mut(col_str) {
            *c += delta;
        } else {
            row_map.insert(col_str.to_string(), delta);
        }
    }

    counts
}

// ---- main ----

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: compute_heatmap_data <input.bed[.gz]> <out.tsv>");
        std::process::exit(1);
    }
    let bed_path = &args[1];
    let out_path = &args[2];

    let counts = load_counts(bed_path);

    let mut row_list: Vec<&String> = counts.keys().collect();
    row_list.sort_unstable();
    let col_set: HashSet<&String> = counts.values().flat_map(HashMap::keys).collect();
    let mut col_list: Vec<&String> = col_set.into_iter().collect();
    col_list.sort_unstable();

    // Write TSV: header row then one row per row-label
    let mut out = std::io::BufWriter::new(
        File::create(out_path).unwrap_or_else(|e| panic!("Cannot create {out_path}: {e}"))
    );

    write!(out, "label").unwrap();
    for col in &col_list {
        write!(out, "\t{col}").unwrap();
    }
    writeln!(out).unwrap();

    for row in &row_list {
        write!(out, "{row}").unwrap();
        for col in &col_list {
            let count = counts[*row].get(*col).copied().unwrap_or(0);
            write!(out, "\t{count}").unwrap();
        }
        writeln!(out).unwrap();
    }

    eprintln!("Saved {out_path}");
}
