use std::{cmp::{max, Reverse}, io::Write, ops::Range};
use jseqio::seq_db::SeqDB;
use crate::{color_storage::SimpleColorStorage, single_colored_kmers::{ColorHierarchy, SingleColoredKmers}};
use crate::traits::*;

pub trait RunWriter: Send {
    fn write_header(&mut self);
    fn write_run(&mut self, seq_id: isize, run_color: Option<usize>, range: Range<usize>);
    fn flush(&mut self);
}

impl<T: RunWriter + ?Sized> RunWriter for &mut T {
    fn write_header(&mut self) { (**self).write_header() }
    fn write_run(&mut self, seq_id: isize, run_color: Option<usize>, range: Range<usize>) { (**self).write_run(seq_id, run_color, range) }
    fn flush(&mut self) { (**self).flush() }
}

pub struct OutputWriter<W: Write> {
    out: W,
    seq_names: Option<Vec<String>>,
    color_names: Option<Vec<String>>,
    root_id: usize,
    report_misses: bool,
    print_header: bool,
}

impl<W: Write> OutputWriter<W> {
    pub fn new(out: W, seq_names: Option<Vec<String>>, color_names: Option<Vec<String>>, root_id: usize, report_misses: bool, print_header: bool) -> Self {
        if let Some(names) = &color_names {
            for name in names {
                if name == "none" {
                    eprintln!("Error: can not use \"none\" as a color name because it is a reserved name");
                }
                if name == "root" {
                    eprintln!("Error: can not use \"root\" as a color name because it is a reserved name");
                }
            }
        }
        Self { out, seq_names, color_names, root_id, report_misses, print_header }
    }

    #[cfg(test)]
    pub fn into_inner(self) -> W {
        self.out
    }
}

impl<W: Write + Send> RunWriter for OutputWriter<W> {
    fn write_header(&mut self) {
        if self.print_header {
            let seq_col = if self.seq_names.is_some() { "query_name" } else { "query_rank" };
            let color_col = if self.color_names.is_some() { "color_name" } else { "color" };
            writeln!(self.out, "{seq_col}\tfrom_kmer\tto_kmer\t{color_col}").unwrap();
        }
    }

    fn write_run(&mut self, seq_id: isize, run_color: Option<usize>, range: Range<usize>) {
        if range.is_empty() { return; }
        if run_color.is_none() && !self.report_misses { return; }

        let from = range.start;
        let to = range.end;

        match &self.seq_names {
            Some(names) => write!(self.out, "{}", &names[seq_id as usize]).unwrap(),
            None => write!(self.out, "{seq_id}").unwrap(),
        }
        write!(self.out, "\t{from}\t{to}\t").unwrap();
        match run_color {
            None => write!(self.out, "{}", if self.color_names.is_some() { "none" } else { "-" }).unwrap(),
            Some(c) => match &self.color_names {
                Some(names) => write!(self.out, "{}", &names[c]).unwrap(),
                None => write!(self.out, "{c}").unwrap(),
            },
        }
        writeln!(self.out).unwrap();
    }

    fn flush(&mut self) {
        self.out.flush().unwrap();
    }
}


struct QueryBatch {
    seqs: SeqDB,

    batch_id: usize,
    sequence_starts: Vec<usize>, // Source sequence changes at these answer indices 
    chars_in_batch: usize,
    kmers_in_batch: usize,
    k: usize,
}

#[derive(Debug)]
struct ProcessedQueryBatch {
    result: Vec<Option<usize>>,

    batch_id: usize,
    sequence_starts: Vec<usize> 
}

impl QueryBatch {

    fn new(batch_id: usize, k: usize) -> Self {
        Self {
            seqs: SeqDB::new(),
            batch_id,
            sequence_starts: vec![],
            chars_in_batch: 0,
            kmers_in_batch: 0,
            k
        }
    }

    fn push(&mut self, seq: &[u8], extend_prev: bool) {
        if !extend_prev {
            self.sequence_starts.push(self.kmers_in_batch);
        }
        self.seqs.push_seq(seq);
        self.kmers_in_batch += kmers_in_n(self.k, seq.len());
        self.chars_in_batch += seq.len();

    }

    fn len(&self) -> usize {
        self.chars_in_batch
    }

    fn run<A: ColoredKmerLookupAlgorithm>(self, kmer_lookup_algo: &A) -> ProcessedQueryBatch {
        let total_query_kmers = self.seqs.iter().fold(0_usize, |acc, rec| 
            acc + kmers_in_n(self.k, rec.seq.len()) 
        );
        let mut color_ids = Vec::<Option<usize>>::with_capacity(total_query_kmers);

        for rec in self.seqs.iter() {
            for color in kmer_lookup_algo.lookup_kmers(rec.seq, self.k) {
                color_ids.push(color);
            }
        }

        assert_eq!(color_ids.len(), total_query_kmers);
        ProcessedQueryBatch{
            result: color_ids,
            batch_id: self.batch_id,
            sequence_starts: self.sequence_starts,
        }
    }
}

impl PartialEq for ProcessedQueryBatch{
    fn eq(&self, other: &Self) -> bool {
        self.batch_id == other.batch_id
    }
}
impl Eq for ProcessedQueryBatch {}

impl PartialOrd for ProcessedQueryBatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other)) // Using the total order
    }
}

impl Ord for ProcessedQueryBatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.batch_id.cmp(&other.batch_id)
    }
}

fn kmers_in_n(k: usize, n: usize) -> usize {
    max(0, n as isize - k as isize + 1) as usize
}

struct OutputState {
    cur_seq_id: isize,
    run_open: Option<usize>,
    run_len: usize,
    run_color: Option<usize>, 
}

fn output_batch_result(batch: &ProcessedQueryBatch, state: &mut OutputState, writer: &mut dyn RunWriter) {

    let cur_seq_id = &mut state.cur_seq_id;
    let run_open = &mut state.run_open;
    let run_len = &mut state.run_len;
    let run_color = &mut state.run_color;

    let mut starts_ptr = batch.sequence_starts.iter().peekable();
    for (i, color) in batch.result.iter().enumerate() {
        while starts_ptr.peek().is_some_and(|&&s| s == i) {
            // New sequence starts. This closes the currently open run, if exists
            if let Some(p) = run_open {
                writer.write_run(*cur_seq_id, *run_color, *p..(*p + *run_len));
                *run_open = None;
                *run_len = 0;
            }
            starts_ptr.next();
            *cur_seq_id += 1;
        }

        match run_open {
            None => {
                // We are at the start of a sequence -> open a new run
                *run_open = Some(0);
                *run_len = 1;
                *run_color = *color;
            },
            Some(p) => {
                // See if we can extend the current run
                if *run_color == *color {
                    // Extend
                    *run_len += 1;
                } else {
                    // Run ends
                    writer.write_run(*cur_seq_id, *run_color, *p .. (*p + *run_len));
                    *run_open = Some(*p + *run_len);
                    *run_len = 1;
                    *run_color = *color;
                }
            }
        }
    }

}

// Returns the total number of k-mers in all the received batches
fn output_thread(query_results: crossbeam::channel::Receiver<ProcessedQueryBatch>, writer: &mut dyn RunWriter) -> usize {
    let mut batch_buffer = std::collections::BinaryHeap::<Reverse<ProcessedQueryBatch>>::new(); // Reverse makes it a min heap
    let mut n_kmers_processed = 0_usize;
    let mut next_batch_id = 0_usize;

    let mut output_state = OutputState {
        cur_seq_id: -1,
        run_open: None, // Run of the same color (None counts as color). Value of None means that no run is active, not even a run of Nones.
        run_len: 0,
        run_color: None,
    };

    writer.write_header();
    while let Ok(batch) = query_results.recv() {
        batch_buffer.push(Reverse(batch)); // Reverse makes this a min heap

        loop { // Print all batches that can now be printed
            let min_batch = batch_buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.batch_id == next_batch_id {
                    output_batch_result(min_batch, &mut output_state, writer);
                    n_kmers_processed += min_batch.result.len();
                    batch_buffer.pop();
                    next_batch_id += 1;
                } else {
                    break; // Not ready to print min_batch yet
                }
            } else {
                break; // Batch buffer is empty -> We are done
            }
        }

    }

    // All batches processed

    // The last run of the last batch remains open (unless it's closed by the start of a
    // sequence that has no k-mers). Let's write it.
    if let Some(p) = output_state.run_open {
        writer.write_run(output_state.cur_seq_id, output_state.run_color, p..(p+output_state.run_len));
    }

    n_kmers_processed
    // Channel is dropped (= closed) here.
}

// Batch size is in nucleotides (= bytes)
pub fn lookup_parallel<A: ColoredKmerLookupAlgorithm + Send + Sync>(n_threads: usize, mut queries: impl sbwt::SeqStream + Send, kmer_lookup_algo: &A, batch_size: usize, k: usize, mut writer: impl RunWriter) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(n_threads); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<ProcessedQueryBatch>(n_threads);

    let query_start = std::time::Instant::now();
    let mut n_bases_processed = 0_usize;

    let n_kmers_processed = std::thread::scope(|s| {

        let reader_handle = s.spawn(|| {

            // Reader thread that pushes batches for workers

            let mut batch = QueryBatch::new(0, k); // Initialize an empty batch
            while let Some(seq) = queries.stream_next() {
                let n = seq.len();
                let b = batch_size;

                crate::util::process_kmers_in_pieces(seq, k, batch_size, |piece_idx, piece: &[u8]|{
                    batch.push(piece, piece_idx > 0);

                    //log::info!("Pushing batch of size {} ({} batches in channel)", batch.len(), batch_send.len());
                    if batch.len() >= b {
                        // Swap the current batch with an empty batch, and send it to processing
                        let next_batch_id = batch.batch_id + 1;
                        let mut batch_to_send = QueryBatch::new(next_batch_id, k); // Empty batch
                        std::mem::swap(&mut batch, &mut batch_to_send);
                        batch_send.send(batch_to_send).unwrap();
                    }
                });

                n_bases_processed += n;
            }

            // Push the last remaining non-full batch. Can be empty but that's ok.
            batch_send.send(batch).unwrap();

            log::info!("All input read");
            drop(batch_send); // Close the channel
        });

        let writer_handle = s.spawn(|| {
            let n_kmers = output_thread(output_recv, &mut writer); // Returns number of k-mers processed
            writer.flush(); // They say this needs to be done because errors during drop are ignored
            n_kmers
        });

        let mut worker_handles = Vec::new();
        for _ in 0..n_threads {
            let output_send_clone = output_send.clone(); // Moved into worker
            let batch_recv_clone = batch_recv.clone(); // Moved into worker
            worker_handles.push(s.spawn(move || {
                while let Ok(batch) = batch_recv_clone.recv() {
                   let processed_batch = batch.run(kmer_lookup_algo);
                   output_send_clone.send(processed_batch).unwrap();
                }
            }));
        }

        // Wait for threads to finish
        reader_handle.join().unwrap(); // All work batches pushed to workers
        worker_handles.into_iter().for_each(|w| w.join().unwrap()); // All batches processed
        drop(output_send); // All output written to channel -> can close the channel
        let n_kmers_processed = writer_handle.join().unwrap(); // All output written out

        #[allow(clippy::let_and_return)] // Is clearer to give it a name
        n_kmers_processed
    });

    let query_duration = query_start.elapsed();
    
    log::info!("{} k-mers in {} base pairs queried in {} seconds (excluding index loading time)", n_kmers_processed, n_bases_processed, query_duration.as_secs());
    log::info!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / n_kmers_processed as f64);
    log::info!("Query time per base pair: {} nanoseconds", query_duration.as_nanos() as f64 / n_bases_processed as f64);

}

#[cfg(test)]
mod tests {
    use rand_chacha::rand_core::{RngCore, SeedableRng};
    use sbwt::{BitPackedKmerSortingMem, SeqStream};

    use crate::{color_storage::SimpleColorStorage, parallel_queries::{OutputWriter, lookup_parallel}, single_colored_kmers::{ColorHierarchy, LcsWrapper, SingleColoredKmers}};

    struct SingleSeqStream {
        seq: Vec<u8>,
        pos: usize,
    }

    impl SingleSeqStream {
        fn new(seq: Vec<u8>) -> Self {
            Self { seq, pos: 0 }
        }
    }

    impl SeqStream for SingleSeqStream {
        fn stream_next(&mut self) -> Option<&[u8]> {
            if self.pos == 0 {
                self.pos = 1;
                Some(&self.seq)
            } else {
                None
            }
        }
    }

    struct MultiSeqStream {
        seqs: Vec<Vec<u8>>,
        pos: usize,
    }

    impl MultiSeqStream {
        fn new(seqs: Vec<Vec<u8>>) -> Self {
            Self { seqs, pos: 0 }
        }
    }

    impl SeqStream for MultiSeqStream {
        fn stream_next(&mut self) -> Option<&[u8]> {
            if self.pos < self.seqs.len() {
                let ret = &self.seqs[self.pos];
                self.pos += 1;
                Some(ret)
            } else {
                None
            }
        }
    }

    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    enum Color {
        NonRoot(usize),
        Root,
    }

    fn random_test(batch_size: usize) {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(125);

        let mut sequences: Vec<Vec<u8>> = Vec::new();

        let k = 7;

        // Generate DNA sequences of random length. 
        // Each k-mer in the sequences must be unique! The following algorithm ensures that,
        // but it may loop forever if we get unlucky. The RNG seed is chosen so that this does not happen.
        let mut kmer_to_color = std::collections::HashMap::<Vec<u8>, Color>::new();
        for color in 0..2 {
            let len = (rng.next_u64() % 5000 + 1) as usize;
            let mut seq: Vec<u8> = Vec::new();
            while seq.len() < len {
                let base = rng.next_u64() % 4;
                let nucleotide = match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
                seq.push(nucleotide);
                if seq.len() >= k {
                    let kmer = &seq[seq.len()-k..];
                    if let Some(old) = kmer_to_color.get(kmer) {
                        if let Color::NonRoot(old_color) = old {
                            if *old_color != color {
                                kmer_to_color.insert(kmer.to_vec(), Color::Root);
                            } // Else the color matches the old color -> keep as is
                        } // Else we have multiple colors -> keep as is
                    } else {
                        // Nothing stored yet -> store single color
                        kmer_to_color.insert(kmer.to_vec(), Color::NonRoot(color));
                    }
                }
            }
            sequences.push(seq);
        }

        // Build SBWT
        // Use in-memory construction
        eprintln!("Building SBWT...");
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false)
            .k(k)
            .build_lcs(true)
            .n_threads(3)
            .precalc_length(3)
            .add_all_dummy_paths(true)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
        .run_from_vecs(&sequences);
        let lcs = lcs.unwrap();

        let seqstreams: Vec<SingleSeqStream> = sequences.iter().map(|s| SingleSeqStream::new(s.clone())).collect();
        eprintln!("Building SingleColoredKmers...");
        let color_names: Vec<String> = (0..sequences.len()).map(|i| format!("{}", i)).collect();
        let sck = SingleColoredKmers::<LcsWrapper, SimpleColorStorage>::new(sbwt, lcs, seqstreams, 3, ColorHierarchy::new_star(color_names));
        eprintln!("SingleColoredKmers built");

        // Generate 1000 random queries of lengths between 1 and 100
        let mut queries: Vec<Vec<u8>> = Vec::new();
        let n_queries = 1000;
        for query_id in 0..n_queries {
            let len = rng.next_u64() % 100 + 1;
            let mut query: Vec<u8> = (0..len).map(|_| {
                let base = rng.next_u64() % 4;
                match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            }).collect();

            if query_id == n_queries - 1 {
                // Make sure the last query contains a known k-mer so that the color last run remains open
                query.extend_from_slice(kmer_to_color.keys().next().unwrap());
            }
            queries.push(query);
        }

        let out_vec = Vec::<u8>::new();
        let out = std::io::Cursor::new(out_vec);
        let mut writer = OutputWriter::new(out, None, None, sck.color_hierarchy().root(), false, true);
        lookup_parallel(2, MultiSeqStream::new(queries.clone()), &sck, batch_size, k, &mut writer);

        // Parse output tsv line by line
        let output_str = String::from_utf8(writer.into_inner().into_inner()).unwrap();
        let output_lines = output_str.lines();
        // For each query, the starting positions and colors of found k-mers
        let mut found_kmers: Vec::<Vec::<(usize,Color)>> = vec![Vec::new(); queries.len()]; 
        for (line_idx, line) in output_lines.enumerate() {
            if line_idx == 0 { // tsv header
                assert_eq!(line, "query_rank\tfrom_kmer\tto_kmer\tcolor");
            } else {
                let mut fields = line.split('\t');
                let seq_id: usize = fields.next().unwrap().parse().unwrap();
                let start: usize = fields.next().unwrap().parse().unwrap();
                let end: usize = fields.next().unwrap().parse().unwrap();
                let color_token = fields.next().unwrap();
                let c = color_token.parse::<usize>().unwrap();
                let color = if c == sck.color_hierarchy().root() { Color::Root } else { Color::NonRoot(c) };
                for i in start..end {
                    found_kmers[seq_id].push((i, color));
                }
            }
        }

        // Verify against known answers
        for query_id in 0..queries.len() {
            let mut true_answer = Vec::<(usize, Color)>::new();
            for (i, kmer) in queries[query_id].windows(k).enumerate() {
                if let Some(color) = kmer_to_color.get(kmer) {
                    true_answer.push((i, *color));
                }
            }
            eprintln!("{:?}", true_answer);
            eprintln!("{:?}", found_kmers[query_id]);
            assert_eq!(true_answer, found_kmers[query_id]);
        }
    }

    #[test]
    fn run_random_tests() {
        random_test(7);
        random_test(8);
        random_test(1000);
    }

    fn random_test_short_queries(batch_size: usize) {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(999);

        let index_k = 7_usize;
        let query_k = 4_usize;
        assert!(query_k < index_k);

        // Build random sequences and index them with index_k.
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        for _ in 0..2 {
            let len = (rng.next_u64() % 5000 + 1) as usize;
            let seq: Vec<u8> = (0..len).map(|_| match rng.next_u64() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            sequences.push(seq);
        }

        // Build ground truth: for each query_k-mer, the color is determined by which
        // input sequences contain that query_k-mer as a substring.
        let mut short_kmer_to_color = std::collections::HashMap::<Vec<u8>, Color>::new();
        for (color, seq) in sequences.iter().enumerate() {
            for short_kmer in seq.windows(query_k) {
                match short_kmer_to_color.get(short_kmer) {
                    None => { short_kmer_to_color.insert(short_kmer.to_vec(), Color::NonRoot(color)); }
                    Some(Color::NonRoot(old)) if *old != color => {
                        short_kmer_to_color.insert(short_kmer.to_vec(), Color::Root);
                    }
                    _ => {}
                }
            }
        }

        // Build SBWT with index_k
        eprintln!("Building SBWT...");
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false)
            .k(index_k)
            .build_lcs(true)
            .n_threads(3)
            .precalc_length(3)
            .add_all_dummy_paths(true) // Required for multi-k support
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
        .run_from_vecs(&sequences);
        let lcs = lcs.unwrap();

        let seqstreams: Vec<SingleSeqStream> = sequences.iter().map(|s| SingleSeqStream::new(s.clone())).collect();
        eprintln!("Building SingleColoredKmers...");
        let color_names: Vec<String> = (0..sequences.len()).map(|i| format!("{}", i)).collect();
        let sck = SingleColoredKmers::<LcsWrapper, SimpleColorStorage>::new(sbwt, lcs, seqstreams, 3, ColorHierarchy::new_star(color_names));
        eprintln!("SingleColoredKmers built");

        // Generate random queries of lengths between 1 and 50
        let n_queries = 1000;
        let queries: Vec<Vec<u8>> = (0..n_queries).map(|_| {
            let len = rng.next_u64() % 50 + 1;
            (0..len).map(|_| match rng.next_u64() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect()
        }).collect();

        let out_vec = Vec::<u8>::new();
        let out = std::io::Cursor::new(out_vec);
        let mut writer = OutputWriter::new(out, None, None, sck.color_hierarchy().root(), false, true);
        lookup_parallel(2, MultiSeqStream::new(queries.clone()), &sck, batch_size, query_k, &mut writer);

        let output_str = String::from_utf8(writer.into_inner().into_inner()).unwrap();
        let output_lines = output_str.lines();
        let mut found_kmers: Vec<Vec<(usize, Color)>> = vec![Vec::new(); queries.len()];
        for (line_idx, line) in output_lines.enumerate() {
            if line_idx == 0 {
                assert_eq!(line, "query_rank\tfrom_kmer\tto_kmer\tcolor");
            } else {
                let mut fields = line.split('\t');
                let seq_id: usize = fields.next().unwrap().parse().unwrap();
                let start: usize = fields.next().unwrap().parse().unwrap();
                let end: usize = fields.next().unwrap().parse().unwrap();
                let color_token = fields.next().unwrap();
                let c = color_token.parse::<usize>().unwrap();
                let color = if c == sck.color_hierarchy().root() { Color::Root } else { Color::NonRoot(c) };
                for i in start..end {
                    found_kmers[seq_id].push((i, color));
                }
            }
        }

        for query_id in 0..queries.len() {
            let mut true_answer = Vec::<(usize, Color)>::new();
            for (i, short_kmer) in queries[query_id].windows(query_k).enumerate() {
                if let Some(color) = short_kmer_to_color.get(short_kmer) {
                    true_answer.push((i, *color));
                }
            }
            assert_eq!(true_answer, found_kmers[query_id], "mismatch at query_id={}", query_id);
        }
    }

    #[test]
    fn run_random_tests_short_queries() {
        random_test_short_queries(4);
        random_test_short_queries(5);
        random_test_short_queries(1000);
    }

    #[test]
    fn serialize_roundtrip() {
        let k = 7_usize;
        let sequences: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGT".to_vec(),
            b"TTTTACGTGGGG".to_vec(),
        ];

        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false)
            .k(k)
            .build_lcs(true)
            .n_threads(1)
            .precalc_length(3)
            .add_all_dummy_paths(true)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
            .run_from_vecs(&sequences);
        let lcs = lcs.unwrap();

        let color_names = vec!["seq0".to_string(), "seq1".to_string()];
        let seqstreams: Vec<SingleSeqStream> = sequences.iter()
            .map(|s| SingleSeqStream::new(s.clone()))
            .collect();
        let original = SingleColoredKmers::<LcsWrapper, SimpleColorStorage>::new(
            sbwt, lcs, seqstreams, 1, ColorHierarchy::new_star(color_names),
        );

        // Serialize
        let mut buf = Vec::<u8>::new();
        original.serialize(&mut buf);

        // Deserialize
        let loaded = SingleColoredKmers::<LcsWrapper, SimpleColorStorage>::load(&mut buf.as_slice());

        // Check structural equality
        assert_eq!(original.k(), loaded.k());
        assert_eq!(original.n_kmers(), loaded.n_kmers());
        assert_eq!(original.color_names(), loaded.color_names());
        assert_eq!(original.color_hierarchy(), loaded.color_hierarchy());

        // Check every SBWT position has the same color
        for i in 0..original.n_kmers() {
            assert_eq!(original.get_color(i), loaded.get_color(i), "color mismatch at position {i}");
        }
    }
}
