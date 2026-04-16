//! Index construction. Consumes input sequences, marks colors using a fold
//! over the LCA (or priority-aware LCA) of co-occurring labels, and assembles
//! an [`HksIndex`]. This module is the single place that knows about
//! priorities; the index type itself is priority-unaware.

use std::ops::Range;
use std::sync::atomic::{AtomicU16, AtomicU32, AtomicU64, AtomicU8};
use std::time::{Duration, Instant};

use bitvec::prelude::*;
use crossbeam::channel::{Receiver, RecvTimeoutError};
use jseqio::seq_db::SeqDB;
use sbwt::{ContractLeft, LcsArray, MatchingStatisticsIterator, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};

use crate::color_storage::SimpleColorStorage;
use crate::lca_tree::LcaTree;
use crate::priority_lca::PriorityLca;
use crate::single_colored_kmers::{ColorHierarchy, FeatureSet, HksIndex};
use crate::traits::*;

/// Build a new index from input sequences. If `priorities` is `Some`, uses
/// priority-aware LCA during the color-merge fold (see [`PriorityLca`]);
/// otherwise uses standard LCA. Priorities are dropped once construction
/// finishes — they are not stored in the returned index.
pub fn build<L, C, T>(
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    input_streams: Vec<T>,
    n_threads: usize,
    hierarchy: ColorHierarchy,
    feature_set_name: &str,
    priorities: Option<Vec<usize>>,
) -> HksIndex<L, C>
where
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess,
    C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>,
    T: SeqStream + Send,
{
    let color_storage = mark_colors_with_priorities::<T, LcsArray>(
        &sbwt, &lcs, input_streams, n_threads, &hierarchy, priorities,
    );

    log::info!("Indexing color id array");
    let color_assignments = C::from(color_storage);
    let fs = FeatureSet { color_assignments, hierarchy, name: feature_set_name.to_owned() };
    HksIndex::<L, C>::new_given_feature_sets(sbwt, lcs, vec![fs])
}

/// Append a new feature set to an existing index.
pub fn add_feature_set<L, C, T>(
    index: &mut HksIndex<L, C>,
    input_streams: Vec<T>,
    n_threads: usize,
    hierarchy: ColorHierarchy,
    feature_set_name: &str,
    priorities: Option<Vec<usize>>,
) -> Result<(), String>
where
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync,
    C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>,
    T: SeqStream + Send,
{
    if index.feature_sets().iter().any(|fs| fs.name == feature_set_name) {
        return Err(format!("Feature set name \"{feature_set_name}\" already exists in index"));
    }

    let color_storage = mark_colors_with_priorities::<T, L>(
        index.sbwt(), index.lcs(), input_streams, n_threads, &hierarchy, priorities,
    );

    log::info!("Indexing color id array");
    let color_assignments = C::from(color_storage);
    index.push_feature_set(FeatureSet { color_assignments, hierarchy, name: feature_set_name.to_owned() });
    Ok(())
}

/// Resolve priorities (or absence thereof) into a merge closure and run
/// `mark_colors_dispatch`. This is the single choke point where priorities
/// exist; everything below sees only an `Fn(usize, usize) -> usize`.
fn mark_colors_with_priorities<T, CL>(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &CL,
    input_streams: Vec<T>,
    n_threads: usize,
    hierarchy: &ColorHierarchy,
    priorities: Option<Vec<usize>>,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    CL: ContractLeft + Sync,
{
    let required_bit_width = SimpleColorStorage::required_bit_width(hierarchy.n_nodes() + 1);
    let tree = hierarchy.tree();

    log::info!("Marking colors");
    match priorities {
        Some(p) => {
            let plca = PriorityLca::new(tree, p)
                .unwrap_or_else(|e| panic!("Invalid node priorities: {e}"));
            mark_colors_dispatch::<T, CL, _>(sbwt, lcs, input_streams, n_threads, required_bit_width, tree, &|a, b| plca.plca(a, b))
        }
        None => {
            mark_colors_dispatch::<T, CL, _>(sbwt, lcs, input_streams, n_threads, required_bit_width, tree, &|a, b| tree.lca(a, b))
        }
    }
}

/// Thin dispatcher that picks the atomic int width for the color vector and
/// calls `mark_colors`. Extracted so the per-closure monomorphization of
/// `mark_colors` lives in one place.
fn mark_colors_dispatch<T, CL, F>(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &CL,
    input_streams: Vec<T>,
    n_threads: usize,
    required_bit_width: usize,
    color_hierarchy: &LcaTree,
    merge: &F,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    CL: ContractLeft + Sync,
    F: Fn(usize, usize) -> usize + Sync,
{
    if required_bit_width <= 8 {
        mark_colors::<T, Vec<AtomicU8>, CL, F>(sbwt, lcs, input_streams, n_threads, color_hierarchy, merge)
    } else if required_bit_width <= 16 {
        mark_colors::<T, Vec<AtomicU16>, CL, F>(sbwt, lcs, input_streams, n_threads, color_hierarchy, merge)
    } else if required_bit_width <= 32 {
        mark_colors::<T, Vec<AtomicU32>, CL, F>(sbwt, lcs, input_streams, n_threads, color_hierarchy, merge)
    } else {
        mark_colors::<T, Vec<AtomicU64>, CL, F>(sbwt, lcs, input_streams, n_threads, color_hierarchy, merge)
    }
}

fn mark_colors<T, A, CL, F>(
    sbwt: &SbwtIndex<SubsetMatrix>,
    lcs: &CL,
    input_streams: Vec<T>,
    n_threads: usize,
    color_hierarchy: &LcaTree,
    merge: &F,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    A: AtomicColorVec + Send + Sync,
    CL: ContractLeft + Sync,
    F: Fn(usize, usize) -> usize + Sync,
{
    let color_ids = A::new(sbwt.n_sets());
    let si = StreamingIndex {
        extend_right: sbwt,
        contract_left: lcs,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };
    let n_colors = color_hierarchy.n_nodes();

    let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let n_bases_processed = AtomicU64::new(0);
    std::thread::scope(|scope| { thread_pool.install(|| {

        let (quit_print_send, quit_print_recv) = crossbeam::channel::unbounded::<bool>();
        let _progress_printer = scope.spawn({
            let n_bases_processed = &n_bases_processed;
            move || { progress_print_thread(n_bases_processed, quit_print_recv); }
        });

        let (batch_send, batch_recv) = crossbeam::channel::bounded::<ColoringBatch>(4);

        let reader_handle = scope.spawn(move || {
            let mut batch = ColoringBatch { dbs: vec![], dummy_mer_dbs: vec![], total_len: 0 };
            let b = 10000;
            let k = sbwt.k();
            for (color, mut stream) in input_streams.into_iter().enumerate() {
                while let Some(seq) = stream.stream_next() {
                    crate::util::for_each_run_with_key(seq, |c| IS_DNA[*c as usize], |mut run_range: Range<usize>| {
                        if !run_range.is_empty() && IS_DNA[seq[run_range.start] as usize] {
                            if run_range.len() >= k {
                                run_range = run_range.start..run_range.start + (k - 1);
                            }
                            let mer = &seq[run_range.clone()];
                            batch.push_dummy_mer(color, mer);
                        }
                    });

                    crate::util::process_kmers_in_pieces(seq, k, b, |_piece_idx, piece: &[u8]| {
                        batch.push(color, piece);
                        if batch.total_len >= b {
                            let mut batch_to_send = ColoringBatch { dbs: vec![], dummy_mer_dbs: vec![], total_len: 0 };
                            std::mem::swap(&mut batch, &mut batch_to_send);
                            batch_send.send(batch_to_send).unwrap();
                        }
                    });
                }
            }
            if batch.total_len > 0 {
                batch_send.send(batch).unwrap();
            }
        });

        let mut worker_handles = Vec::new();
        for _ in 0..n_threads {
            let batch_recv_clone = batch_recv.clone();
            let si_ref = &si;
            let n_bases_processed_ref = &n_bases_processed;
            let color_ids_ref = &color_ids;
            worker_handles.push(scope.spawn(move || {
                while let Ok(batch) = batch_recv_clone.recv() {
                    batch.run(si_ref, color_ids_ref, n_bases_processed_ref, merge);
                }
            }));
        }

        reader_handle.join().unwrap();
        for w in worker_handles {
            w.join().unwrap();
        }
        quit_print_send.send(true).unwrap();
    })});

    log::info!("Bitpacking color id array");
    let mut compressed_colors = SimpleColorStorage::new(sbwt.n_sets(), n_colors);
    let mut total_some_count = 0_usize;
    let mut total_none_count = 0_usize;
    for i in 0..sbwt.n_sets() {
        let cv = color_ids.read(i);
        match cv {
            Some(_) => total_some_count += 1,
            None => total_none_count += 1,
        }
        compressed_colors.set_color(i, cv);
    }

    log::info!("Colored {total_some_count} sbwt positions");
    log::info!("{total_none_count} sbwt positions left uncolored");

    compressed_colors
}

fn progress_print_thread(n_bases_processed: &AtomicU64, quit_signal: Receiver<bool>) {
    log::info!("Processing up to 2n bases, where n is the number of bases in the input.");
    let print_interval = 10; // seconds
    let mut last_print_time = Instant::now();
    let mut last_count = 0_u64;
    loop {
        match quit_signal.recv_timeout(Duration::from_secs(print_interval)) {
            Ok(_) => return,
            Err(RecvTimeoutError::Timeout) => {
                let count = n_bases_processed.load(std::sync::atomic::Ordering::Relaxed);
                let dcount = count - last_count;
                let elapsed = last_print_time.elapsed().as_secs_f64();
                log::info!("{dcount} bases processed (total {count}) ({:.2} Mbases/sec)", dcount as f64 / elapsed / 1e6);
                last_count = count;
                last_print_time = Instant::now();
            }
            Err(RecvTimeoutError::Disconnected) => return,
        }
    }
}

// This bit vector of length 256 marks the ascii values of acgtACGT.
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

impl<T: AtomicUint> AtomicColorVec for Vec<T> {
    fn new(len: usize) -> Self {
        (0..len).map(|_| T::new(T::max_value())).collect()
    }

    fn update<F: Fn(usize, usize) -> usize>(&self, i: usize, x: usize, merge: &F) {
        assert!(x != Self::none_sentinel(), "x must not be the none sentinel");
        self[i].fetch_update(|cur| {
            if cur == Self::none_sentinel() {
                x
            } else {
                merge(cur, x)
            }
        });
    }

    fn read(&self, i: usize) -> Option<usize> {
        let x = self[i].load(std::sync::atomic::Ordering::Relaxed);
        if x == Self::none_sentinel() { None } else { Some(x) }
    }

    fn none_sentinel() -> usize {
        T::max_value()
    }
}

struct ColoringBatch {
    dbs: Vec<(usize, SeqDB)>,
    dummy_mer_dbs: Vec<(usize, SeqDB)>,
    total_len: usize,
}

impl ColoringBatch {
    fn push(&mut self, color: usize, seq: &[u8]) {
        let mut extended = false;
        if let Some((last_color, last_db)) = self.dbs.last_mut() {
            if *last_color == color {
                last_db.push_seq(seq);
                extended = true;
            }
        }
        if !extended {
            let mut db = SeqDB::new();
            db.push_seq(seq);
            self.dbs.push((color, db));
        }
        self.total_len += seq.len();
    }

    fn push_dummy_mer(&mut self, color: usize, mer: &[u8]) {
        let mut extended = false;
        if let Some((last_color, last_db)) = self.dummy_mer_dbs.last_mut() {
            if *last_color == color {
                last_db.push_seq(mer);
                extended = true;
            }
        }
        if !extended {
            let mut db = SeqDB::new();
            db.push_seq(mer);
            self.dummy_mer_dbs.push((color, db));
        }
        self.total_len += mer.len();
    }

    fn run<V, CL, F>(&self, si: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, CL>, color_ids: &V, progress_counter: &AtomicU64, merge: &F)
    where
        V: AtomicColorVec,
        CL: ContractLeft,
        F: Fn(usize, usize) -> usize,
    {
        let k = si.k;
        let mut thread_progress = 0_usize;

        for (color, db) in self.dbs.iter() {
            for rec in db.iter() {
                let seq = rec.seq;
                let ms = si.matching_statistics_iter(seq);
                ms.enumerate().for_each(|(i, (len, range))| {
                    if len == k {
                        debug_assert!(range.len() == 1);
                        color_ids.update(range.start, *color, merge);
                    } else if cfg!(debug_assertions) && i >= k - 1 {
                        let kmer = &seq[i - (k - 1)..=i];
                        let all_acgt = kmer.iter().all(|c| IS_DNA[*c as usize]);
                        if all_acgt {
                            panic!("Error: k-mer {} not found in sbwt", String::from_utf8_lossy(kmer));
                        }
                    }
                    thread_progress += 1;
                    if thread_progress == 10000 {
                        progress_counter.fetch_add(10000, std::sync::atomic::Ordering::Relaxed);
                        thread_progress = 0;
                    }
                });
            }
        }

        for (color, db) in self.dummy_mer_dbs.iter() {
            for rec in db.iter() {
                let mer = rec.seq;
                let ms = si.matching_statistics_iter(mer);
                ms.for_each(|(len, range)| {
                    assert!(range.len() > 0);
                    assert!(len < k);
                    let colex = range.start;
                    color_ids.update(colex, *color, merge);
                });
            }
        }

        progress_counter.fetch_add(thread_progress as u64, std::sync::atomic::Ordering::Relaxed);
    }
}
