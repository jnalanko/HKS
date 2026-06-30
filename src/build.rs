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
use crate::single_colored_kmers::{ColorHierarchy, HksBase, HksIndex, Labeling};
use crate::traits::*;

/// Build a new labeling from an existing index (sbwt + lcs) and input streams.
/// The result can be serialized to a standalone labeling file.
pub fn build_labeling<L, C, T>(
    base: &HksBase<L>,
    input_streams: Vec<T>,
    n_threads: usize,
    hierarchy: ColorHierarchy,
    labeling_name: &str,
    priorities: Option<Vec<usize>>,
) -> Labeling<C>
where
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync,
    C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>,
    T: SeqStream + Send,
{
    let color_storage = mark_colors_with_priorities::<T, C, L>(
        base, input_streams, n_threads, &hierarchy, priorities,
    );

    log::info!("Indexing color id array");
    let color_assignments = C::from(color_storage);
    Labeling { color_assignments, hierarchy, name: labeling_name.to_owned() }
}

/// Resolve priorities (or absence thereof) into a merge closure and run
/// `mark_colors_dispatch`. This is the single choke point where priorities
/// exist; everything below sees only an `Fn(usize, usize) -> usize`.
fn mark_colors_with_priorities<T, C, L>(
    base: &HksBase<L>,
    input_streams: Vec<T>,
    n_threads: usize,
    hierarchy: &ColorHierarchy,
    priorities: Option<Vec<usize>>,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>,
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync,
{
    let required_bit_width = SimpleColorStorage::required_bit_width(hierarchy.n_nodes() + 1);
    let tree = hierarchy.tree();

    log::info!("Marking colors");
    match priorities {
        Some(p) => {
            let plca = PriorityLca::new(tree, p)
                .unwrap_or_else(|e| panic!("Invalid node priorities: {e}"));
            let plca_override = |a,b| Some(plca.plca(a, b));
            mark_colors_dispatch::<T, L, _>(base, input_streams, n_threads, required_bit_width, tree, &plca_override)
        }
        None => {
            let no_override = |_: usize, _: usize| None;
            mark_colors_dispatch::<T, L, _>(base, input_streams, n_threads, required_bit_width, tree, &no_override)
        }
    }
}

/// Thin dispatcher that picks the atomic int width for the color vector and
/// calls `mark_colors`. Extracted so the per-closure monomorphization of
/// `mark_colors` lives in one place. The lca_override function takes a pair
/// of nodes and returns Some(node) if we want to override the LCA with that
/// node instead, otherwise None.
fn mark_colors_dispatch<T, L, F>(
    base: &HksBase<L>,
    input_streams: Vec<T>,
    n_threads: usize,
    required_bit_width: usize,
    color_hierarchy: &LcaTree,
    lca_override: &F,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync,
    F: Fn(usize, usize) -> Option<usize> + Sync,
{
    if required_bit_width <= 8 {
        mark_colors::<T, Atomic64BitAlignedColorBuf<AtomicU8>, L, F>(base, input_streams, n_threads, color_hierarchy, lca_override)
    } else if required_bit_width <= 16 {
        mark_colors::<T, Atomic64BitAlignedColorBuf<AtomicU16>, L, F>(base, input_streams, n_threads, color_hierarchy, lca_override)
    } else if required_bit_width <= 32 {
        mark_colors::<T, Atomic64BitAlignedColorBuf<AtomicU32>, L, F>(base, input_streams, n_threads, color_hierarchy, lca_override)
    } else {
        mark_colors::<T, Atomic64BitAlignedColorBuf<AtomicU64>, L, F>(base, input_streams, n_threads, color_hierarchy, lca_override)
    }
}

fn mark_colors<T, A, L, F>(
    base: &HksBase<L>,
    input_streams: Vec<T>,
    n_threads: usize,
    color_hierarchy: &LcaTree,
    lca_override: &F,
) -> SimpleColorStorage
where
    T: SeqStream + Send,
    A: AtomicColorVec + Send + Sync,
    L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess + Sync,
    F: Fn(usize, usize) -> Option<usize> + Sync,
{

    let sbwt = base.sbwt();
    let lcs = base.lcs();
    let color_ids = A::new(sbwt.n_sets());
    let si = StreamingIndex {
        extend_right: sbwt,
        contract_left: lcs,
        n: sbwt.n_sets(),
        k: sbwt.k(),
    };
    let n_colors = color_hierarchy.n_nodes();


    log::info!("Computing dummy node marks");
    let dummy_marks = si.extend_right.compute_dummy_node_marks();
    let dummy_marks_ref = &dummy_marks; // To capture by value into a closure

    log::info!("Coloring");

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
                            if run_range.len() >= k-1 {
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
                    batch.run(si_ref, color_ids_ref, n_bases_processed_ref, color_hierarchy, lca_override, dummy_marks_ref);
                }
            }));
        }

        reader_handle.join().unwrap();
        for w in worker_handles {
            w.join().unwrap();
        }
        quit_print_send.send(true).unwrap();
    })});

    log::info!("Bitpacking color id array in place");
    let (compressed_colors, total_some_count, total_none_count) =
        color_ids.into_color_storage(n_colors);

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

/// An atomic color vector backed by a `Vec<AtomicU64>` so that the same allocation can
/// be reinterpreted as a packed `BitVec<u64, Lsb0>` at the end without copying. The
/// backing `AtomicU64`s give the storage the `UnsafeCell` semantics that interior
/// mutation through any `&AtomicU*` sub-view requires. Workers see a `&[T]` view
/// (where `T` is `AtomicU8/16/32/64`); after they join, the buffer is converted to
/// `Vec<u64>` via `AtomicU64::into_inner` (zero-realloc thanks to the std `Map`/
/// `Collect` specialization for same-size/align element types) for the in-place pack.
pub struct Atomic64BitAlignedColorBuf<T: AtomicUint> {
    buf: Vec<AtomicU64>, // sized to ceil(n * sizeof(T) / 8) words; each initialized to u64::MAX
    n: usize,            // number of logical elements
    _marker: std::marker::PhantomData<T>,
}

impl<T: AtomicUint> Atomic64BitAlignedColorBuf<T> {
    fn atomic_slice(&self) -> &[T] {
        // SAFETY: T is AtomicU8/16/32/64, which has the same in-memory representation
        // as the corresponding integer. The buffer is a `Vec<AtomicU64>` whose base
        // pointer is aligned to 8 (>= align_of::<T>()), and `AtomicU64` contains
        // `UnsafeCell<u64>` so the storage permits interior mutation through any
        // `&AtomicU*` re-borrow. The first `n * size_of::<T>()` bytes are initialized
        // (each `AtomicU64` is u64::MAX, so every sub-word is T::max_value()). During
        // the worker phase we only access through this `&[T]` view, so no mixed-size
        // atomic accesses happen concurrently.
        unsafe { std::slice::from_raw_parts(self.buf.as_ptr() as *const T, self.n) }
    }
}

impl<T: AtomicUint> AtomicColorVec for Atomic64BitAlignedColorBuf<T> {
    fn new(len: usize) -> Self {
        let bytes = len.checked_mul(std::mem::size_of::<T>()).expect("color buffer size overflow");
        let words = bytes.div_ceil(8);
        // u64::MAX = all-1s, so every sub-word of width sizeof(T) reads as T::max_value(),
        // i.e. the "no color assigned" sentinel.
        let buf: Vec<AtomicU64> = (0..words).map(|_| AtomicU64::new(u64::MAX)).collect();
        Self {
            buf,
            n: len,
            _marker: std::marker::PhantomData,
        }
    }

    fn update<F: Fn(usize, usize) -> Option<usize>>(&self, i: usize, x: usize, hierarchy: &LcaTree, lca_override: &F) {
        assert!(x != Self::none_sentinel(), "x must not be the none sentinel");
        self.atomic_slice()[i].fetch_update(|cur| {
            if cur == Self::none_sentinel() {
                x
            } else if let Some(y) = lca_override(cur, x) {
                y
            } else {
                hierarchy.lca(cur, x)
            }
        });
    }

    fn read(&self, i: usize) -> Option<usize> {
        let x = self.atomic_slice()[i].load(std::sync::atomic::Ordering::Relaxed);
        if x == Self::none_sentinel() { None } else { Some(x) }
    }

    fn none_sentinel() -> usize {
        T::max_value()
    }

    fn into_color_storage(self, n_colors: usize) -> (SimpleColorStorage, usize, usize) {
        let n = self.n;
        let bits_per_color = SimpleColorStorage::required_bit_width(n_colors);
        let none_packed: usize = (1usize << bits_per_color) - 1;
        let elem_bits = std::mem::size_of::<T>() * 8;
        let none_unpacked = T::max_value();

        // Caller (the dispatch above) picks the atomic width so this holds.
        assert!(bits_per_color <= elem_bits);

        // Drop the atomic wrapper. `AtomicU64` and `u64` have identical size and align,
        // so std's in-place collect specialization reuses the same allocation.
        let mut buf: Vec<u64> = self.buf.into_iter().map(AtomicU64::into_inner).collect();

        let mut n_colored = 0_usize;
        let mut n_uncolored = 0_usize;

        // Forward in-place pack. For each i:
        //   read   bits [i*elem_bits .. (i+1)*elem_bits]   (the unpacked value)
        //   write  bits [i*bits_per_color .. (i+1)*bits_per_color]
        // Since bits_per_color <= elem_bits, the write range never overlaps any element
        // we haven't read yet.
        if n > 0 {
            let slice = BitSlice::<u64, Lsb0>::from_slice_mut(&mut buf);
            for i in 0..n {
                let r = i * elem_bits;
                let raw: usize = slice[r..r + elem_bits].load_le();
                let packed = if raw == none_unpacked {
                    n_uncolored += 1;
                    none_packed
                } else {
                    n_colored += 1;
                    raw
                };
                let w = i * bits_per_color;
                slice[w..w + bits_per_color].store_le(packed);
            }
        }

        // Trim the backing Vec to just the words we need, then release the tail.
        let total_bits = n * bits_per_color;
        let total_words = total_bits.div_ceil(64);
        buf.truncate(total_words);
        buf.shrink_to_fit();

        let mut colors = BitVec::<u64, Lsb0>::from_vec(buf);
        colors.truncate(total_bits);

        // Zero the dead bits in the last word so the backing storage is a pure
        // function of the logical contents (matters for byte-identical serialization).
        colors.set_uninitialized(false);

        (SimpleColorStorage::from_packed(colors, n_colors), n_colored, n_uncolored)
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

    fn run<V, CL, F>(&self, si: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, CL>, color_ids: &V, progress_counter: &AtomicU64, color_hierarchy: &LcaTree, lca_override: &F, dummy_marks: &BitSlice)
    where
        V: AtomicColorVec,
        CL: ContractLeft,
        F: Fn(usize, usize) -> Option<usize>,
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
                        color_ids.update(range.start, *color, color_hierarchy, lca_override);
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
                    assert!(dummy_marks[colex]); // This must be dummy (otherwise not all dummies are present in the SBWT)
                    color_ids.update(colex, *color, color_hierarchy, lca_override);
                });
            }
        }

        progress_counter.fetch_add(thread_progress as u64, std::sync::atomic::Ordering::Relaxed);
    }
}
