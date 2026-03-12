use std::ops::Range;

// Splits `seq` into pieces of length at most `max_piece_len`, such that the pieces overlap
// by k-1 characters. Calls f on pairs (piece_idx, piece). 
// If seq has length less than k, calls f with a single empty piece.
pub fn process_kmers_in_pieces<F: FnMut(usize, &[u8])>(seq: &[u8], k: usize, max_piece_len: usize, mut f: F) {
    let n = seq.len();
    let b = max_piece_len;

    if n < k {
        // No full k-mers in this sequence.
        f(0, b"");
        return;
    }

    // Let b be batch size and n be the length of the sequence.
    // Split the sequence into m pieces of length b except for the
    // last sequence that can have a shorter length, such that the
    // pieces overlap by k-1 characters and cover the whole sequence.
    // How many pieces will be there be? That is, what is the smallest
    // m such that b + (m-1)*(b-(k-1)) >= n? We must have b-k+1 > 0, or 
    // otherwise the inequality flips the wrong way around.
    // Assuming b-k+1 > 0, the solution is: m >= (n-k+1) / (b-k+1).
    // So we take the ceil of the righthand side.

    assert!(b as isize - k as isize + 1 > 0); // b-k+1 > 0
    let m = (n-k+1).div_ceil(b-k+1);
    assert!(m > 0); // This should be true since we checked for n < k earlier

    for piece_idx in 0..m {
        let pieces_before = piece_idx;
        let start = b*piece_idx - (k-1)*pieces_before;
        let piece = if piece_idx < m-1 {
            // Not the last piece: has full length b
            &seq[start..start+b]
        } else {
            &seq[start..] // Until the end (can have length shorter than b)
        };

        f(piece_idx, piece);
    }
}

pub fn for_each_run_with_key<T: Eq, KeyType: Eq, F1: Fn(&T) -> KeyType, F2: FnMut(Range<usize>)>(items: &[T], key_fn: F1, mut callback: F2) {
    if items.is_empty() { return }

    let mut run_start = 0;
    let n = items.len();
    for i in 1..n {
        if key_fn(&items[i]) != key_fn(&items[i-1]) {
            callback(run_start..i);
            run_start = i;
        }
    }
    // Final run
    callback(run_start..n);
}

/// Returns a mutable slice for every region in `regions`.
///
/// ## Preconditions  (guaranteed by the caller)
/// * All ranges are inside `0..v.len()`.
/// * The ranges are sorted by `start` and do **not** overlap.
pub(crate) fn split_to_mut_regions<'a>(
    v: &'a mut [u64],
    regions: &[Range<usize>],
) -> Vec<&'a mut [u64]> {
    let mut result = Vec::with_capacity(regions.len());
    let mut tail: &mut [u64] = v;
    let mut consumed = 0; // absolute index we have reached in `v`

    for r in regions {
        // translate the absolute `start` to an index inside `tail`
        let rel_start = r.start - consumed;
        let len       = r.end   - r.start;

        // First split off everything before the wanted region …
        let (_, after_start)    = tail.split_at_mut(rel_start);
        // ... then split that remainder into the desired region and the rest.
        let (region, after_end) = after_start.split_at_mut(len);

        result.push(region); // keep the region
        tail = after_end; // keep working with the suffix
        consumed = r.end; // advance the absolute cursor
    }

    result
}