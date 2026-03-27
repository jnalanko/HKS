use std::{io::{Read, Write}, ops::Range};
use bitvec::prelude::*;

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

pub fn serialize_bitvec_u64(bv: &BitVec::<u64, Lsb0>, mut out: impl Write) {
    // Serialize using the same format as used by serde + bincode
    let type_id = b"bitvec::order::Lsb0";
    bincode::serialize_into(&mut out, &type_id.len()).unwrap();
    out.write_all(type_id).unwrap();
    let word_size: u8 = 64;
    bincode::serialize_into(&mut out, &word_size).unwrap();
    let what_is_this: u8 = 0; // I don't know what this byte does, but it's there in the bincode format.
    bincode::serialize_into(&mut out, &what_is_this).unwrap();

    let n_real_bits = bv.len();

    bincode::serialize_into(&mut out, &n_real_bits).unwrap();
    let n_words = n_real_bits.div_ceil(word_size as usize);
    bincode::serialize_into(&mut out, &n_words).unwrap();
    let words = bv.as_raw_slice();
    let words_as_bytes: &[u8] = bytemuck::cast_slice(words);
    out.write_all(words_as_bytes).unwrap();
}

pub fn load_bitvec_u64(mut input: impl Read) -> BitVec::<u64, Lsb0> {
    let type_id_len: u64 = bincode::deserialize_from(&mut input).unwrap();
    let mut type_id_bytes = vec![0u8; type_id_len as usize];
    input.read_exact(&mut type_id_bytes).unwrap();
    let type_id_str = std::str::from_utf8(&type_id_bytes).unwrap();
    assert_eq!(type_id_str, "bitvec::order::Lsb0");

    let word_size: u8 = bincode::deserialize_from(&mut input).unwrap();
    assert_eq!(word_size, 64);
    let what_is_this: u8 = bincode::deserialize_from(&mut input).unwrap();
    assert_eq!(what_is_this, 0); // I don't know what this is, but it's zero in the bincode format

    let n_real_bits: usize = bincode::deserialize_from(&mut input).unwrap();
    let n_words: usize = bincode::deserialize_from(&mut input).unwrap();

    // Deserialize words
    let mut words = vec![0u64; n_words];
    let words_as_bytes = bytemuck::cast_slice_mut(&mut words);
    input.read_exact(words_as_bytes).unwrap();

    let mut bv = BitVec::<u64, Lsb0>::from_vec(words);
    bv.truncate(n_real_bits); 

    bv
}

#[cfg(test)]
mod tests {

    use super::*;
    use rand::{SeedableRng, Rng};

    #[test]
    fn bitvec_custom_serialize_and_store() {
        // Check that our custom serialize and load match for the serde + bincode implementation.
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let n_bits = 1100;
        let bv: BitVec::<u64, Lsb0> = (0..n_bits).map(|_| rng.gen_bool(0.5)).collect();

        // Round trip with our functions
        let mut serialized = Vec::new();
        serialize_bitvec_u64(&bv, &mut serialized);
        // Print the serialized data as hex for debugging
        for byte in &serialized {
            eprint!("{:02x} ", byte);
        }
        eprintln!();

        let deserialized = load_bitvec_u64(&serialized[..]);
        assert_eq!(bv, deserialized);

        // Check that the custom serialize format matches the serde + bincode format.
        let mut bincode_serialized = Vec::new();
        bincode::serialize_into(&mut bincode_serialized, &bv).unwrap();
        assert_eq!(serialized, bincode_serialized);
        
    }
}