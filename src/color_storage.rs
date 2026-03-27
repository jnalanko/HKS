use crate::lca_tree::LcaTree;
use crate::traits::*;

use bitvec_sds::traits::RandomAccessU32;
use rayon::iter::{ParallelBridge, ParallelIterator};
use serde::{Serialize, Deserialize};
use bitvec::prelude::*;
use std::io::{Read, Write};
use std::ops::Range;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleColorStorage {
    colors: BitVec::<u64, Lsb0>,
    bits_per_color: usize,
    n_colors: usize, // Number of colors, not including the special "none" color.
}

impl bitvec_sds::traits::RandomAccessU32 for SimpleColorStorage {
    fn len(&self) -> usize {
        self.colors.len() / self.bits_per_color
    }

    fn get(&self, idx: usize) -> u32 {
        self.colors[idx*self.bits_per_color .. (idx+1)*self.bits_per_color].load_le()
    }
}

impl MySerialize for SimpleColorStorage {
    fn serialize(&self, mut out: &mut impl Write) {
        bincode::serialize_into(&mut out, &self.n_colors).unwrap();
        bincode::serialize_into(&mut out, &self.bits_per_color).unwrap();

        // Use our own serialization for this because the bincode deserializer
        // uses a doubling buffer which can have 2x overhead.
        crate::util::serialize_bitvec_u64(&self.colors, &mut out);
    }

    fn load(mut input: &mut impl Read) -> Box<Self> {
        let n_colors: usize = bincode::deserialize_from(&mut input).unwrap();
        let bits_per_color: usize = bincode::deserialize_from(&mut input).unwrap();

        // Use our own deserialization for this because the bincode deserializer
        // uses a doubling buffer which can have 2x overhead.
        let colors: BitVec<u64, Lsb0> = crate::util::load_bitvec_u64(&mut input);

        Box::new(SimpleColorStorage { n_colors, colors, bits_per_color })
    }
}

impl ColorStorage for SimpleColorStorage {

    fn get_color(&self, colex: usize) -> Option<usize> {
        Self::get_color_from_slice(&self.colors, self.bits_per_color, colex)
    }

    fn set_color(&mut self, colex: usize, value: Option<usize>) {
        Self::set_color_in_slice(&mut self.colors, self.bits_per_color, colex, value);
    }

    fn get_color_of_range(&self, range: Range<usize>, color_hierarchy: &LcaTree) -> Option<usize> {
        // This is O(|range| in the worst case)

        if range.is_empty() { return None }

        let mut lca: Option<usize> = None;
        for colex in range {
            lca = color_hierarchy.lca_options(lca, self.get_color(colex));
            if lca == Some(color_hierarchy.root()) {
                // Already at root -> can early exit here
                return lca
            }
        }
        lca
    }

    fn substitute_lca_for_s_mer_ranges<L: LcsAccess + Send + Sync>(&mut self, s: usize, hierarchy: &LcaTree, lcs: &L, n_threads: usize) {
        let n = self.len(); // Number of elements
        let n_bits = n * self.bits_per_color;
        let total_words = n_bits.next_multiple_of(64) / 64;
        let block_size_bits = n_bits.div_ceil(n_threads).next_multiple_of(64*self.bits_per_color);
        let block_size_words = block_size_bits / 64;

        let mut word_ranges = Vec::<Range<usize>>::new();
        for b in 0..n_threads {
            let start = b * block_size_words;
            if start >= total_words { break; }
            let end = ((b+1) * block_size_words).min(total_words);
            word_ranges.push(start..end);
        }
        let raw_data = self.colors.as_raw_mut_slice();
        assert!(raw_data.len() == total_words);
        let mut color_slices = crate::util::split_to_mut_regions(raw_data, &word_ranges);

        // Compute and fill in the LCA of all colors in the given colex range,
        // addressed relative to start_element_colex within bv.
        let fill_lca_range = |bv: &mut BitSlice<u64, Lsb0>, bits_per_color: usize, start_element_colex: usize, run: Range<usize>| {
            let mut combined: Option<usize> = None;
            for colex in run.clone() {
                let color = SimpleColorStorage::get_color_from_slice(bv, bits_per_color, colex - start_element_colex);
                combined = hierarchy.lca_options(combined, color);
            }
            for colex in run {
                SimpleColorStorage::set_color_in_slice(bv, bits_per_color, colex - start_element_colex, combined);
            }
        };

        color_slices.iter_mut().enumerate().par_bridge().for_each(|(slice_idx, slice)| {
            let bv = bitvec::slice::BitSlice::from_slice_mut(slice);
            let start_element_colex = word_ranges[slice_idx].start * 64 / self.bits_per_color;
            // The last block may be padded to a full word, so cap at the true element count.
            let n_elements = (bv.len() / self.bits_per_color).min(n - start_element_colex);

            // Sweep through every maximal run of positions whose consecutive LCS >= s
            // (i.e. all k-mers in the run share a common s-mer). Compute the LCA of all
            // colors in the run and write it back to every position in the run.
            let mut run_colex_start = start_element_colex;
            for rel_element in 1..=n_elements {
                let run_colex_end = start_element_colex + rel_element;
                let run_continues = rel_element < n_elements && lcs.get_lcs(run_colex_end) >= s;
                if !run_continues {
                    if run_colex_end - run_colex_start > 1 { // Avoid wasted work: only need to do LCA for runs longer than 1
                        fill_lca_range(bv, self.bits_per_color, start_element_colex, run_colex_start..run_colex_end);
                    }
                    run_colex_start = run_colex_end;
                }
            }
        });

        // Handle runs that cross block split points sequentially.
        // The parallel phase wrote partial LCAs on each side; since LCA is associative
        // we can re-read those values, find the full run extent, take the LCAs, and write back.
        let bv = BitSlice::from_slice_mut(self.colors.as_raw_mut_slice());
        for word_range in word_ranges {
            let start_element = word_range.start * 64 / self.bits_per_color;

            // Find full extent of the cross-boundary run
            let mut run_start = start_element;
            while run_start > 0 && lcs.get_lcs(run_start) >= s {
                run_start -= 1;
            }
            let mut run_end = start_element + 1;
            while run_end < n && lcs.get_lcs(run_end) >= s {
                run_end += 1;
            }

            fill_lca_range(bv, self.bits_per_color, 0, run_start..run_end);
        }

    }
}

impl SimpleColorStorage {

    fn get_color_from_slice(slice: &BitSlice<u64, Lsb0>, bits_per_color: usize, colex: usize) -> Option<usize> {
        let x: usize = slice[colex*bits_per_color .. (colex+1)*bits_per_color].load_le();
        if x == (1 << bits_per_color) - 1 { // Max value is reserved for None
            None
        } else {
            Some(x)
        }
    }

    fn set_color_in_slice(slice: &mut BitSlice<u64, Lsb0>, bits_per_color: usize, colex: usize, value: Option<usize>) {
        let x = match value {
            None => (1 << bits_per_color) - 1,
            Some(x) => {
                assert!(x < (1 << bits_per_color) - 1);
                x
            }
        };

        slice[colex*bits_per_color .. (colex+1)*bits_per_color].store_le(x);
    }

    pub fn new(len: usize, n_colors: usize) -> Self {
        let bits_per_color = Self::required_bit_width(n_colors);
        SimpleColorStorage {
            n_colors,
            colors: bitvec![u64, Lsb0; 0; len * bits_per_color],
            bits_per_color,
        }
    }

    pub fn required_bit_width(n_colors: usize) -> usize {
        log2_ceil(n_colors + 1) // +1 is for the special "none" value 
    }

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }

}

fn log2_ceil(x: usize) -> usize {
    assert!(x > 0);
    let mut v = 1_usize;
    let mut bits = 0_usize;
    while v < x {
        v <<= 1;
        bits += 1;
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log2_ceil(){
        assert_eq!(log2_ceil(1), 0);
        assert_eq!(log2_ceil(2), 1);
        assert_eq!(log2_ceil(3), 2);
        assert_eq!(log2_ceil(4), 2);
        assert_eq!(log2_ceil(5), 3);
        assert_eq!(log2_ceil(6), 3);
        assert_eq!(log2_ceil(7), 3);
        assert_eq!(log2_ceil(8), 3);
    }
}