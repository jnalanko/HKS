use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use jseqio::{reader::DynamicFastXReader, reverse_complement_in_place};
use sbwt::SeqStream;

pub struct ChainedInputStream{
    paths: Vec<PathBuf>,
    cur_file: Option<DynamicFastXReader>,
    seq_buf: Vec<u8>, // Local buffer from which we can borrow (can not use the buffer of cur_file for lifetime reasons)
    cur_file_idx : usize, // Index of the db currently being iterated over
}

impl ChainedInputStream {
    pub fn new(filenames: Vec<PathBuf>) -> Self {
        let first_file = filenames.first().map(|f| DynamicFastXReader::from_file(f)
            .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", f.display())));
        Self {paths: filenames, cur_file: first_file, seq_buf: vec![], cur_file_idx: 0}
    }
}

/// Opens the input file when the first sequence is read
pub struct LazyFileSeqStream { 
    path: PathBuf,
    stream: Option<DynamicFastXReader>,
    seq_buf: Vec<u8>, // Local buffer from which we can borrow (can not use the buffer of cur_file for lifetime reasons)
    revcomps_enabled: bool,
    revcomp_next: bool,
}

impl LazyFileSeqStream {
    pub fn new(filename: PathBuf, revcomps_enabled: bool) -> Self {
        Self{path: filename, stream: None, seq_buf: vec![], revcomps_enabled, revcomp_next: false}
    }
}

impl SeqStream for ChainedInputStream {

    fn stream_next(&mut self) -> Option<&[u8]> {
        self.seq_buf.clear();
        if let Some(f) = self.cur_file.as_mut() {
            if let Some(rec) = f.read_next().unwrap() {
                self.seq_buf.extend_from_slice(rec.seq);
                Some(&self.seq_buf)
            } else {
                // File is finished -> open the next file
                self.cur_file_idx += 1;
                self.cur_file = if self.cur_file_idx == self.paths.len() {
                    None // All files procesed
                } else {
                    let new_file = DynamicFastXReader::from_file(&self.paths[self.cur_file_idx])
                        .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", self.paths[self.cur_file_idx].display()));
                    Some(new_file)
                };

                self.stream_next()
            }
        } else {
            None    
        }
    }
}

impl SeqStream for LazyFileSeqStream {

    fn stream_next(&mut self) -> Option<&[u8]> {
        if self.revcomps_enabled {
            if self.revcomp_next {
                reverse_complement_in_place(&mut self.seq_buf);
                self.revcomp_next = !self.revcomp_next;
                return Some(&self.seq_buf);
            } else {
                self.revcomp_next = !self.revcomp_next;
            }
        }


        if self.stream.is_none() {
            self.stream = Some(DynamicFastXReader::from_file(&self.path)
                .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", self.path.display())));
        }

        let s = self.stream.as_mut().unwrap();
        s.read_next().unwrap().map(|r| {
            self.seq_buf.clear();
            self.seq_buf.extend_from_slice(r.seq);
            self.seq_buf.as_slice()
        })
    }
}

/// Lazily reads one sequence from a shared reader on first access, then yields
/// that sequence (and optionally its reverse complement) as a color stream.
pub struct SingleSeqStream {
    reader: Arc<Mutex<DynamicFastXReader>>,
    seq: Vec<u8>,
    revcomps_enabled: bool,
    phase: u8, // 0=unread, 1=revcomp pending, 2=done
}

impl SingleSeqStream {
    pub fn new(reader: Arc<Mutex<DynamicFastXReader>>, revcomps_enabled: bool) -> Self {
        Self { reader, seq: Vec::new(), revcomps_enabled, phase: 0 }
    }
}

impl SeqStream for SingleSeqStream {
    fn stream_next(&mut self) -> Option<&[u8]> {
        match self.phase {
            0 => {
                {
                    let mut guard = self.reader.lock().unwrap();
                    let rec = guard.read_next().unwrap()
                        .expect("SingleSeqStream: shared reader exhausted prematurely");
                    self.seq.clear();
                    self.seq.extend_from_slice(rec.seq);
                } // guard and rec dropped here
                self.phase = if self.revcomps_enabled { 1 } else { 2 };
                Some(&self.seq)
            }
            1 => {
                reverse_complement_in_place(&mut self.seq);
                self.phase = 2;
                Some(&self.seq)
            }
            _ => {
                self.seq = Vec::new(); // free the allocation
                None
            }
        }
    }
}