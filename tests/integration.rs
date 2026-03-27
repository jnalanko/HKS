use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

/*
Currently, these test just that all the commands run succesfully.
Outputs are not checked.
*/

// These are set by Cargo when compiling
const BIN: &str = env!("CARGO_BIN_EXE_hks");
const PROJECT_DIR: &str = env!("CARGO_MANIFEST_DIR");

static COUNTER: AtomicU64 = AtomicU64::new(0);

fn tmp_dir() -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::SeqCst);
    let dir = PathBuf::from(PROJECT_DIR)
        .join("target")
        .join("test-tmp")
        .join(format!("{}", id));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn hks() -> Command {
    let mut cmd = Command::new(BIN);
    cmd.current_dir(PROJECT_DIR);
    cmd
}

fn build_basic_index(out: &Path) {
    let status = hks()
        .args(["build", "-s", "10", "--label-by-file", "example/file_of_files.txt", "-o"])
        .arg(out)
        .status()
        .unwrap();
    assert!(status.success(), "basic index build failed");
}

// --- build ---

#[test]
fn build_label_by_file() {
    let dir = tmp_dir();
    let status = hks()
        .args(["build", "-s", "10", "--label-by-file", "example/file_of_files.txt", "-o"])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_label_by_seq() {
    let dir = tmp_dir();
    // A single FASTA file where each sequence becomes one label
    let combined = dir.join("combined.fna");
    std::fs::write(
        &combined,
        ">seqA\nACGTACGTGTCGTA\n>seqB\nACGTGCTGAGCA\n",
    )
    .unwrap();
    let status = hks()
        .args(["build", "-s", "10", "--label-by-seq"])
        .arg(&combined)
        .args(["-o"])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_with_hierarchy() {
    let dir = tmp_dir();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "--hierarchy",
            "example/hierarchy.txt",
            "-o",
        ])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_with_custom_labels() {
    let dir = tmp_dir();
    let labels_file = dir.join("labels.txt");
    std::fs::write(&labels_file, "labelA\nlabelB\nlabelC\nlabelD\n").unwrap();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "--labels",
        ])
        .arg(&labels_file)
        .args(["-o"])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_with_unitigs() {
    let dir = tmp_dir();
    // The unitigs file must contain all k-mers from all input files
    let unitigs = dir.join("unitigs.fna");
    let content = std::fs::read_to_string(
        PathBuf::from(PROJECT_DIR).join("example/A.fna"),
    )
    .unwrap()
        + &std::fs::read_to_string(PathBuf::from(PROJECT_DIR).join("example/B.fna")).unwrap()
        + &std::fs::read_to_string(PathBuf::from(PROJECT_DIR).join("example/C.fna")).unwrap()
        + &std::fs::read_to_string(PathBuf::from(PROJECT_DIR).join("example/D.fna")).unwrap();
    std::fs::write(&unitigs, content).unwrap();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "--unitigs",
        ])
        .arg(&unitigs)
        .args(["-o"])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_forward_only() {
    let dir = tmp_dir();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "--forward-only",
            "-o",
        ])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_n_threads() {
    let dir = tmp_dir();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "-t",
            "2",
            "-o",
        ])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn build_external_memory() {
    let dir = tmp_dir();
    let tmp_work = dir.join("tmp");
    std::fs::create_dir_all(&tmp_work).unwrap();
    let status = hks()
        .args([
            "build",
            "-s",
            "10",
            "--label-by-file",
            "example/file_of_files.txt",
            "--external-memory",
        ])
        .arg(&tmp_work)
        .args(["-o"])
        .arg(dir.join("index.hks"))
        .status()
        .unwrap();
    assert!(status.success());
}

// --- lookup ---

#[test]
fn lookup_basic() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_with_k() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["-k", "5"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_report_label_ids() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["--report-label-ids"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_report_query_names() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["--report-query-names"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_report_misses() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["--report-misses"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_no_header() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["--no-header"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_n_threads() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["-t", "2"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn lookup_batch_size() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["lookup", "-q", "example/query.fasta", "-i"])
        .arg(&index)
        .args(["--batch-size", "100"])
        .status()
        .unwrap();
    assert!(status.success());
}

// --- stats ---

#[test]
fn stats_basic() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["stats", "-i"])
        .arg(&index)
        .status()
        .unwrap();
    assert!(status.success());
}

// --- node-stats ---

#[test]
fn node_stats_basic() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["node-stats", "--index"])
        .arg(&index)
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn node_stats_report_label_names() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["node-stats", "--index"])
        .arg(&index)
        .args(["--report-label-ids"])
        .status()
        .unwrap();
    assert!(status.success());
}

#[test]
fn node_stats_n_threads() {
    let dir = tmp_dir();
    let index = dir.join("index.hks");
    build_basic_index(&index);
    let status = hks()
        .args(["node-stats", "--index"])
        .arg(&index)
        .args(["-t", "2"])
        .status()
        .unwrap();
    assert!(status.success());
}

