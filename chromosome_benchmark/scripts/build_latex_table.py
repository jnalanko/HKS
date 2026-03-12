#!/usr/bin/env python3
import sys
import re

def parse_stderr(lines):
    rss_re = re.compile(r"Maximum resident set size.*:\s*(\d+)")

    elapsed_seconds = None
    max_rss_bytes = None

    for line in lines:
        # Match elapsed time — supports h:mm:ss, m:ss, or s.s formats
        if "Elapsed (wall clock) time" in line:
            # Two formats based on whether the time is more than a hour:
            # Elapsed (wall clock) time (h:mm:ss or m:ss): 1:13:01
            # Elapsed (wall clock) time (h:mm:ss or m:ss): 51:00.54
            s = line.split(" ")[-1].strip()
            if s.count(":") == 2: # h:mm:ss
                hours = int(s.split(":")[0])
                mins = int(s.split(":")[1])
                secs = int(s.split(":")[2])
                elapsed_seconds = hours * 60*60 + mins * 60 + secs
            else: # m:ss
                mins = int(s.split(":")[0])
                secs = float(s.split(":")[1])
                elapsed_seconds = mins * 60 + secs

        # Match max RSS (in kilobytes)
        elif "Maximum resident set size" in line:
            match = rss_re.search(line)
            if match:
                max_rss_bytes = int(match.group(1)) * 1024  # convert KB to bytes


    return {"elapsed_seconds": elapsed_seconds, "max_rss_bytes": max_rss_bytes}

k_list = [15,31,47,63]
m_list = [15,22,31]
t = 32
query_total_bases = 3117292070 
dks_index_size = 11150448137

log_base_dir = "./final_logs"

# Parse index sizes from index_sizes.txt
# Format: -rw-r--r-- 1 user group SIZE month day time ./kraken_k{k}_m{m}/hash.k2d
index_sizes = {}  # (k, m) -> bytes
path_re = re.compile(r'kraken_k(\d+)_m(\d+)/hash\.k2d')
try:
    for line in open(f"{log_base_dir}/index_sizes.txt"):
        parts = line.split()
        if len(parts) < 9:
            continue
        path_match = path_re.search(parts[8])
        if path_match:
            k, m = int(path_match.group(1)), int(path_match.group(2))
            index_sizes[(k, m)] = int(parts[4])
except Exception as e:
    sys.stderr.write(str(e) + "\n")

data = {}
for k in k_list:
    for m in m_list:
        if m > k: continue

        query_filename = f"{log_base_dir}/logs/kraken-query-k{k}-m{m}-t{t}.log"
        try:
            res = parse_stderr(open(query_filename).readlines())
            data[(k, m)] = res
        except Exception as e:
            sys.stderr.write(str(e) + "\n")

build_data = {}
for k in k_list:
    for m in m_list:
        if m > k: continue
        build_filename = f"{log_base_dir}/logs/kraken-build-k{k}-m{m}-t{t}.log"
        try:
            res = parse_stderr(open(build_filename).readlines())
            build_data[(k, m)] = res
        except Exception as e:
            sys.stderr.write(str(e) + "\n")

dks_build_time = None
dks_build_mem = None
dks_build_phase2_time = None
dks_build_phase2_mem = None
try:
    r1 = parse_stderr(open(f"{log_base_dir}/logs/sbwt-in-mem-k63-t{t}.log").readlines())
    r2 = parse_stderr(open(f"{log_base_dir}/logs/dks-build-k63-t{t}.log").readlines())
    if r1['elapsed_seconds'] is not None and r2['elapsed_seconds'] is not None:
        dks_build_time = r1['elapsed_seconds'] + r2['elapsed_seconds']
    if r1['max_rss_bytes'] is not None and r2['max_rss_bytes'] is not None:
        dks_build_mem = max(r1['max_rss_bytes'], r2['max_rss_bytes'])
    dks_build_phase2_time = r2['elapsed_seconds']
    dks_build_phase2_mem = r2['max_rss_bytes']
except Exception as e:
    sys.stderr.write(str(e) + "\n")

ns_per_bp_re = re.compile(r"Query time per base pair:\s*([\d.]+)\s*nanoseconds")

dks_data = {}
for k in k_list:
    query_filename = f"{log_base_dir}/logs/dks-query-k{k}-t{t}.log"
    try:
        ns_per_bp = None
        for line in open(query_filename):
            match = ns_per_bp_re.search(line)
            if match:
                ns_per_bp = float(match.group(1))
                break
        dks_data[k] = ns_per_bp
    except Exception as e:
        sys.stderr.write(str(e) + "\n")

valid_ms = [m for m in m_list]

def fmt_time(k, m):
    if m > k or (k, m) not in data:
        return "NA"
    r = data[(k, m)]
    if r['elapsed_seconds'] is None:
        return "?"
    throughput_mbps = query_total_bases / r['elapsed_seconds'] / 2**20
    return f"{throughput_mbps:.1f}"

def fmt_size(k, m):
    if m > k:
        return "NA"
    return f"{index_sizes[(k, m)] / 2**30:.1f}" if (k, m) in index_sizes else "?"

def fmt_build_time(k, m):
    if m > k or (k, m) not in build_data:
        return "NA"
    r = build_data[(k, m)]
    return f"{r['elapsed_seconds'] / 60:.1f}" if r['elapsed_seconds'] is not None else "?"

def fmt_build_mem(k, m):
    if m > k or (k, m) not in build_data:
        return "NA"
    r = build_data[(k, m)]
    return f"{r['max_rss_bytes'] / 2**30:.1f}" if r['max_rss_bytes'] is not None else "?"

# Build LaTeX table: columns = k values (each with two subcolumns), rows = m values
print(r"\begin{table}[h]")
print(r"\centering")
print(r"\small")
print(r"\setlength{\tabcolsep}{6pt}")
print(r"\renewcommand{\arraystretch}{1.2}")
col_spec = "|l|" + "cccc|" * len(k_list)
print(r"\begin{tabular}{" + col_spec + r"}")
print(r"\hline")

# Top header: k values spanning four subcolumns each
top = ""
for i, k in enumerate(k_list):
    fmt = r"c|"
    top += f" & " + r"\multicolumn{4}{" + fmt + r"}{$k=" + str(k) + r"$}"
print(top + r" \\")

# Sub-header: build time / build mem / index size / query throughput
sub = ""
for k in k_list:
    sub += r" & $B_t$ & $B_m$ & $I$ & $Q_t$"
print(sub + r" \\")
print(r"\hline")

# Data rows
for m in valid_ms:
    row = f"Kraken-m{m}"
    for k in k_list:
        row += f" & {fmt_build_time(k, m)} & {fmt_build_mem(k, m)} & {fmt_size(k, m)} & {fmt_time(k, m)}"
    row += r" \\"
    print(row)

print(r"\hline")

# DKS row
dks_row = r"\textsc{dks}"
for k in k_list:
    ns_per_bp = dks_data.get(k)
    if ns_per_bp is not None:
        throughput_mbps = 1e9 / ns_per_bp / 2**20  # ns/bp -> Mibp/s
        dks_row += f" & {dks_build_time / 60:.1f} & {dks_build_mem / 2**30:.1f} & {dks_index_size / 2**30:.1f} & {throughput_mbps:.1f}" if dks_build_time is not None else f" & ? & ? & {dks_index_size / 2**30:.1f} & {throughput_mbps:.1f}"
    else:
        dks_row += f" & ? & ? & {dks_index_size / 2**30:.1f} & ?"
dks_row += r" \\"
print(dks_row)

print(r"\hline")
print(r"\end{tabular}")
print(r"\caption{Comparison of Kraken (with minimizer length $m \leq k$, written as Kraken-m) and \textsc{dks} at $t=32$ threads and varying k-mer length $k$. All \textsc{dks} results use the same index built for maximum $k = 63$. Columns: build time ($B_t$, minutes), peak memory during build ($B_m$, GiB), index size ($I$, GiB), and query throughput ($Q_t$, MiB/s). \textsc{dks} construction proceeds in two phases: first, an SBWT index is built in memory, and second, the \textsc{dks} index is constructed from the SBWT. The $B_t$ and $B_m$ values reported for \textsc{dks} are for the whole two-stage pipeline. The second phase took 189 seconds and used 15.6 GiB of RAM. A slower disk-based option exists for SBWT construction, taking 37 minutes, using 15.9 GiB of RAM and 92.9 GiB of disk. \textcolor{red}{Need to subtract index loading time from Kraken (though it should be close to negligible since the indexes are so small). The result on DKS k = 15 is slow probably because it becomes IO-bottlenecked by the Lurstre file system. Could be avoided by writing to /dev/null (it's not a real limitation, a fast SSD would also fix it)}}")
print(r"\label{tab:performance}")
print(r"\end{table}")

