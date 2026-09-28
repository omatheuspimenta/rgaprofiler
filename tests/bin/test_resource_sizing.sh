#!/usr/bin/env bash
# Runs the REAL conf/base.config (through tests/probe/) under mocked hosts and asserts on the
# cpus/memory the per-chunk labels resolve to. Needs `nextflow` on PATH; no containers, no GPU.
# Run: bash tests/bin/test_resource_sizing.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROBE="$HERE/../probe"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fails=0
check() { if [[ "$2" == "$3" ]]; then echo "ok   - $1"; else echo "FAIL - $1: expected '$2', got '$3'"; fails=$((fails+1)); fi; }
le()    { if (( $2 <= $3 )); then echo "ok   - $1 ($2 <= $3)"; else echo "FAIL - $1: $2 > $3"; fails=$((fails+1)); fi; }

# Mocked hosts: RGAPROFILER_* env vars are read by bin/detect_host_resources.sh.
printf 'MemTotal:       527433728 kB\n' > "$TMP/meminfo_big"     # 503 GiB
printf 'MemTotal:         8388608 kB\n' > "$TMP/meminfo_tiny"    # 8 GiB
head -c $((300 * 1024)) /dev/zero | tr '\0' 'A' > "$TMP/chunk_small.fa"        # ~300 KB  (~ a 300-sequence chunk)
head -c $((200 * 1024 * 1024)) /dev/zero | tr '\0' 'A' > "$TMP/chunk_huge.fa"  # 200 MB   (an unchunked whole proteome)

probe() { # nproc meminfo chunk [extra nextflow args...] -> prints probe lines
    local nproc=$1 meminfo=$2 chunk=$3; shift 3
    ( cd "$TMP" && RGAPROFILER_NPROC="$nproc" RGAPROFILER_MEMINFO_FILE="$meminfo" \
        nextflow run "$PROBE/main.nf" -c "$PROBE/nextflow.config" --chunk "$chunk" -work-dir "$TMP/work" "$@" 2>&1 | grep -E '^(CHUNK_|TOOL )' )
}
val() { printf '%s\n' "$1" | grep "^$2 " | sed -n "s/.*$3=\([0-9]*\).*/\1/p"; }
tval() { val "$1" "TOOL $2" "$3"; }
# Real-looking chunks for the per-tool models: 300 proteins of 450 residues, and a single
# 35,000-residue protein (titin-sized; DeepTMHMM's VRAM need grows with the longest one).
awk 'BEGIN { s = sprintf("%450s", ""); gsub(/ /, "M", s); for (i = 1; i <= 300; i++) printf(">p%d\n%s\n", i, s) }' > "$TMP/chunk_300x450.fa"
awk 'BEGIN { s = sprintf("%35000s", ""); gsub(/ /, "M", s); printf(">titin\n%s\n", s) }' > "$TMP/chunk_35k.fa"

# 1. Large host, small (typical --num_blocks 1000) chunk: modest, chunk-sized requests --
#    NOT the old 150GB-per-task whole-proteome value that starved every other process.
out=$(probe 256 "$TMP/meminfo_big" "$TMP/chunk_small.fa")
med_mem=$(val "$out" CHUNK_MEDIUM mem_gb); high_mem=$(val "$out" CHUNK_HIGH mem_gb)
le "256cpu/503GB, small chunk: GPU-tool memory is chunk-sized (<= 32GB, was 150GB)" "$med_mem" 32
le "256cpu/503GB, small chunk: InterProScan memory is chunk-sized (<= 32GB, was 100GB)" "$high_mem" 32
check "256cpu/503GB: process_medium_chunk cpus" 6 "$(val "$out" CHUNK_MEDIUM cpus)"   # use_gpu=false => CPU fallback sizing
check "256cpu/503GB: process_high_chunk cpus"   4 "$(val "$out" CHUNK_HIGH cpus)"
# => many chunks fit at once: 503GB / that per-task memory
echo "info - concurrent InterProScan chunks that fit in RAM alone: $(( 503 / high_mem ))"

# 2. Memory scales with chunk size (unchunked/whole-proteome path stays safe on a big host)
out=$(probe 256 "$TMP/meminfo_big" "$TMP/chunk_huge.fa")
big_mem=$(val "$out" CHUNK_MEDIUM mem_gb)
if (( big_mem > med_mem )); then echo "ok   - bigger chunk asks for more memory ($med_mem -> $big_mem GB)"; else echo "FAIL - memory did not scale with chunk size ($med_mem -> $big_mem)"; fails=$((fails+1)); fi
le "huge chunk never exceeds 90% of host RAM" "$big_mem" $(( 503 * 9 / 10 ))

# 3. Tiny host: a task must NEVER ask for more than the machine can provide (silent deadlock)
out=$(probe 2 "$TMP/meminfo_tiny" "$TMP/chunk_huge.fa")
le "2cpu/8GB: GPU-tool cpus <= host"    "$(val "$out" CHUNK_MEDIUM cpus)"    2
le "2cpu/8GB: GPU-tool memory <= host"  "$(val "$out" CHUNK_MEDIUM mem_gb)"  8
le "2cpu/8GB: InterProScan cpus <= host"   "$(val "$out" CHUNK_HIGH cpus)"   2
le "2cpu/8GB: InterProScan memory <= host" "$(val "$out" CHUNK_HIGH mem_gb)" 8

# 4. Explicit user limits still win over auto-scaling (process.resourceLimits)
cat > "$TMP/limits.config" <<'L'
process.resourceLimits = [cpus: 3, memory: 10.GB, time: 1.h]
L
out=$(probe 256 "$TMP/meminfo_big" "$TMP/chunk_huge.fa" -c "$TMP/limits.config")
le "resourceLimits caps cpus"   "$(val "$out" CHUNK_HIGH cpus)"   3
le "resourceLimits caps memory" "$(val "$out" CHUNK_HIGH mem_gb)" 10

# 5. Unusable/failed detection degrades to fixed defaults instead of crashing
out=$(probe 4 /nonexistent/meminfo "$TMP/chunk_small.fa")
check "meminfo unreadable: pipeline still runs and sizes both tasks" 2 "$(printf '%s\n' "$out" | grep -c '^CHUNK_')"

# 6. Per-tool models (CPU mode: the probe runs with use_gpu=false) on a 24-CPU/503GB host
out=$(probe 24 "$TMP/meminfo_big" "$TMP/chunk_300x450.fa")
check "SignalP6 reserves the 8 threads it really runs" 8 "$(tval "$out" SIGNALP6 cpus)"
le    "SignalP6 memory is its own (~3.5GB measured), not the old shared 14GB+" "$(tval "$out" SIGNALP6 mem_gb)" 6
check "DeepTMHMM on CPU (local) reserves the whole host (its threads = host cores)" 24 "$(tval "$out" DEEPTMHMM cpus)"
check "DeepLoc2 keeps 6 CPUs on CPU" 6 "$(tval "$out" DEEPLOC2 cpus)"
check "DeepCoil2 keeps 6 CPUs on CPU (its -n_cpu)" 6 "$(tval "$out" DEEPCOIL2 cpus)"
check "typical chunk keeps the 8h time floor" 8 "$(tval "$out" SIGNALP6 time_h)"
out=$(probe 24 "$TMP/meminfo_big" "$TMP/chunk_300x450.fa" --deeptmhmm_cpu_threads 6)
check "--deeptmhmm_cpu_threads 6: DeepTMHMM reserves just 6" 6 "$(tval "$out" DEEPTMHMM cpus)"

# 7. DeepTMHMM's GPU-lock need follows the chunk's longest protein (11.7GB measured at 35k)
out_typ=$(probe 24 "$TMP/meminfo_big" "$TMP/chunk_300x450.fa"); out_long=$(probe 24 "$TMP/meminfo_big" "$TMP/chunk_35k.fa")
check "stats: longest protein read from the chunk" 35000 "$(tval "$out_long" DEEPTMHMM longest)"
le    "DeepTMHMM VRAM need, typical chunk" "$(tval "$out_typ" DEEPTMHMM vram_mb)" 5000
le    "DeepTMHMM VRAM need, 35k protein, covers the measured 11662MiB" 11662 "$(tval "$out_long" DEEPTMHMM vram_mb)"

# 8. Time scales with the chunk: a whole proteome in one chunk gets far more than 8h on CPU
out=$(probe 24 "$TMP/meminfo_big" "$TMP/chunk_huge.fa")
le "SignalP6 on CPU, 200MB chunk: time limit grows past the 8h floor" 9 "$(tval "$out" SIGNALP6 time_h)"

if (( fails )); then echo "$fails failure(s)"; exit 1; fi
echo "all resource-sizing tests passed"
