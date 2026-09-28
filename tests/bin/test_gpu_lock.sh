#!/usr/bin/env bash
# Tests for bin/gpu_lock.sh: mutual exclusion, slot count, device pinning, fail-open.
# Needs `flock` (util-linux). No GPU required. Run: bash tests/bin/test_gpu_lock.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fails=0
check() { if [[ "$2" == "$3" ]]; then echo "ok   - $1"; else echo "FAIL - $1: expected '$2', got '$3'"; fails=$((fails+1)); fi; }

# Mimics the real call site: a `set -e` wrapper that SOURCES the lock, then does "GPU work".
cat > "$TMP/holder.sh" <<'H'
#!/bin/bash
set -e
cd "$(mktemp -d)"
source "$REPO/bin/gpu_lock.sh" "$SLOTS" ${NEED:-}
echo "$(date +%s.%N) start ${CUDA_VISIBLE_DEVICES:-none} ${NEED:-}" >> "$EVENTS"
sleep 1
echo "$(date +%s.%N) end ${NEED:-}" >> "$EVENTS"
H
peak() { sort -n "$1" | awk '/start/{c++; if(c>m)m=c} /end/{c--} END{print m+0}'; }
run_holders() { # n slots [extra env...]
    local n=$1 slots=$2; shift 2; : > "$TMP/ev"
    for _ in $(seq "$n"); do env REPO="$REPO" SLOTS="$slots" EVENTS="$TMP/ev" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks" "$@" bash "$TMP/holder.sh" & done
    wait
}

# fake single-GPU host so detection is deterministic (no real nvidia-smi involved)
mkdir -p "$TMP/bin"; for t in bash awk sed grep tr head sort wc nproc cat dirname date sleep mktemp id mkdir flock seq env; do ln -sf "$(command -v $t)" "$TMP/bin/$t"; done
printf '#!/bin/sh\necho 20470\n' > "$TMP/bin/nvidia-smi"; chmod +x "$TMP/bin/nvidia-smi"

run_holders 5 1 PATH="$TMP/bin"
check "1 slot: never more than 1 concurrent holder" 1 "$(peak "$TMP/ev")"
check "1 slot: all 5 tasks eventually ran"          5 "$(grep -c start "$TMP/ev")"

run_holders 6 2 PATH="$TMP/bin"
check "2 slots: exactly 2 concurrent holders"       2 "$(peak "$TMP/ev")"

# a killed holder frees its slot immediately (kernel releases the flock)
: > "$TMP/ev"
( cd "$TMP" && env REPO="$REPO" SLOTS=1 EVENTS="$TMP/ev2" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_kill" PATH="$TMP/bin" bash -c 'source "$REPO/bin/gpu_lock.sh" 1; exec sleep 30' ) & victim=$!
sleep 1
timeout 2 env REPO="$REPO" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_kill" PATH="$TMP/bin:$PATH" bash -c 'cd "$(mktemp -d)"; source "$REPO/bin/gpu_lock.sh" 1' ; rc=$?
check "lock is really held while the holder lives (2nd acquirer blocks)" 124 "$rc"
kill -9 "$victim" 2>/dev/null; wait "$victim" 2>/dev/null; sleep 0.5
start=$(date +%s)
( cd "$TMP" && env REPO="$REPO" SLOTS=1 RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_kill" PATH="$TMP/bin" bash -c 'source "$REPO/bin/gpu_lock.sh" 1; true' )
check "lock freed after holder is killed (next task not blocked)" 1 "$(( $(date +%s) - start < 5 ))"

# device pinning: with 2 GPUs restricted via CUDA_VISIBLE_DEVICES, each task gets one of them
run_holders 4 1 PATH="$TMP/bin" CUDA_VISIBLE_DEVICES=3,5
check "2 devices x 1 slot: 2 concurrent" 2 "$(peak "$TMP/ev")"
check "tasks pinned to the allowed devices only" "" "$(grep start "$TMP/ev" | awk '{print $3}' | grep -v -E '^(3|5)$' || true)"

# 0 disables locking entirely (unbounded concurrency)
run_holders 4 0 PATH="$TMP/bin"
check "slots=0: no locking" 4 "$(peak "$TMP/ev")"

# fail-open: unwritable lock dir must not fail (or hang) the task
: > "$TMP/ev"
( cd "$TMP" && env REPO="$REPO" SLOTS=1 EVENTS="$TMP/ev" RGAPROFILER_GPU_LOCK_DIR=/proc/nonexistent/locks PATH="$TMP/bin" bash "$TMP/holder.sh" ) 2>/dev/null
check "unwritable lock dir: task still runs (fail-open)" 1 "$(grep -c start "$TMP/ev")"

# --- VRAM-weighted mode ('auto' + a per-task VRAM need). The fake card has 20470 MiB ->
# 20470*85%/512 = 33 units; 4600 MiB = 9 units, 2200 MiB = 5, 'all' = 33.
# peak_units: max over time of the summed units of concurrently running holders.
peak_units() { sort -n "$1" | awk '
    function u(n) { return n == "all" ? 33 : int((n + 511) / 512) }
    /start/ { c += u($4); if (c > m) m = c } /end/ { c -= u($3) } END { print m + 0 }'; }
run_mixed() { # run_mixed need... : one holder per need, all started at once
    : > "$TMP/ev"
    for need in "$@"; do env REPO="$REPO" SLOTS=auto NEED="$need" EVENTS="$TMP/ev" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_w" PATH="$TMP/bin" bash "$TMP/holder.sh" & done
    wait
}

run_mixed 4600 4600 4600 4600 4600 4600
check "weighted: 4600MiB tasks -> 3 fit on the card at once" 3 "$(peak "$TMP/ev")"
check "weighted: all 6 ran" 6 "$(grep -c start "$TMP/ev")"

run_mixed all all all
check "weighted: 'all' is exclusive" 1 "$(peak "$TMP/ev")"

run_mixed all 4600 2200 all 4600 2200 2200 4600 all 2200
check "weighted mix: summed units never exceed the card's 33" 1 "$(( $(peak_units "$TMP/ev") <= 33 ))"
check "weighted mix: light tools really share the card" 1 "$(( $(peak "$TMP/ev") > 1 ))"
check "weighted mix: no deadlock, all 10 ran" 10 "$(grep -c start "$TMP/ev")"

# an explicit per-GPU count still means N tasks of any size (need is ignored)
: > "$TMP/ev"
for _ in 1 2 3 4; do env REPO="$REPO" SLOTS=2 NEED=all EVENTS="$TMP/ev" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_n" PATH="$TMP/bin" bash "$TMP/holder.sh" & done; wait
check "explicit --gpu_concurrency 2 ignores VRAM need" 2 "$(peak "$TMP/ev")"

# unreadable VRAM (nvidia-smi prints garbage) -> weighted request falls back to 1 per GPU
mkdir -p "$TMP/bin_novram"; cp -a "$TMP/bin/." "$TMP/bin_novram/"
printf '#!/bin/sh\necho "ERR!"\n' > "$TMP/bin_novram/nvidia-smi"; chmod +x "$TMP/bin_novram/nvidia-smi"
: > "$TMP/ev"
for _ in 1 2 3; do env REPO="$REPO" SLOTS=auto NEED=2200 EVENTS="$TMP/ev" RGAPROFILER_GPU_LOCK_DIR="$TMP/locks_nv" PATH="$TMP/bin_novram" bash "$TMP/holder.sh" & done; wait
check "unknown VRAM: weighted request is serialized" 1 "$(peak "$TMP/ev")"

if (( fails )); then echo "$fails failure(s)"; exit 1; fi
echo "all gpu-lock tests passed"
