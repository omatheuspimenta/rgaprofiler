#!/usr/bin/env bash
# Tests for bin/detect_host_resources.sh against mocked hosts (no real GPU needed).
# Run: bash tests/bin/test_host_resources.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$HERE/../../bin/detect_host_resources.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fails=0

check() { # name expected actual
    if [[ "$2" == "$3" ]]; then echo "ok   - $1"; else echo "FAIL - $1: expected '$2', got '$3'"; fails=$((fails+1)); fi
}
get() { printf '%s\n' "$1" | sed -n "s/^$2=//p"; }

# PATH containing only the tools the script needs (no nvidia-smi), plus optional fake ones.
mkdir -p "$TMP/basebin"
for t in bash awk sed grep tr head sort wc nproc cat dirname; do ln -sf "$(command -v $t)" "$TMP/basebin/$t"; done
fake_smi() { # dir, vram lines...
    mkdir -p "$1"; local d="$1"; shift
    printf '#!/bin/sh\nprintf "%%s\\n" %s\n' "$(printf "'%s' " "$@")" > "$d/nvidia-smi"; chmod +x "$d/nvidia-smi"
}
printf 'MemTotal:       527433728 kB\n' > "$TMP/meminfo_big"     # ~503 GiB
printf 'MemTotal:        16384000 kB\n' > "$TMP/meminfo_small"   # ~15 GiB

# 1. large host, one 20GB GPU
fake_smi "$TMP/smi1" 20470
out=$(PATH="$TMP/smi1:$TMP/basebin" RGAPROFILER_NPROC=256 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big" bash "$SCRIPT")
check "large host: cpus"     256   "$(get "$out" cpus)"
check "large host: mem_gb"   503   "$(get "$out" mem_gb)"
check "large host: gpus"     1     "$(get "$out" gpus)"
check "large host: vram"     20470 "$(get "$out" gpu_vram_mb)"

# 2. small host, no GPU at all
out=$(PATH="$TMP/basebin" RGAPROFILER_NPROC=4 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_small" bash "$SCRIPT")
check "small host: cpus"     4 "$(get "$out" cpus)"
check "small host: mem_gb"   15 "$(get "$out" mem_gb)"
check "small host: gpus"     0 "$(get "$out" gpus)"
check "small host: vram"     0 "$(get "$out" gpu_vram_mb)"
check "no gpu: slots >= 1"   1 "$(get "$out" gpu_slots_per_gpu)"

# 3. cgroup limit below MemTotal wins; "max" (unlimited) is ignored
echo $((8*1024*1024*1024)) > "$TMP/cg_8g"; echo max > "$TMP/cg_max"
out=$(PATH="$TMP/basebin" RGAPROFILER_NPROC=4 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big" RGAPROFILER_CGROUP_MEM_FILE="$TMP/cg_8g" bash "$SCRIPT")
check "cgroup limit caps mem_gb" 8 "$(get "$out" mem_gb)"
out=$(PATH="$TMP/basebin" RGAPROFILER_NPROC=4 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big" RGAPROFILER_CGROUP_MEM_FILE="$TMP/cg_max" bash "$SCRIPT")
check "cgroup 'max' ignored"     503 "$(get "$out" mem_gb)"

# 4. two GPUs of different size: gpus=2, vram = the smaller one
fake_smi "$TMP/smi2" 24576 12288
out=$(PATH="$TMP/smi2:$TMP/basebin" RGAPROFILER_NPROC=8 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big" bash "$SCRIPT")
check "2 gpus: count"        2     "$(get "$out" gpus)"
check "2 gpus: min vram"     12288 "$(get "$out" gpu_vram_mb)"

# 5. CUDA_VISIBLE_DEVICES restricts the GPU count
out=$(CUDA_VISIBLE_DEVICES=1 PATH="$TMP/smi2:$TMP/basebin" RGAPROFILER_NPROC=8 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big" bash "$SCRIPT")
check "CUDA_VISIBLE_DEVICES=1: gpus" 1 "$(get "$out" gpus)"

# 6. a wedged driver printing text instead of numbers -> treated as no GPU, not a crash
fake_smi "$TMP/smi3" "NVIDIA-SMI has failed because it couldn't communicate with the NVIDIA driver."
out=$(PATH="$TMP/smi3:$TMP/basebin" RGAPROFILER_NPROC=4 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_small" bash "$SCRIPT")
check "garbage nvidia-smi output: gpus" 0 "$(get "$out" gpus)"

# 7. --gpu-slots: explicit N per GPU x GPUs, never below 1, 'auto' is safe
env_base=(PATH="$TMP/smi2:$TMP/basebin" RGAPROFILER_NPROC=8 RGAPROFILER_MEMINFO_FILE="$TMP/meminfo_big")
check "--gpu-slots 3 on 2 GPUs" 6 "$(env "${env_base[@]}" bash "$SCRIPT" --gpu-slots 3)"
check "--gpu-slots auto (unknown footprint => 1/GPU)" 2 "$(env "${env_base[@]}" RGAPROFILER_GPU_TOOL_PEAK_VRAM_MB=unknown bash "$SCRIPT" --gpu-slots auto)"
check "--gpu-slots on GPU-less host >= 1" 1 "$(PATH="$TMP/basebin" bash "$SCRIPT" --gpu-slots auto)"
# with a known footprint, slots follow VRAM: 12288MiB * 85% / 4000MiB = 2 per GPU
check "--gpu-slots auto from VRAM" 4 "$(env "${env_base[@]}" RGAPROFILER_GPU_TOOL_PEAK_VRAM_MB=4000 bash "$SCRIPT" --gpu-slots auto)"
check "--gpu-slots auto: footprint > VRAM still 1" 2 "$(env "${env_base[@]}" RGAPROFILER_GPU_TOOL_PEAK_VRAM_MB=99999 bash "$SCRIPT" --gpu-slots auto)"

# 8. the shipped VRAM constant: strict serialization on the 20GB dev card, more only on big cards
for spec in "20470 1" "40960 1" "49152 2" "81920 3"; do
    set -- $spec; fake_smi "$TMP/smi_v$1" "$1"
    check "shipped constant: ${1}MiB card => $2 slot(s)/GPU" "$2" "$(PATH="$TMP/smi_v$1:$TMP/basebin" bash "$SCRIPT" | sed -n 's/^gpu_slots_per_gpu=//p')"
done

if (( fails )); then echo "$fails failure(s)"; exit 1; fi
echo "all host-resource tests passed"
