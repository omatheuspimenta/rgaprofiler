#!/usr/bin/env bash
set -euo pipefail

# detect_host_resources.sh — report the CPUs, RAM and GPU(s) of the host this runs on.
#
# Used lazily, per task, by conf/base.config (through the `ext.host` closure) to size
# per-task resource requests, and once at config-parse time to size the GPU throttle
# (`maxForks`). Safe to run standalone:
#   ./bin/detect_host_resources.sh
#   ./bin/detect_host_resources.sh --gpu-slots auto
#
# Default output, one `key=value` per line:
#   cpus=<int>               CPUs usable by this process (nproc: respects cpusets/affinity)
#   mem_gb=<int>             total RAM in whole GB, capped at the cgroup memory limit if any
#   gpus=<int>               visible NVIDIA GPUs (respects CUDA_VISIBLE_DEVICES); 0 if none
#   gpu_vram_mb=<int>        total VRAM of the smallest visible GPU, in MiB; 0 if none
#   gpu_slots_per_gpu=<int>  how many GPU tasks may safely share ONE GPU at once (>= 1),
#                            derived from gpu_vram_mb -- see GPU_TOOL_PEAK_VRAM_MB below
#
# `--gpu-slots [N|auto]`: print only the total number of concurrent GPU-task slots
# (per-GPU slots x visible GPUs, never below 1). N is an explicit per-GPU value
# (--gpu_concurrency); 'auto' (the default) uses gpu_slots_per_gpu above.
#
# Anything that can't be determined confidently falls back to the SAFE side: unknown
# VRAM/footprint means 1 slot per GPU (strict serialization), never more.
#
# NOTE: like detect_gpu.sh, this describes the host it runs on. Under executors where
# the Nextflow launch host differs from the compute node (most cluster schedulers),
# conf/base.config only uses it for the `local` executor.
#
# Test hooks (env vars, unset in normal use) so a large/small/GPU-less host can be
# mocked without touching /proc or a real GPU -- see tests/bin/:
#   RGAPROFILER_NPROC             overrides `nproc`
#   RGAPROFILER_MEMINFO_FILE      overrides /proc/meminfo
#   RGAPROFILER_CGROUP_MEM_FILE   overrides the cgroup memory limit file (default:
#                                 cgroup v2 memory.max, else v1 memory.limit_in_bytes)
#   nvidia-smi is looked up on PATH, so prepend a directory with a fake one.

# Worst-case VRAM (MiB) used by any ONE of the GPU-capable tools processing one chunk,
# MEASURED on an RTX A4500 (20GB) as peak nvidia-smi memory.used while each task held the
# GPU alone, for a ~300-sequence chunk (real --num_blocks 1000 size) / a 4913-sequence chunk
# (default --fasta_qc_chunk_size 5000 size):
#   DeepCoil2  10.9GB / 19.2GB   <- grows with chunk size (TensorFlow allocator), the worst case
#   DeepTMHMM   3.7GB /  4.0GB
#   DeepLoc2    3.1GB /  3.1GB
#   SignalP6    1.9GB /  1.9GB
# gpu_slots_per_gpu assumes every concurrent task needs the worst of these, i.e. it is
# deliberately conservative: it yields 1 slot (strict serialization) on cards up to 40GB,
# 2 on 48GB, 3 on 80GB. It is tool-agnostic on purpose; a VRAM-weighted lock (light
# tools sharing a card with each other) would be a natural extension.
GPU_TOOL_PEAK_VRAM_MB="${RGAPROFILER_GPU_TOOL_PEAK_VRAM_MB:-19200}"
# Fraction of a GPU's VRAM that concurrent tasks may plan to use (leaves headroom for
# the CUDA context/fragmentation and other users of the card).
GPU_VRAM_USABLE_FRACTION_PCT=85

# --- CPUs ---------------------------------------------------------------------------
cpus="${RGAPROFILER_NPROC:-$(nproc 2>/dev/null || echo 1)}"

# --- RAM ----------------------------------------------------------------------------
meminfo="${RGAPROFILER_MEMINFO_FILE:-/proc/meminfo}"
mem_kb="$(awk '/^MemTotal:/ {print $2; exit}' "$meminfo" 2>/dev/null || true)"
mem_kb="${mem_kb:-0}"
mem_gb=$(( mem_kb / 1024 / 1024 ))

# A cgroup limit (docker/Slurm/systemd) is what the kernel will actually enforce, and is
# often far below MemTotal. Only honoured if it is a real number (v2 says "max" when unlimited).
cg_file="${RGAPROFILER_CGROUP_MEM_FILE:-}"
if [[ -z "$cg_file" ]]; then
    for f in /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory/memory.limit_in_bytes; do
        [[ -r "$f" ]] && { cg_file="$f"; break; }
    done
fi
if [[ -n "$cg_file" && -r "$cg_file" ]]; then
    cg_bytes="$(head -n1 "$cg_file" 2>/dev/null || true)"
    if [[ "$cg_bytes" =~ ^[0-9]+$ ]]; then
        cg_gb=$(( cg_bytes / 1024 / 1024 / 1024 ))
        if (( cg_gb > 0 && (mem_gb == 0 || cg_gb < mem_gb) )); then mem_gb=$cg_gb; fi
    fi
fi

# --- GPUs ---------------------------------------------------------------------------
gpus=0
gpu_vram_mb=0
if command -v nvidia-smi >/dev/null 2>&1; then
    vram_list="$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null || true)"
    # keep only plain integers (a wedged driver prints error text, not numbers)
    vram_list="$(printf '%s\n' "$vram_list" | tr -d ' \r' | grep -E '^[0-9]+$' || true)"
    if [[ -n "$vram_list" ]]; then
        gpus="$(printf '%s\n' "$vram_list" | wc -l)"
        gpu_vram_mb="$(printf '%s\n' "$vram_list" | sort -n | head -n1)"
    fi
fi
# CUDA_VISIBLE_DEVICES restricts which GPUs this pipeline may use.
cvd="${CUDA_VISIBLE_DEVICES:-}"
if [[ -n "$cvd" && "$cvd" != "NoDevFiles" ]] && (( gpus > 0 )); then
    n_cvd="$(printf '%s' "$cvd" | tr ',' '\n' | grep -c . || true)"
    (( n_cvd < gpus )) && gpus=$n_cvd
fi

# --- GPU slots ----------------------------------------------------------------------
slots_per_gpu=1
if [[ "$GPU_TOOL_PEAK_VRAM_MB" =~ ^[0-9]+$ ]] && (( GPU_TOOL_PEAK_VRAM_MB > 0 && gpu_vram_mb > 0 )); then
    slots_per_gpu=$(( gpu_vram_mb * GPU_VRAM_USABLE_FRACTION_PCT / 100 / GPU_TOOL_PEAK_VRAM_MB ))
    (( slots_per_gpu < 1 )) && slots_per_gpu=1
fi

if [[ "${1:-}" == "--gpu-slots" ]]; then
    want="${2:-auto}"
    if [[ "$want" =~ ^[0-9]+$ && "$want" -ge 1 ]]; then per_gpu=$want; else per_gpu=$slots_per_gpu; fi
    n_gpus=$gpus
    (( n_gpus < 1 )) && n_gpus=1
    echo $(( per_gpu * n_gpus ))
    exit 0
fi

echo "cpus=${cpus}"
echo "mem_gb=${mem_gb}"
echo "gpus=${gpus}"
echo "gpu_vram_mb=${gpu_vram_mb}"
echo "gpu_slots_per_gpu=${slots_per_gpu}"
