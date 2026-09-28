#!/usr/bin/env bash
set -euo pipefail

# gpu_max_forks.sh — the `maxForks` value for a process_gpu-labelled process
# (conf/base.config), evaluated once when the config is parsed.
#
# Usage: gpu_max_forks.sh <use_gpu> <gpu_concurrency> [vram_mb|all]
#
# Prints nothing (= no cap) when GPU tasks won't take the GPU lock at all: --use_gpu
# false, --gpu_concurrency 0, or --use_gpu auto finding no usable GPU. Otherwise prints
# how many of this process's tasks can hold the lock at once (bin/detect_host_resources.sh
# --gpu-slots), so no more are admitted -- each reserving CPUs and RAM -- than can run.
# [vram_mb|all] is the process's own peak VRAM need, for the VRAM-weighted lock
# (bin/gpu_lock.sh); without it every task counts as one slot.

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
use_gpu="${1:-auto}"
concurrency="${2:-}"
need="${3:-}"

[[ "$concurrency" == "0" || "$use_gpu" == "false" ]] && exit 0
if [[ "$use_gpu" == "auto" ]] && ! bash "$here/detect_gpu.sh"; then exit 0; fi

# --gpu_concurrency unset reaches here as "null" (or empty) -> 'auto'
[[ "$concurrency" =~ ^[0-9]+$ ]] || concurrency=auto
bash "$here/detect_host_resources.sh" --gpu-slots "$concurrency" ${need:+"$need"}
