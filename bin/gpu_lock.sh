#!/usr/bin/env bash
# gpu_lock.sh — host-level mutual exclusion for GPU-labelled tasks. SOURCE it, don't run it:
#
#   source bin/gpu_lock.sh <slots_per_gpu|auto|0>
#
# It is wired in as the `beforeScript` of the `process_gpu` label (conf/base.config).
# Nextflow runs `beforeScript` on the HOST, inside the task's `.command.run` wrapper
# and before the container is launched (verified against Nextflow 26.04), so a flock
# taken here is held on the host for the whole container run and is released by the
# kernel the moment the wrapper (and its children) exit -- normal end, failure, or kill.
# No bind-mounted lock directory and no `flock` inside the container images are needed,
# and it behaves the same for docker, podman, singularity/apptainer and conda.
#
# Why this exists: Nextflow's `local` executor has no notion of a GPU slot, so without
# it DeepCoil2/DeepLoc2/SignalP6/DeepTMHMM chunk tasks all hit the same card at once
# (real CUDA OOMs). Exclusivity here is independent of the CPU/memory admission
# budget, so unrelated CPU-bound processes (InterProScan chunks) can neither starve
# nor be starved by it.
#
# Semantics
#   - One lock file per (GPU, slot) under $RGAPROFILER_GPU_LOCK_DIR
#     (default: ${TMPDIR:-/tmp}/rgaprofiler-gpu-locks-<uid>). Lock files are host-local
#     by design: flock is not reliable on network filesystems.
#   - <slots_per_gpu>: max concurrent tasks per GPU. 'auto' derives it from VRAM
#     (bin/detect_host_resources.sh; 1 whenever that can't be determined confidently).
#     0 disables locking entirely (for schedulers that already allocate GPUs).
#   - On acquiring a slot, exports CUDA_VISIBLE_DEVICES/CUDA_DEVICE_ORDER so the task
#     lands on the GPU whose slot it holds (forwarded into docker by containerOptions
#     `-e`; singularity/apptainer inherit the host environment). If the caller already
#     restricted CUDA_VISIBLE_DEVICES, only those devices are used.
#   - Fails OPEN: if locking is impossible (no `flock`, unwritable lock dir, no GPU
#     found) it prints a warning and lets the task run unlocked -- i.e. the old
#     un-serialized behaviour -- rather than failing the whole run over a lock.
#   - Wait/acquire timing is appended to `.gpu_lock.log` in the task directory. NOTE:
#     Nextflow starts the task's realtime clock before this runs, so trace.txt's
#     `realtime` for a GPU task includes any time it spent waiting here.

_rgaprofiler_gpu_lock() {
    local per_gpu="${1:-auto}"
    local here
    here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    [[ "$per_gpu" == "0" ]] && return 0

    if ! command -v flock >/dev/null 2>&1; then
        echo "[rgaprofiler] WARNING: 'flock' not found on host; GPU tasks will NOT be serialized." | tee -a .gpu_lock.log >&2
        return 0
    fi

    local lockdir="${RGAPROFILER_GPU_LOCK_DIR:-${TMPDIR:-/tmp}/rgaprofiler-gpu-locks-$(id -u)}"
    if ! mkdir -p "$lockdir" 2>/dev/null || [[ ! -w "$lockdir" ]]; then
        echo "[rgaprofiler] WARNING: GPU lock dir '$lockdir' not writable; GPU tasks will NOT be serialized." | tee -a .gpu_lock.log >&2
        return 0
    fi

    # Which GPUs, and how many slots each.
    local info gpus slots
    info="$(bash "$here/detect_host_resources.sh" 2>/dev/null || true)"
    gpus="$(printf '%s\n' "$info" | sed -n 's/^gpus=//p')"
    slots="$(printf '%s\n' "$info" | sed -n 's/^gpu_slots_per_gpu=//p')"
    [[ "$per_gpu" =~ ^[0-9]+$ ]] && slots="$per_gpu"
    [[ "$slots" =~ ^[0-9]+$ && "$slots" -ge 1 ]] || slots=1
    # A forced `--use_gpu true` on a host where detection sees no GPU still gets one lock.
    [[ "$gpus" =~ ^[0-9]+$ && "$gpus" -ge 1 ]] || gpus=1

    local devs=() d
    if [[ -n "${CUDA_VISIBLE_DEVICES:-}" && "${CUDA_VISIBLE_DEVICES}" != "NoDevFiles" ]]; then
        IFS=',' read -r -a devs <<< "${CUDA_VISIBLE_DEVICES}"
    else
        for (( d = 0; d < gpus; d++ )); do devs+=("$d"); done
    fi

    local total=$(( ${#devs[@]} * slots ))
    local t0 dev s fd name got
    t0=$(date +%s)

    while :; do
        for dev in "${devs[@]}"; do
            name="${dev//[^A-Za-z0-9_-]/_}"
            for (( s = 0; s < slots; s++ )); do
                exec {fd}>>"$lockdir/gpu${name}.slot${s}" || continue
                got=0
                # A single slot overall: block on the kernel's wait queue instead of polling.
                if (( total == 1 )); then flock "$fd" && got=1; else flock -n "$fd" && got=1; fi
                if (( got )); then
                    # Deliberately NOT closed: the fd (inherited by the container client)
                    # is what holds the lock until the whole task wrapper exits.
                    export CUDA_VISIBLE_DEVICES="$dev"
                    export CUDA_DEVICE_ORDER=PCI_BUS_ID
                    echo "acquired gpu${name}.slot${s} (of ${total}) after $(( $(date +%s) - t0 ))s at epoch $(date +%s.%N)" >> .gpu_lock.log
                    return 0
                fi
                exec {fd}>&-
            done
        done
        sleep "$(( 1 + RANDOM % 3 ))"   # jittered poll; only reached when total > 1 and all slots busy
    done
}

_rgaprofiler_gpu_lock "$@" || true
