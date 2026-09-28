#!/usr/bin/env bash
# gpu_lock.sh — host-level mutual exclusion for GPU-labelled tasks. SOURCE it, don't run it:
#
#   source bin/gpu_lock.sh <slots_per_gpu|auto|0> [vram_mb|all]
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
#   - [vram_mb|all] (only with 'auto'): VRAM-weighted sharing instead. Each GPU's usable
#     VRAM (GPU_VRAM_USABLE_FRACTION_PCT of the smallest card, as in
#     detect_host_resources.sh) is split into GPU_LOCK_UNIT_MB units, one lock file each,
#     and the task holds ceil(vram_mb / unit) of them on ONE GPU -- so light tools share a
#     card while their summed peak VRAM still fits. 'all' takes every unit (exclusive).
#     Several units are grabbed atomically under a short per-GPU mutex, so two tasks can
#     never deadlock each holding part of what they need; the units themselves are still
#     plain flocks, released by the kernel whenever the holder exits. Falls back to one
#     exclusive slot per GPU when the VRAM can't be read.
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
    local per_gpu="${1:-auto}" need_mb="${2:-}"
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

    # Which GPUs, how many lock slots each, and how many of them this task needs.
    local info gpus slots vram need
    info="$(bash "$here/detect_host_resources.sh" 2>/dev/null || true)"
    gpus="$(printf '%s\n' "$info" | sed -n 's/^gpus=//p')"
    vram="$(printf '%s\n' "$info" | sed -n 's/^gpu_vram_mb=//p')"
    slots="$(printf '%s\n' "$info" | sed -n 's/^gpu_slots_per_gpu=//p')"
    need=1
    [[ "$per_gpu" =~ ^[0-9]+$ ]] && slots="$per_gpu"
    [[ "$slots" =~ ^[0-9]+$ && "$slots" -ge 1 ]] || slots=1
    # VRAM-weighted mode: slots become VRAM units, and the task needs several of them.
    if [[ "$per_gpu" == "auto" && -n "$need_mb" && "$vram" =~ ^[0-9]+$ && "$vram" -gt 0 ]]; then
        slots="$(bash "$here/detect_host_resources.sh" --gpu-units)"
        [[ "$slots" =~ ^[0-9]+$ && "$slots" -ge 1 ]] || slots=1
        need="$(bash "$here/detect_host_resources.sh" --gpu-units "$need_mb")"
        [[ "$need" =~ ^[0-9]+$ && "$need" -ge 1 && "$need" -le "$slots" ]] || need=$slots
    fi
    # A forced `--use_gpu true` on a host where detection sees no GPU still gets one lock.
    [[ "$gpus" =~ ^[0-9]+$ && "$gpus" -ge 1 ]] || gpus=1

    local devs=() d
    if [[ -n "${CUDA_VISIBLE_DEVICES:-}" && "${CUDA_VISIBLE_DEVICES}" != "NoDevFiles" ]]; then
        IFS=',' read -r -a devs <<< "${CUDA_VISIBLE_DEVICES}"
    else
        for (( d = 0; d < gpus; d++ )); do devs+=("$d"); done
    fi

    local total=$(( ${#devs[@]} * slots ))
    local t0 dev s fd mfd name
    local -a held
    t0=$(date +%s)

    while :; do
        for dev in "${devs[@]}"; do
            name="${dev//[^A-Za-z0-9_-]/_}"
            held=()
            if (( total == 1 )); then
                # A single slot overall: block on the kernel's wait queue instead of polling.
                exec {fd}>>"$lockdir/gpu${name}.slot0" || continue
                flock "$fd" && held+=("$fd")
            else
                # Grab `need` free slots of this GPU atomically w.r.t. other acquirers.
                exec {mfd}>>"$lockdir/gpu${name}.mutex" || continue
                flock "$mfd"
                for (( s = 0; s < slots && ${#held[@]} < need; s++ )); do
                    exec {fd}>>"$lockdir/gpu${name}.slot${s}" || continue
                    if flock -n "$fd"; then held+=("$fd"); else exec {fd}>&-; fi
                done
                exec {mfd}>&-
            fi
            if (( ${#held[@]} == need )); then
                # Deliberately NOT closed: these fds (inherited by the container client)
                # are what hold the lock until the whole task wrapper exits.
                export CUDA_VISIBLE_DEVICES="$dev"
                export CUDA_DEVICE_ORDER=PCI_BUS_ID
                echo "acquired ${need} of ${slots} slot(s) on gpu${name} (${total} in all) after $(( $(date +%s) - t0 ))s at epoch $(date +%s.%N)" >> .gpu_lock.log
                return 0
            fi
            for fd in "${held[@]}"; do exec {fd}>&-; done
        done
        sleep "$(( 1 + RANDOM % 3 ))"   # jittered poll; only reached when total > 1 and too few slots are free
    done
}

_rgaprofiler_gpu_lock "$@" || true
