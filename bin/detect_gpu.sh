#!/usr/bin/env bash
set -euo pipefail

# detect_gpu.sh — exit 0 if a usable NVIDIA GPU is visible on this host, exit 1 otherwise.
#
# Used by nextflow.config to resolve `--use_gpu auto`. Safe to run standalone:
#   ./bin/detect_gpu.sh && echo "GPU available" || echo "No GPU"
#
# NOTE: this only tells you whether a GPU is visible on the host where it runs.
# Under executors where the Nextflow launch host differs from the host(s) that
# actually execute GPU-labelled processes (e.g. most cluster/HPC schedulers),
# 'auto' detection here reflects the launch host, not necessarily the compute
# node — force --use_gpu true/false explicitly in that case instead of 'auto'.

command -v nvidia-smi >/dev/null 2>&1 || exit 1
nvidia-smi -L 2>/dev/null | grep -q '^GPU' || exit 1
exit 0
