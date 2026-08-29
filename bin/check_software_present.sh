#!/usr/bin/env bash
set -euo pipefail

# check_software_present.sh — pre-flight validator for user-supplied,
# license-gated third-party software required by rgaprofiler.
#
# Usage:
#   check_software_present.sh <tool> [base_dir]
#
# <tool>      one of: interproscan, deeptmhmm, signalp6, deeploc2
# [base_dir]  directory to check under. Defaults to $RGAPROFILER_SOFTWARES_DIR
#             or ./softwares. For 'interproscan' this should be the InterProScan
#             installation directory itself (i.e. params.interproscan_db), since
#             interproscan_db is an independent pipeline parameter, not a
#             subdirectory of softwares_dir.
#
# Exits 0 if every REQUIRED file/dir for the tool is present (missing OPTIONAL
# files are reported as warnings on stderr but do not fail the check).
# Exits 1 with an actionable message (what's missing, where to get it, where
# to put it) if any REQUIRED file/dir is missing.
# Exits 2 for a usage error (unknown tool name).
#
# See docs/software-setup.md for the full per-tool download/placement guide.

tool="${1:?Usage: $(basename "$0") <tool> [base_dir]}"
base_dir="${2:-${RGAPROFILER_SOFTWARES_DIR:-softwares}}"

missing_required=()
missing_optional=()

require() { [ -e "$1" ] || missing_required+=("$1"); }
optional() { [ -e "$1" ] || missing_optional+=("$1"); }

help_msg=""

case "$tool" in
    interproscan)
        # base_dir = the InterProScan install dir itself (params.interproscan_db)
        require "${base_dir}/interproscan.sh"
        require "${base_dir}/data"
        optional "${base_dir}/bin/phobius/1.01/phobius.pl"
        optional "${base_dir}/bin/signalp/4.1/signalp"
        optional "${base_dir}/bin/tmhmm/2.0c/decodeanhmm"
        help_msg="Run ./setup_interproscandb.sh from the pipeline root, or download/unpack
InterProScan yourself from:
  https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
then pass its directory via --interproscan_db /path/to/interproscan-5.78-109.0

The licensed sub-analyses (Phobius, SignalP 4.1, TMHMM 2.0c) are optional —
InterProScan runs fine without them, just skipping those specific analyses.
To enable them, see:
  https://github.com/ebi-pf-team/interproscan/wiki/ActivatingLicensedAnalyses"
        ;;
    deeptmhmm)
        tool_dir="${base_dir}/DeepTMHMM/DeepTMHMM-Academic-License-v1.0"
        for f in deeptmhmm_cv_0.model deeptmhmm_cv_1.model deeptmhmm_cv_2.model \
                 deeptmhmm_cv_3.model deeptmhmm_cv_4.model esm_model_alphabet.pt \
                 esm_model_args.pt esm_model_state_dict.pt; do
            require "${tool_dir}/${f}"
        done
        help_msg="DeepTMHMM is distributed under a CC BY-NC-SA 4.0 (non-commercial) academic
license and its model weights cannot be redistributed with this pipeline.
Obtain the DeepTMHMM Academic License release yourself and place the 5 model
checkpoints (deeptmhmm_cv_0..4.model) and the 3 ESM1b weight files
(esm_model_alphabet.pt, esm_model_args.pt, esm_model_state_dict.pt) under:
  ${tool_dir}/"
        ;;
    signalp6)
        tool_dir="${base_dir}/SignalP6/signalp-6-package/models"
        gpu_dir="${base_dir}/SignalP6/signalp-6-package/models_gpu"
        # 3rd arg: whether GPU execution was requested/detected for this run (see
        # workflows/rgaprofiler.nf). SignalP6 has no runtime --device flag -- GPU vs CPU
        # is baked into the weight files themselves -- so a *separate*, GPU-converted
        # copy of the weights is only required when GPU is actually in play.
        use_gpu="${3:-false}"
        # This pipeline defaults to --mode slow-sequential (see modules/local/signalp6),
        # since it's the only mode whose weights are actually distributed/installed
        # in this repo's softwares/ -- so its weight set is required, not optional.
        # The 'fast' (distilled_model_signalp6.pt) and 'slow' (ensemble_model_signalp6.pt)
        # modes are separate downloads (see docs/software-setup.md) and stay optional.
        require "${tool_dir}/sequential_models_signalp6"
        optional "${tool_dir}/distilled_model_signalp6.pt"
        optional "${tool_dir}/ensemble_model_signalp6.pt"
        help_msg="SignalP 6.0 model weights are distributed separately by DTU Health Tech
under an academic license and cannot be redistributed with this pipeline.
Download the SignalP 6.0 package from:
  https://services.healthtech.dtu.dk/services/SignalP-6.0/
and place its models/ directory contents (at least the 'slow-sequential' mode's
sequential_models_signalp6/ weights) under:
  ${tool_dir}/"
        if [ "${use_gpu}" = "true" ]; then
            require "${gpu_dir}/sequential_models_signalp6"
            help_msg="${help_msg}

GPU execution is requested/detected for this run. SignalP6 has no runtime --device
flag -- GPU vs CPU is baked into the weight files themselves -- so it needs a
*separate*, GPU-converted copy of the weights above, at ${gpu_dir}/.
One-time setup (run once, from the signalp6 image so signalp6_convert_models and
its Python dependencies are available):
  docker run --rm -v \"\$(pwd)/${base_dir}/SignalP6\":/data quay.io/signalp6:baseline bash -c '
      cp -r /data/signalp-6-package/models /data/signalp-6-package/models_gpu &&
      signalp6_convert_models gpu /data/signalp-6-package/models_gpu
  '
This copies the CPU weights first so the original models/ is left untouched, then
converts the copy in place. See docs/software-setup.md. (Or pass --use_gpu false to
run on CPU instead -- slower, but needs no extra setup.)"
        fi
        ;;
    deeploc2)
        tool_dir="${base_dir}/DeepLoc2/DeepLoc2/models"
        require "${tool_dir}/models_esm1b"
        require "${tool_dir}/models_prott5"
        require "${tool_dir}/ESM1b_alphabet.pkl"
        require "${tool_dir}/ProtT5_alphabet.pkl"
        # The default 'Fast' model additionally needs the underlying 650M-parameter
        # ESM1b transformer (a torch-hub weight, distinct from the DTU classifier
        # heads above) -- not DTU-restricted, but too large (~7.3GB) to bake into
        # the baseline image, so it's supplied the same way. 'Accurate' mode's
        # ProtT5-XL encoder (Hugging Face, ~11GB) is intentionally not checked here:
        # it's optional (not the pipeline's default model) and downloads/caches
        # itself into the same torch_cache mount on first use if HF_HOME is set there.
        torch_cache_dir="${base_dir}/DeepLoc2/torch_cache"
        require "${torch_cache_dir}/hub/checkpoints/esm1b_t33_650M_UR50S.pt"
        require "${torch_cache_dir}/hub/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt"
        help_msg="DeepLoc 2 model weights are distributed separately by DTU Health Tech
under an academic license and cannot be redistributed with this pipeline.
Download DeepLoc 2 from:
  https://services.healthtech.dtu.dk/services/DeepLoc-2.1/
and place its models directory under:
  ${tool_dir}/

DeepLoc 2's default 'Fast' model also needs the ESM1b base encoder (fair-esm,
freely downloadable but ~7.3GB, so it's not baked into the image either). Fetch
it once (e.g. by running fair-esm's pretrained.load_model_and_alphabet(\"esm1b_t33_650M_UR50S\")
or letting a real DeepLoc2 install download it) and place the two resulting
files under:
  ${torch_cache_dir}/hub/checkpoints/"
        ;;
    *)
        echo "ERROR: unknown tool '${tool}'. Known tools: interproscan, deeptmhmm, signalp6, deeploc2" >&2
        exit 2
        ;;
esac

if [ "${#missing_optional[@]}" -gt 0 ]; then
    echo "WARNING: [${tool}] optional file(s) not found under '${base_dir}' (some analyses will be skipped):" >&2
    printf '  - %s\n' "${missing_optional[@]}" >&2
fi

if [ "${#missing_required[@]}" -gt 0 ]; then
    echo "ERROR: [${tool}] required file(s)/dir(s) not found under '${base_dir}':" >&2
    printf '  - %s\n' "${missing_required[@]}" >&2
    echo >&2
    echo "${help_msg}" >&2
    exit 1
fi

exit 0
