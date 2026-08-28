# DeepTMHMM 1.0 (DTU Health Tech / BioLib) container -- "baseline" image.
#
# ACADEMIC / NON-COMMERCIAL LICENSE (CC BY-NC-SA 4.0): this is the most
# restrictive license among this pipeline's tools. Do not use DeepTMHMM
# outputs for commercial purposes without separately obtaining rights to do
# so -- see LICENSE.txt (embedded below) and docs/software-setup.md. The
# package also embeds pytorch-crf (pytorchcrf/LICENSE.txt) and the ESM1b
# model architecture (ESM-License.txt); both are carried into this image
# for attribution.
#
# DeepTMHMM ships as bare scripts, not a pip package (no setup.py at all) --
# predict.py resolves its 8 weight files (5 CV checkpoints + 3 ESM1b files)
# via bare relative filenames in the current working directory, with no
# --model-dir-style override. So this image vendors only the *code*
# (docker/deeptmhmm/src/, confirmed weight-free) into a fixed, world-writable
# /opt/deeptmhmm/, and the module symlinks the real weights there at runtime
# (see modules/local/deeptmhmm/main.nf) -- same reasoning as deeploc2's
# models/ symlink trick, forced here by predict.py having no CLI override at all.
#
# The upstream README installs torch==1.5.0+cu92 (2020-era, CUDA 9.2), which
# this repo confirmed does NOT run on modern Ampere/Ada GPUs ("RuntimeError:
# cublas runtime error"). softwares/DeepTMHMM/.../README_py310.txt documents
# a validated fix (Python 3.10, a modern torch build, and patching every
# torch.load() call for PyTorch 2.6+'s weights_only default change) that was
# already confirmed end-to-end on this host's real RTX A4500 -- this image
# reproduces that working configuration instead of the broken original one.
# The two torch.load() calls that upstream's fix left unpatched
# (docker/deeptmhmm/src/util.py) are patched here too, for consistency;
# experiments/tmhmm3/tm_models.py's remaining unpatched torch.load() calls
# are dead code for inference (training/distributed-eval only, never called
# by predict.py) and are deliberately left as-is.

FROM python:3.10-slim-bookworm

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000

# procps (for `ps`) is required by Nextflow itself to collect task resource
# metrics inside the container.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

# Pinned to softwares/DeepTMHMM/.../README_py310.txt's validated modern
# replacement stack, confirmed working against this host's real GPU.
RUN pip install --upgrade pip==24.0 setuptools==69.5.1 wheel && \
    pip install \
        torch==2.13.0 \
        biopython==1.79 \
        fair-esm==0.4.0 \
        h5py==3.7.0 \
        matplotlib==3.5.2 \
        numpy==1.23.4 \
        PeptideBuilder==1.1.0 \
        tqdm==4.64.1

# Vendored DeepTMHMM source (code only, no weights -- see header comment above).
COPY src /opt/deeptmhmm

RUN groupadd --gid 1000 deeptmhmm \
    && useradd --uid 1000 --gid deeptmhmm --create-home --shell /bin/bash deeptmhmm \
    && chmod -R a+rwX /opt/deeptmhmm

# Nextflow's docker profile runs containers as the host UID, which has no
# passwd entry in the image and thus no write access via ownership -- /opt/deeptmhmm
# must stay writable by an arbitrary UID (the module symlinks the staged weight
# files there at task runtime, predict.py resolves them by bare relative filename).
USER deeptmhmm
WORKDIR /data
