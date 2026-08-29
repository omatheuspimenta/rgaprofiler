# SignalP 6.0 (DTU Health Tech) container -- "baseline" image.
#
# Unlike DeepLoc2, SignalP6's own package_data ('signalp/model_weights/*') is
# genuinely empty in the distributed source (just a README -- the real ~9.2GB
# weight checkpoints live in a *sibling* 'models/' directory outside the
# Python package and are copied in by hand post-install, per the package's
# own README). So the vendored source here (docker/signalp6/src/, code only,
# confirmed weight-free) can be pip-installed directly at build time without
# risking baking in any DTU-licensed weights.
#
# At runtime, the module points --model_dir at the mounted
# softwares_dir/SignalP6/signalp-6-package/models/ (same pattern as
# interproscan_db) -- SignalP6 exposes this as a first-class CLI flag, so no
# symlink trick is needed (unlike deeploc2). Only the 'sequential_models_signalp6'
# ("slow-sequential" mode) weight set is present in this repo's softwares/, so
# that's this pipeline's default --mode (see modules/local/signalp6/main.nf).

FROM python:3.8-slim-bookworm

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000

# procps (for `ps`) is required by Nextflow itself to collect task resource
# metrics inside the container. zip/unzip are required by the vendored
# signalp6_convert_models script (see docs/software-setup.md's GPU conversion
# step) -- it unzips a checkpoint's TorchScript archive, patches it, then
# re-zips it; without these it fails with "command not found" partway through.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        zip \
        unzip \
    && rm -rf /var/lib/apt/lists/*

# Pinned to the real, working `signalp6` conda env on this host (confirmed to
# reproduce the committed R570 ground truth, softwares/SignalP6/signalp-6-package/r570/).
RUN pip install --upgrade pip==24.0 setuptools==69.5.1 wheel && \
    pip install \
        numpy==1.24.4 \
        matplotlib==3.7.5 \
        tqdm==4.70.0 \
        torch==1.13.1

# Vendored SignalP6 source (code only -- signalp/model_weights/ is empty by
# upstream design, see header comment above). No cmdclass/atexit hooks in its
# setup.py, so a plain `pip install .` is safe and gives us the real `signalp6`
# console script (unlike deeploc2, which needed a hand-written shim).
COPY src/signalp-6-package /opt/signalp6-src
RUN pip install /opt/signalp6-src --no-deps && rm -rf /opt/signalp6-src

RUN groupadd --gid 1000 signalp6 \
    && useradd --uid 1000 --gid signalp6 --create-home --shell /bin/bash signalp6

# Nextflow's docker profile runs containers as the host UID, not this image's
# own uid-1000 user -- USER here just keeps the image sane for anyone running
# it directly (matches docker/phobius, docker/interproscan).
USER signalp6
WORKDIR /data
