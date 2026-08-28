# DeepLoc 2.1 (DTU Health Tech) container -- "baseline" image.
#
# DeepLoc2's own package (github.com/Bio-Ontology-Research-Group/DeepLoc-2.0 /
# the DTU-distributed source tarball) is not on public PyPI, and its setup.py
# bundles the DTU model checkpoints as package_data *and* auto-downloads the
# 650M-parameter ESM1b encoder at install time via an atexit hook. Neither of
# those is redistributable/bakeable here (DTU academic-use license), so this
# image vendors only DeepLoc2's Python source (src/DeepLoc2/*.py, code only,
# no weights -- see docker/deeploc2/src/) and installs its dependencies
# pinned to the versions confirmed working against this repo's real R570
# reference run (softwares/DeepLoc2/env.yml).
#
# At runtime, the module bind-mounts (via Nextflow `path` staging, same
# pattern as the interproscan module):
#   - the DTU model checkpoints  -> DeepLoc2/models (pkg_resources-resolved)
#   - the ESM1b torch-hub cache  -> $TORCH_HOME
# See modules/local/deeploc2/main.nf and docs/software-setup.md.

FROM python:3.8-slim-bookworm

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

# Pinned to softwares/DeepLoc2/env.yml, the lockfile from a real working
# DeepLoc 2.1 install on this host (conda env `deeploc20`) -- confirmed to
# reproduce the committed R570 ground truth (softwares/DeepLoc2/r570/).
RUN pip install --upgrade pip==24.0 setuptools==69.5.1 wheel && \
    pip install \
        numpy==1.24.3 \
        pandas==1.5.3 \
        scipy==1.10.1 \
        scikit-learn==1.3.0 \
        matplotlib==3.7.2 \
        biopython==1.83 \
        h5py==3.11.0 \
        sentencepiece==0.2.0 \
        fair-esm==0.4.0 \
        onnxruntime==1.19.2 \
        torch==2.4.1 \
        pytorch-lightning==2.4.0 \
        transformers==4.46.3

# Vendored DeepLoc2 source (code only -- no models/, see header comment above).
COPY src/DeepLoc2 /opt/deeploc2/DeepLoc2
COPY bin/deeploc2 /usr/local/bin/deeploc2
RUN chmod +x /usr/local/bin/deeploc2

ENV PYTHONPATH="/opt/deeploc2"

RUN groupadd --gid 1000 deeploc2 \
    && useradd --uid 1000 --gid deeploc2 --create-home --shell /bin/bash deeploc2 \
    && chmod -R a+rwX /opt/deeploc2

# Nextflow's docker profile runs containers as the host UID (`-u $(id -u):$(id -g)`,
# see nextflow.config's `gpu`/general docker profile), not this image's own uid-1000
# user -- so /opt/deeploc2/DeepLoc2 must stay writable by an arbitrary UID (the module
# symlinks the staged model weights to DeepLoc2/models at task runtime; deliberately
# NOT pre-created here -- `ln -sfn X models` nests inside an existing real directory
# instead of replacing it, so `models` must not exist until the module creates it).
USER deeploc2
WORKDIR /data
