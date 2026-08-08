# DeepCoil (labstructbioinf/DeepCoil) container
#
# DeepCoil's PyPI release pins a 2020-era ML stack.
# The repo requires python>=3.7,<3.9 -> we use 3.7 to fetch pre-compiled wheels
# for older packages like spacy 2.1.9 and avoid legacy build-from-source errors.

FROM python:3.7-slim-bullseye AS base

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_DEFAULT_TIMEOUT=1000

# Build tools are needed to compile jsonnet (allennlp dependency) from source,
# plus git for fetching example data from the repo. libopenblas speeds up numpy.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        cmake \
        git \
        curl \
        libopenblas-dev \
        libssl-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/deepcoil

# --- Step 1: numpy + cython pinned first.
RUN pip install --upgrade pip==23.3.2 setuptools==65.5.1 wheel && \
    pip install \
        numpy==1.18.5 \
        "cython<3.0"

# --- Step 2: TensorFlow, pinned to the version deepcoil expects.
RUN pip install tensorflow==2.3.0

# --- Step 3: torch, pinned within deepcoil's torch<2.0,>=1.2 range.
RUN pip install torch==1.7.1

# --- Step 3.5: jsonnet pinned to an older version (2020 era) so allennlp
# doesn't fetch a modern version that requires Python 3.8+ build tools.
RUN pip install jsonnet==0.15.0

# --- Step 3.8: Pin spaCy's C-extensions to their 2019 versions. This prevents
# pip from fetching modern source distributions that require Cython >= 3.1.
RUN pip install \
        "murmurhash==1.0.2" \
        "cymem==2.0.2" \
        "preshed==2.0.1" \
        "blis==0.2.2" \
        "thinc==7.0.8"

# --- Step 4: allennlp 0.9.0 and its overrides dependency.
# Modern 'overrides' (v7+) does strict runtime type-checking that breaks
# allennlp 0.9.0. Pinning to a 3.x version fixes the TypeError.
RUN pip install overrides==3.1.0 allennlp==0.9.0

# --- Step 5: the remaining direct deepcoil dependencies.
RUN pip install \
        "biopython>=1.77" \
        "pandas==1.1.5" \
        "seaborn==0.11.2" \
        "matplotlib==3.3.4"

# --- Step 6: deepcoil itself.
RUN pip install deepcoil==2.0.1 --no-deps

# Fetch the repo purely for its bundled example FASTA/plot and the tests
RUN git clone --depth 1 https://github.com/labstructbioinf/DeepCoil.git /opt/deepcoil/repo

# Download adtional file
RUN apt-get update && apt-get install -y unzip && \
    mkdir -p /usr/local/lib/python3.7/site-packages/deepcoil/weights/seqvec && \
    curl -k -L -o seqvec.zip https://lbs.cent.uw.edu.pl/static/files/seqvec.zip && \
    unzip seqvec.zip -d /usr/local/lib/python3.7/site-packages/deepcoil/weights/seqvec && \
    rm seqvec.zip && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data