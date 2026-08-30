# omatheuspimenta/rgaprofiler

[![GitHub Actions CI Status](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/nf-test.yml/badge.svg)](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/linting.yml/badge.svg)](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

**omatheuspimenta/rgaprofiler** is a bioinformatics pipeline that predicts RGAs (Resistance Gene Analogs) in plant proteomes. Given one or more protein FASTA files, it cleans and deduplicates the input, runs six independent protein-prediction tools in parallel (DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6, DeepTMHMM), and combines their outputs into per-protein RGA family/subclass calls using the classification logic from [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor), plus a self-contained HTML summary report. See [`docs/output.md`](docs/output.md) for the full output structure.

1. Input QC: deduplicate, strip trailing stop codons, split into chunks ([`FASTA_QC`](docs/output.md#fasta_qc))
2. Run six prediction tools in parallel: coiled-coil domains ([`DeepCoil2`](docs/output.md#deepcoil2)), signal peptides + TM topology ([`Phobius`](docs/output.md#phobius)), domain/functional annotation ([`InterProScan`](docs/output.md#interproscan)), subcellular localization ([`DeepLoc2`](docs/output.md#deeploc2)), signal peptides ([`SignalP6`](docs/output.md#signalp6)), transmembrane helices ([`DeepTMHMM`](docs/output.md#deeptmhmm))
3. Classify each protein as an RGA (family/subclass) from the combined evidence ([`RGA_CLASSIFY`](docs/output.md#rga-classification))
4. Render a self-contained HTML summary report ([`RGA_REPORT`](docs/output.md#summary-report))

DeepCoil2, DeepLoc2, SignalP6, and DeepTMHMM can run on a GPU (`--use_gpu`); Phobius, InterProScan, and the RGA classification/report steps are CPU-only. Several of the underlying tools (InterProScan's database, DeepTMHMM/SignalP6/DeepLoc2's model weights) are license-gated and must be downloaded separately by the user — see [`docs/software-setup.md`](docs/software-setup.md).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow.

> [!IMPORTANT]
> Every one of the steps below is required before **any** run of this pipeline succeeds — including the small bundled test in step 5. There is no "just try it first" shortcut: the pipeline checks for all four license-gated software/database sets (InterProScan, DeepTMHMM, SignalP6, DeepLoc2) before it starts any work, test data or not, and refuses to run if any of them is missing.

Follow these steps in order, on the machine that will actually run the pipeline:

### 1. Get the pipeline

```bash
git clone https://github.com/omatheuspimenta/rgaprofiler.git
cd rgaprofiler
```

Clone it rather than letting `nextflow run omatheuspimenta/rgaprofiler` fetch it for you on first use. The next steps download several tens of GB of license-gated software into a `softwares/` folder inside the pipeline directory and run a setup script "from the pipeline root" — cloning gives you an obvious, visible folder to do that in. (Once everything below is downloaded, you can switch to running `nextflow run omatheuspimenta/rgaprofiler` from anywhere and point `--softwares_dir`/`--interproscan_db` at this cloned folder if you prefer; that's an optional convenience, not required.)

### 2. Install Nextflow

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (`>=25.10.4`) and make sure it's on your `PATH`:

```bash
nextflow -version
```

### 3. Install Docker

Install [Docker](https://docs.docker.com/get-docker/) and make sure it's running:

```bash
docker info
```

This pipeline's tools ship as ready-to-use public Docker images on GHCR — Docker pulls them automatically the first time each is needed, so there is nothing to build yourself. `-profile docker` is the only container profile this pipeline is built and tested with (see [`docs/usage.md`](docs/usage.md#-profile) for why the others aren't recommended).

### 4. Download the license-gated software (all four, every time)

A few of the tools this pipeline runs need model weights or databases that cannot legally be bundled into the Docker images, so you download them yourself, once, into a `softwares/` folder that stays local to your machine (git-ignored, never uploaded anywhere):

| Tool | What | Required even for the test data? |
| --- | --- | --- |
| InterProScan | Its full release + member-database data (tens of GB) | Yes |
| DeepTMHMM | 5 model checkpoints + 3 ESM1b weight files (academic license) | Yes |
| SignalP 6.0 | Model weights for one run mode (academic license) | Yes |
| DeepLoc 2 | Classifier checkpoints + ESM1b base encoder (academic license) | Yes |

Go to [`docs/software-setup.md`](docs/software-setup.md) now and work through all four "Per-tool setup" sections there — it tells you exactly where to go, what to download, and what command to run to unpack/place each one. Come back here once every tool's files are in place.

Then verify nothing is missing before you burn time on a run that will just fail at the preflight check:

```bash
./bin/check_software_present.sh interproscan   <path/to/interproscan-5.XX-YY.0>
./bin/check_software_present.sh deeptmhmm      softwares
./bin/check_software_present.sh signalp6       softwares
./bin/check_software_present.sh deeploc2       softwares
```

Each prints `ERROR: ...` and exits non-zero if something required is missing, naming exactly what and where to get it. Don't move on until all four exit cleanly.

### 5. Run the bundled test data

Confirm your Nextflow/Docker/software setup works end to end using the small dataset bundled with the pipeline (a handful of real sequences, finishes in a few minutes; the first run additionally downloads ~7 Docker images, a few GB total):

```bash
nextflow run . \
   -profile docker,test \
   --interproscan_db /path/to/interproscan-5.XX-YY.0 \
   --outdir results_test
```

`--interproscan_db` is the path `docs/software-setup.md`'s InterProScan setup step gave you; the other three tools' software is picked up automatically from `./softwares` (step 4). If this finishes without errors, your setup is ready for real data.

### 6. Prepare a samplesheet for your own data

```csv title="samplesheet.csv"
sample,fasta
sample1,/path/to/sample1.protein.fasta
```

Each row is one protein FASTA to profile. `sample` is a free-form identifier; `fasta` must be an **absolute path** that exists and ends in `.fa`/`.fasta` (optionally gzipped). See [`docs/usage.md`](docs/usage.md#samplesheet-input) for the full explanation, including why a relative path here is risky.

### 7. Run the pipeline on your own data

```bash
nextflow run . \
   -profile docker \
   --input samplesheet.csv \
   --interproscan_db /path/to/interproscan-5.XX-YY.0 \
   --outdir <OUTDIR>
```

If your institution provides its own [nf-core/configs](https://github.com/nf-core/configs) profile, you can add it alongside, e.g. `-profile docker,<institute>`. See [`docs/usage.md`](docs/usage.md) for the full list of parameters, profiles (including `long_running` for a full-scale proteome), and GPU options (`--use_gpu`).

For a large proteome, add `--num_blocks <N>` (e.g. `--num_blocks 1000`) to split each sample's input into that many sequence blocks — DeepCoil2, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM then each run once per block instead of once on the whole proteome, letting Nextflow schedule more independent tasks in parallel (still bounded by your executor/resource configuration). This matters most for DeepCoil2, which can fail or become impractical on a very large single-task input. See [`docs/output.md`](docs/output.md) for how per-block outputs are merged back together.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

> [!TIP]
> A `-params-file` (YAML/JSON) is the more robust option of the two, especially for
> typed parameters like `--num_blocks`/`--use_gpu` — see [`docs/usage.md`](docs/usage.md#running-the-pipeline)
> and the filled-in example at [`assets/params.example.yml`](assets/params.example.yml).

## Credits

omatheuspimenta/rgaprofiler was originally written by Pimenta-Zanon, M. H..

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- Add a Zenodo DOI citation here after the pipeline's first tagged release, and a matching "Cite with Zenodo" badge at the top of this file. -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
