# omatheuspimenta/rgaprofiler

[![GitHub Actions CI Status](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/nf-test.yml/badge.svg)](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/linting.yml/badge.svg)](https://github.com/omatheuspimenta/rgaprofiler/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/omatheuspimenta/rgaprofiler)

## Introduction

**omatheuspimenta/rgaprofiler** is a bioinformatics pipeline that predicts RGAs (Resistance Gene Analogs) in plant proteomes. Given one or more protein FASTA files, it cleans and deduplicates the input, runs six independent protein-prediction tools in parallel (DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6, DeepTMHMM), and combines their outputs into per-protein RGA family/subclass calls using the classification logic from [`SugarcaneTranscriptomics`](https://github.com/omatheuspimenta/SugarcaneTranscriptomics), plus a self-contained HTML summary report. See [`docs/output.md`](docs/output.md) for the full output structure.

1. Input QC: deduplicate, strip trailing stop codons, split into chunks ([`FASTA_QC`](docs/output.md#fasta_qc))
2. Run six prediction tools in parallel: coiled-coil domains ([`DeepCoil2`](docs/output.md#deepcoil2)), signal peptides + TM topology ([`Phobius`](docs/output.md#phobius)), domain/functional annotation ([`InterProScan`](docs/output.md#interproscan)), subcellular localization ([`DeepLoc2`](docs/output.md#deeploc2)), signal peptides ([`SignalP6`](docs/output.md#signalp6)), transmembrane helices ([`DeepTMHMM`](docs/output.md#deeptmhmm))
3. Classify each protein as an RGA (family/subclass) from the combined evidence ([`RGA_CLASSIFY`](docs/output.md#rga-classification))
4. Render a self-contained HTML summary report ([`RGA_REPORT`](docs/output.md#summary-report))

DeepCoil2, DeepLoc2, SignalP6, and DeepTMHMM can run on a GPU (`--use_gpu`); Phobius, InterProScan, and the RGA classification/report steps are CPU-only. Several of the underlying tools (InterProScan's database, DeepTMHMM/SignalP6/DeepLoc2's model weights) are license-gated and must be downloaded separately by the user — see [`docs/software-setup.md`](docs/software-setup.md).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fasta
sample1,/path/to/sample1.protein.fasta
```

Each row is one protein FASTA to profile. `sample` is a free-form identifier; `fasta` must exist and end in `.fa`/`.fasta` (optionally gzipped).

You'll also need the license-gated software/databases each tool expects under `--softwares_dir` (default `./softwares`) — see [`docs/software-setup.md`](docs/software-setup.md) for exactly what to download and where to put it, and `--interproscan_db` pointed at your InterProScan database directory (no default, since it varies per install).

Now, you can run the pipeline using:

```bash
nextflow run omatheuspimenta/rgaprofiler \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --interproscan_db /path/to/interproscan-5.XX-YY.0 \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

omatheuspimenta/rgaprofiler was originally written by Pimenta-Zanon, M. H..

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- Add a Zenodo DOI citation here after the pipeline's first tagged release; update the badge at the top of this file too. -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
