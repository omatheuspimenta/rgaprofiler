# omatheuspimenta/rgaprofiler: Usage

> _Every parameter, with its description and default, is listed by `nextflow run . --help`
> (`--help <parameter>` for the details of one, `--show_hidden` to include advanced ones)._

## Introduction

This pipeline predicts RGAs (Resistance Gene Analogs) in plant proteomes — it takes one
or more protein FASTA files, not sequencing reads. See the main [`README.md`](../README.md)
for a description of what it does and [`docs/output.md`](output.md) for what it produces.

> [!TIP]
> First time running this pipeline? Follow the numbered walkthrough in the
> [main `README.md`](../README.md#usage) instead of this page — it takes you from a fresh
> clone to a finished test run in order. This page is a reference for everything beyond that:
> every parameter, profile, and advanced option.

## Samplesheet input

You will need to create a samplesheet listing the protein FASTA file(s) you'd like to
analyse before running the pipeline. It's a comma-separated file with exactly 2 columns
and a header row, as shown below.

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
sample,fasta
R570,/absolute/path/to/R570.protein.fasta
another_sample,/absolute/path/to/another_sample.protein.fasta
```

| Column   | Description                                                                                                                             |
| -------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `sample` | A name for this proteome. Each tool's results for it go in a folder with this name (e.g. `rga/<sample>/`), and file names start with it. Spaces are automatically converted to underscores. |
| `fasta`  | Full path to a protein FASTA file for this sample. Must exist and end in `.fa`/`.fasta` (optionally gzipped, e.g. `.fasta.gz`).         |

Unlike read-based nf-core pipelines, there's no concept of "multiple runs of the same
sample" here (no lanes to concatenate) — one row is one FASTA to profile. If you have
multiple FASTA files that should be treated as a single proteome, concatenate them
yourself before listing the result as one row.

**Important**: use an **absolute path** in the `fasta` column, not a relative one.
Nextflow resolves a relative path in the samplesheet against the directory you _launch_
`nextflow run` from, not against wherever `samplesheet.csv` itself lives — a relative
path that looks correct can silently fail to resolve if you run the pipeline from a
different directory than the one you wrote the samplesheet in.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline, from a clone of this repository (see
[`README.md`](../README.md#1-get-the-pipeline) for why cloning is recommended over letting
Nextflow fetch the pipeline itself), is as follows:

```bash
nextflow run . \
    --input ./samplesheet.csv \
    --interproscan_db /path/to/interproscan-5.XX-YY.0 \
    --outdir ./results \
    -profile docker
```

If you've already populated `softwares/` on this machine and would rather not keep a local
clone around, you can instead let Nextflow pull and cache the pipeline itself, and point
`--softwares_dir`/`--interproscan_db` at wherever that software lives:

```bash
nextflow run omatheuspimenta/rgaprofiler \
    --input ./samplesheet.csv \
    --interproscan_db /path/to/interproscan-5.XX-YY.0 \
    --softwares_dir /path/to/softwares \
    --outdir ./results \
    -profile docker
```

`--interproscan_db` has no default (it varies per install) and is always required. Several
tools also need license-gated model weights/databases you supply yourself under
`--softwares_dir` (default `./softwares`) — see [`docs/software-setup.md`](software-setup.md).

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

> [!TIP]
> Both ways work: command-line flags (`--num_blocks 1000`) and a params file. A params
> file is the easiest way to keep and share the exact settings of a run.

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run omatheuspimenta/rgaprofiler -profile docker,long_running -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
interproscan_db: '/path/to/interproscan-5.78-109.0'
num_blocks: 'auto'      # or a number, e.g. 1000
use_gpu: 'auto'         # 'auto' (default), 'true' or 'false'
```

A filled-in copy of this is committed at [`assets/params.example.yml`](../assets/params.example.yml)
— copy it and edit the paths rather than starting from scratch.

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

If you're running from a clone (`nextflow run .`, the recommended path above), update it the
normal git way:

```bash
git pull
```

If instead you're running the pipeline by name (`nextflow run omatheuspimenta/rgaprofiler`),
Nextflow automatically pulls the pipeline code from GitHub on first use and stores it as a
cached version. When running the pipeline after this, it will always use the cached version
if available - even if the pipeline has been updated since. To make sure that you're running
the latest version of the pipeline, make sure that you regularly update the cached version of
the pipeline:

```bash
nextflow pull omatheuspimenta/rgaprofiler
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [omatheuspimenta/rgaprofiler releases page](https://github.com/omatheuspimenta/rgaprofiler/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

Keep the **chunking** fixed too (`--num_blocks` / `--fasta_qc_chunk_size`, see
[Chunking and reproducibility](#chunking-and-reproducibility)): which proteins share a
chunk slightly affects a few tools' numbers. Every run's parameters are saved to
`<outdir>/pipeline_info/params_<timestamp>.json`.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Which setup fits you?

The defaults work on any machine; these additions make the most of yours:

| Your situation | What to add |
|---|---|
| A computer with an NVIDIA GPU | nothing — the GPU is detected and used automatically (`--use_gpu auto`); just do SignalP6's one-time [GPU step](software-setup.md#signalp-60) during setup |
| A computer **without** a GPU | nothing is required; see [Running on a machine without a GPU](#running-on-a-machine-without-a-gpu) for what to expect and how to speed it up |
| A whole proteome (tens of thousands of proteins or more) | `-profile docker,long_running` and `--num_blocks auto` |
| A small protein set (hundreds to a few thousand proteins) | nothing; `--num_blocks auto` lets it use all your CPUs |
| Several proteomes at once | one row per proteome in the samplesheet — each gets its own result folders |
| A cluster (Slurm, SGE, …) | your institution's profile (`-profile docker,<institute>`), an explicit `--num_blocks N`, and [`resourceLimits`](#limiting-what-the-pipeline-may-use) |

## Sequence batching (`--num_blocks`)

To run fast and to cope with very large proteomes, the pipeline splits each cleaned
proteome into **chunks** and runs the six prediction tools once per chunk, in parallel.
Each tool's results are then merged back into one result per sample, so the files in
`--outdir` look the same however many chunks were used.

You choose the chunking with **one** of:

- **`--num_blocks auto`** — the pipeline picks the number of chunks from the number of
  proteins and the CPUs of your machine. The easiest choice on a single machine.
- **`--num_blocks <N>`** — exactly `N` chunks (e.g. `--num_blocks 1000`), balanced by
  protein count. Use this on a cluster, or to reproduce an earlier run exactly.
- **`--fasta_qc_chunk_size <N>`** (default `5000`, used when `--num_blocks` isn't set) — a
  fixed number of proteins per chunk instead.

```bash
# A whole proteome on one machine
nextflow run . \
    -profile docker,long_running \
    --input samplesheet.csv \
    --interproscan_db /path/to/interproscan-5.XX-YY.0 \
    --num_blocks auto \
    --outdir results
```

More chunks means more, smaller tasks that *can* run in parallel; how many actually run
at once is decided by your machine (or cluster) and the pipeline's resource settings.

**Checking what happened:** the log prints `Sequence batching:` lines — which setting is in
effect and, once the input is cleaned, how many chunks were made. With `--num_blocks auto`
the exact number chosen is in the FASTA_QC task's log (`num_blocks auto: ... rerun with
--num_blocks N to reproduce exactly`). If the default chunking would leave most of your
CPUs idle, the pipeline warns and suggests `--num_blocks auto`. The chunks themselves are
saved in `<outdir>/fasta/<sample>/<sample>_clean_chunks/`.

### Choosing a chunk count

| Situation | Suggestion |
|---|---|
| Any run on a single machine | `--num_blocks auto` |
| Whole proteome, explicit number | about one chunk per 300–500 proteins (e.g. 300,000 proteins → `--num_blocks 1000`) |
| Small protein set | at least ~100 proteins per chunk — every chunk pays a fixed start-up cost to load the models |
| GPU with 20 GB or less | at most ~5,000 proteins per chunk (DeepCoil2's GPU memory grows with chunk size) |
| Cluster | an explicit `--num_blocks`; more chunks = more, shorter jobs |

<details markdown="1">
<summary>Technical details</summary>

- `--num_blocks auto` makes chunks of 100–500 proteins: as many as needed for about two
  rounds of InterProScan tasks (the slowest step, 4 CPUs each) across the machine's CPUs,
  so no core idles at the end. The floor of 100 keeps the per-task model loading
  (~15–25 s per GPU tool, far more on CPU) from dominating; the cap of 500 bounds
  DeepCoil2's GPU memory and each task's run time. Examples: 2,000 proteins → 4 chunks on
  an 8-CPU machine, 20 on a 256-CPU one; a 300,000-protein proteome → 600. It reads the
  CPUs of the machine that runs Nextflow, which is why it is meant for single machines.
- Chunks are balanced by protein count (`seqkit split2 --by-part`); a proteome with fewer
  proteins than `N` gets one chunk per protein.
- There is no upper limit on the number of chunks: above 1,000 chunks per sample the merge
  steps automatically work in two stages (tested with 60,000 chunks).
- Phobius is single-threaded, which is why it is chunked too: on the full sugarcane R570
  proteome it took 8.8 h as a single task and ~7 min across 1,000 chunks, with an identical
  result.

</details>

### Chunking and reproducibility

Every protein is predicted exactly once, whatever the chunking. But a few tools process
proteins in small groups internally, so **which proteins share a chunk can very slightly
change their numbers**. Measured on the same 100 proteins in two different chunks:

| Tool | Effect of changing the chunking |
|---|---|
| Phobius, DeepLoc2, DeepTMHMM | none |
| SignalP6 | probabilities change by at most 0.000007; no prediction changed |
| DeepCoil2 | per-residue scores change by at most 0.003 on GPU (0.19 on CPU); no residue crossed the 0.5 call threshold |
| InterProScan | Gene3D domain boundaries/e-values differed for 1 of the 100 proteins |

That is why the pipeline never changes your chunking by itself. To compare runs exactly,
use the same `--num_blocks` / `--fasta_qc_chunk_size` (every run's parameters are saved in
`<outdir>/pipeline_info/params_<timestamp>.json`).

Even with identical settings, two things vary slightly from run to run (in any version of
the pipeline): the row order of InterProScan's table, and — on GPU — an occasional
last-digit change in a DeepCoil2 score (e.g. 0.247 vs 0.248), which comes from DeepCoil2's
GPU arithmetic. Neither changed an RGA call in any comparison we made. Results from a CPU
run and a GPU run also differ at this floating-point level.

## Running on a machine without a GPU

Nothing special is needed: with the default `--use_gpu auto`, the pipeline checks for an
NVIDIA GPU and, if there is none, runs DeepCoil2, DeepLoc2, SignalP6 and DeepTMHMM on the
CPU. (Use `--use_gpu false` to force CPU mode on a machine that does have a GPU.) SignalP6
does not need the GPU-converted weights (`models_gpu/`) in this case.

```bash
nextflow run . \
    -profile docker,long_running \
    --input samplesheet.csv \
    --interproscan_db /path/to/interproscan-5.XX-YY.0 \
    --num_blocks auto \
    --outdir results
```

**What to expect:**

- **Results** are equivalent to a GPU run; numbers differ only at floating-point level
  (e.g. SignalP6 probabilities in the 5th decimal). The RGA calls were identical on the
  test data.
- **Speed:** a CPU run is much slower. For a chunk of 300 proteins, each GPU tool takes about
  a minute on a GPU, but on 6–8 CPU cores: DeepCoil2 ~9–15 min, DeepLoc2 ~11–15 min,
  SignalP6 ~19–25 min. Use `-profile long_running` for anything beyond small protein sets.
- **DeepTMHMM** always uses one thread per CPU core of the whole machine (this is built
  into DeepTMHMM). So that several copies don't fight over the CPUs, on a single machine
  the pipeline runs DeepTMHMM one chunk at a time on CPU, each with the whole machine to
  itself (the other tools' tasks run before and after it). On machines with very many cores (100+) DeepTMHMM becomes very
  slow with that many threads, and you can opt in to fewer:

  ```bash
  --deeptmhmm_cpu_threads 6
  ```

  With this, DeepTMHMM runs with 6 threads and many chunks side by side (on a 256-core
  machine: ~14 min per chunk instead of ~3.75 h). The trade-off: its intermediate
  per-protein files (`deeptmhmm/<sample>/embeddings/`) then differ at floating-point
  level — the predicted topologies were identical in testing. Leave it unset to get
  DeepTMHMM's output exactly as DeepTMHMM itself produces it.

If you pass `--use_gpu true` (or `-profile gpu`) on a machine without a usable GPU, the
pipeline warns at start-up: every GPU task would fail to start its container.

## Resource sizing and GPU sharing

**You normally don't need to set any resources.** The pipeline works out what each task
needs:

- Every chunk task asks for memory and time based on what is actually in its chunk (how
  many proteins, how long they are) and on the tool — measured for each tool, not guessed.
- On a single machine (the `local` executor), it reads the machine's CPUs, memory and GPUs
  and never asks for more than the machine has. On a cluster it uses fixed per-chunk
  defaults instead, since the machine that launches Nextflow says nothing about the
  compute nodes.
- A task that runs out of memory is retried once with double the memory (and time).
- GPU tools share the GPU by memory (see [GPU sharing](#gpu-sharing---gpu_concurrency)).

### Limiting what the pipeline may use

On a shared machine or a cluster, cap every task with Nextflow's
[`resourceLimits`](https://www.nextflow.io/docs/latest/reference/process.html#resourcelimits)
in a small config file passed with `-c`:

```groovy title="limits.config"
process.resourceLimits = [ cpus: 64, memory: 200.GB, time: 48.h ]
```

```bash
nextflow run . -profile docker --input samplesheet.csv ... -c limits.config
```

On a cluster, set `time` to your queue's maximum. If a chunk task still runs out of
memory, use more chunks (a larger `--num_blocks`).

<details markdown="1">
<summary>How requests are computed (technical details)</summary>

DeepCoil2, DeepLoc2, SignalP6, DeepTMHMM (label `process_medium_chunk`) and InterProScan
(`process_high_chunk`) run once per chunk. At submission time `conf/base.config` reads the
chunk once (`task.ext.chunk`: number of proteins, total residues, longest protein) and sizes
the task from it. None of this changes any output — only what each task reserves and how
long it may run.

**Memory** = a per-tool floor + a term proportional to the chunk's FASTA size, × the
attempt number. Peak memory was measured for every task of a real 1,000-chunk R570 run
(joined to each chunk's composition) and for single proteins of 1,000–35,000 residues:

| Tool | 300-protein chunk | 4,913-protein chunk | one 35,000-residue protein | Requested (300-protein chunk) |
|---|---|---|---|---|
| SignalP6 | 3.1–3.6 GB | 4.3 GB | 3.5 GB (1.9 GB on CPU) | 5 GB |
| DeepTMHMM | 5.4–5.5 GB (6.5 GB on CPU) | 5.5 GB | 5.4 GB | 9 GB |
| DeepLoc2 | 9.9–10.3 GB | 10.3 GB | 10.0 GB | 13 GB |
| DeepCoil2 | 5.7–6.8 GB (4.0 GB on CPU) | 16.0 GB | 5.0 GB | 9 GB |
| InterProScan | 4.5–7.8 GB | ~13–14 GB (extrapolated) | 6.0 GB | 11 GB |

Memory is set mostly by the model, not by protein length (these tools truncate or window
long proteins); only DeepCoil2 and InterProScan grow with the amount in the chunk.

**Time** = the larger of the previous fixed limit (8 h; 16 h for InterProScan) and 5× the
measured throughput for the chunk's residues on that device (GPU or CPU), × the attempt
number — so a large chunk on CPU is not killed at a fixed limit. `-profile long_running`
replaces these with its own fixed multi-day limits.

**CPUs** match what each tool really uses: InterProScan 4 (it used ~3.1 cores at `-cpu 4`
or `6`); GPU tools 4; on CPU, DeepCoil2 and DeepLoc2 6, SignalP6 8 (it always runs 8
threads), and DeepTMHMM the whole machine (see [above](#running-on-a-machine-without-a-gpu)).

</details>

### GPU sharing (`--gpu_concurrency`)

On a single machine, GPU tasks share the GPU **by memory**: small tools (SignalP6,
DeepLoc2, DeepTMHMM on typical proteins) run side by side while they fit on the card, and
DeepCoil2 always gets the GPU to itself. This does not change any result, and it makes the
GPU part of a run about 1.5× faster (35.7 instead of 55.6 min on 2,400 proteins).

- **Default (unset):** share by memory, as above. If the GPU's memory can't be read, one
  task at a time.
- **`--gpu_concurrency N`:** exactly `N` GPU tasks at a time per GPU, whatever the tool
  (`1` = strictly one at a time).
- **`--gpu_concurrency 0`:** no limit — for clusters whose scheduler already hands out GPUs
  per job.

<details markdown="1">
<summary>How GPU sharing works (technical details)</summary>

Nextflow's `local` executor has no notion of a GPU, so the pipeline uses a host-level lock
(`bin/gpu_lock.sh`, `flock`-based, held for the whole task and released automatically if
the task dies). Each GPU's usable memory (85%) is divided into 512 MiB units, and a task
holds as many units as its tool's measured peak GPU memory (`ext.gpu_vram_mb` in
`conf/base.config`, with ~15% headroom):

- SignalP6 ~2.2 GB and DeepLoc2 ~3.6 GB — flat, they truncate long proteins;
- DeepTMHMM — computed per chunk from its longest protein (measured 3.7 GB at 2,500
  residues, 5.6 GB at 10,000, 11.7 GB at 35,000; ~4.1 GB for a typical plant chunk), so a
  chunk with a titin-sized protein runs almost alone on the card;
- DeepCoil2 — the whole GPU: its memory grows with both protein length and chunk size
  (~10 GB for 300 proteins, ~19 GB for ~5,000), so no fixed estimate is safe.

On a 20 GB card that allows up to 3 DeepTMHMM, 4 DeepLoc2 or 6 SignalP6 tasks at once, or
a mix. These tools leave the GPU mostly idle while loading models and doing CPU-side work
(21–44% average use when running alone), which is where the speed-up comes from. The same
chunk run alone and while sharing gave identical results for all four tools. Extremely long
proteins can still exceed a small card with a single tool (DeepTMHMM ~11.7 GB at 35,000
residues).

`maxForks` also limits how many tasks per tool wait for the GPU at once, so thousands of
chunks don't sit reserving CPUs and memory while they wait. Lock files live in
`${TMPDIR:-/tmp}/rgaprofiler-gpu-locks-<uid>` (change with the `RGAPROFILER_GPU_LOCK_DIR`
environment variable; it must be on a local disk, not NFS). If locking is impossible, the
task warns and runs without the lock. Each GPU task's `.gpu_lock.log` records what it held
and how long it waited (that wait counts in the trace's `realtime`). On multi-GPU machines
each task is pinned to the GPU it holds (`CUDA_VISIBLE_DEVICES` is honoured).

</details>

## Running on another organism, or with different classification parameters

This pipeline is **organism-agnostic**: it takes a protein FASTA in, not a genome or a
species-specific reference. There is no `--organism`/`--species`/`--genome` parameter to
set, and (unlike genome-based nf-core pipelines) nothing in `nextflow_schema.json` needs
touching to point the pipeline at a different species.

- **`FASTA_QC` and all six prediction tools** (DeepCoil2, Phobius, InterProScan, DeepLoc2,
  SignalP6, DeepTMHMM) are general-purpose protein predictors — they run identically
  regardless of which organism the input proteins came from. To run on a different
  organism, list its protein FASTA(s) in the samplesheet (see
  [Samplesheet input](#samplesheet-input) above) and run the pipeline exactly as
  documented — nothing else changes.
- **`--interproscan_db`** is InterProScan's member-database data (Pfam, PANTHER, Gene3D,
  …). It is not organism- or clade-specific — the same install works for every taxon
  InterProScan supports — so it never needs to change per organism either.

The **one** organism-specific piece of the whole pipeline is the RGA classification
ruleset itself: which InterProScan accessions count as which feature (NB-ARC, TIR, LRR,
…), the coiled-coil score threshold, the consensus policies, and the ordered
family/subclass rules built on top of them. All of it lives in one file, vendored from
[`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor):

```
docker/rga_classify/src/code/rgas/config/rga_config.yaml
```

That file's own header states the intent directly: _"Everything that is organism-,
database- or threshold-specific lives here. The Python code contains NO accessions and
NO magic numbers."_ — `rgas_prediction.py` (the script `RGA_CLASSIFY` runs) never hard-codes
a domain accession or a cut-off; it only reads this file. Every default in it was
calibrated against the pipeline's own reference dataset (the R570 sugarcane proteome),
so it works out of the box for plant RGA calling in general, but individual values may
need revisiting for a different plant lineage — or for a non-plant genome. `rga_config.yaml`
is heavily commented with the rationale for every threshold and rule; skim it before
changing anything, since several of the choices documented there (e.g. the CC threshold
below) are deliberate judgement calls, not tool-recommended defaults.

> [!NOTE]
> The file's header also points to `docs/rga/README.md` for a fuller "how to adapt to
> another organism" walkthrough. That document is part of the upstream
> [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor) project, not of this
> pipeline's own vendored copy or its `docs/` folder — it is not available here. The
> sections below are this pipeline's own self-contained equivalent.

### What to check when adapting to a new organism

Open `rga_config.yaml` and look at, in likely order of relevance:

- **`interproscan_features`** — the Pfam/InterPro/SMART/PROSITE/PRINTS/Gene3D accessions
  that count as evidence for each feature (`NB-ARC`, `TIR`, `RPW8`, `LRR`, `STTK`, `LysM`,
  `CC`). These are curated domain models, not R570-specific calls, so most of them apply
  to any plant proteome unchanged. The file gives a worked example of when this is
  **not** true: NACHT (`PF05729`/`IPR007111`) is currently kept out of `NB-ARC` and only
  reported via `watch_accessions`, because it never fires in sugarcane — but NACHT NLRs
  are the norm in fungal and animal genomes, so classifying those would mean moving it
  into `interproscan_features.NB-ARC` first.
- **`coiled_coil.threshold`** (default `0.5`) — explicitly documented as "a deliberate
  midpoint choice", not a value DeepCoil2 itself recommends. If you have an
  organism-appropriate reference (the file uses the Rx N-terminal domain, `PF18052`, as
  its own precision/recall check), re-tuning this threshold against it is worth doing
  before trusting CNL/CN/TM-CC calls on a new species.
- **`rules`** — the ordered, mutually-exclusive family/subclass rules (`CNL`, `TNL`,
  `RNL`, `LRR-RLK`, …). These encode standard NLR/RLK/RLP nomenclature and are reusable
  across plants as-is, but two spots are called out in-file as organism/method choices
  rather than fixed logic: the `other-RLK` rule (priority 13) ships **commented out**
  because the pipeline's sugarcane-specific default (Rody et al. 2019) excludes it —
  uncomment it for RGAugury-style scope instead; and the `any_of: [[TM, SP]]` extension
  on the `LRR-RLP`/`LysM-RLP` rules is documented as an addition beyond the published
  rule, with the one-line change to revert it noted right above the rules.
- **`ids.locus_regex`** — matches R570's own protein-ID convention
  (`SoffiXsponR570.7os1g018900.1.p` → locus `SoffiXsponR570.7os1g018900`) to build the
  locus-level summary (`rga_predictions_by_locus.tsv`). This **will** need updating to
  match your organism's own annotation ID format, or locus-level collapsing will be
  silently wrong (or empty) for anything but R570-style IDs.
- **`ectodomain_features`** and the `deeploc` localisation lists — extend these if your
  organism's receptors rely on ectodomains or subcellular compartments the defaults
  don't cover.

### How to apply a change

Three options, in increasing order of how much of the ruleset you're changing:

**1. Tune a handful of values, no file edit at all.** Everything under
`consensus_group` in `rgas_prediction.py --help` is a plain CLI flag with no file/path
involved, so it can be passed straight through via `task.ext.args` in a custom config
(`-c custom.config`, never via `-params-file`/`--flag`, see the
[warning above](#running-the-pipeline)):

```groovy title="custom.config"
process {
    withName: 'RGA_CLASSIFY' {
        ext.args = '--cc-threshold 0.6 --tm-policy union --min-lrr-copies 2'
    }
}
```

```bash
nextflow run . -profile docker --input samplesheet.csv \
    --interproscan_db /path/to/interproscan-5.XX-YY.0 \
    --outdir results -c custom.config
```

Available flags: `--tm-policy`, `--sp-policy`, `--cc-policy`, `--cc-threshold`,
`--cc-min-length`, `--cc-max-gap`, `--min-lrr-copies`, `--rga-only`/`--keep-non-rga`.

**2. Point at your own full `rga_config.yaml`, no image rebuild.** Copy
`docker/rga_classify/src/code/rgas/config/rga_config.yaml`, edit your copy, and pass
`--config` — but unlike the paths this pipeline resolves through its own Nextflow
`path` process inputs (e.g. `--interproscan_db`, `--softwares_dir`), a raw host path
inside `ext.args` is **not** automatically bind-mounted into the container, so add it to
`containerOptions` explicitly:

```groovy title="custom.config"
process {
    withName: 'RGA_CLASSIFY' {
        ext.args = '--config /abs/path/to/my_organism_config.yaml'
        containerOptions = '-v /abs/path/to:/abs/path/to'
    }
}
```

**3. Make it the pipeline's own default.** Edit
`docker/rga_classify/src/code/rgas/config/rga_config.yaml` directly in this repo, then
build the image under a new tag
(`docker build -t <registry>/rga_classify:<new_tag> -f docker/rga_classify/Dockerfile docker/rga_classify`),
push it to a registry your machines can pull from, and point the `container` line of both
`modules/local/rga_classify/main.nf` and `modules/local/rga_report/main.nf` (the report
reuses the same image) at the new tag. This is the right choice
once a new organism's config is settled and you want every future run to use it without
passing `-c`/`ext.args` at all.

### Different pipeline-level parameters

The above is about the RGA *classification* ruleset. For tuning the wrapper pipeline
itself instead — chunking (`--num_blocks`/`--fasta_qc_chunk_size`), GPU usage
(`--use_gpu`), or where license-gated software lives (`--softwares_dir`) — see
[Sequence batching](#sequence-batching---num_blocks) above and the full parameter
reference generated from `nextflow_schema.json`; none of those are organism-specific
either.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> **This pipeline is built for and tested with `-profile docker`.** Every tool this pipeline
> wraps (DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6, DeepTMHMM, `rga_classify`) ships
> as a public, pre-built Docker image on GHCR (`ghcr.io/omatheuspimenta/...`) — Docker pulls
> these automatically, so there is nothing to build yourself. Conda is **not** supported for
> this pipeline: only `FASTA_QC` declares a conda environment, so `-profile conda` will fail on
> every other step. Singularity/Podman/Apptainer/Charliecloud can in principle pull the same
> public images (they all understand plain Docker/OCI images), but that path has not been
> validated for this pipeline — Docker is the one to reach for.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

> [!NOTE]
> The list below is every profile actually defined in this pipeline's `nextflow.config` —
> checked against the file directly, not copied from the generic nf-core template. If you
> spot one here that doesn't match, `grep -A2 "^profiles {" nextflow.config` is the source
> of truth.

- `test`
  - A profile with a complete configuration for automated testing against the small (4-sequence) bundled dataset
  - Includes links to test data so needs no other parameters (except `--interproscan_db`, always required)
- `test_full`
  - Like `test`, but points at the real, full-scale R570 proteome test dataset instead of the small subsample. Pair with `long_running` (see below) — it's too large for the default resource/time budgets.
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above.
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `mamba`
  - Like `conda`, but installs the environment with [Mamba](https://mamba.readthedocs.io/) instead of Conda's own resolver.
- `debug`
  - Developer-oriented profile: dumps environment hashes, keeps `work/` after the run (`cleanup = false`), and prints the hostname before each process. Not needed for normal use.
- `arm64`
  - Targets an ARM64 host by building/pulling ARM64 container variants via [Wave](https://seqera.io/wave/) rather than this pipeline's own (x86-64) published images.
- `emulate_amd64`
  - The reverse case: forces Docker to run this pipeline's normal (x86-64) images under QEMU emulation on a non-x86-64 host (`--platform=linux/amd64`), instead of `arm64`'s native-build approach.
- `gpu`
  - This pipeline-specific profile (not part of the standard nf-core set) is a convenience shortcut for `--use_gpu true` — equivalent to passing that flag directly. Prefer `--use_gpu auto|true|false` (default `auto`; see the main [`README.md`](../README.md)) for anything beyond a quick default-on toggle.
- `long_running`
  - This pipeline-specific profile (not part of the standard nf-core set) raises every process's time budget to 1–10 days, sized for a real, full-scale proteome rather than the small test dataset. Slow steps (InterProScan, SignalP6, …) can otherwise be killed for running "too long" on real data. Add it alongside your other profiles, e.g. `-profile docker,long_running` — no config file editing needed. See `conf/long_running.config` if you need to raise the numbers even further.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

The pipeline sizes its own requests (see [Resource sizing](#resource-sizing-and-gpu-sharing)),
so you rarely need to change them. If a task fails because it ran out of memory or time,
it is automatically retried once with double the memory and time; if it fails again, the
run stops and tells you which task failed. To cap everything, use
[`resourceLimits`](#limiting-what-the-pipeline-may-use).

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

Unlike most nf-core pipelines, this one does **not** pull from
[biocontainers](https://biocontainers.pro/)/[bioconda](https://bioconda.github.io/) —
none of DeepCoil2, Phobius, InterProScan, DeepLoc2, SignalP6, DeepTMHMM or `rga_classify`
are packaged there. Every module's `container` directive instead points at this
pipeline's own image on GHCR (`ghcr.io/omatheuspimenta/<tool>:<tag>`), built from the
Dockerfiles under `docker/`; conda is not a supported alternative either (see
[`-profile`](#-profile) above — only `FASTA_QC` declares a conda environment).

To use a different image tag for a tool, override that module's `container` directive
via a custom config (`-c custom.config`), e.g.
`process { withName: 'DEEPTMHMM' { container = 'ghcr.io/omatheuspimenta/deeptmhmm:my-tag' } }`.
To change what an image actually contains, edit the corresponding `docker/<tool>/Dockerfile`,
build it under a new tag (`docker build -t <registry>/<tool>:<new_tag> -f docker/<tool>/Dockerfile docker/<tool>`,
run from the pipeline root), push it, and point the module's `container` line at it.

### Custom Tool Arguments

Every module accepts extra CLI flags for its underlying tool via `task.ext.args`, set
per-process in a custom config (`-c custom.config`, not `-params-file`/`--flag` — see the
[warning above](#running-the-pipeline)):

```groovy title="custom.config"
process {
    withName: 'SIGNALP6' {
        ext.args = '--mode fast'
    }
}
```

See [`conf/modules.config`](../conf/modules.config) for this pipeline's own `ext.args`
usage, and the
["Running on another organism"](#running-on-another-organism-or-with-different-classification-parameters)
section above for a worked example on `RGA_CLASSIFY` specifically.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
