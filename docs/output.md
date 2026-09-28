# omatheuspimenta/rgaprofiler: Output

## Introduction

This page describes what the pipeline writes to your results directory (`--outdir`) and
how to read it. All paths below are relative to that directory.

**Start here:** the main result for each sample is
`rga/<sample>/rga_predictions.tsv` (one row per protein, with its RGA call), and the
quickest overview is `summary_report/<sample>/report.html` (open it in a browser).

## How the results are organised

Every tool has its own folder, and inside it **one folder per sample** (the `sample`
column of your samplesheet), so the samples of a multi-sample run never overwrite each
other:

```
<outdir>/
├── fasta/<sample>/             cleaned input FASTA (+ the chunks it was split into)
├── deepcoil2/<sample>/         coiled-coil predictions
├── phobius/<sample>/           signal peptides + transmembrane topology
├── interproscan/<sample>/      protein domains
├── deeploc2/<sample>/          subcellular localisation
├── signalp6/<sample>/          signal peptides
├── deeptmhmm/<sample>/         transmembrane topology
├── rga/<sample>/               ★ RGA classification — the main result
├── summary_report/<sample>/    ★ one-page HTML summary
└── pipeline_info/              run reports, parameters and software versions
```

## Pipeline overview

For each sample, the pipeline:

1. Cleans the protein FASTA and splits it into chunks ([FASTA_QC](#fasta_qc)).
2. Runs six prediction tools on every chunk and merges each tool's per-chunk results back
   into one result per sample: [DeepCoil2](#deepcoil2) (coiled coils),
   [Phobius](#phobius) (signal peptides + transmembrane topology),
   [InterProScan](#interproscan) (protein domains), [DeepLoc2](#deeploc2) (subcellular
   localisation), [SignalP6](#signalp6) (signal peptides) and [DeepTMHMM](#deeptmhmm)
   (transmembrane topology).
3. Combines all six into a per-protein RGA (Resistance Gene Analog) call
   ([RGA classification](#rga-classification)), using the classification logic of
   [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor).
4. Summarises the calls in a one-page HTML report ([Summary report](#summary-report)).

```
                              ┌─ DeepCoil2 ─────┐
                              ├─ Phobius ───────┤
input.fasta ──▶ FASTA_QC ──▶  ├─ InterProScan ──┼──▶ RGA classification ──▶ Summary report
  (per sample)  (clean, chunk)├─ DeepLoc2 ──────┤
                              ├─ SignalP6 ──────┤
                              └─ DeepTMHMM ─────┘
                               (each: one task per chunk, then merged per sample)
```

Chunking is invisible in the results: every protein is predicted exactly once and each
tool's files contain every protein of the sample. (Which proteins share a chunk can nudge
a few tools' numbers very slightly — see
[Chunking and reproducibility](usage.md#chunking-and-reproducibility).)

### FASTA_QC

<details markdown="1">
<summary>Output files</summary>

- `fasta/<sample>/`
  - `<sample>_clean.fasta`: your input after removing duplicate protein IDs, trailing
    stop codons (`*`) and stray characters, and uppercasing. This is what every tool
    actually analyses.
  - `<sample>_clean_chunks/<sample>_clean.part_NNN.fasta.fasta`: the same sequences split
    into chunks (whole FASTA records only). How many chunks is set by `--num_blocks` or
    `--fasta_qc_chunk_size` — see [Sequence batching](usage.md#sequence-batching---num_blocks).

</details>

Every sequence that was changed or dropped is listed as a `WARNING` in the FASTA_QC task's log (`.command.err` in its work directory).

### DeepCoil2

<details markdown="1">
<summary>Output files</summary>

- `deepcoil2/<sample>/`
  - One `<protein_id>.out` per protein (the ID with punctuation removed): a per-residue
    table with columns `aa`, `cc`, `raw_cc`, `prob_a`, `prob_d` — the coiled-coil
    probability and the heptad-position (`a`/`d`) probabilities at every residue.
    Proteins shorter than 20 residues are skipped by DeepCoil2 and have no file.

</details>

[DeepCoil2](https://github.com/labstructbioinf/DeepCoil) predicts coiled-coil domains.
Runs on GPU when available.

### Phobius

<details markdown="1">
<summary>Output files</summary>

- `phobius/<sample>/<sample>_phobius.tsv`: Phobius short format — one row per protein, in
  the order of your input, with the number of transmembrane helices (`TM`), whether a
  signal peptide was found (`SP`) and the predicted topology string.

</details>

[Phobius](https://software.sbc.su.se/cgi-bin/request.cgi?project=phobius) predicts signal
peptides and transmembrane topology together. CPU only.

### InterProScan

<details markdown="1">
<summary>Output files</summary>

- `interproscan/<sample>/<sample>_interpro.tsv`: the standard InterProScan TSV (no header
  row) — one row per domain/site hit, from every member database (Pfam, PANTHER, Gene3D,
  PROSITE, CDD, …). This is where the NB-ARC, LRR, TIR, RPW8 and coiled-coil domain
  evidence for the RGA calls comes from. Row order is not meaningful (InterProScan itself
  writes rows in a different order from run to run).

</details>

[InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) annotates protein
domains and functional sites. CPU only. Needs its database downloaded once
(`--interproscan_db`, see [software setup](software-setup.md)). InterProScan's own
licensed Phobius/SignalP 4.1/TMHMM 2.0c analyses are not needed: those signals come from
this pipeline's dedicated Phobius, SignalP6 and DeepTMHMM steps.

### DeepLoc2

<details markdown="1">
<summary>Output files</summary>

- `deeploc2/<sample>/<sample>_deeploc2.csv`: one row per protein — predicted subcellular
  localisation(s), the probability of each compartment, sorting signals and membrane type.

</details>

[DeepLoc2](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) predicts subcellular
localisation (using its "Fast" model). Runs on GPU when available. Localisation never
decides an RGA class; it only adjusts the reported confidence.

### SignalP6

<details markdown="1">
<summary>Output files</summary>

- `signalp6/<sample>/`
  - `<sample>_signalp6_predictions.txt`: one row per protein — predicted signal-peptide
    type (Sec/SPI, Sec/SPII, Tat/SPI, …, or `OTHER` for none) with the probability of each
    type and the cleavage site.
  - `<sample>_signalp6.gff3`, `region_output.gff3`: the predicted signal peptides (and
    their regions) in GFF3 format.
  - `processed_entries.fasta`: the sequences SignalP6 scored.
  - `chunk_N_output.json`: SignalP6's full raw output, one file per chunk.

</details>

[SignalP 6.0](https://github.com/fteufel/signalp-6.0) predicts all five types of signal
peptide. Runs on GPU when available (with a GPU-converted copy of its weights, see
[software setup](software-setup.md)).

### DeepTMHMM

<details markdown="1">
<summary>Output files</summary>

- `deeptmhmm/<sample>/`
  - `<sample>_deeptmhmm.gff3`: per-protein topology regions (signal peptide, inside,
    outside, transmembrane helix, beta strand).
  - `<sample>_predicted_topologies.3line`: the same calls in DeepTMHMM's three-lines-per-
    protein format (header with the predicted type, sequence, topology string).
  - `embeddings/`, `probabilities/`: DeepTMHMM's intermediate per-protein files (named by
    a hash of the sequence), kept for completeness.
  - `summaries/chunk_N_deeptmhmm_results.md`: DeepTMHMM's short run summary, one per chunk.

</details>

[DeepTMHMM](https://dtu.biolib.com/DeepTMHMM) predicts alpha-helical and beta-barrel
transmembrane topology. Runs on GPU when available.

### RGA classification

<details markdown="1">
<summary>Output files</summary>

- `rga/<sample>/`
  - `rga_predictions.tsv`: **the main result** — every protein, one row each: the
    evidence gathered from all six tools (`sp_signalp`, `sp_phobius`, `n_tm_phobius`,
    `n_tm_deeptmhmm`, `cc_deepcoil`, `cc_coils`, `cc_rx_domain`,
    `predicted_localization`, `features_found`, …) and the call: `is_rga`, `rga_family`
    (e.g. `NLR`) and `rga_subclass` (e.g. `CNL`).
  - `rga_predictions_rga_only.tsv`: the same table, RGA candidates only.
  - `rga_predictions_by_locus.tsv`: one row per gene locus (collapses the transcript
    models of the same gene) — use this for gene-level counts in polyploid genomes.
  - `rga_summary_counts.tsv`: how many proteins fall in each RGA family/subclass.
  - `rga_domain_evidence_long.tsv`: every individual piece of evidence behind the calls,
    one row per hit.
  - `report.html` / `report.md`: the classifier's own detailed report (methods, rules
    applied, confidence, warnings, top candidates).
  - `run_metadata.json`: exact command, settings and input checksums, for reproducibility.
  - `accession_audit.tsv`, `unmatched_ids_report.tsv`: how protein IDs were matched across
    the six tools' files — check these if counts look off.
  - `cc_policy_sensitivity.tsv`, `cc_segment_sensitivity.tsv`: how the coiled-coil-based
    classes would change under the other coiled-coil settings.
  - `cache/`, `logs/`: the classifier's cache and log.

</details>

The classification logic comes from [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor)
(see [`CITATIONS.md`](../CITATIONS.md)). Its rules and thresholds are configurable — see
[Running on another organism](usage.md#running-on-another-organism-or-with-different-classification-parameters).

### Summary report

<details markdown="1">
<summary>Output files</summary>

- `summary_report/<sample>/report.html`: a single self-contained page (open it in any
  browser, no internet needed) with the number of proteins and RGA candidates, tables of
  RGA families and subclasses, how many proteins each tool found evidence for, links to
  this sample's detailed results in the other folders, and the software versions used.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`, `execution_timeline_*.html`, `pipeline_dag_*.html`:
    Nextflow's reports on the run (run time, CPU and memory used by every task).
  - `params_*.json`: every parameter the run used — keep it to repeat a run exactly.
  - `rgaprofiler_software_versions.yml`: the version of every tool that ran.

</details>

These are generated by [Nextflow](https://docs.seqera.io/platform-cloud/reports/overview)
and are the first place to look when troubleshooting a run or checking how long each step
took.
