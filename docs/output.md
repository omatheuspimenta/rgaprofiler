# omatheuspimenta/rgaprofiler: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below are created in the results directory (`--outdir`) after the pipeline has finished. All paths are relative to the top-level results directory.

Every prediction-tool directory is named after the process's first underscore-token (e.g. `RGA_CLASSIFY` publishes to `rga/`, `RGA_REPORT` to `summary_report/`) — a pre-existing, shared `publishDir` convention (`conf/modules.config`), not something specific to any one tool.

## Pipeline overview

The pipeline takes a protein FASTA per sample and:

1. Cleans/deduplicates it and splits it into sequence blocks/chunks ([FASTA_QC](#fasta_qc)).
2. Runs six independent protein-prediction tools against the cleaned FASTA: [DeepCoil2](#deepcoil2) (coiled-coil domains), [Phobius](#phobius) (signal peptides + TM topology), [InterProScan](#interproscan) (protein domains/functional annotation), [DeepLoc2](#deeploc2) (subcellular localization), [SignalP6](#signalp6) (signal peptides), [DeepTMHMM](#deeptmhmm) (transmembrane helices). Five of these six — every one except Phobius — run once per FASTA_QC chunk rather than once on the whole proteome, then have their own `*_MERGE` step reassemble the per-chunk outputs into one result per sample (see [FASTA_QC](#fasta_qc) below).
3. Combines all six tools' outputs into per-protein RGA (Resistance Gene Analog) classifications ([RGA classification](#rga-classification)), using the classification logic from [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor).
4. Renders a self-contained HTML summary of those classifications ([Summary report](#summary-report)).

```
                              ┌─ DeepCoil2 ────▶ DeepCoil2_MERGE ─────┐
                              ├─ Phobius ─────────────────────────────┤
input.fasta ──▶ FASTA_QC ──▶ ├─ InterProScan ─▶ InterProScan_MERGE ──┼──▶ RGA_CLASSIFY ──▶ RGA_REPORT
                              ├─ DeepLoc2 ─────▶ DeepLoc2_MERGE ──────┤
                              ├─ SignalP6 ─────▶ SignalP6_MERGE ──────┤
                              └─ DeepTMHMM ────▶ DeepTMHMM_MERGE ─────┘
```

### FASTA_QC

<details markdown="1">
<summary>Output files</summary>

- `fasta/`
  - `<sample>_clean.fasta`: the input FASTA after deduplication, stripping the trailing stop codon (`*`) some proteomes ship, and uppercasing — this cleaned file is what every downstream tool actually receives (directly, for Phobius; split into chunks, for the other five tools below).
  - `<sample>_clean_chunks/<sample>_clean.part_NNN.fasta.fasta`: the same cleaned sequences split into chunks — complete FASTA records only, never an arbitrary line split. DeepCoil2, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM each run once per chunk (their own `*_MERGE` process then reassembles the per-chunk results into one file/directory per sample, published under each tool's own output directory below — chunking is otherwise invisible downstream).

  Chunk count/size is controlled by one of two mutually exclusive parameters:
  - `--num_blocks <N>` (e.g. `--num_blocks 1000`): split into (up to) `N` chunks, balanced as evenly as possible by sequence count. Not hard-coded — set as high as needed for a very large proteome (this is what keeps DeepCoil2 in particular from ever having to process an entire proteome as a single task). A higher value gives Nextflow's executor more independent, smaller tasks to schedule in parallel; it does not itself force that many tasks to run concurrently — that remains governed by your `-profile`/executor/resource configuration.
  - `--fasta_qc_chunk_size <N>` (default `5000`, used only when `--num_blocks` is unset): a fixed number of sequences per chunk instead, so the chunk *count* scales with input size.

</details>

### DeepCoil2

<details markdown="1">
<summary>Output files</summary>

- `deepcoil2/results/`
  - One `<protein_id>.out` per input sequence: a per-residue TSV (`aa`, `cc`, `raw_cc`, `prob_a`, `prob_d`) giving the predicted coiled-coil probability and heptad-register (`a`/`d` core position) probabilities at every residue.

</details>

[DeepCoil2](https://github.com/labstructbioinf/DeepCoil) predicts coiled-coil domains. GPU-capable (`--use_gpu`, see [usage docs](../README.md#usage)). Runs once per FASTA_QC chunk (`--num_blocks`/`--fasta_qc_chunk_size`, see [FASTA_QC](#fasta_qc)) — every input sequence appears in exactly one chunk's `.out` file, so `DeepCoil2_MERGE` just collects them into the single `results/` directory published here.

### Phobius

<details markdown="1">
<summary>Output files</summary>

- `phobius/<sample>_phobius.tsv`: short-format Phobius output — one row per protein with predicted transmembrane helix count and signal-peptide call.

</details>

[Phobius](https://software.sbc.su.se/cgi-bin/request.cgi?project=phobius) predicts signal peptides and transmembrane topology jointly. CPU-only.

### InterProScan

<details markdown="1">
<summary>Output files</summary>

- `interproscan/<sample>_interpro.tsv`: standard 14/15-column InterProScan TSV — one row per domain/site hit per protein, across every member database InterProScan runs (Pfam, PANTHER, Gene3D, PROSITE, HAMAP, CDD, …). This is the primary source of the NB-ARC, LRR, TIR, RPW8, and coiled-coil domain calls the RGA classification step relies on.

</details>

[InterProScan](https://www.ebi.ac.uk/interpro/about/interproscan/) does protein domain/functional-site annotation. CPU-only. Requires a pre-downloaded database (`--interproscan_db`, see [`docs/software-setup.md`](software-setup.md)); by design this pipeline's committed reference runs did **not** enable the licensed Phobius/SignalP-4.1/TMHMM-2.0c sub-analyses within InterProScan itself — those signals come from this pipeline's own dedicated Phobius/SignalP6/DeepTMHMM modules instead. Like DeepCoil2, runs once per FASTA_QC chunk; `InterProScan_MERGE` concatenates the per-chunk TSVs (no header row) into the single file published here.

### DeepLoc2

<details markdown="1">
<summary>Output files</summary>

- `deeploc2/results/<sample>_deeploc2.csv`: one row per protein — predicted subcellular localization(s), per-class probabilities, predicted sorting signals, and membrane-protein type.

</details>

[DeepLoc2](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) predicts subcellular localization. GPU-capable (`--use_gpu`); runs the "Fast" model by default. Runs once per FASTA_QC chunk; `DeepLoc2_MERGE` keeps the first chunk's CSV header and concatenates every chunk's data rows into the single file published here.

### SignalP6

<details markdown="1">
<summary>Output files</summary>

- `signalp6/results/`
  - `<sample>_signalp6_predictions.txt`: per-protein predicted signal-peptide type (Sec/SPI, Sec/SPII, Tat/SPI, …) and cleavage-site position/probability.
  - `<sample>_signalp6.gff3` / `region_output.gff3`: the same calls in GFF3 form.
  - `chunk_N_output.json`: full per-residue probability output, one file per FASTA_QC chunk (a per-chunk JSON object, not a per-protein list, so it has no lossless line-level merge across chunks — every chunk's copy is kept individually rather than dropped or naively concatenated into invalid JSON).
  - `processed_entries.fasta`: the exact sequences SignalP6 scored (post its own internal filtering), across every chunk.

</details>

[SignalP6](https://github.com/fteufel/signalp-6.0) predicts signal peptides (all five types). Runs in `slow-sequential` mode by default (the only weight set this pipeline's reference `softwares/SignalP6/` install ships). GPU-capability is a property of which weight files are staged, not a CLI flag — see [`docs/software-setup.md`](software-setup.md). Runs once per FASTA_QC chunk; `SignalP6_MERGE` reassembles `_predictions.txt`/`.gff3`/`region_output.gff3` (keeping one shared header, then every chunk's data rows/blocks) and concatenates `processed_entries.fasta`, into the files published here.

### DeepTMHMM

<details markdown="1">
<summary>Output files</summary>

- `deeptmhmm/results/`
  - `<sample>_deeptmhmm.gff3`: per-protein predicted region boundaries (signal peptide / inside / outside / transmembrane helix / beta-barrel strand).
  - `<sample>_predicted_topologies.3line`: the same topology calls in DeepTMHMM's compact three-line-per-protein format.
  - `summaries/chunk_N_deeptmhmm_results.md`: a short run summary, one per FASTA_QC chunk.
  - `embeddings/`, `probabilities/`: intermediate per-protein ESM1b embeddings and per-residue class probabilities, collected across every chunk.

</details>

[DeepTMHMM](https://dtu.biolib.com/DeepTMHMM) predicts alpha/beta transmembrane topology. GPU-capable — auto-detects `torch.cuda.is_available()` with no CLI flag needed, driven by the same `--use_gpu` setting as the other GPU-capable tools. Runs once per FASTA_QC chunk; `DeepTMHMM_MERGE` keeps one shared GFF3 header, concatenates the per-chunk blocks/records, and collects the embeddings/probabilities into the files published here.

### RGA classification

<details markdown="1">
<summary>Output files</summary>

- `rga/rga_out/`
  - `rga_predictions.tsv`: the main per-protein result — every input protein, one row each, with the harmonised evidence columns pulled from all six upstream tools (`sp_signalp`, `sp_phobius`, `n_tm_phobius`, `n_tm_deeptmhmm`, `cc_deepcoil`, `cc_coils`, `cc_rx_domain`, `predicted_localization`, `features_found`, …) plus the final call: `is_rga`, `rga_family` (e.g. `NLR`), `rga_subclass` (e.g. `CNL`).
  - `rga_predictions_rga_only.tsv` / `rga_predictions_by_locus.tsv`: the same predictions filtered to RGA-positive proteins only, and collapsed to one row per locus (multiple transcript models of the same gene).
  - `rga_summary_counts.tsv`: aggregate counts per RGA family/subclass.
  - `rga_domain_evidence_long.tsv`: the individual domain hits behind each classification, one row per hit (long format).
  - `accession_audit.tsv`, `unmatched_ids_report.tsv`: bookkeeping on how protein IDs were parsed/matched across the six input files.
  - `cc_policy_sensitivity.tsv`, `cc_segment_sensitivity.tsv`: sensitivity of the coiled-coil call to the classifier's internal thresholding policy.
  - `report.html` / `report.md`, `run_metadata.json`: `rga_classify`'s own run report and metadata (distinct from this pipeline's own [summary report](#summary-report) below, which is generated separately from these same TSVs).
  - `cache/`, `logs/`: intermediate parser cache and the run log.

</details>

Ports the RGA classification logic from [`rgapredictor`](https://github.com/omatheuspimenta/rgapredictor) (vendored unmodified, see [`CITATIONS.md`](../CITATIONS.md)). CPU-only.

### Summary report

<details markdown="1">
<summary>Output files</summary>

- `summary_report/report_out/report.html`: a single self-contained HTML page (inline CSS, no external assets) built directly from `rga_classify`'s `rga_predictions.tsv`/`rga_summary_counts.tsv` — protein/RGA-candidate counts, an RGA family table, an RGA subclass table, a per-tool evidence-contribution table, links to every tool's detailed output, and a software-versions table.

</details>

This pipeline does not use MultiQC (see [`README.md`](../README.md)); this lightweight report is the summary output instead.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_*.html`, `execution_timeline_*.html`, `pipeline_dag_*.html`.
  - `rgaprofiler_software_versions.yml`: collated tool/software versions for every process that ran, including `RGA_REPORT` itself.
  - `params_*.json`: the parameters used for the run.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
