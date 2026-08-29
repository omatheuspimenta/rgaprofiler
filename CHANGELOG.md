# omatheuspimenta/rgaprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of omatheuspimenta/rgaprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--num_blocks`: split each sample's input FASTA into a fixed number of sequence blocks (instead of a fixed sequences-per-chunk size via `--fasta_qc_chunk_size`), which DeepCoil2, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM now all run once per chunk over (previously only InterProScan was chunked) — each with its own `*_MERGE` process reassembling the per-chunk outputs back into one result per sample. Prevents DeepCoil2 in particular from ever being forced to process an entire proteome as a single task.

### `Fixed`

### `Dependencies`

### `Deprecated`
