# omatheuspimenta/rgaprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of omatheuspimenta/rgaprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--num_blocks`: split each sample's input FASTA into a fixed number of sequence blocks (instead of a fixed sequences-per-chunk size via `--fasta_qc_chunk_size`), which DeepCoil2, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM now all run once per chunk over (previously only InterProScan was chunked) — each with its own `*_MERGE` process reassembling the per-chunk outputs back into one result per sample. Prevents DeepCoil2 in particular from ever being forced to process an entire proteome as a single task.

### `Fixed`

- CLI-provided integer/boolean parameters (e.g. `--num_blocks 1000`) failing pipeline parameter validation (`Value is [string] but should be [integer]`) under Nextflow 26.04's v2 syntax parser, which always parses CLI flags as strings. `validation.lenientMode` (a prior attempt at this fix) does not perform this cast and was a red herring; the actual fix is nf-schema's `cast_cli_params` option, which requires `nf-schema>=2.7.2` (bumped from `2.5.1`) and is now enabled via `cli_typecast: true` in `utils_nfcore_rgaprofiler_pipeline`'s call to `UTILS_NFSCHEMA_PLUGIN`. See [nf-core/blog: Why parameters are strings all of a sudden](https://nf-co.re/blog/2026/parameter-types). A params file (`-params-file`, see `assets/params.example.yml`) remains the more robust option for typed parameters regardless, since the cast only affects validation, not the runtime `params` values themselves.

### `Dependencies`

- Bumped the `nf-schema` plugin pin from `2.5.1` to `2.7.2` (see `Fixed`, above).

### `Deprecated`
