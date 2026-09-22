# omatheuspimenta/rgaprofiler: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of omatheuspimenta/rgaprofiler, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- `--num_blocks`: split each sample's input FASTA into a fixed number of sequence blocks (instead of a fixed sequences-per-chunk size via `--fasta_qc_chunk_size`), which DeepCoil2, InterProScan, DeepLoc2, SignalP6 and DeepTMHMM now all run once per chunk over (previously only InterProScan was chunked) — each with its own `*_MERGE` process reassembling the per-chunk outputs back into one result per sample. Prevents DeepCoil2 in particular from ever being forced to process an entire proteome as a single task.
- `--gpu_concurrency`: how many GPU tasks may share one GPU at once on the `local` executor (default: derived from the GPU's VRAM, which is 1 — strict serialization — on cards up to 40GB or when VRAM can't be determined; `0` disables the lock and throttle). See `docs/usage.md`, "Resource sizing and GPU sharing".
- `bin/detect_host_resources.sh`: reports a host's CPUs, RAM (cgroup-aware), GPU count and VRAM. Used lazily per task (`task.ext.host` in `conf/base.config`, memoized for the run) to size chunk-task requests, never above what the host can provide, and only on the `local` executor. Explicit `process.resourceLimits` still take precedence.
- `bin/gpu_lock.sh`: a real host-level `flock` semaphore for GPU tasks (see `Fixed`). Each GPU task's `.gpu_lock.log` records the slot it took and its wait.
- New labels `process_medium_chunk` (DeepCoil2, DeepLoc2, SignalP6, DeepTMHMM) and `process_high_chunk` (InterProScan), sized for one chunk rather than a whole proteome.
- Tests: `tests/bin/test_host_resources.sh`, `test_gpu_lock.sh` and `test_resource_sizing.sh` (mocked 256-CPU/503GB, 2-CPU/8GB and GPU-less hosts; lock mutual exclusion, slot count, kill-release and fail-open) and `tests/resources.nf.test` (chunk labels under `-profile test`'s CI-sized `resourceLimits`).

### `Fixed`

- CLI-provided integer/boolean parameters (e.g. `--num_blocks 1000`) failing pipeline parameter validation (`Value is [string] but should be [integer]`) under Nextflow 26.04's v2 syntax parser, which always parses CLI flags as strings. `validation.lenientMode` (a prior attempt at this fix) does not perform this cast and was a red herring; the actual fix is nf-schema's `cast_cli_params` option, which requires `nf-schema>=2.7.2` (bumped from `2.5.1`) and is now enabled via `cli_typecast: true` in `utils_nfcore_rgaprofiler_pipeline`'s call to `UTILS_NFSCHEMA_PLUGIN`. See [nf-core/blog: Why parameters are strings all of a sudden](https://nf-co.re/blog/2026/parameter-types). A params file (`-params-file`, see `assets/params.example.yml`) remains the more robust option for typed parameters regardless, since the cast only affects validation, not the runtime `params` values themselves.
- **Chunked tools starving each other at high `--num_blocks`** (observed on a full-scale proteome with `-profile docker,long_running --num_blocks 1000` on a 256-CPU/503GB/20GB-GPU host: InterProScan kept progressing while DeepCoil2, DeepLoc2, SignalP6 and DeepTMHMM stalled). Two compounding causes:
  - **Memory sized for the wrong execution model.** `long_running`'s `process_medium` (150GB) was calibrated for these four tools running once against a whole proteome, but they run once per chunk. Measured peak RSS per chunk task on that run was 3.8–10.5GB (InterProScan 8.8GB, which was likewise reserved at 100GB), so only ~3 `process_medium` tasks and ~5 InterProScan tasks ever fit in RAM, competing for the same 503GB. The five chunked tools now carry their own labels whose memory is a floor plus a term proportional to the chunk's FASTA size (times `task.attempt`), capped at 90% of the host on the `local` executor (roughly 14–15GB per GPU-tool chunk and 11–12GB per InterProScan chunk at ~300 sequences, versus 150GB/100GB before); `long_running` only raises their time ceiling. `process_medium` under `long_running` (now used only by `RGA_CLASSIFY`, measured ~2.8GB) drops from 150GB to 16GB.
  - **GPU exclusivity faked through CPU admission control.** `process_gpu` reserved 75% of the host's CPUs per task so Nextflow would run one at a time; InterProScan's continuous stream of small CPU tasks kept refilling freed capacity, so such a task could starve indefinitely while the GPU sat idle, and it ignored VRAM entirely. Replaced by a real host-level lock (`bin/gpu_lock.sh`, sourced from `beforeScript`, held for the whole task, released by the kernel if the task dies) plus a per-process `maxForks` cap so `--num_blocks 1000` doesn't park hundreds of idle tasks on the lock. `process_gpu` no longer inflates `cpus`.
- Measured GPU footprint, for reference (RTX A4500, one chunk per task, ~300 / ~4900 sequences): peak VRAM DeepCoil2 10.9 / 19.2GB (grows with chunk size — keep chunks around the default 5000 sequences or fewer on a 20GB card), DeepTMHMM 3.7 / 4.0GB, DeepLoc2 3.1 / 3.1GB, SignalP6 1.9 / 1.9GB; peak RSS DeepCoil2 6.4 / 16GB, DeepLoc2 10.3 / 10.3GB, DeepTMHMM 5.4 / 5.5GB, SignalP6 3.5 / 4.3GB, InterProScan 6.5GB / 7.6GB (at ~300 / ~1000 sequences).
- InterProScan chunk tasks now request 6 CPUs instead of 12. `-cpu` is only an upper bound on its threads: measured use on a real 1000-chunk run was median 3.3 cores (max 4.0), and one chunk ran in 16m24s at `-cpu 6` vs 16m45s at `-cpu 12`. Twice as many chunks can now be admitted at once.

### `Dependencies`

- Bumped the `nf-schema` plugin pin from `2.5.1` to `2.7.2` (see `Fixed`, above).

### `Deprecated`
