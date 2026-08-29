# Third-party software setup

Several tools this pipeline runs cannot legally be redistributed as part of the pipeline's Docker images — either because they carry an academic/non-commercial license (DeepTMHMM), or because their model weights/databases are distributed separately by their authors (SignalP 6.0, DeepLoc 2, InterProScan's licensed sub-analyses).

To handle this, every tool's Docker image is a **"baseline" image**: it bakes in only redistributable content (OS, runtime, non-restricted code/dependencies) and never bakes in license-gated weights, binaries, or databases. You supply those yourself, the pipeline validates they're present before running the corresponding step, and they're mounted into the container at runtime. Because the images themselves never embed restricted content, they're safe to publish and pull publicly — only the software _you_ place under `softwares/` is subject to each tool's own license terms.

`softwares/` is gitignored and untracked — nothing you place there is ever committed to this repository.

> [!IMPORTANT]
> All four sections below (InterProScan, DeepTMHMM, SignalP 6.0, DeepLoc 2) are **mandatory**, not pick-and-choose — the pipeline checks for all four before it runs anything at all, including the small bundled `-profile test` dataset. There is no flag to skip a tool. DeepCoil2 and Phobius need nothing from you (see the bottom of this page) — only the four above involve a manual download.

## Setup checklist

Work through these in order, from the pipeline's root directory (`cd rgaprofiler` after cloning — see the main [README](../README.md#1-get-the-pipeline)):

1. [ ] **InterProScan** — run `./setup_interproscandb.sh` (downloads + unpacks automatically). → [jump to section](#interproscan)
2. [ ] **DeepTMHMM** — request the Academic License release, place 8 files under `softwares/DeepTMHMM/DeepTMHMM-Academic-License-v1.0/`. → [jump to section](#deeptmhmm)
3. [ ] **SignalP 6.0** — download from DTU, place weights under `softwares/SignalP6/signalp-6-package/models/`. → [jump to section](#signalp-60)
4. [ ] **DeepLoc 2** — download from DTU, place weights under `softwares/DeepLoc2/DeepLoc2/models/` and the ESM1b encoder under `softwares/DeepLoc2/torch_cache/`. → [jump to section](#deeploc-2)
5. [ ] **Verify everything** with `bin/check_software_present.sh` (one command per tool — see [below](#how-the-check-works)).

## How the check works

Before running a process that needs license-gated software, the pipeline runs [`bin/check_software_present.sh`](../bin/check_software_present.sh) against the relevant directory. If anything required is missing, the pipeline fails immediately with a message listing exactly what's missing and how to get it — you'll never get partway through a run only to have a tool fail deep inside a container with a cryptic error.

You can also run the check yourself at any time, e.g.:

```bash
./bin/check_software_present.sh interproscan softwares/InterProScan/interproscan-5.78-109.0
./bin/check_software_present.sh deeptmhmm softwares
./bin/check_software_present.sh signalp6 softwares
./bin/check_software_present.sh deeploc2 softwares
```

By default, tools other than InterProScan look under `${params.softwares_dir}` (default `<pipeline root>/softwares`, overridable with `--softwares_dir`).

## Per-tool setup

### InterProScan

InterProScan is not bundled at all — its release + member-database data are tens of GB. Point `--interproscan_db` at your own installation.

**Step 1 — download and unpack.** From the pipeline root, run:

```bash
./setup_interproscandb.sh
```

This downloads InterProScan 5.78-109.0 directly from EBI (`https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz`, no account/registration needed), unpacks it under `databases/interproscan/interproscan-5.78-109.0`, and runs its own `setup.py` configuration step. Expect this to take a while (large download).

**Step 2 — note the path it prints.** The script ends with a line like:

```
Done! In Nextflow, use: --interproscan_db /full/path/to/databases/interproscan/interproscan-5.78-109.0
```

Use that exact path as `--interproscan_db` in every pipeline command from now on.

If you already have an InterProScan installation elsewhere (e.g. this repo's own `softwares/InterProScan/interproscan-5.78-109.0/`, from a prior manual run), point `--interproscan_db` at that instead — no need to re-download.

**Required**: `interproscan.sh` and the `data/` directory inside that install.

**Optional — licensed sub-analyses**: InterProScan ships placeholder files instead of the Phobius, SignalP 4.1, and TMHMM 2.0c binaries it can optionally call, because those are separately licensed. InterProScan runs fine without them (just skipping those specific analyses). To enable them, download each tool yourself and follow: https://github.com/ebi-pf-team/interproscan/wiki/ActivatingLicensedAnalyses

### DeepTMHMM

**License**: CC BY-NC-SA 4.0 (non-commercial, academic use only) — this is the most restrictive license among this pipeline's tools. Do not use DeepTMHMM outputs for commercial purposes without separately obtaining rights to do so.

**Step 1 — get the release.** Go to https://dtu.biolib.com/DeepTMHMM and request/obtain the "DeepTMHMM Academic License v1.0" release from its authors (DTU Health Tech / BioLib). If you have trouble locating it, BioLib's own licensing contact is `licensing@biolib.com` (from the license text bundled in the release itself). You'll receive an archive that unpacks to a folder named `DeepTMHMM-Academic-License-v1.0/` containing a `README.txt`, a `LICENSE.txt`, and the model files below.

**Step 2 — place the files.** Copy (or move) that whole folder so its contents end up directly under:

```
softwares/DeepTMHMM/DeepTMHMM-Academic-License-v1.0/
```

**Required files** at that path:

- `deeptmhmm_cv_0.model` … `deeptmhmm_cv_4.model` (5 cross-validation checkpoints)
- `esm_model_alphabet.pt`, `esm_model_args.pt`, `esm_model_state_dict.pt` (bundled ESM1b weights)

**Note on the CUDA stack**: upstream's own installation instructions pin `torch==1.5.0+cu92` (2020-era, CUDA 9.2), which does not run on modern GPUs (`RuntimeError: cublas runtime error` on Ampere/Ada cards). `docker/deeptmhmm/Dockerfile` uses a modern replacement stack (Python 3.10, a current PyTorch build) instead — validated end-to-end on this host's real GPU — with the small number of `torch.load()` calls on DeepTMHMM's inference path patched for PyTorch 2.6+'s `weights_only` default change. This is a packaging fix only; it doesn't change model behavior or outputs.

### SignalP 6.0

**License**: DTU Health Tech academic-use license; model weights are distributed separately from the pipeline.

SignalP 6.0 ships three alternative run modes (`fast`, `slow`, `slow-sequential`), each needing its own separately-downloaded weight file(s) — "your download only included the one you picked" (DTU's own words). This pipeline defaults to **`slow-sequential`**, since that's the mode whose weights (`sequential_models_signalp6/`) are actually present in this repo's own reference install; it takes ~6x longer than `fast` but needs no more RAM.

**Step 1 — download.** Go to https://services.healthtech.dtu.dk/services/SignalP-6.0/, accept DTU's academic license/terms on that page, and download the Linux package. DTU distributes each run mode's weights as separate downloads/variants of the package — make sure whichever one you get includes the **`slow-sequential`** weights, since that's the mode this pipeline uses by default (if in doubt, re-check the download page for a `slow_sequential`/`slow-sequential` labelled option, or contact DTU support if only `fast` is offered). You'll get a `signalp-6*.tar.gz`-style archive that unpacks to a `signalp-6-package/` folder containing (among other things) a `models/` subdirectory.

**Step 2 — place the files.** Copy that package's `models/` directory contents so you end up with:

```
softwares/SignalP6/signalp-6-package/models/sequential_models_signalp6/
```

i.e. the whole `signalp-6-package/` folder you downloaded should live under `softwares/SignalP6/`.

**Required**, under `softwares_dir/SignalP6/signalp-6-package/models/`:

- `sequential_models_signalp6/` (the `slow-sequential` mode's weights — this pipeline's default)

**Optional** (only needed if you override `--mode` via `task.ext.args` in `conf/modules.config`):

- `distilled_model_signalp6.pt` (`fast` mode)
- `ensemble_model_signalp6.pt` (`slow` mode)

**GPU note**: unlike DeepLoc2, SignalP6 has no `--device`/`-d` runtime flag — whether it uses a GPU is a property of the weight _files themselves_. The pipeline handles this automatically: whenever GPU is requested/detected (`--use_gpu auto`, the default, or `--use_gpu true`), `workflows/rgaprofiler.nf` looks for a **separate**, GPU-converted copy of the weights at `softwares_dir/SignalP6/signalp-6-package/models_gpu/` and points SignalP6 at that directory instead of the normal CPU one — you just need to have produced that directory once, ahead of time (SignalP6 has no way to convert its own weights at run time, and doing so on every run would be needlessly slow anyway). One-time setup, run from the signalp6 image so `signalp6_convert_models` and its Python dependencies are available (needs the same real GPU + NVIDIA Container Toolkit as running the pipeline itself — this actually rewrites tensors onto a CUDA device, it isn't a pure format conversion):

```bash
docker run --rm --gpus all -u $(id -u):$(id -g) -v "$(pwd)/softwares/SignalP6":/data ghcr.io/omatheuspimenta/signalp6:6.0h bash -c '
    cp -r /data/signalp-6-package/models /data/signalp-6-package/models_gpu &&
    signalp6_convert_models gpu /data/signalp-6-package/models_gpu
'
```

This copies the CPU weights first, so the original `models/` directory is left untouched — `--use_gpu false` still works afterwards, unaffected. `-u $(id -u):$(id -g)` matters: the container needs to write into your bind-mounted `softwares/`, which it can't do running as its own internal user. If `models_gpu/` doesn't exist yet and GPU is requested, the pipeline's pre-flight check (`bin/check_software_present.sh`) fails fast with this exact command rather than letting SignalP6 itself fail deep inside a container. Pass `--use_gpu false` to skip all of this and run on CPU (slower, but needs no extra setup).

### DeepLoc 2

**License**: DTU Health Tech academic-use license (the classifier checkpoints below); the ESM1b base encoder is a separate, freely-downloadable fair-esm weight, kept out of the image only because of its size (~7.3GB), not its license.

**Step 1 — download the DTU package.** Go to https://services.healthtech.dtu.dk/services/DeepLoc-2.1/, accept DTU's academic license/terms on that page, and download the package. You'll get a `DeepLoc-2.*`-style archive/repo that unpacks to a folder containing a `models/` subdirectory with the classifier checkpoints.

**Step 2 — place the classifier checkpoints.** Copy that package's `models/` directory contents so you end up with:

```
softwares/DeepLoc2/DeepLoc2/models/
```

**Required files** at that path:

- `models_esm1b/`, `models_prott5/` (per-fold checkpoint directories)
- `ESM1b_alphabet.pkl`, `ProtT5_alphabet.pkl`

**Step 3 — get the ESM1b base encoder** (a separate, freely-downloadable weight — not DTU-restricted, just too large, ~7.3GB, to bake into the Docker image). The pipeline's default "Fast" model needs this in addition to step 2. The easiest way to obtain it is to let a real, standalone DeepLoc 2 install download it once:

```bash
pip install .              # from the DTU package you downloaded in step 1
deeploc2 -f <any.fasta>     # first run downloads and caches the weights automatically
```

That populates `~/.cache/torch/hub/checkpoints/`. Copy (or symlink) that `hub/checkpoints/` directory so you end up with:

```
softwares/DeepLoc2/torch_cache/hub/checkpoints/esm1b_t33_650M_UR50S.pt
softwares/DeepLoc2/torch_cache/hub/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt
```

`torch_cache/` is laid out exactly like a `$TORCH_HOME` directory (i.e. it _is_ one — you can point `TORCH_HOME` at it directly outside the pipeline too). The "Accurate" model (ProtT5-XL, not used by default) needs a further ~11GB Hugging Face download and is not covered by the pre-flight check.

## Tools that need no manual setup

DeepCoil2 and Phobius bake in everything they need at Docker build time (DeepCoil2's weights come from its own PyPI package; Phobius's binary tarball is vendored under `docker/phobius/`) — nothing to place under `softwares/` for these two.

`rga_classify` (the final RGA-calling step, vendored from [omatheuspimenta/rgapredictor](https://github.com/omatheuspimenta/rgapredictor)) needs nothing under `softwares/` either — it's pure Python/pandas/PyYAML classification logic, no model weights, no license-gated database, no network access at runtime.
