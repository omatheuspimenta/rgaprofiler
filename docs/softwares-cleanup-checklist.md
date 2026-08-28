# `softwares/` cleanup checklist

`softwares/` is gitignored and was never tracked by git, so nothing below risks losing committed history — but it currently holds ~600GB of local data, a meaningful chunk of which looks unintentional (build artifacts, runtime scratch directories, duplicate copies). This is a manual checklist for the maintainer to action at their convenience; **nothing here is deleted automatically** by the pipeline or by any agent working on this repo.

Verify a size/necessity before deleting anything — this list reflects what was observed during the initial repository inventory (2026-08-27) and may be stale by the time you read it.

## DeepTMHMM (`softwares/DeepTMHMM/`, ~524G total)

- [ ] `DeepTMHMM-venv/` (~1.3G) and `DeepTMHMM-venv-new/` (~4.8G) — two fully-materialized Python virtualenvs checked into the tree. Once `docker/deeptmhmm/Dockerfile` exists (Stage 4c) and reproduces this environment inside the container, these local venvs are redundant for pipeline purposes (they may still be useful for your own manual/interactive runs — check before deleting).
- [ ] `pytorch-1.5/torch-1.5.0+cu92-cp38-cp38-linux_x86_64.whl` (576M) — the exact CUDA 9.2 torch wheel; keep a copy somewhere if `docker/deeptmhmm/Dockerfile`'s build needs to fetch it from a stable location, but it doesn't need to live inside this working tree.
- [ ] `__MACOSX/` — leftover from a macOS zip extraction; safe to delete.
- [ ] `DeepTMHMM-Academic-License-v1.0/results_r570/embeddings/` (bulk of the 524G) — per-protein embedding zip archives from the R570 reference run. **Do not delete** — Stage 5/Stage 6 use this run's other outputs (`TMRs.gff3`, `predicted_topologies.3line`) as ground truth, and this directory may still be referenced. Confirm nothing needs it before removing.

## InterProScan (`softwares/InterProScan/`, ~77G total)

- [ ] `interproscan-5.78-109.0-64-bit.tar.gz` — the original download tarball, sitting next to its own already-unpacked contents. Redundant once unpacked; safe to delete if you're confident you can re-download it if needed.
- [ ] `interproscan-5.78-109.0/temp/` (~20G) and `interproscan-5.78-109.0/work/` (~134M) — runtime scratch directories left over from a prior manual run. Safe to delete; InterProScan recreates them as needed.

## DeepLoc2 (`softwares/DeepLoc2/`)

- [ ] `build/lib/DeepLoc2/` — a `setuptools` build artifact that duplicates the entire `DeepLoc2/models/` weight tree a second time. Safe to delete once you've confirmed the top-level `DeepLoc2/DeepLoc2/models/` (the one `bin/check_software_present.sh deeploc2` validates against) is intact.

## SignalP6 (`softwares/SignalP6/`)

- [ ] `signalp-6-package/build/lib/signalp/model_weights/` — same pattern as DeepLoc2: a duplicated build artifact. Safe to delete once you've confirmed `signalp-6-package/models/` is intact.

## DeepCoil (`softwares/DeepCoil/`)

- [ ] `Dockerfile` — near-duplicate of the canonical `docker/deepcoil2/Dockerfile`, but with a divergent `ENTRYPOINT ["deepcoil"]` that would actually break Nextflow's container execution model (Nextflow needs to run its own command inside the container; an image whose entrypoint is fixed to `deepcoil` would try to pass Nextflow's own launch command as arguments to the `deepcoil` CLI instead of executing it). Delete this file — `docker/deepcoil2/Dockerfile` is the one actually used by the `modules/local/deepcoil2` module and should remain the single source of truth for this image.

## General

- [ ] Re-run `du -sh softwares/*` after acting on the above to confirm the expected reclaimed space, and re-run `./bin/check_software_present.sh <tool>` for every tool afterward to confirm nothing required was accidentally removed.
