# Publishing the pipeline's Docker images

This pipeline builds seven "baseline" images (§4.1 of `docs/software-setup.md`'s
license-gating rationale, decision #1 in `PLAN.md`): every one of them bakes in only
redistributable content — OS, runtime, non-restricted code/dependencies — and never bakes
in license-gated weights, binaries, or databases. Those live under `softwares/<Tool>/`,
supplied by the user, and are bind-mounted in at runtime. **Because of that, every image
listed here is safe to publish publicly** — there is nothing in any of them that requires a
license the pipeline doesn't already have a right to redistribute.

This repository is published on GitHub only (not submitted to nf-core), so these images are
meant for **GHCR under your own GitHub account/org** (`ghcr.io/omatheuspimenta/...`), not
`quay.io` or any nf-core-managed registry. `quay.io` already appears in `nextflow.config`
(`docker.registry = 'quay.io'`) and in a couple of modules' `container` directives — that's
just this pipeline's *default registry prefix* for resolving bare image names against a
locally-built tag (see `PLAN.md` Stage 4a's bug notes on this), not an indication that
anything is actually hosted on quay.io. Nothing in this pipeline is currently published
anywhere; every image referenced by a module is built and tagged locally for development
and testing.

## The seven images

| Image | Dockerfile | Suggested tag | Size (dev build) |
|---|---|---|---|
| `deepcoil2` | `docker/deepcoil2/Dockerfile` | `2.0.2` (DeepCoil version) | ~11.8GB |
| `phobius` | `docker/phobius/Dockerfile` | `1.01` (Phobius version) | ~190MB |
| `interproscan` | `docker/interproscan/Dockerfile` | `5.78-109.0` (InterProScan release) | ~620MB |
| `deeploc2` | `docker/deeploc2/Dockerfile` | `1.0.0` | ~9.1GB |
| `signalp6` | `docker/signalp6/Dockerfile` | `6.0h` | ~5.6GB |
| `deeptmhmm` | `docker/deeptmhmm/Dockerfile` | `1.0` | ~8.1GB |
| `rga_classify` | `docker/rga_classify/Dockerfile` | `0.0.1` (vendored `sugarcane-rga` version) | ~1.1GB |

`rga_report` (`modules/local/rga_report/`) has no image of its own — it deliberately reuses
`rga_classify`'s image (see `PLAN.md` Stage 7), so there is nothing separate to publish for it.

Five of the seven modules already carry a commented-out `ghcr.io/omatheuspimenta/...`
`container` line right next to the active local one (`deepcoil2`, `deeploc2`, `signalp6`,
`deeptmhmm`, `rga_classify`) — the tags above match those exactly. `phobius` and
`interproscan` don't have one yet; the table above is the convention to extend to them
whenever they're actually published (see "Switching the pipeline over" below).

## One-time setup: authenticate to GHCR

```bash
# A GitHub personal access token (classic or fine-grained) with at least
# `write:packages` scope (and `read:packages` to pull it back down later).
export CR_PAT=<your token>
echo "$CR_PAT" | docker login ghcr.io -u omatheuspimenta --password-stdin
```

GHCR packages are **private by default** on first push. If these are meant to be pulled by
Nextflow's `docker` profile without authentication (the normal case for a public pipeline),
make each package public afterwards: on GitHub, go to the package's page
(`github.com/omatheuspimenta?tab=packages`) → **Package settings** → **Change visibility** →
**Public**. Do this deliberately, one package at a time — don't default every future package
to public without checking what's actually safe to publish first, even though everything in
the table above already is.

## Build, tag, and push

All seven `docker build` invocations use the pattern already established throughout this
repo (see `PLAN.md` Stages 3–4, 6): run from the pipeline root, with the Dockerfile's own
directory as the build context (so a Dockerfile's `COPY src/...` lines resolve correctly and
the build never has to reach outside `docker/<tool>/`, in particular never into the
gitignored, ~600GB `softwares/`).

```bash
cd /path/to/omatheuspimenta-rgaprofiler

for tool_tag in \
    "deepcoil2:2.0.2" \
    "phobius:1.01" \
    "interproscan:5.78-109.0" \
    "deeploc2:1.0.0" \
    "signalp6:6.0h" \
    "deeptmhmm:1.0" \
    "rga_classify:0.0.1"
do
    tool="${tool_tag%%:*}"
    tag="${tool_tag##*:}"
    docker build -t "ghcr.io/omatheuspimenta/${tool}:${tag}" -f "docker/${tool}/Dockerfile" "docker/${tool}"
    docker push "ghcr.io/omatheuspimenta/${tool}:${tag}"
done
```

Build one at a time and check each `docker build` output if you're publishing for the first
time — some of these (`deeploc2`, `deeptmhmm`, `signalp6`) take several minutes and multiple
gigabytes each; see each Dockerfile's own header comment for what it does and why.

Also push a `latest` tag per image if you want `ghcr.io/omatheuspimenta/<tool>:latest` to
resolve to whatever you just built:

```bash
docker tag "ghcr.io/omatheuspimenta/${tool}:${tag}" "ghcr.io/omatheuspimenta/${tool}:latest"
docker push "ghcr.io/omatheuspimenta/${tool}:latest"
```

## Switching the pipeline over to the published images

Right now every module's `container` directive points at a **local** tag
(`deepcoil:nextflow`, `signalp6:baseline`, `rga_classify:baseline`, …), which is what this
whole engagement has built and validated against (`PLAN.md` Stages 3–7). Publishing images
does **not** by itself change what the pipeline runs — nothing in `modules/local/*/main.nf`
is touched by the steps above.

To actually make a module pull from GHCR instead of requiring a local build, swap which
`container` line is active in that module's `main.nf`:

```groovy
// container 'deeptmhmm:baseline'                          // was active
container 'ghcr.io/omatheuspimenta/deeptmhmm:1.0'           // now active
```

Do this one module at a time and re-run that module's `nf-test` (real container, `-profile
docker`) before moving to the next — the whole point of validating everything against local
builds first was to catch problems (like the several found during Stages 3–7, e.g. the
`quay.io` registry-prefix gotcha, the `ln -sfn` directory-vs-file collisions) before they'd
also have to be debugged through a slow push/pull cycle. `phobius` and `interproscan` need
their own `ghcr.io/omatheuspimenta/...` build+push (using the tags in the table above) before
they have anything to switch to; the other five already do once published.

## Rebuilding and republishing

Every image here is expected to change over time (a Dockerfile fix, a dependency bump). There
is no CI automation for this (deliberately — see `PLAN.md` Stage 8, which calls a publish
GitHub Actions workflow optional and maintainer-gated); publishing is a manual step you run
yourself with the loop above. Bump the tag when the *contents* of an image change in a way
that matters (a real fix, not just a rebuild with identical layers) so that a pinned
`container` directive in a module continues to mean something specific — don't silently
overwrite a version tag that something else might already be depending on; push a new one
and update the `container` line instead.
