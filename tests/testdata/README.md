# Test data

`r570_subsample.fasta` — 4 real protein sequences hand-picked (not random) from the full R570
proteome (`softwares/data_test/SofficinarumxspontaneumR570_771_v2.1.protein.fa`, 299,731
sequences), chosen so `-profile test` exercises a genuine positive case for every tool in the
pipeline, not just "does it run." Each was selected by grepping the real, pre-existing R570
ground-truth output already committed under `softwares/*/r570*`/`results_r570/` — see
`docs/software-setup.md` and `PLAN.md` §6 Stage 5 for how each was found.

| Sequence ID | Length | Why it's here |
|---|---|---|
| `SoffiXsponR570.7os1g018900.1.p` | 76 aa | Signal peptide positive — SignalP6 predicts `SP`, cleavage site 30–31 (`SignalP6/.../r570/prediction_results.txt`); also DeepTMHMM's `signal`/`outside` (0 TMRs) case. Used throughout Stages 3–4's manual confirmation runs. |
| `SoffiXsponR570.7os1g005900.1.p` | 128 aa | Transmembrane positive — DeepTMHMM predicts 3 TM helices (`DeepTMHMM/.../results_r570/TMRs.gff3`). |
| `SoffiXsponR570.09Cg101200.1.p` | 234 aa | Coiled-coil positive — DeepCoil2 predicts a strong, sustained coiled-coil region (174/234 residues with `cc` probability > 0.7 in `DeepCoil/r570/results/`). |
| `SoffiXsponR570.10Eg022100.1.p` | 1160 aa | **RGA positive** — a real CC-NB-LRR (CNL) resistance gene: InterProScan calls an `Rx N-terminal domain` (coiled-coil), two `NB-ARC domain` (Pfam PF00931) hits, a `Leucine-rich repeat region`, `PANTHER:PTHR23155 DISEASE RESISTANCE PROTEIN RP`, and `Coils` hits — the canonical domain architecture this whole pipeline exists to detect (`InterProScan/.../results/r570_interpro.tsv`). This is the most directly on-topic validation case in the test set. |

Lengths are 1 residue shorter than the real R570 proteome's own sequences: the trailing stop
codon (`*`) that real proteome FASTAs commonly carry is already stripped here. `FASTA_QC`
normally does this before any tool sees the data, but per-module `nf-test`s call each module
directly (bypassing `FASTA_QC`) — Phobius and InterProScan both fail outright on a trailing `*`
(confirmed while building this test set; see PLAN.md Stage 3's bug log for the original
discovery), so the committed FASTA here is pre-cleaned rather than relying on every test to
clean it itself.

This is deliberately small (4 sequences, ~1.6kb total) to keep `-profile test` fast — DeepTMHMM's
CPU-bound topology decoding is the slowest step per sequence (see nf-test timings under
`modules/local/deeptmhmm/tests/`), so the count is kept minimal rather than the sequences
themselves being trimmed shorter (all 4 are used at their real, full length — length itself is
part of what's being validated for the TM and RGA-domain cases).

`samplesheet_full.csv` references the real, full R570 proteome directly from
`softwares/data_test/` (not copied — that file is 200MB and gitignored). It only works once
`softwares/data_test/` is populated per `docs/software-setup.md`; `-profile test_full` is an
opt-in, full-scale validation run, not part of routine/CI testing.
