# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# GenomicLayers → Drosophila dm6

## What this is
Extending GenomicLayers (Dave Gerrard's R package) to model chromatin state across the Drosophila dm6 genome (~140Mb, ~780,000 nucleosomes) with ML-optimised parameters. Chat is used for design/architecture; this repo is where it gets built.

## Core architecture: 4 marks as ONE coupled reader-writer network
NOT four parallel tracks. GenomicLayers factors read both sequence AND existing layer states, so marks are wired together through mark-mark interactions.

Marks (functional alphabet):
- K4me3 — active promoter
- K36me3 — active gene body
- K9me3 — constitutive silencing
- K27me3 — Polycomb repression

Couplings:
1. Transcription as shared layer: K4me3 at promoter licenses an "active" layer across the gene body; K36me3 writer reads that layer (emergent targeting, not bespoke).
2. Antagonism: PRC2 inhibited by K4me3/K36me3 (Schmitges/Yuan 2011), so K27 writer skips active regions; K9/K27 mutually exclusive. Domain boundaries emerge from this.
3. Within-mark feedback (spreading): K27->EED->K27, K9->HP1/SUV39H->K9.

Central engineering risk: rate-balancing the coupled feedback loops (oscillation / latching to wrong stable state / runaway). Most ML effort goes here, not "picking marks".

## Sequence targeting
Small set of sequence-anchored inputs propagate through the network: TATA (K4), satellite repeats (K9), PRE motifs (K27).

## Phases
1. Bayesian opt (scikit-optimize) on sacCer3, sub-saturation promotion for all 9 states
2. Port to dm6 with PRE/satellite/TATA motifs
3. ML optimisation on dm6
4. K9 spreading (HP1/SUV39H)
5. Gene expression readout vs modENCODE RNA-seq (key payoff)
6. Multi-cell population simulation (N stochastic "cells")

Currently: Phase 1 (in progress — see Phase 1 progress below).

## Stack
R + GenomicLayers, Python + scikit-optimize, Rcpp, Snakemake/Nextflow.

R deps (GenomicLayers, Biostrings, GenomicRanges, BSgenome.Scerevisiae.UCSC.sacCer3,
TxDb.Scerevisiae.UCSC.sacCer3.sgdGene) are already installed locally — no setup needed.
Python: `.venv/` at repo root has scikit-optimize installed (`python/requirements.txt`).

## Phase 1 progress

Ported the most mature prior sacCer3 model (`Model 5 tf_active +DNA seq.R` from the
OneDrive research folder) into `R/model5_tf_active_dna_seq.R`, with hardcoded
OneDrive/Windows paths replaced by relative paths (run from repo root).

The script now takes the rate-balancing knobs as CLI args instead of hardcoded values,
so a Python optimizer can drive it:
`--meUp_sampler_K4/K9/K27`, `--meDown_K4/K9/K27` (sampling rates per iteration),
`--K27_spread_abu` (PRC2/EED self-recruitment rate), `--tf_K9_target_frac`,
`--tf_K27_target_frac` (tf_active pool sizes), plus `--n_iter`, `--verbose`, `--seed`,
`--sim_tag`. It writes a machine-readable `output/Model_5_Rule3_only/<sim_tag>.objective.json`
(RMSE_all, RMSE_me3, per-mark coverage) after each run.

`python/optimize_phase1.py` drives this with `skopt.gp_minimize` over the 9 knobs above,
shelling out to Rscript per trial, logging every trial to `output/skopt_trials.csv` and
the best params to `output/skopt_best.json`.

**Timing finding (important for planning runs):** `verbose=TRUE` (the original script's
default) makes GenomicLayers print a line per binding-factor match check, and this
printing — not the underlying computation — was the dominant cost: a 100-iteration
run took ~90 minutes with verbose on. With `verbose=FALSE`, a 50-iteration run took
~24.5 minutes and reached RMSE_all=3.16 (close to the original 100-iteration run's
fit), so 50 iterations is a reasonable proxy for optimization. Per-iteration cost is
NOT constant — it grows as marks spread across the genome (more GRanges → costlier
overlap checks each cycle), so runtime doesn't scale linearly with `--n_iter`.

Practical implication: a real optimization campaign (15-20 skopt trials at
`--n_iter=50`) is a ~6-8 hour job. Run it unattended/in the background between
sessions rather than blocking an interactive session — that run hasn't been done yet.

## Repo / GitHub
Connected to `Rolo101761/Rules-Based-Simulation-Project` on GitHub (public) — this was
an existing near-empty repo from prior sacCer3 research, now the home for this port.
Local git identity for this repo only: `Amin <amin.miah@icloud.com>` (not global).
No Co-Authored-By lines in commits — user's name only.

## Lessons carried from prior sacCer3 work
- tf_active pools must balance coverage vs spatial specificity — random fill dilutes spatial signal.
- K27 spreading must decouple promotion from tf_active expansion or it runs away.
- Saturation promotion prevents independent me1/me2 tuning — use sub-saturation.
- Empirical k-mer scan across all 65,536 8-mers is the best motif-ID approach.

## Conventions
- GenomicLayers installs from GitHub, not CRAN.
- Local repo, not cloud-synced.
- Flag partial/stubbed work as partial — don't present a stub as finished.
