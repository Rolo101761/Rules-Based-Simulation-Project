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

## Tech + refs

### Objective function (Phase 1+)
Score simulated mark abundances against target abundances using RMSE + R² (as in prior MSc sacCer3 work). skopt minimises RMSE across all 9 states, sub-saturation promotion params are the search space.

Note: Coupling 3 (K27→EED→K27) is the self-recruitment/positive-feedback design (MSc "Rule 3A") — not the earlier two-speed-eraser variant (MSc "Rule 3") gated on K9 presence. Don't conflate the two; 3A was the one that worked.

## Phase 1 progress

`R/model5_tf_active_dna_seq.R` (filename kept for continuity, content corrected —
see below) is the sacCer3 DNA-motif-seeded coupled model, ported from the OneDrive
research folder with hardcoded OneDrive/Windows paths replaced by relative paths
(run from repo root).

**Correction applied:** initially ported from `Model 5 tf_active +DNA seq.R`, which
per the MSc 5-model framework (M2 = all rules active, M3/M4/M5 = one rule each) is
the Rule-3-only ablation — its header even says so, but stray Rule 1/2 comments left
in the binding-factor code made it look like all three rules were wired when they
weren't (confirmed by reading the actual `profile.layers`, not the comments — the
antagonist mark was simply missing from `bf_promotion_2_3`/`bf_K27_promotion_2_3`/
`bf_meUp`/`bf_K27_meUp`). Fixed by porting the genuinely-all-rules-active logic from
`Model 2 tf_active + DNA sequence.R` (M2 in the MSc framework): antagonism is now
enforced at both the initial targeting stage AND the me2→me3 promotion step, plus a
per-iteration cleanup pass that strips any promotion flag that slipped through a
same-iteration ordering race. Verified: RULE 1 and RULE 2 self-validation both
genuinely PASS now (K4me3+K27me3 bivalency 0.00%, K9me3+K4me3 co-occurrence 0.35%,
vs. 13.18%/coincidental-2.87% before the fix). Output dir renamed
`Model_5_Rule3_only` → `sacCer3_coupled_model` to stop implying it's an ablation.
Lesson: when porting from prior research scripts, verify binding-factor
`profile.layers`/`profile.marks` directly — comments in this codebase have drifted
from the code they describe at least once already.

The script now takes the rate-balancing knobs as CLI args instead of hardcoded values,
so a Python optimizer can drive it:
`--meUp_sampler_K4/K9/K27`, `--meDown_K4/K9/K27` (sampling rates per iteration),
`--K27_spread_abu` (PRC2/EED self-recruitment rate), `--tf_K9_target_frac`,
`--tf_K27_target_frac` (tf_active pool sizes), plus `--n_iter`, `--verbose`, `--seed`,
`--sim_tag`. It writes a machine-readable `output/sacCer3_coupled_model/<sim_tag>.objective.json`
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

## Future expansions
- Ablate one coupling at a time (transcription-sharing, antagonism, spreading) against a no-couplings baseline and an all-couplings full model.
- MSc 5-model framework (for ablation experiments): M1 = no rules (baseline) · M2 = all rules active · M3/M4/M5 = one rule active each (isolates that rule's contribution). Same logic applies here.
