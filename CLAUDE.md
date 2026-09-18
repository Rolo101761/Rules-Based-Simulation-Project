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

Currently: Phase 1.

## Stack
R + GenomicLayers, Python + scikit-optimize, Rcpp, Snakemake/Nextflow.

## Lessons carried from prior sacCer3 work
- tf_active pools must balance coverage vs spatial specificity — random fill dilutes spatial signal.
- K27 spreading must decouple promotion from tf_active expansion or it runs away.
- Saturation promotion prevents independent me1/me2 tuning — use sub-saturation.
- Empirical k-mer scan across all 65,536 8-mers is the best motif-ID approach.

## Conventions
- GenomicLayers installs from GitHub, not CRAN.
- Local repo, not cloud-synced.
- Flag partial/stubbed work as partial — don't present a stub as finished.
