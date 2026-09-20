# ============================================================
# dm6 coupled 4-mark model — Phase 2 SKELETON (development slice: chr4 by default)
# ============================================================
# PARTIAL / NOT PRODUCTION. What this is: the sacCer3 coupled reader-writer network (Rules 1-3, from
# R/coupled_rules.R) ported onto dm6, with the three sequence-anchored inputs seeded from the Phase 2
# motif scan (R/dm6_motif_scan.R, R/dm6_tata_variants.R):
#   K4  <- TATA variants that are actually TSS-specific in dm6 (TATAAAAG 23.9x, TATATAAG 13.5x at TSS -30..-16;
#          degenerate TATAWAWR is ~noise: TATATATA alone is 40% of its hits at 1.5x)
#   K9  <- satellite repeat units: (AATAT)3, (AAGAG)3, dodeca CCCGTACTCGGT  (+ random fill, as in yeast)
#   K27 <- PRE motifs: GAGAGAGAG (GAF) and GCCATHWT (Pho)                   (+ random fill)
#
# !! TARGETS ARE PLACEHOLDERS !!  No dm6 empirical coverage targets exist in this repo yet (ENCODE has no usable
# fly histone set; modENCODE is on GEO in dm3 coordinates and the liftover is parked). The numbers in
# `placeholder_targets` are rough round figures, NOT data. Fit numbers from this script are meaningless as
# biology until they are replaced.
#
# !! The objective here is NOT comparable to the sacCer3 9-state RMSE !!  Only the states in `scored_states`
# are scored (default: H3K4me3, H3K9me3, H3K27me3 = 3 states; add H3K4me1/2 with --score_K4me12=TRUE once real
# targets exist). RMSE_scored is an RMS over those N states only. The JSON records n_scored and the flag.
#
# Known omissions vs the sacCer3 script: no dwell-time / transition-matrix tracking, no per-iteration plots,
# no yeast-specific subtelomere/HML/HMR spatial tests. Added instead: seed-motif report and a TSS-promoter
# K4me3 enrichment check. chr4 has ~0 satellite seed motifs (0 dodeca, 0 (AAGAG)3, 8 (AATAT)3) so on the chr4
# slice K9 is driven almost entirely by random fill: the slice exercises TATA->K4 and PRE->K27 seeding, not
# satellite->K9 seeding.
#
# Rates (--meUp_sampler_*, --meDown_*, --K27_spread_abu) are given in "sacCer3-equivalent" units and are
# multiplied by (nucleosomes on this slice / --ref_nucs) so Phase-1-tuned values transfer across genome sizes.
#
# Run from the repo root:
#   Rscript R/dm6_coupled_model.R --chroms=chr4 --n_iter=50 --sim_tag=dev01
# Whole-slice options: --chroms=chr4,chrY (comma-separated seqnames)

.default_params <- list(
  chroms             = "chr4",   # comma-separated dm6 seqnames to simulate (development slice)
  meUp_sampler_K4    = 12000,    # sacCer3-equivalent per-iteration rates (sacCer3 script defaults, untuned for dm6)
  meDown_K4          = 23000,
  meUp_sampler_K9    = 12000,
  meDown_K9          = 13000,
  meUp_sampler_K27   = 12000,
  meDown_K27         = 19000,
  K27_spread_abu     = 5000,
  ref_nucs           = 64991,    # nucleosome count of the sacCer3 runs the rates above were tuned on
  tf_K4_stateWidth   = 600,      # bp of tf_K4_active laid down around each TATA hit (untuned)
  tf_K9_stateWidth   = 3000,     # ... around each satellite-unit hit (untuned)
  tf_K27_stateWidth  = 1500,     # ... around each PRE-motif hit (untuned)
  tf_K9_target_frac  = 0.25,     # fraction of nucleosomes in tf_K9_active after random fill (untuned; yeast used 0.78)
  tf_K27_target_frac = 0.25,     # fraction of nucleosomes in tf_K27_active after random fill (untuned; yeast used 0.40)
  score_K4me12       = FALSE,    # also score H3K4me1/me2 (no confirmed dm6 targets -> off)
  n_iter             = 50,
  verbose            = FALSE,
  seed               = NA,
  sim_tag            = ""
)

.parse_cli_args <- function(defaults) {
  params <- defaults
  for (a in commandArgs(trailingOnly = TRUE)) {
    if (!startsWith(a, "--")) next
    kv <- sub("^--", "", a)
    parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    if (!key %in% names(defaults)) stop(sprintf("Unknown parameter: --%s", key))
    value <- paste(parts[-1], collapse = "=")
    default_val <- defaults[[key]]
    params[[key]] <- if (is.character(default_val)) value else if (is.logical(default_val)) as.logical(value) else as.numeric(value)
  }
  params
}
params <- .parse_cli_args(.default_params)
if (!is.na(params$seed)) set.seed(params$seed)

simName <- "dm6_coupled_model"
outputDir <- paste0("output/", simName, "/")
if (!file.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
subSimName <- if (nzchar(params$sim_tag)) params$sim_tag else paste0("dm6_", gsub(",", "-", params$chroms))

# PLACEHOLDER targets (% of slice covered). NA = no target at all (state is reported, never scored).
placeholder_targets <- c(H3K4me1 = 8.0, H3K4me2 = 4.0, H3K4me3 = 3.0,
                         H3K9me1 = NA,  H3K9me2 = NA,  H3K9me3 = 10.0,
                         H3K27me1 = NA, H3K27me2 = NA, H3K27me3 = 8.0)
scored_states <- c("H3K4me3", "H3K9me3", "H3K27me3", if (params$score_K4me12) c("H3K4me1", "H3K4me2"))
TARGETS_ARE_PLACEHOLDER <- TRUE
cat("\n########################################################################\n")
cat("# dm6 model: coverage TARGETS ARE PLACEHOLDERS (not empirical data).\n")
cat(sprintf("# Scored states (N=%d): %s\n", length(scored_states), paste(scored_states, collapse = ", ")))
cat("# RMSE_scored is NOT comparable to the sacCer3 9-state RMSE.\n")
cat("########################################################################\n\n")


# LOAD LIBRARIES ----------------------------------------------------------

library(GenomicLayers)
source("R/removeShortGRanges.R")
source("R/removeGRangesBySize.R")
suppressMessages({
  require(Biostrings)
  require(BSgenome.Dmelanogaster.UCSC.dm6)
  require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

slice_chroms <- strsplit(params$chroms, ",", fixed = TRUE)[[1]]
genome <- keepBSgenomeSequences(BSgenome.Dmelanogaster.UCSC.dm6, slice_chroms)
genome_size <- sum(as.numeric(seqlengths(genome)))
nucleosomeWidth <- 147
saturationAbundance <- 1000000
cat(sprintf("Slice: %s (%.2f Mb)\n", params$chroms, genome_size / 1e6))


# SEQUENCE-ANCHORED SEED FACTORS (K4 / K9 / K27 tf_active pools) ----------

motif_bf <- function(name, pattern, layer, width) {
  createBindingFactor.DNA_consensus(name = name, patternString = pattern, mod.layers = layer, mod.marks = 1, stateWidth = width)
}
motif_factors <- list(
  # K4: TSS-specific TATA variants (see R/dm6_tata_variants.R)
  motif_bf("bf_TF_K4_TATAAAAG", "TATAAAAG", "tf_K4_active", params$tf_K4_stateWidth),
  motif_bf("bf_TF_K4_TATATAAG", "TATATAAG", "tf_K4_active", params$tf_K4_stateWidth),
  # K9: satellite units
  motif_bf("bf_TF_K9_AATAT3", "AATATAATATAATAT", "tf_K9_active", params$tf_K9_stateWidth),
  motif_bf("bf_TF_K9_AAGAG3", "AAGAGAAGAGAAGAG", "tf_K9_active", params$tf_K9_stateWidth),
  motif_bf("bf_TF_K9_dodeca", "CCCGTACTCGGT", "tf_K9_active", params$tf_K9_stateWidth),
  # K27: PRE motifs
  motif_bf("bf_TF_K27_GAF9", "GAGAGAGAG", "tf_K27_active", params$tf_K27_stateWidth),
  motif_bf("bf_TF_K27_Pho", "GCCATHWT", "tf_K27_active", params$tf_K27_stateWidth)
)
names(motif_factors) <- sapply(motif_factors, function(b) b$name)


# RULE SET (binding factors + lists; VERBATIM from the sacCer3 model) -----

source("R/coupled_rules.R")


# ABUNDANCES ---------------------------------------------------------------
# Built after the nucleosome layout is known, because rates are rescaled by slice size (see below).


# CREATE LAYERSET ----------------------------------------------------------

layer_names <- c("sampled_meUp", "sampled_meDown", "sampled_K9_meUp", "sampled_K27_meUp", "nucleosome",
                 "tf_K4_active", "tf_K9_active", "tf_K27_active",
                 "H3K4me_promotion", "H3K4me_any", "H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3",
                 "H3K9me_promotion", "H3K9me_any", "H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3",
                 "H3K27me_promotion", "H3K27me_any", "H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3")
scLayerSetNuc <- createLayerSet.BSgenome(genome = genome, layer.names = layer_names, n.layers = length(layer_names), verbose = FALSE)

scLayerSetNuc$layerSet[["nucleosome"]] <- randGrangesBigGenome(genome = genome,
                                                               sizeFunc = function(value, n) rep(x = value, times = n),
                                                               gapFunc = function(n, value) rpois(n = n, lambda = value),
                                                               argsSizeFunc = list(value = 147, n = 5),
                                                               argsGapFunc = list(value = 40, n = 7))
total_nucs <- length(scLayerSetNuc$layerSet$nucleosome)
rate_scale <- total_nucs / params$ref_nucs
cat(sprintf("Nucleosomes on slice: %d (rate scale vs sacCer3 reference = %.4f)\n", total_nucs, rate_scale))
sc_rate <- function(x) max(1, round(x * rate_scale))


# SEED tf_active LAYERS ----------------------------------------------------

cat("\n=== Seeding tf_active pools from sequence motifs ===\n")
seed_report <- data.frame(factor = names(motif_factors), layer = sapply(motif_factors, function(b) names(b$mods)),
                          regions_added = NA_integer_, stringsAsFactors = FALSE)
for (i in seq_along(motif_factors)) {
  lay <- names(motif_factors[[i]]$mods)
  before <- length(scLayerSetNuc$layerSet[[lay]])
  res <- tryCatch(runLayerBinding.BSgenome(layerList = scLayerSetNuc, factorSet = motif_factors[i],
                                           bf.abundances = setNames(saturationAbundance, names(motif_factors)[i]),
                                           verbose = params$verbose),
                  error = function(e) { cat(sprintf("  %-20s FAILED/no hits: %s\n", names(motif_factors)[i], conditionMessage(e))); NULL })
  if (!is.null(res)) scLayerSetNuc <- res
  seed_report$regions_added[i] <- length(scLayerSetNuc$layerSet[[lay]]) - before
  cat(sprintf("  %-20s -> %-13s +%d regions\n", names(motif_factors)[i], lay, seed_report$regions_added[i]))
}

pool_cov <- function(layer) 100 * sum(width(reduce(scLayerSetNuc$layerSet[[layer]]))) / genome_size

# Random fill up to target fraction of nucleosomes (as in the yeast model); K4 has no fill (motif-only pool)
fill_pool <- function(layer, target_frac) {
  in_pool <- unique(queryHits(findOverlaps(scLayerSetNuc$layerSet$nucleosome, scLayerSetNuc$layerSet[[layer]])))
  need <- round(target_frac * total_nucs) - length(in_pool)
  cat(sprintf("  %s: %d nucleosomes (%.1f%%) from motifs; target %.0f%%", layer, length(in_pool), 100 * length(in_pool) / total_nucs, 100 * target_frac))
  if (need > 0) {
    pick <- sample(setdiff(seq_len(total_nucs), in_pool), min(need, total_nucs - length(in_pool)))
    scLayerSetNuc$layerSet[[layer]] <<- c(scLayerSetNuc$layerSet[[layer]], scLayerSetNuc$layerSet$nucleosome[pick])
    cat(sprintf(" -> random fill +%d\n", length(pick)))
  } else cat(" -> motifs already exceed target, no fill\n")
  invisible(length(in_pool))
}
motif_nucs_K9 <- fill_pool("tf_K9_active", params$tf_K9_target_frac)
motif_nucs_K27 <- fill_pool("tf_K27_active", params$tf_K27_target_frac)
cat(sprintf("  tf_K4_active: %d regions, %.2f%% of slice (motif-only pool)\n", length(scLayerSetNuc$layerSet$tf_K4_active), pool_cov("tf_K4_active")))
cat(sprintf("Pool coverage: K4 %.2f%% | K9 %.2f%% | K27 %.2f%%\n", pool_cov("tf_K4_active"), pool_cov("tf_K9_active"), pool_cov("tf_K27_active")))
if (length(scLayerSetNuc$layerSet$tf_K4_active) == 0) cat("WARNING: tf_K4_active is empty (no TATA hit on this slice) -> no K4 marks can form.\n")
seed_report <- rbind(seed_report, data.frame(factor = c("K9 motif-seeded nucleosomes", "K27 motif-seeded nucleosomes"),
                                             layer = c("tf_K9_active", "tf_K27_active"),
                                             regions_added = c(motif_nucs_K9, motif_nucs_K27)))
write.csv(seed_report, paste0(outputDir, subSimName, ".seed_report.csv"), row.names = FALSE)


# STATE READOUT ------------------------------------------------------------

get_all_states_fast <- function(layerSet) {
  nuc <- layerSet$nucleosome
  st <- lapply(c(K4 = "H3K4", K9 = "H3K9", K27 = "H3K27"), function(m) {
    s <- rep("me0", length(nuc))
    for (lv in c(3, 2, 1)) {
      lay <- layerSet[[paste0(m, "me", lv)]]
      if (length(lay) > 0) { h <- unique(queryHits(findOverlaps(nuc, lay))); s[h[s[h] == "me0"]] <- paste0("me", lv) }
    }
    s
  })
  data.frame(K4 = st$K4, K9 = st$K9, K27 = st$K27, stringsAsFactors = FALSE)
}


# SIMULATION LOOP (same ordering / cleanup as the sacCer3 model) ------------

bfAbund_methylate <- list(
  K4 = c(sc_rate(params$meUp_sampler_K4), rep(saturationAbundance, 5)),
  K9 = c(sc_rate(params$meUp_sampler_K9), rep(saturationAbundance, 5)),
  K27 = c(sc_rate(params$meUp_sampler_K27), saturationAbundance, sc_rate(params$K27_spread_abu), rep(saturationAbundance, 4)))
names(bfAbund_methylate$K4) <- names(bfList_H3K4_methylate)
names(bfAbund_methylate$K9) <- names(bfList_H3K9_methylate)
names(bfAbund_methylate$K27) <- names(bfList_H3K27_methylate)
mk_demeth <- function(down_rate, lst) { a <- c(sc_rate(down_rate), rep(saturationAbundance, 4)); names(a) <- names(lst); a }
bf_demeth_K4 <- mk_demeth(params$meDown_K4, bfList_H3K4_demethyl)
bf_demeth_K9 <- mk_demeth(params$meDown_K9, bfList_H3K9_demethyl)
bf_demeth_K27 <- mk_demeth(params$meDown_K27, bfList_H3K27_demethyl)

bfList_eFirst <- c(bfList_H3K4_demethyl, bfList_H3K4_methylate, bfList_H3K9_demethyl, bfList_H3K9_methylate,
                   bfList_H3K27_demethyl, bfList_H3K27_methylate)
bfAbund_eFirst <- c(bf_demeth_K4, bfAbund_methylate$K4, bf_demeth_K9, bfAbund_methylate$K9, bf_demeth_K27, bfAbund_methylate$K27)
bfList_eFirst <- bfList_eFirst[names(bfAbund_eFirst)]
# The samplers' meDown abundance sits in slot 1 of each demethylation vector (rate given above); keep explicit for clarity:
bfAbund_eFirst["bf_meDown"] <- sc_rate(params$meDown_K4)
bfAbund_eFirst["bf_K9_meDown"] <- sc_rate(params$meDown_K9)
bfAbund_eFirst["bf_K27_meDown"] <- sc_rate(params$meDown_K27)

layersToClean <- c("sampled_meUp", "sampled_K9_meUp", "sampled_K27_meUp",
                   "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any", "H3K4me_demethylate",
                   "H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any", "H3K9me_demethylate",
                   "H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any", "H3K27me_demethylate")

cat(sprintf("\nRunning %d iterations (%s)\n", params$n_iter, subSimName))
scLayerSetBothAbund <- scLayerSetNuc
t0 <- Sys.time()
for (i in 1:params$n_iter) {
  scLayerSetBothAbund <- runLayerBinding.BSgenome(layerList = scLayerSetBothAbund, factorSet = bfList_eFirst, bf.abundances = bfAbund_eFirst,
                                                  verbose = params$verbose, collect.stats = T, keep.stats = T)
  for (thisLayer in layersToClean) {
    scLayerSetBothAbund$layerSet[[thisLayer]] <- removeShortGRanges(x = scLayerSetBothAbund$layerSet[[thisLayer]], minSize = nucleosomeWidth)
  }
  # Rule enforcement safety net (same as sacCer3): strip promotion flags that slipped past antagonism within an iteration
  ls <- scLayerSetBothAbund$layerSet
  if (length(ls$H3K4me_promotion) > 0 && length(ls$H3K9me3) > 0) {
    b <- findOverlaps(ls$H3K4me_promotion, ls$H3K9me3)
    if (length(b) > 0) scLayerSetBothAbund$layerSet$H3K4me_promotion <- ls$H3K4me_promotion[-unique(queryHits(b))]
  }
  ls <- scLayerSetBothAbund$layerSet
  if (length(ls$H3K27me_promotion) > 0 && length(ls$H3K4me3) > 0) {
    b <- findOverlaps(ls$H3K27me_promotion, ls$H3K4me3)
    if (length(b) > 0) scLayerSetBothAbund$layerSet$H3K27me_promotion <- ls$H3K27me_promotion[-unique(queryHits(b))]
  }
  scLayerSetBothAbund$layerSet[["nucleosome"]] <- removeGRangesBySize(x = scLayerSetBothAbund$layerSet[["nucleosome"]], verbose = FALSE, minSize = nucleosomeWidth, maxSize = nucleosomeWidth)
  for (s in c("sampled_meUp", "sampled_K9_meUp", "sampled_K27_meUp")) {
    scLayerSetBothAbund$layerSet[[s]] <- removeGRangesBySize(x = scLayerSetBothAbund$layerSet[[s]], verbose = FALSE, maxSize = 0)
  }
  if (i %% 10 == 0) cat(sprintf("  iter %d  (%.1f min elapsed)\n", i, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat(sprintf("Simulation wall time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))


# STATISTICS ---------------------------------------------------------------

final <- scLayerSetBothAbund$layerSet
current_states <- get_all_states_fast(final)
cov_of <- function(layer) if (length(final[[layer]]) == 0) 0 else 100 * sum(width(reduce(final[[layer]]))) / genome_size
mark_names <- names(placeholder_targets)
reached <- sapply(mark_names, cov_of)

cat("\n=== Genome coverage on slice (targets are PLACEHOLDERS; NA = no target) ===\n")
for (m in mark_names) cat(sprintf("%-9s %6.2f%%   target %s %s\n", m, reached[m],
                                  ifelse(is.na(placeholder_targets[m]), "  n/a", sprintf("%5.1f%%", placeholder_targets[m])),
                                  ifelse(m %in% scored_states, "[scored]", "")))
sc <- scored_states
RMSE_scored <- sqrt(mean((reached[sc] - placeholder_targets[sc])^2))
cat(sprintf("\nRMSE_scored over N=%d states (%s): %.3f  [PLACEHOLDER targets; NOT comparable to sacCer3 9-state RMSE]\n",
            length(sc), paste(sc, collapse = ","), RMSE_scored))

calc_overlap <- function(l1, l2) {
  if (length(final[[l1]]) == 0 || length(final[[l2]]) == 0) return(0)
  100 * length(unique(queryHits(findOverlaps(final[[l1]], final[[l2]])))) / length(final[[l1]])
}
K4_K27 <- calc_overlap("H3K4me3", "H3K27me3"); K9_K4 <- calc_overlap("H3K9me3", "H3K4me3"); K9_K27 <- calc_overlap("H3K9me3", "H3K27me3")
cat("\n=== Rule checks ===\n")
cat(sprintf("Rule 1  K4me3+K27me3 co-occurrence: %.2f%% (%s)\n", K4_K27, ifelse(K4_K27 < 5, "PASS <5%", "FAIL")))
cat(sprintf("Rule 2  K9me3+K4me3  co-occurrence: %.2f%% (%s)\n", K9_K4, ifelse(K9_K4 < 5, "PASS <5%", "FAIL")))
cat(sprintf("        K9me3+K27me3 co-occurrence: %.2f%% (informational)\n", K9_K27))
if (length(final$H3K27me3) > 0) {
  blocks <- reduce(final$H3K27me3, min.gapwidth = 100)
  cat(sprintf("Rule 3A K27me3 blocks: %d, mean %.0f bp (%.1f nucleosomes) (%s)\n", length(blocks), mean(width(blocks)), mean(width(blocks)) / nucleosomeWidth,
              ifelse(mean(width(blocks)) > 2 * nucleosomeWidth, "domains forming", "no domain formation")))
} else cat("Rule 3A: no K27me3 formed\n")

# dm6-specific spatial check: K4me3 at annotated promoters (TSS -300..+100) vs slice-wide
tx <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
tx <- tx[as.character(seqnames(tx)) %in% slice_chroms]
prom <- suppressWarnings(trim(unique(promoters(tx, upstream = 300, downstream = 100))))
tss_K4_frac <- NA; slice_K4_frac <- NA
if (length(prom) > 0 && length(final$H3K4me3) > 0) {
  tss_K4_frac <- 100 * mean(countOverlaps(prom, final$H3K4me3) > 0)
  slice_K4_frac <- 100 * length(unique(queryHits(findOverlaps(final$nucleosome, final$H3K4me3)))) / length(final$nucleosome)
  cat(sprintf("\nK4me3 at annotated promoters: %.1f%% of %d promoters vs %.1f%% of all nucleosomes (%.1fx)\n",
              tss_K4_frac, length(prom), slice_K4_frac, tss_K4_frac / max(slice_K4_frac, 1e-9)))
  cat("  (K4 seeds come from TSS-specific TATA variants only, so <1x here means TATA seeds miss most promoters)\n")
}


# SAVE ---------------------------------------------------------------------

write.csv(as.data.frame.table(table(K4 = current_states$K4, K9 = current_states$K9, K27 = current_states$K27)),
          file = paste0(outputDir, subSimName, ".state_counts.csv"), row.names = FALSE)
write.csv(data.frame(mark = mark_names, coverage_pct = as.numeric(reached), placeholder_target_pct = as.numeric(placeholder_targets),
                     scored = mark_names %in% scored_states), file = paste0(outputDir, subSimName, ".coverage.csv"), row.names = FALSE)
writeLines(c(
  "{",
  sprintf('  "RMSE_scored": %.6f,', RMSE_scored),
  sprintf('  "n_scored": %d,', length(sc)),
  sprintf('  "scored_states": [%s],', paste(sprintf('"%s"', sc), collapse = ", ")),
  '  "targets_are_placeholder": true,',
  '  "comparable_to_sacCer3_9state_RMSE": false,',
  sprintf('  "slice": "%s",', params$chroms),
  sprintf('  "nucleosomes": %d,', total_nucs),
  sprintf('  "rule1_K4me3_K27me3_pct": %.4f,', K4_K27),
  sprintf('  "rule2_K9me3_K4me3_pct": %.4f,', K9_K4),
  sprintf('  "coverage": {%s}', paste(sprintf('"%s": %.6f', mark_names, reached), collapse = ", ")),
  "}"), paste0(outputDir, subSimName, ".objective.json"))
cat(sprintf("\nOBJECTIVE_JSON: %sobjective.json\n", paste0(outputDir, subSimName, ".")))
cat("DONE (skeleton run; placeholder targets)\n")
