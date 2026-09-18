# ============================================================
# MODEL 5 - DNA motif-seeded tf_active per mark RULE 3 ONLY
# ============================================================
# Extension of Dave's architecture with biologically-grounded targeting:
#   K4: TATA box (TATATAAA) marks active/regulated promoters
#        - ~20% of yeast genes (Basehoar et al. 2004, Cell)
#   #   K9: GTGTGGGT (8-mer subtelomere-enriched TG-repeat) — 217 hits in sacCer3,
#       29.6x enriched at subtelomeric regions (63.6% of hits within 20kb of chr ends).
#       stateWidth=17000 → ~24% motif coverage + random fill to 66%.
#       Identified by empirical k-mer enrichment scan of subtelomeric vs internal regions.
#   K27: Random subset (yeast lacks PRC2/H3K27me3 biologically)
#
# Architecture matches Dave's original exactly:
# - Non-persistent promotion and demethylation flags
# - Saturation abundances for all promotion/demethylation steps
# - Genome-wide sampler, tf_active filtering at bf_meUp stage
# - bf_meDown targets all nucleosomes
#
# Rules:
#   Rule 1: K4me3 blocks K27me3 (at bf_meUp and me2→me3)
#   Rule 2: K9me3 blocks K4me3 (at bf_meUp and me2→me3)
#   Rule 3: K27me3 self-recruitment (spreading)
# ============================================================

# Run from the project root (e.g. `Rscript R/model5_tf_active_dna_seq.R`).
#
# Rate-balancing knobs are overridable via --key=value CLI args, e.g.:
#   Rscript R/model5_tf_active_dna_seq.R --meUp_sampler_K4=9000 --n_iter=30 --verbose=FALSE --sim_tag=trial01
# This is the interface the Phase 1 scikit-optimize loop (python/optimize_phase1.py) drives.

.default_params <- list(
  meUp_sampler_K4    = 12000,  # H3K4 nucleosome sampling rate per iteration
  meDown_K4          = 23000,  # H3K4 demethylation sampling rate per iteration
  meUp_sampler_K9    = 12000,  # H3K9 nucleosome sampling rate per iteration
  meDown_K9          = 13000,  # H3K9 demethylation sampling rate per iteration
  meUp_sampler_K27   = 12000,  # H3K27 nucleosome sampling rate per iteration
  meDown_K27         = 19000,  # H3K27 demethylation sampling rate per iteration
  K27_spread_abu     = 5000,   # bf_K27_meUp_spread abundance (PRC2/EED self-recruitment rate)
  tf_K9_target_frac  = 0.78,   # target fraction of nucleosomes in tf_K9_active pool (motif-seeded + random fill)
  tf_K27_target_frac = 0.40,   # target fraction of nucleosomes in tf_K27_active pool (random, no PRC2 motif in yeast)
  n_iter             = 100,    # number of runLayerBinding.BSgenome iterations
  verbose            = TRUE,   # GenomicLayers per-binding-factor match logging
  seed               = NA,     # RNG seed; NA leaves runs stochastic (the script default before parameterization)
  sim_tag            = ""      # appended to output filenames so parallel/repeated trials don't clobber each other
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
    params[[key]] <- if (is.character(default_val)) {
      value
    } else if (is.logical(default_val)) {
      as.logical(value)
    } else {
      as.numeric(value)
    }
  }
  params
}

params <- .parse_cli_args(.default_params)
if (!is.na(params$seed)) set.seed(params$seed)

simName <- "Model_5_Rule3_only"
outputDir <- paste0("output/", simName, "/")
if(!file.exists(outputDir)) dir.create(outputDir, recursive = TRUE)


# LOAD LIBRARIES ----------------------------------------------------------

library(GenomicLayers)
source("R/removeShortGRanges.R")
source("R/removeGRangesBySize.R")
source("R/plotLayers.R")

require(Biostrings)
require(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3
nucleosomeWidth <- 147


# ========== DNA MOTIF BINDING FACTORS (run once to seed tf_active) ==========

# K4: TATA box at regulated yeast promoters (Basehoar et al. 2004, Cell)
# stateWidth controls spread from each motif hit
# Adjust stateWidth after checking coverage from initial scan
# TATATAAA found 1,012 hits in yeast genome.
# target K4 coverage -24% empirical 19.7% and 20% buffer to account for stochatism.
# # statewidth (0.24 * 12,157,105) / 1,012 = 270bp
# verified with actual coverage percentage 27.5%
bf_TF_K4 <- createBindingFactor.DNA_consensus(name="bf_TF_K4",
                                              patternString = "TATATAAA",
                                              mod.layers = "tf_K4_active",
                                              mod.marks = 1,
                                              stateWidth = 1700)

bf_TF_K9 <- createBindingFactor.DNA_consensus(name="bf_TF_K9",
                                              patternString = "GTGTGGGT",
                                              mod.layers = "tf_K9_active",
                                              mod.marks = 1,
                                              stateWidth = 17000)


# ========== H3K4 METHYLATION FACTORS ==========

bf_meUp_sampler <- createBindingFactor.layer_region(name="bf_meUp_sampler", 
                                                    profile.layers="nucleosome", 
                                                    patternLength = nucleosomeWidth,
                                                    profile.marks = 1,
                                                    mod.layers="sampled_meUp",
                                                    mod.marks=1,
                                                    stateWidth= nucleosomeWidth)

# Requires tf_K4_active AND sampled AND NOT K9me3 (Rule 2)
bf_meUp <- createBindingFactor.layer_region(name="bf_meUp", 
                                            profile.layers=c("tf_K4_active", "sampled_meUp"), 
                                            patternLength = 6,
                                            profile.marks = c(1, 1),
                                            mod.layers="H3K4me_promotion",
                                            mod.marks=1,
                                            stateWidth= 6 + (2*nucleosomeWidth))

# Non-persistent promotion flags (Dave's pattern)
bf_promotion_0_1 <- createBindingFactor.layer_region(name="bf_promotion_0_1", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1,0,0,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 1,0,0,1))

bf_promotion_1_2 <- createBindingFactor.layer_region(name="bf_promotion_1_2", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1,1,0,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,1,0,1))

# Rule 2: K9me3 blocks K4 at me2→me3
bf_promotion_2_3 <- createBindingFactor.layer_region(name="bf_promotion_2_3", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1,0,1,0), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,0,1,1))

bf_promotion_3_3 <- createBindingFactor.layer_region(name="bf_promotion_3_3", 
                                                     patternLength = nucleosomeWidth, 
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1,0,0,1), 
                                                     mod.layers = c("H3K4me_promotion", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                     mod.marks = c(0, 0,0,1,1))


# ========== H3K4 DEMETHYLATION (Dave's pattern — targets all nucleosomes) ==========

bf_meDown <- createBindingFactor.layer_region(name="bf_meDown", 
                                              profile.layers="nucleosome",
                                              patternLength = nucleosomeWidth,
                                              profile.marks = 1,
                                              mod.layers="H3K4me_demethylate",
                                              mod.marks=1,
                                              stateWidth=nucleosomeWidth)

# Non-persistent demethylation (Dave's pattern)
bf_demethylate_3_2 <- createBindingFactor.layer_region(name="bf_demethylate_3_2", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,0,1), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 0,1,0,1))

bf_demethylate_2_1 <- createBindingFactor.layer_region(name="bf_demethylate_2_1", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,1,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 1,0,0,1))

bf_demethylate_1_0 <- createBindingFactor.layer_region(name="bf_demethylate_1_0", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 1,0,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 0,0,0,0))

bf_demethylate_0_0 <- createBindingFactor.layer_region(name="bf_demethylate_0_0", 
                                                       patternLength = nucleosomeWidth, 
                                                       profile.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                       profile.marks = c(1, 0,0,0), 
                                                       mod.layers = c("H3K4me_demethylate", "H3K4me1", "H3K4me2", "H3K4me3", "H3K4me_any"),
                                                       mod.marks = c(0, 0,0,0,0))


# ========== H3K9 METHYLATION FACTORS ==========

bf_K9_meUp_sampler <- createBindingFactor.layer_region(name="bf_K9_meUp_sampler", 
                                                       profile.layers="nucleosome", 
                                                       patternLength = nucleosomeWidth,
                                                       profile.marks = 1,
                                                       mod.layers="sampled_K9_meUp",
                                                       mod.marks=1,
                                                       stateWidth= nucleosomeWidth)

# K9 requires tf_K9_active only (no antagonism — dominance hierarchy)
bf_K9_meUp <- createBindingFactor.layer_region(name="bf_K9_meUp", 
                                               profile.layers=c("tf_K9_active", "sampled_K9_meUp"), 
                                               patternLength = 6,
                                               profile.marks = c(1, 1),
                                               mod.layers="H3K9me_promotion",
                                               mod.marks=1,
                                               stateWidth= 6 + (2*nucleosomeWidth))

# Non-persistent promotion (Dave's pattern)
bf_K9_promotion_0_1 <- createBindingFactor.layer_region(name="bf_K9_promotion_0_1", 
                                                        patternLength = nucleosomeWidth, 
                                                        profile.layers = c("nucleosome", "H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                        profile.marks = c(1,1,0,0,0), 
                                                        mod.layers = c("H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                        mod.marks = c(0, 1,0,0,1))

bf_K9_promotion_1_2 <- createBindingFactor.layer_region(name="bf_K9_promotion_1_2", 
                                                        patternLength = nucleosomeWidth, 
                                                        profile.layers = c("nucleosome", "H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                        profile.marks = c(1,1,1,0,0), 
                                                        mod.layers = c("H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                        mod.marks = c(0, 0,1,0,1))

bf_K9_promotion_2_3 <- createBindingFactor.layer_region(name="bf_K9_promotion_2_3", 
                                                        patternLength = nucleosomeWidth, 
                                                        profile.layers = c("nucleosome", "H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                        profile.marks = c(1,1,0,1,0), 
                                                        mod.layers = c("H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                        mod.marks = c(0, 0,0,1,1))

bf_K9_promotion_3_3 <- createBindingFactor.layer_region(name="bf_K9_promotion_3_3", 
                                                        patternLength = nucleosomeWidth, 
                                                        profile.layers = c("nucleosome", "H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                        profile.marks = c(1,1,0,0,1), 
                                                        mod.layers = c("H3K9me_promotion", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                        mod.marks = c(0, 0,0,1,1))


# ========== H3K9 DEMETHYLATION (Dave's pattern) ==========

bf_K9_meDown <- createBindingFactor.layer_region(name="bf_K9_meDown", 
                                                 profile.layers="nucleosome",
                                                 patternLength = nucleosomeWidth,
                                                 profile.marks = 1,
                                                 mod.layers="H3K9me_demethylate",
                                                 mod.marks=1,
                                                 stateWidth=nucleosomeWidth)

bf_K9_demethylate_3_2 <- createBindingFactor.layer_region(name="bf_K9_demethylate_3_2", 
                                                          patternLength = nucleosomeWidth, 
                                                          profile.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                          profile.marks = c(1, 0,0,1), 
                                                          mod.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                          mod.marks = c(0, 0,1,0,1))

bf_K9_demethylate_2_1 <- createBindingFactor.layer_region(name="bf_K9_demethylate_2_1", 
                                                          patternLength = nucleosomeWidth, 
                                                          profile.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                          profile.marks = c(1, 0,1,0), 
                                                          mod.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                          mod.marks = c(0, 1,0,0,1))

bf_K9_demethylate_1_0 <- createBindingFactor.layer_region(name="bf_K9_demethylate_1_0", 
                                                          patternLength = nucleosomeWidth, 
                                                          profile.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                          profile.marks = c(1, 1,0,0), 
                                                          mod.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                          mod.marks = c(0, 0,0,0,0))

bf_K9_demethylate_0_0 <- createBindingFactor.layer_region(name="bf_K9_demethylate_0_0", 
                                                          patternLength = nucleosomeWidth, 
                                                          profile.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3"),
                                                          profile.marks = c(1, 0,0,0), 
                                                          mod.layers = c("H3K9me_demethylate", "H3K9me1", "H3K9me2", "H3K9me3", "H3K9me_any"),
                                                          mod.marks = c(0, 0,0,0,0))


# ========== H3K27 METHYLATION FACTORS ==========

bf_K27_meUp_sampler <- createBindingFactor.layer_region(name="bf_K27_meUp_sampler", 
                                                        profile.layers="nucleosome", 
                                                        patternLength = nucleosomeWidth,
                                                        profile.marks = 1,
                                                        mod.layers="sampled_K27_meUp",
                                                        mod.marks=1,
                                                        stateWidth= nucleosomeWidth)

# Rule 1: K4me3 blocks K27
bf_K27_meUp <- createBindingFactor.layer_region(name="bf_K27_meUp", 
                                                profile.layers=c("tf_K27_active", "sampled_K27_meUp"), 
                                                patternLength = 6,
                                                profile.marks = c(1, 1),
                                                mod.layers="H3K27me_promotion",
                                                mod.marks=1,
                                                stateWidth= 6 + (2*nucleosomeWidth))

# Rule 3: K27me3 spreading via PRC2 EED
bf_K27_meUp_spread <- createBindingFactor.layer_region(name="bf_K27_meUp_spread", 
                                                       profile.layers=c("nucleosome", "H3K27me_any", "H3K4me3"),
                                                       patternLength = nucleosomeWidth,
                                                       profile.marks = c(1, 1, 0),
                                                       mod.layers=c("H3K27me_promotion"),
                                                       mod.marks=c(1),
                                                       stateWidth= 3*nucleosomeWidth)

# Non-persistent promotion (Dave's pattern)
bf_K27_promotion_0_1 <- createBindingFactor.layer_region(name="bf_K27_promotion_0_1", 
                                                         patternLength = nucleosomeWidth, 
                                                         profile.layers = c("nucleosome", "H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                         profile.marks = c(1,1,0,0,0), 
                                                         mod.layers = c("H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                         mod.marks = c(0, 1,0,0,1))

bf_K27_promotion_1_2 <- createBindingFactor.layer_region(name="bf_K27_promotion_1_2", 
                                                         patternLength = nucleosomeWidth, 
                                                         profile.layers = c("nucleosome", "H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                         profile.marks = c(1,1,1,0,0), 
                                                         mod.layers = c("H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                         mod.marks = c(0, 0,1,0,1))

# Rule 1: K4me3 blocks K27 at me2→me3
bf_K27_promotion_2_3 <- createBindingFactor.layer_region(name="bf_K27_promotion_2_3", 
                                                         patternLength = nucleosomeWidth, 
                                                         profile.layers = c("nucleosome", "H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                         profile.marks = c(1,1,0,1,0), 
                                                         mod.layers = c("H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                         mod.marks = c(0, 0,0,1,1))

bf_K27_promotion_3_3 <- createBindingFactor.layer_region(name="bf_K27_promotion_3_3", 
                                                         patternLength = nucleosomeWidth, 
                                                         profile.layers = c("nucleosome", "H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                         profile.marks = c(1,1,0,0,1), 
                                                         mod.layers = c("H3K27me_promotion", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                         mod.marks = c(0, 0,0,1,1))


# ========== H3K27 DEMETHYLATION ==========

bf_K27_meDown <- createBindingFactor.layer_region(name="bf_K27_meDown", 
                                                  profile.layers="nucleosome",
                                                  patternLength = nucleosomeWidth,
                                                  profile.marks = 1,
                                                  mod.layers="H3K27me_demethylate",
                                                  mod.marks=1,
                                                  stateWidth=nucleosomeWidth)

bf_K27_demethylate_3_2 <- createBindingFactor.layer_region(name="bf_K27_demethylate_3_2", 
                                                           patternLength = nucleosomeWidth, 
                                                           profile.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                           profile.marks = c(1, 0,0,1),
                                                           mod.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                           mod.marks = c(0, 0,1,0,1))

bf_K27_demethylate_2_1 <- createBindingFactor.layer_region(name="bf_K27_demethylate_2_1", 
                                                           patternLength = nucleosomeWidth, 
                                                           profile.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                           profile.marks = c(1, 0,1,0),
                                                           mod.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                           mod.marks = c(0, 1,0,0,1))

bf_K27_demethylate_1_0 <- createBindingFactor.layer_region(name="bf_K27_demethylate_1_0", 
                                                           patternLength = nucleosomeWidth, 
                                                           profile.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                           profile.marks = c(1, 1,0,0),
                                                           mod.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                           mod.marks = c(0, 0,0,0,0))

bf_K27_demethylate_0_0 <- createBindingFactor.layer_region(name="bf_K27_demethylate_0_0", 
                                                           patternLength = nucleosomeWidth, 
                                                           profile.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                           profile.marks = c(1, 0,0,0), 
                                                           mod.layers = c("H3K27me_demethylate", "H3K27me1", "H3K27me2", "H3K27me3", "H3K27me_any"),
                                                           mod.marks = c(0, 0,0,0,0))


# FACTOR LISTS ------------------------------------------------------------

bfList_H3K4_methylate <- list(bf_meUp_sampler=bf_meUp_sampler, bf_meUp=bf_meUp, 
                              bf_promotion_0_1=bf_promotion_0_1, bf_promotion_1_2=bf_promotion_1_2, 
                              bf_promotion_2_3=bf_promotion_2_3, bf_promotion_3_3=bf_promotion_3_3)

bfList_H3K4_demethyl <- list(bf_meDown=bf_meDown, bf_demethylate_3_2=bf_demethylate_3_2, 
                             bf_demethylate_2_1=bf_demethylate_2_1, bf_demethylate_1_0=bf_demethylate_1_0, 
                             bf_demethylate_0=bf_demethylate_0_0)

bfList_H3K9_methylate <- list(bf_K9_meUp_sampler=bf_K9_meUp_sampler, bf_K9_meUp=bf_K9_meUp, 
                              bf_K9_promotion_0_1=bf_K9_promotion_0_1, bf_K9_promotion_1_2=bf_K9_promotion_1_2, 
                              bf_K9_promotion_2_3=bf_K9_promotion_2_3, bf_K9_promotion_3_3=bf_K9_promotion_3_3)

bfList_H3K9_demethyl <- list(bf_K9_meDown=bf_K9_meDown, bf_K9_demethylate_3_2=bf_K9_demethylate_3_2, 
                             bf_K9_demethylate_2_1=bf_K9_demethylate_2_1, bf_K9_demethylate_1_0=bf_K9_demethylate_1_0, 
                             bf_K9_demethylate_0=bf_K9_demethylate_0_0)

bfList_H3K27_methylate <- list(bf_K27_meUp_sampler=bf_K27_meUp_sampler, bf_K27_meUp=bf_K27_meUp, 
                               bf_K27_meUp_spread=bf_K27_meUp_spread,
                               bf_K27_promotion_0_1=bf_K27_promotion_0_1, bf_K27_promotion_1_2=bf_K27_promotion_1_2, 
                               bf_K27_promotion_2_3=bf_K27_promotion_2_3, bf_K27_promotion_3_3=bf_K27_promotion_3_3)

bfList_H3K27_demethyl <- list(bf_K27_meDown=bf_K27_meDown, bf_K27_demethylate_3_2=bf_K27_demethylate_3_2,
                              bf_K27_demethylate_2_1=bf_K27_demethylate_2_1, bf_K27_demethylate_1_0=bf_K27_demethylate_1_0,
                              bf_K27_demethylate_0=bf_K27_demethylate_0_0)


# ABUNDANCES (Dave's pattern) ---------------------------------------------

saturationAbundance <- 1000000

bf_methylate_abu <- c(12000, rep(saturationAbundance, 5))
names(bf_methylate_abu) <- names(bfList_H3K4_methylate)

bf_K9_methylate_abu <- c(12000, rep(saturationAbundance, 5))
names(bf_K9_methylate_abu) <- names(bfList_H3K9_methylate)

bf_K27_methylate_abu <- c(params$meUp_sampler_K27, saturationAbundance, params$K27_spread_abu, rep(saturationAbundance, 4))
names(bf_K27_methylate_abu) <- names(bfList_H3K27_methylate)

bf_demethylate_abu <- c(0, rep(saturationAbundance, 4))
names(bf_demethylate_abu) <- names(bfList_H3K4_demethyl)

bf_K9_demethylate_abu <- c(0, rep(saturationAbundance, 4))
names(bf_K9_demethylate_abu) <- names(bfList_H3K9_demethyl)

bf_K27_demethylate_abu <- c(0, rep(saturationAbundance, 4))
names(bf_K27_demethylate_abu) <- names(bfList_H3K27_demethyl)


# CREATE LAYERSET ---------------------------------------------------------

scLayerSetNuc <- createLayerSet.BSgenome(genome=genome, 
                                         layer.names=c("sampled_meUp", "sampled_meDown",
                                                       "sampled_K9_meUp", "sampled_K27_meUp",
                                                       "nucleosome", 
                                                       "tf_K4_active", "tf_K9_active", "tf_K27_active",
                                                       "H3K4me_promotion", "H3K4me_any", "H3K4me_demethylate",
                                                       "H3K4me1", "H3K4me2", "H3K4me3",
                                                       "H3K9me_promotion", "H3K9me_any", "H3K9me_demethylate",
                                                       "H3K9me1", "H3K9me2", "H3K9me3",
                                                       "H3K27me_promotion", "H3K27me_any", "H3K27me_demethylate",
                                                       "H3K27me1", "H3K27me2", "H3K27me3"), 
                                         n.layers=26, verbose=TRUE)


# POSITION NUCLEOSOMES ----------------------------------------------------

scLayerSetNuc$layerSet[["nucleosome"]] <- randGrangesBigGenome(genome=genome, 
                                                               sizeFunc = function(value, n) rep(x=value, times=n), 
                                                               gapFunc = function(n, value) rpois(n = n, lambda = value), 
                                                               argsSizeFunc = list(value=147, n=5),
                                                               argsGapFunc = list(value=40, n=7))

save(scLayerSetNuc, file = paste0(outputDir, "scLayerSetNuc.object.Rdata"))
total_nucs <- length(scLayerSetNuc$layerSet$nucleosome)
cat(sprintf("Total nucleosomes: %d\n", total_nucs))


# SEED tf_active LAYERS WITH DNA MOTIFS -----------------------------------

# K4: Scan genome for TATA boxes, mark surrounding regions as tf_K4_active
cat("\n=== Scanning for TATA box (TATATAAA) ===\n")
tf_K4_factorSet <- list(bf_TF_K4=bf_TF_K4)
tf_K4_abu <- c(saturationAbundance)
names(tf_K4_abu) <- names(tf_K4_factorSet)

scLayerSetNuc <- runLayerBinding.BSgenome(layerList = scLayerSetNuc,
                                          factorSet = tf_K4_factorSet,
                                          bf.abundances = tf_K4_abu,
                                          verbose = params$verbose)

K4_cov <- sum(width(reduce(scLayerSetNuc$layerSet$tf_K4_active))) / sum(as.numeric(seqlengths(genome))) * 100
cat(sprintf("tf_K4_active coverage: %.2f%%\n", K4_cov))
cat(sprintf("tf_K4_active regions: %d\n", length(scLayerSetNuc$layerSet$tf_K4_active)))


cat("\n=== Scanning for subtelomere-enriched TG-repeat (GTGTGGGT) ===\n")
tf_K9_factorSet <- list(bf_TF_K9=bf_TF_K9)
tf_K9_abu <- c(saturationAbundance)
names(tf_K9_abu) <- names(tf_K9_factorSet)

scLayerSetNuc <- runLayerBinding.BSgenome(layerList = scLayerSetNuc,
                                          factorSet = tf_K9_factorSet,
                                          bf.abundances = tf_K9_abu,
                                          verbose = params$verbose)

K9_cov <- sum(width(reduce(scLayerSetNuc$layerSet$tf_K9_active))) / sum(as.numeric(seqlengths(genome))) * 100
cat(sprintf("tf_K9_active motif coverage: %.2f%%\n", K9_cov))
cat(sprintf("tf_K9_active regions: %d\n", length(scLayerSetNuc$layerSet$tf_K9_active)))

# Random fill to reach 62% target
cat("\n--- K9: Random fill to reach target pool size ---\n")
nuc_in_K9 <- unique(queryHits(findOverlaps(scLayerSetNuc$layerSet$nucleosome, 
                                           scLayerSetNuc$layerSet$tf_K9_active)))
cat(sprintf("Nucleosomes in motif-seeded K9 regions: %d\n", length(nuc_in_K9)))

target_K9_nucs <- round(params$tf_K9_target_frac * total_nucs)
remaining_needed <- target_K9_nucs - length(nuc_in_K9)

if(remaining_needed > 0) {
  nuc_not_in_K9 <- setdiff(1:total_nucs, nuc_in_K9)
  random_fill <- sample(nuc_not_in_K9, min(remaining_needed, length(nuc_not_in_K9)), replace=FALSE)
  random_K9_regions <- scLayerSetNuc$layerSet$nucleosome[random_fill]
  scLayerSetNuc$layerSet$tf_K9_active <- c(scLayerSetNuc$layerSet$tf_K9_active, random_K9_regions)
  cat(sprintf("Random fill nucleosomes added: %d\n", length(random_fill)))
}

K9_cov <- sum(width(reduce(scLayerSetNuc$layerSet$tf_K9_active))) / sum(as.numeric(seqlengths(genome))) * 100
cat(sprintf("Final tf_K9_active coverage: %.2f%% (target 62%%)\n", K9_cov))


# K27: Random subset (yeast lacks PRC2 — no biological motif available)
cat("\n=== K27: Random targeting (no PRC2 in yeast) ===\n")
nK27 <- round(params$tf_K27_target_frac * total_nucs)
K27_nucs <- sample(1:total_nucs, nK27, replace=FALSE)
scLayerSetNuc$layerSet$tf_K27_active <- scLayerSetNuc$layerSet$nucleosome[K27_nucs]

K27_cov <- sum(width(reduce(scLayerSetNuc$layerSet$tf_K27_active))) / sum(as.numeric(seqlengths(genome))) * 100
cat(sprintf("tf_K27_active coverage: %.2f%%\n", K27_cov))
cat(sprintf("tf_K27_active sites: %d\n", length(scLayerSetNuc$layerSet$tf_K27_active)))

# Report overlap between tf pools
cat("\n=== tf_active pool overlaps ===\n")
K4_K9_overlap <- length(findOverlaps(scLayerSetNuc$layerSet$tf_K4_active, 
                                     scLayerSetNuc$layerSet$tf_K9_active))
cat(sprintf("tf_K4 ∩ tf_K9 overlapping regions: %d (Rule 2 test sites)\n", K4_K9_overlap))

K4_K27_overlap <- length(findOverlaps(scLayerSetNuc$layerSet$tf_K4_active, 
                                      scLayerSetNuc$layerSet$tf_K27_active))
cat(sprintf("tf_K4 ∩ tf_K27 overlapping regions: %d (Rule 1 test sites)\n", K4_K27_overlap))

K9_K27_overlap <- length(findOverlaps(scLayerSetNuc$layerSet$tf_K9_active, 
                                      scLayerSetNuc$layerSet$tf_K27_active))
cat(sprintf("tf_K9 ∩ tf_K27 overlapping regions: %d (co-occurrence sites)\n", K9_K27_overlap))

# IMPORTANT: Check if K4/K9 coverages match expected pool sizes
# If TATA gives too few hits, increase stateWidth in bf_TF_K4
# If telomeric gives too few hits, increase stateWidth in bf_TF_K9 or use shorter motif
# Target: K4 ~24%, K9 ~62%, K27 ~31%
cat(sprintf("\n=== COVERAGE CHECK ===\n"))
cat(sprintf("K4 target: ~24%%, actual: %.1f%% %s\n", K4_cov, ifelse(K4_cov < 15, "⚠ TOO LOW - increase stateWidth", ifelse(K4_cov > 40, "⚠ TOO HIGH", "✓"))))
cat(sprintf("K9 target: ~62%%, actual: %.1f%% %s\n", K9_cov, ifelse(K9_cov < 40, "⚠ TOO LOW - increase stateWidth or use shorter motif", ifelse(K9_cov > 80, "⚠ TOO HIGH", "✓"))))
cat(sprintf("K27 target: ~31%%, actual: %.1f%% %s\n", K27_cov, ifelse(K27_cov < 20, "⚠ TOO LOW", ifelse(K27_cov > 45, "⚠ TOO HIGH", "✓"))))


# TRACKING

total_nucleosomes <- total_nucs

prev_states <- data.frame(K4 = rep("me0", total_nucleosomes),
                          K9 = rep("me0", total_nucleosomes),
                          K27 = rep("me0", total_nucleosomes),
                          stringsAsFactors = FALSE)
current_states <- prev_states

dwell_K4 <- rep(0, total_nucleosomes)
dwell_K9 <- rep(0, total_nucleosomes)
dwell_K27 <- rep(0, total_nucleosomes)

collected_dwell_K4 <- list(me0=c(), me1=c(), me2=c(), me3=c())
collected_dwell_K9 <- list(me0=c(), me1=c(), me2=c(), me3=c())
collected_dwell_K27 <- list(me0=c(), me1=c(), me2=c(), me3=c())

transition_K4 <- matrix(0, nrow=4, ncol=4, dimnames=list(c("me0","me1","me2","me3"), c("me0","me1","me2","me3")))
transition_K9 <- matrix(0, nrow=4, ncol=4, dimnames=list(c("me0","me1","me2","me3"), c("me0","me1","me2","me3")))
transition_K27 <- matrix(0, nrow=4, ncol=4, dimnames=list(c("me0","me1","me2","me3"), c("me0","me1","me2","me3")))

get_all_states_fast <- function(layerSet) {
  K4_states <- rep("me0", total_nucleosomes)
  K9_states <- rep("me0", total_nucleosomes)
  K27_states <- rep("me0", total_nucleosomes)
  nuc_ranges <- layerSet$nucleosome
  if(length(layerSet$H3K4me3) > 0) K4_states[unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K4me3)))] <- "me3"
  if(length(layerSet$H3K4me2) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K4me2))); K4_states[h[K4_states[h]=="me0"]] <- "me2" }
  if(length(layerSet$H3K4me1) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K4me1))); K4_states[h[K4_states[h]=="me0"]] <- "me1" }
  if(length(layerSet$H3K9me3) > 0) K9_states[unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K9me3)))] <- "me3"
  if(length(layerSet$H3K9me2) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K9me2))); K9_states[h[K9_states[h]=="me0"]] <- "me2" }
  if(length(layerSet$H3K9me1) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K9me1))); K9_states[h[K9_states[h]=="me0"]] <- "me1" }
  if(length(layerSet$H3K27me3) > 0) K27_states[unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K27me3)))] <- "me3"
  if(length(layerSet$H3K27me2) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K27me2))); K27_states[h[K27_states[h]=="me0"]] <- "me2" }
  if(length(layerSet$H3K27me1) > 0) { h <- unique(queryHits(findOverlaps(nuc_ranges, layerSet$H3K27me1))); K27_states[h[K27_states[h]=="me0"]] <- "me1" }
  return(data.frame(K4=K4_states, K9=K9_states, K27=K27_states, stringsAsFactors=FALSE))
}

update_tracking_fast <- function(iteration) {
  K4_changed <- which(current_states$K4 != prev_states$K4)
  K9_changed <- which(current_states$K9 != prev_states$K9)
  K27_changed <- which(current_states$K27 != prev_states$K27)
  if(length(K4_changed) > 0) { for(idx in K4_changed) { from <- prev_states$K4[idx]; to <- current_states$K4[idx]; transition_K4[from, to] <<- transition_K4[from, to] + 1; collected_dwell_K4[[from]] <<- c(collected_dwell_K4[[from]], dwell_K4[idx]); dwell_K4[idx] <<- 0 } }
  if(length(K9_changed) > 0) { for(idx in K9_changed) { from <- prev_states$K9[idx]; to <- current_states$K9[idx]; transition_K9[from, to] <<- transition_K9[from, to] + 1; collected_dwell_K9[[from]] <<- c(collected_dwell_K9[[from]], dwell_K9[idx]); dwell_K9[idx] <<- 0 } }
  if(length(K27_changed) > 0) { for(idx in K27_changed) { from <- prev_states$K27[idx]; to <- current_states$K27[idx]; transition_K27[from, to] <<- transition_K27[from, to] + 1; collected_dwell_K27[[from]] <<- c(collected_dwell_K27[[from]], dwell_K27[idx]); dwell_K27[idx] <<- 0 } }
  dwell_K4 <<- dwell_K4 + 1; dwell_K9 <<- dwell_K9 + 1; dwell_K27 <<- dwell_K27 + 1
  prev_states <<- current_states
}

# SIMULATION LOOP ---------------------------------------------------------

scLayerSetBoth <- scLayerSetNuc

bfLIst_both <- c(bfList_H3K4_methylate, bfList_H3K4_demethyl, 
                 bfList_H3K9_methylate, bfList_H3K9_demethyl,
                 bfList_H3K27_methylate, bfList_H3K27_demethyl)
bfAbund_both <- c(bf_methylate_abu, bf_demethylate_abu, 
                  bf_K9_methylate_abu, bf_K9_demethylate_abu,
                  bf_K27_methylate_abu, bf_K27_demethylate_abu)

# cleanup pattern
layersToClean <- c("sampled_meUp", "sampled_K9_meUp", "sampled_K27_meUp",
                   'H3K4me_promotion', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K4me_any', 'H3K4me_demethylate',
                   'H3K9me_promotion', 'H3K9me1', 'H3K9me2', 'H3K9me3', 'H3K9me_any', 'H3K9me_demethylate',
                   'H3K27me_promotion', 'H3K27me1', 'H3K27me2', 'H3K27me3', 'H3K27me_any', 'H3K27me_demethylate')

clean <- TRUE
abundSpread <- c(1000)

bfAbund_both_eFirst <- c(bf_demethylate_abu, bf_methylate_abu,
                         bf_K9_demethylate_abu, bf_K9_methylate_abu,
                         bf_K27_demethylate_abu, bf_K27_methylate_abu)
bfLIst_both_eFirst <- bfLIst_both[names(bfAbund_both_eFirst)]

for(thisAbund in abundSpread) {
  scLayerSetBothAbund <- scLayerSetBoth
  
  bfAbund_both_eFirst["bf_meUp_sampler"] <- params$meUp_sampler_K4
  bfAbund_both_eFirst["bf_meDown"] <- params$meDown_K4
  bfAbund_both_eFirst["bf_K9_meUp_sampler"] <- params$meUp_sampler_K9
  bfAbund_both_eFirst["bf_K9_meDown"] <- params$meDown_K9
  bfAbund_both_eFirst["bf_K27_meUp_sampler"] <- params$meUp_sampler_K27
  bfAbund_both_eFirst["bf_K27_meDown"] <- params$meDown_K27

  subSimName <- if (nzchar(params$sim_tag)) params$sim_tag else paste0("Model2B_motif.w", bfAbund_both_eFirst["bf_meUp_sampler"], ".e", thisAbund)
  print(subSimName)

  for(i in 1:params$n_iter) {
    if (params$verbose) print(i)
    scLayerSetBothAbund <- runLayerBinding.BSgenome(layerList = scLayerSetBothAbund,
                                                    factorSet = bfLIst_both_eFirst,
                                                    bf.abundances = bfAbund_both_eFirst,
                                                    verbose = params$verbose, collect.stats = T, keep.stats = T)
    
    current_states <- get_all_states_fast(scLayerSetBothAbund$layerSet)
    if(i > 1) update_tracking_fast(i)
    
    if(clean) {
      for(thisLayer in layersToClean) {
        scLayerSetBothAbund$layerSet[[thisLayer]] <- removeShortGRanges(x=scLayerSetBothAbund$layerSet[[thisLayer]], minSize = nucleosomeWidth)
      }
      
      # no rule enforcement needed for rule 3 script only.
      
      scLayerSetBothAbund$layerSet[["nucleosome"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["nucleosome"]], verbose=T, minSize=nucleosomeWidth, maxSize=nucleosomeWidth)
      scLayerSetBothAbund$layerSet[["sampled_meUp"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["sampled_meUp"]], verbose=T, maxSize=0)
      scLayerSetBothAbund$layerSet[["sampled_K9_meUp"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["sampled_K9_meUp"]], verbose=T, maxSize=0)
      scLayerSetBothAbund$layerSet[["sampled_K27_meUp"]] <- removeGRangesBySize(x=scLayerSetBothAbund$layerSet[["sampled_K27_meUp"]], verbose=T, maxSize=0)
    }
    
    if(params$verbose && i %% 10 == 0) {
      ymax <- max(
        scLayerSetBothAbund$history$nBlocks.H3K4me3,
        scLayerSetBothAbund$history$nBlocks.H3K4me2,
        scLayerSetBothAbund$history$nBlocks.H3K4me1,
        scLayerSetBothAbund$history$nBlocks.H3K9me3,
        scLayerSetBothAbund$history$nBlocks.H3K9me2,
        scLayerSetBothAbund$history$nBlocks.H3K9me1,
        scLayerSetBothAbund$history$nBlocks.H3K27me3,
        scLayerSetBothAbund$history$nBlocks.H3K27me2,
        scLayerSetBothAbund$history$nBlocks.H3K27me1,
        na.rm=TRUE) * 1.15
      
      plot(scLayerSetBothAbund$history$nBlocks.H3K4me3, pch=6, col="black", ylab="nucleosomes", xlab="iteration", main=paste("iteration", i), ylim=c(0,ymax))
      points(scLayerSetBothAbund$history$nBlocks.H3K4me2, pch=4, col="black")
      points(scLayerSetBothAbund$history$nBlocks.H3K4me1, pch=1, col="black")
      points(scLayerSetBothAbund$history$nBlocks.H3K9me3, pch=6, col="red")
      points(scLayerSetBothAbund$history$nBlocks.H3K9me2, pch=4, col="red")
      points(scLayerSetBothAbund$history$nBlocks.H3K9me1, pch=1, col="red")
      points(scLayerSetBothAbund$history$nBlocks.H3K27me3, pch=6, col="blue")
      points(scLayerSetBothAbund$history$nBlocks.H3K27me2, pch=4, col="blue")
      points(scLayerSetBothAbund$history$nBlocks.H3K27me1, pch=1, col="blue")
      legend("topright", legend=c("H3K4me3","H3K4me2","H3K4me1","H3K9me3","H3K9me2","H3K9me1","H3K27me3","H3K27me2","H3K27me1"), 
             pch=c(6,4,1,6,4,1,6,4,1), col=c("black","black","black","red","red","red","blue","blue","blue"), ncol=3, cex=0.5, bty="n")
    }
  } # end iteration loop
  
  # STATISTICS 
  
  cat("\n=== FINAL STATISTICS FOR", subSimName, "===\n")
  
  # 1. STATE COUNTING
  final_states <- table(K4 = current_states$K4, K9 = current_states$K9, K27 = current_states$K27)
  cat("\nFinal State Distribution:\n")
  print(final_states)
  
  bivalent <- sum(current_states$K4 == "me3" & current_states$K27 == "me3")
  K9_K27_coop <- sum(current_states$K9 == "me3" & current_states$K27 == "me3")
  triple_marked <- sum(current_states$K4 == "me3" & current_states$K9 == "me3" & current_states$K27 == "me3")
  
  cat(sprintf("\nBivalent (K4me3+K27me3): %d (%.2f%%)\n", bivalent, 100*bivalent/total_nucleosomes))
  cat(sprintf("K9me3+K27me3 cooperation: %d (%.2f%%)\n", K9_K27_coop, 100*K9_K27_coop/total_nucleosomes))
  cat(sprintf("Triple marked: %d (%.2f%%)\n\n", triple_marked, 100*triple_marked/total_nucleosomes))
  
  # 2. DWELL TIMES
  cat("=== Mean Dwell Times (iterations) ===\n")
  if(length(collected_dwell_K4$me1) > 0) {
    cat(sprintf("H3K4me1: %.2f, me2: %.2f, me3: %.2f\n", 
                mean(collected_dwell_K4$me1), mean(collected_dwell_K4$me2), mean(collected_dwell_K4$me3)))
  }
  if(length(collected_dwell_K9$me1) > 0) {
    cat(sprintf("H3K9me1: %.2f, me2: %.2f, me3: %.2f\n", 
                mean(collected_dwell_K9$me1), mean(collected_dwell_K9$me2), mean(collected_dwell_K9$me3)))
  }
  if(length(collected_dwell_K27$me1) > 0) {
    cat(sprintf("H3K27me1: %.2f, me2: %.2f, me3: %.2f\n", 
                mean(collected_dwell_K27$me1), mean(collected_dwell_K27$me2), mean(collected_dwell_K27$me3)))
  }
  
  # 3. GENOME COVERAGE
  genome_size <- sum(as.numeric(seqlengths(genome)))
  calc_coverage <- function(layer_name) {
    if(length(scLayerSetBothAbund$layerSet[[layer_name]]) == 0) return(0)
    return(100 * sum(width(reduce(scLayerSetBothAbund$layerSet[[layer_name]]))) / genome_size)
  }
  
  cat("\n=== Genome Coverage ===\n")
  K4me3_cov <- calc_coverage("H3K4me3")
  K9me3_cov <- calc_coverage("H3K9me3")
  K27me3_cov <- calc_coverage("H3K27me3")
  cat(sprintf("H3K4me3: %.2f%% (target: 1.2%%)\n", K4me3_cov))
  cat(sprintf("H3K9me3: %.2f%% (target: 21.7%%)\n", K9me3_cov))
  cat(sprintf("H3K27me3: %.2f%% (target: 14.4%%)\n", K27me3_cov))
  K4me1_cov <- calc_coverage("H3K4me1")
  K4me2_cov <- calc_coverage("H3K4me2")
  K9me1_cov <- calc_coverage("H3K9me1")
  K9me2_cov <- calc_coverage("H3K9me2")
  K27me1_cov <- calc_coverage("H3K27me1")
  K27me2_cov <- calc_coverage("H3K27me2")
  
  cat(sprintf("H3K4me1: %.2f%% (target: 12.8%%)\n", K4me1_cov))
  cat(sprintf("H3K4me2: %.2f%% (target: 5.7%%)\n", K4me2_cov))
  cat(sprintf("H3K9me1: %.2f%% (target: 11.1%%)\n", K9me1_cov))
  cat(sprintf("H3K9me2: %.2f%% (target: 18.7%%)\n", K9me2_cov))
  cat(sprintf("H3K27me1: %.2f%% (target: 4.0%%)\n", K27me1_cov))
  cat(sprintf("H3K27me2: %.2f%% (target: 7.5%%)\n", K27me2_cov))
  
  # RMSE calculation
  
  targets <- c(12.8, 5.7, 1.2, 11.1, 18.7, 21.7, 4.0, 7.5, 14.4)
  reached <- c(K4me1_cov, K4me2_cov, K4me3_cov, K9me1_cov, K9me2_cov, K9me3_cov, K27me1_cov, K27me2_cov, K27me3_cov)
  RMSE_all <- sqrt(mean((reached - targets)^2))
  RMSE_me3 <-sqrt(mean((reached[c(3,6,9)]- targets[c(3,6,9)])^2))
  cat(sprintf("\nRMSE (all 9 states): %.3f\n", RMSE_all))
  cat(sprintf("RMSE (me3 only): %.3f\n", RMSE_me3))

  # Machine-readable objective for the Phase 1 scikit-optimize loop
  # (hand-built JSON to avoid a jsonlite dependency — structure is fixed and flat)
  mark_names <- c("H3K4me1","H3K4me2","H3K4me3","H3K9me1","H3K9me2","H3K9me3","H3K27me1","H3K27me2","H3K27me3")
  objective_path <- paste0(outputDir, subSimName, ".objective.json")
  writeLines(c(
    "{",
    sprintf('  "RMSE_all": %.6f,', RMSE_all),
    sprintf('  "RMSE_me3": %.6f,', RMSE_me3),
    sprintf('  "coverage": {%s}', paste(sprintf('"%s": %.6f', mark_names, reached), collapse=", ")),
    "}"
  ), objective_path)
  cat(sprintf("\nOBJECTIVE_JSON: %s\n", objective_path))

  # 4. CO-OCCURRENCE
  calc_overlap <- function(layer1, layer2) {
    if(length(scLayerSetBothAbund$layerSet[[layer1]]) == 0) return(0)
    if(length(scLayerSetBothAbund$layerSet[[layer2]]) == 0) return(0)
    overlaps <- findOverlaps(scLayerSetBothAbund$layerSet[[layer1]], 
                             scLayerSetBothAbund$layerSet[[layer2]])
    return(100 * length(unique(queryHits(overlaps))) / length(scLayerSetBothAbund$layerSet[[layer1]]))
  }
  
  cat("\n=== Peak Co-occurrence ===\n")
  K4_K27_overlap <- calc_overlap("H3K4me3", "H3K27me3")
  K9_K27_overlap <- calc_overlap("H3K9me3", "H3K27me3")
  cat(sprintf("H3K4me3 + H3K27me3: %.2f%% (RULE 1: should be <5%%)\n", K4_K27_overlap))
  cat(sprintf("H3K9me3 + H3K27me3: %.2f%% (empirical: 31.7%%)\n", K9_K27_overlap))
  
  # 5. RULE VALIDATIONS
  
  # RULE 2 VALIDATION: K9me3-K4me3 exclusion
  if(length(scLayerSetBothAbund$layerSet$H3K4me3) > 0 && length(scLayerSetBothAbund$layerSet$H3K9me3) > 0) {
    K9_K4_overlap <- calc_overlap("H3K9me3", "H3K4me3")
    
    cat(sprintf("\n=== RULE 2 TEST: K9me3-K4me3 Antagonism ===\n"))
    cat(sprintf("H3K9me3 + H3K4me3 co-occurrence: %.2f%% (should be <5%%)\n", K9_K4_overlap))
    
    if(K9_K4_overlap < 5) {
      cat("✓ RULE 2 WORKING - K9me3 excludes K4me3!\n")
    } else if(K9_K4_overlap < 10) {
      cat("~ RULE 2 PARTIAL - some exclusion but incomplete\n")
    } else {
      cat("✗ RULE 2 NOT WORKING - heterochromatin not blocking active marks\n")
    }
  }
  
  # RULE 3 VALIDATION: K27me3 domain formation via spreading
  if(length(scLayerSetBothAbund$layerSet$H3K27me3) > 0) {
    K27_blocks <- reduce(scLayerSetBothAbund$layerSet$H3K27me3, min.gapwidth = 100)
    K27_block_sizes <- width(K27_blocks)
    mean_block_size <- mean(K27_block_sizes)
    large_domains <- sum(K27_block_sizes > 3*nucleosomeWidth)
    
    cat(sprintf("\n=== RULE 3A TEST: K27me3 Domain Formation ===\n"))
    cat(sprintf("Total K27me3 blocks: %d\n", length(K27_blocks)))
    cat(sprintf("Mean block size: %.0f bp (%.1f nucleosomes)\n", mean_block_size, mean_block_size/nucleosomeWidth))
    cat(sprintf("Large domains (>3 nucleosomes): %d (%.1f%%)\n", large_domains, 100*large_domains/length(K27_blocks)))
    
    if(mean_block_size > 2*nucleosomeWidth) {
      cat("✓ RULE 3A WORKING - K27me3 forms multi-nucleosome domains!\n")
    } else if(mean_block_size > 1.5*nucleosomeWidth) {
      cat("~ RULE 3A PARTIAL - some spreading occurring\n")
    } else {
      cat("✗ RULE 3A NOT WORKING - no domain formation (isolated nucleosomes)\n")
    }
  }
  
  # 6. TRANSITION PROBABILITIES
  cat("\n=== Transition Probability Matrices ===\n")
  cat("H3K4me:\n")
  print(round(prop.table(transition_K4, margin=1), 3))
  cat("\nH3K9me:\n")
  print(round(prop.table(transition_K9, margin=1), 3))
  cat("\nH3K27me:\n")
  print(round(prop.table(transition_K27, margin=1), 3))
  
  # 7. RULE VALIDATION SUMMARY
  cat("\n=== RULE VALIDATION SUMMARY ===\n")
  if(K4_K27_overlap < 5) {
    cat("✓ RULE 1 PASSED: K4me3-K27me3 bivalency is low\n")
  } else {
    cat("✗ RULE 1 FAILED: bivalency too high\n")
  }
  
  if(exists("K9_K4_overlap") && K9_K4_overlap < 5) {
    cat("✓ RULE 2 PASSED: K9me3-K4me3 exclusion is strong\n")
  } else if(exists("K9_K4_overlap")) {
    cat("✗ RULE 2 FAILED: heterochromatin not excluding active chromatin\n")
  }
  
  cat(sprintf("\nK9me3-K27me3 co-occurrence: %.1f%% (empirical: 31.7%%)\n", K9_K27_overlap))
  if(K9_K27_overlap > 20 && K9_K27_overlap < 45) {
    cat("✓ Independent targeting produces realistic overlap\n")
  } else if(K9_K27_overlap < 10) {
    cat("⚠ Overlap lower than empirical - marks too segregated\n")
  } else {
    cat("⚠ Overlap higher than empirical - check samplers\n")
  }
  
  # 8. SAVE EVERYTHING
  write.csv(as.data.frame.table(final_states), file=paste0(outputDir, subSimName, ".state_counts.csv"))
  write.csv(current_states, file=paste0(outputDir, subSimName, ".final_states.csv"), row.names=FALSE)
  save(collected_dwell_K4, collected_dwell_K9, collected_dwell_K27, 
       file=paste0(outputDir, subSimName, ".dwell_times.Rdata"))
  save(transition_K4, transition_K9, transition_K27, 
       file=paste0(outputDir, subSimName, ".transitions.Rdata"))
  
  # Export coverage data
  coverage_data <- data.frame(
    Mark = c("H3K4me1", "H3K4me2", "H3K4me3", 
             "H3K9me1", "H3K9me2", "H3K9me3", 
             "H3K27me1", "H3K27me2", "H3K27me3"),
    Coverage_Percent = c(K4me1_cov, K4me2_cov, K4me3_cov, 
                         K9me1_cov, K9me2_cov, K9me3_cov, 
                         K27me1_cov, K27me2_cov, K27me3_cov),
    Target_Percent = c(12.8, 5.7, 1.2, 
                       11.1, 18.7, 21.7, 
                       4.0, 7.5, 14.4)
  )
  write.csv(coverage_data, file=paste0(outputDir, subSimName, ".coverage.csv"), row.names=FALSE)
  
  # Export co-occurrence data
  cooccurrence_data <- data.frame(
    Comparison = c("K4me3+K27me3", "K9me3+K4me3", "K9me3+K27me3"),
    Overlap_Percent = c(K4_K27_overlap, 
                        ifelse(exists("K9_K4_overlap"), K9_K4_overlap, NA), 
                        K9_K27_overlap),
    Target_Percent = c(5, 5, 31.7),
    Rule_Test = c("Rule 1: K4me3 blocks K27me3", 
                  "Rule 2: K9me3 blocks K4me3",
                  "Rule 3: K27 self-recruitment")
  )
  write.csv(cooccurrence_data, file=paste0(outputDir, subSimName, ".cooccurrence.csv"), row.names=FALSE)
  
  cat("\n=== All statistics saved ===\n")
  
  # ============================================================
  # SPATIAL VALIDATION - Model 2B DNA motif
  # Run this AFTER the simulation finishes
  # Requires: scLayerSetBothAbund (the final simulation output)
  # ============================================================
  
  # ============================================================
  # SPATIAL VALIDATION 
  # ============================================================
  
  library(GenomicRanges)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
  
  genome <- BSgenome.Scerevisiae.UCSC.sacCer3
  genome_size <- sum(as.numeric(seqlengths(genome)))
  nuc_ranges <- scLayerSetBothAbund$layerSet$nucleosome
  total_nucs <- length(nuc_ranges)
  
  cat("\n============================================================\n")
  cat("SPATIAL VALIDATION ANALYSIS\n")
  cat("============================================================\n")
  
  # ========== K4 VALIDATION 1: K4me3 depletion at K9 overlap sites ==========
  cat("\n=== K4 TEST 1: K4me3 depletion at tf_K4/tf_K9 overlap sites ===\n")
  cat("(Rule 2 emergent — tests whether K9me3 excludes K4me3 spatially)\n\n")
  
  tf_K4 <- scLayerSetBothAbund$layerSet$tf_K4_active
  tf_K9 <- scLayerSetBothAbund$layerSet$tf_K9_active
  K4me3_regions <- scLayerSetBothAbund$layerSet$H3K4me3
  
  overlap_hits <- findOverlaps(tf_K4, tf_K9)
  if(length(overlap_hits) > 0 && length(K4me3_regions) > 0) {
    overlap_regions <- tf_K4[unique(queryHits(overlap_hits))]
    non_overlap_K4 <- tf_K4[-unique(queryHits(overlap_hits))]
    
    K4me3_at_overlap <- length(unique(queryHits(findOverlaps(K4me3_regions, overlap_regions))))
    K4me3_at_K4only <- length(unique(queryHits(findOverlaps(K4me3_regions, non_overlap_K4))))
    
    K9me3_at_overlap <- length(unique(queryHits(findOverlaps(scLayerSetBothAbund$layerSet$H3K9me3, overlap_regions))))
    
    total_overlap <- length(overlap_regions)
    total_K4only <- length(non_overlap_K4)
    
    K4me3_rate_overlap <- K4me3_at_overlap / total_overlap * 100
    K4me3_rate_K4only <- K4me3_at_K4only / total_K4only * 100
    K9me3_rate_overlap <- K9me3_at_overlap / total_overlap * 100
    
    cat(sprintf("tf_K4/tf_K9 overlapping regions: %d\n", total_overlap))
    cat(sprintf("tf_K4-only regions: %d\n", total_K4only))
    cat(sprintf("K4me3 at overlap sites: %.1f%%\n", K4me3_rate_overlap))
    cat(sprintf("K4me3 at K4-only sites: %.1f%%\n", K4me3_rate_K4only))
    cat(sprintf("K9me3 at overlap sites: %.1f%%\n", K9me3_rate_overlap))
    
    if(total_K4only > 0 && K4me3_rate_K4only > 0 && K4me3_rate_overlap < K4me3_rate_K4only) {
      depletion_fold <- K4me3_rate_K4only / max(K4me3_rate_overlap, 0.01)
      cat(sprintf("✓ K4me3 DEPLETED at K9 overlap sites (%.1fx reduction)\n", depletion_fold))
      cat("  Rule 2 (K9me3 blocks K4me3) validated spatially\n")
    } else if(total_K4only == 0) {
      cat("No K4-only regions available — K9 pool covers all K4 sites\n")
      cat("K4/K9 overlap test not applicable at current pool sizes\n")
    } else {
      cat("✗ K4me3 not depleted at K9 overlap sites\n")
    }
  } else {
    cat("Insufficient data for overlap analysis\n")
  }
  
  # ========== K4 VALIDATION 2: K4me3 depletion at subtelomeric regions ==========
  cat("\n\n=== K4 TEST 2: K4me3 depletion at subtelomeric regions ===\n")
  cat("(Rule 2 emergent — K4me3 should be excluded from heterochromatic regions)\n\n")
  
  subtelo <- GRanges()
  for(chr in seqnames(genome)) {
    if(chr == "chrM") next
    chr_len <- seqlengths(genome)[chr]
    subtelo <- c(subtelo,
                 GRanges(chr, IRanges(1, min(20000, chr_len))),
                 GRanges(chr, IRanges(max(1, chr_len - 20000), chr_len)))
  }
  
  subtelo_coverage <- sum(width(reduce(subtelo))) / genome_size * 100
  
  if(length(K4me3_regions) > 0) {
    K4me3_at_subtelo <- length(unique(queryHits(findOverlaps(K4me3_regions, subtelo))))
    total_K4me3 <- length(K4me3_regions)
    
    observed_pct <- K4me3_at_subtelo / total_K4me3 * 100
    expected_pct <- subtelo_coverage
    fold_K4_subtelo <- observed_pct / expected_pct
    
    cat(sprintf("Subtelomeric regions: %.1f%% of genome\n", subtelo_coverage))
    cat(sprintf("K4me3 at subtelomeres: %d / %d (%.1f%%)\n", K4me3_at_subtelo, total_K4me3, observed_pct))
    cat(sprintf("Expected by chance: %.1f%%\n", expected_pct))
    cat(sprintf("Fold enrichment: %.2fx\n", fold_K4_subtelo))
    
    if(fold_K4_subtelo < 0.7) {
      cat("✓ K4me3 DEPLETED at subtelomeres (Rule 2 blocks K4 at heterochromatin)\n")
    } else if(fold_K4_subtelo < 0.95) {
      cat("~ K4me3 slightly depleted at subtelomeres\n")
    } else {
      cat("✗ K4me3 not depleted at subtelomeres\n")
    }
  } else {
    cat("No K4me3 regions found\n")
  }
  
  # ========== K4 VALIDATION 3: K4me3 at internal vs subtelomeric genes ==========
  cat("\n\n=== K4 TEST 3: K4me3 at internal vs subtelomeric genes ===\n")
  cat("(Emergent — model doesn't know expression levels)\n\n")
  
  txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
  genes <- genes(txdb)
  prom <- promoters(genes, upstream=500, downstream=100)
  prom <- trim(prom)
  
  chr_lens <- seqlengths(genome)
  gene_chrs <- as.character(seqnames(genes))
  gene_starts <- start(genes)
  gene_ends <- end(genes)
  
  dist_to_left <- gene_starts
  dist_to_right <- chr_lens[gene_chrs] - gene_ends
  gene_dist_to_end <- pmin(dist_to_left, dist_to_right)
  
  subtelo_genes <- which(gene_dist_to_end < 20000 & gene_chrs != "chrM")
  internal_genes <- which(gene_dist_to_end > 50000 & gene_chrs != "chrM")
  
  if(length(K4me3_regions) > 0) {
    subtelo_prom <- prom[subtelo_genes]
    internal_prom <- prom[internal_genes]
    
    K4me3_subtelo_genes <- length(unique(queryHits(findOverlaps(subtelo_prom, K4me3_regions))))
    K4me3_internal_genes <- length(unique(queryHits(findOverlaps(internal_prom, K4me3_regions))))
    
    subtelo_rate <- K4me3_subtelo_genes / length(subtelo_genes) * 100
    internal_rate <- K4me3_internal_genes / length(internal_genes) * 100
    
    cat(sprintf("Subtelomeric genes (<20kb from end): %d\n", length(subtelo_genes)))
    cat(sprintf("Internal genes (>50kb from end): %d\n", length(internal_genes)))
    cat(sprintf("K4me3 at subtelomeric gene promoters: %d (%.1f%%)\n", K4me3_subtelo_genes, subtelo_rate))
    cat(sprintf("K4me3 at internal gene promoters: %d (%.1f%%)\n", K4me3_internal_genes, internal_rate))
    
    if(internal_rate > subtelo_rate && subtelo_rate > 0) {
      cat(sprintf("✓ K4me3 enriched at internal genes (%.1fx vs subtelomeric)\n", internal_rate / subtelo_rate))
      cat("  Consistent with expression correlation — euchromatic genes retain K4me3\n")
    } else if(internal_rate > subtelo_rate) {
      cat("✓ K4me3 enriched at internal genes (subtelomeric rate = 0)\n")
    } else {
      cat("✗ No K4me3 enrichment at internal vs subtelomeric genes\n")
    }
    
    fisher_mat <- matrix(c(K4me3_subtelo_genes, length(subtelo_genes) - K4me3_subtelo_genes,
                           K4me3_internal_genes, length(internal_genes) - K4me3_internal_genes),
                         nrow=2)
    fisher_res <- fisher.test(fisher_mat)
    cat(sprintf("Fisher's exact test p-value: %e\n", fisher_res$p.value))
  } else {
    cat("No K4me3 regions found\n")
  }
  
  # ========== K9 VALIDATION 1: K9me3 enrichment at subtelomeric regions ==========
  cat("\n\n=== K9 TEST 1: K9me3 enrichment at subtelomeric regions ===\n\n")
  
  K9me3_regions <- scLayerSetBothAbund$layerSet$H3K9me3
  total_K9me3 <- length(K9me3_regions)
  
  if(total_K9me3 > 0) {
    K9me3_at_subtelo <- length(unique(queryHits(findOverlaps(K9me3_regions, subtelo))))
    
    observed_pct <- K9me3_at_subtelo / total_K9me3 * 100
    expected_pct <- subtelo_coverage
    fold_K9_subtelo <- observed_pct / expected_pct
    
    K9me3_not_telo <- total_K9me3 - K9me3_at_subtelo
    total_bins <- round(genome_size / nucleosomeWidth)
    telo_bins <- round(sum(width(reduce(subtelo))) / nucleosomeWidth)
    non_telo_bins <- total_bins - telo_bins
    
    fisher_matrix <- matrix(c(K9me3_at_subtelo, K9me3_not_telo,
                              telo_bins - K9me3_at_subtelo, non_telo_bins - K9me3_not_telo),
                            nrow=2)
    fisher_result <- fisher.test(fisher_matrix)
    
    cat(sprintf("K9me3 regions total: %d\n", total_K9me3))
    cat(sprintf("K9me3 at subtelomeres: %d (%.1f%%)\n", K9me3_at_subtelo, observed_pct))
    cat(sprintf("Expected by chance: %.1f%%\n", expected_pct))
    cat(sprintf("Fold enrichment: %.2fx\n", fold_K9_subtelo))
    cat(sprintf("Fisher's exact test p-value: %e\n", fisher_result$p.value))
    
    if(fold_K9_subtelo > 1.5) {
      cat("✓ K9me3 ENRICHED at subtelomeric regions\n")
    } else if(fold_K9_subtelo > 1.1) {
      cat("~ K9me3 slightly enriched at subtelomeric regions\n")
    } else {
      cat("✗ K9me3 not enriched at subtelomeres\n")
    }
  } else {
    cat("No K9me3 regions found\n")
  }
  
  # ========== K9 VALIDATION 2: K9me3 gradient from chromosome ends ==========
  cat("\n\n=== K9 TEST 2: K9me3 gradient from chromosome ends ===\n\n")
  
  K9me3_status <- rep(0, total_nucs)
  if(total_K9me3 > 0) {
    K9me3_nucs <- unique(queryHits(findOverlaps(nuc_ranges, K9me3_regions)))
    K9me3_status[K9me3_nucs] <- 1
  }
  
  chr_lens <- seqlengths(genome)
  nuc_chrs <- as.character(seqnames(nuc_ranges))
  nuc_starts <- start(nuc_ranges)
  nuc_ends <- end(nuc_ranges)
  
  dist_to_left <- nuc_starts
  dist_to_right <- chr_lens[nuc_chrs] - nuc_ends
  dist_to_end <- pmin(dist_to_left, dist_to_right)
  
  nuclear <- nuc_chrs != "chrM"
  
  breaks <- seq(0, 200000, by=10000)
  bins <- cut(dist_to_end[nuclear], breaks=breaks, include.lowest=TRUE)
  
  K9me3_by_dist <- tapply(K9me3_status[nuclear], bins, mean, na.rm=TRUE)
  nuc_count_by_dist <- tapply(K9me3_status[nuclear], bins, length)
  
  cat("Distance from chr end | K9me3 fraction | Nucleosomes\n")
  cat("--------------------------------------------------\n")
  for(i in 1:length(K9me3_by_dist)) {
    if(!is.na(K9me3_by_dist[i])) {
      cat(sprintf("  %6d - %6d kb  |     %.3f      |   %d\n", 
                  breaks[i]/1000, breaks[i+1]/1000, 
                  K9me3_by_dist[i], nuc_count_by_dist[i]))
    }
  }
  
  mid_points <- (breaks[-length(breaks)] + breaks[-1]) / 2
  valid <- !is.na(K9me3_by_dist)
  if(sum(valid) > 3) {
    if(sd(K9me3_by_dist[valid]) > 0) {
      cor_test <- cor.test(mid_points[valid], K9me3_by_dist[valid], method="spearman")
      cat(sprintf("\nSpearman correlation (distance vs K9me3): rho = %.3f, p = %e\n", 
                  cor_test$estimate, cor_test$p.value))
      
      if(cor_test$estimate < -0.3 && cor_test$p.value < 0.05) {
        cat("✓ K9me3 shows significant DECLINING gradient from telomeres\n")
      } else if(cor_test$estimate < 0) {
        cat("~ K9me3 shows weak declining trend from telomeres\n")
      } else {
        cat("✗ K9me3 does not show telomeric gradient\n")
      }
    } else {
      cat("\nK9me3 fraction is uniform across all bins — no gradient detected\n")
    }
  }
  
  gradient_data <- data.frame(
    distance_kb = mid_points[valid] / 1000,
    K9me3_fraction = K9me3_by_dist[valid],
    nucleosome_count = nuc_count_by_dist[valid]
  )
  write.csv(gradient_data, file=paste0(outputDir, subSimName, ".K9me3_gradient.csv"), row.names=FALSE)
  
  png(filename=paste0(outputDir, subSimName, ".K9me3_gradient.png"), width=2400, height=1600, res=300)
  plot(gradient_data$distance_kb, gradient_data$K9me3_fraction, 
       type="b", pch=16, col="red",
       xlab="Distance from chromosome end (kb)", 
       ylab="K9me3 fraction",
       main="K9me3 enrichment vs distance from telomere")
  abline(h=mean(K9me3_status[nuclear]), lty=2, col="grey")
  legend("topright", legend=c("K9me3 fraction", "Genome average"), 
         col=c("red", "grey"), lty=c(1,2), pch=c(16, NA))
  dev.off()
  cat("\nGradient plot saved\n")
  
  # ========== K9 VALIDATION 3: K9me3 at known heterochromatic loci ==========
  cat("\n\n=== K9 TEST 3: K9me3 at known yeast heterochromatic loci ===\n\n")
  
  known_het <- GRanges(
    seqnames = c("chrIII", "chrIII"),
    ranges = IRanges(start = c(1, 290000), end = c(20000, 300000))
  )
  names(known_het) <- c("HML", "HMR")
  
  if(total_K9me3 > 0) {
    for(i in 1:length(known_het)) {
      locus <- known_het[i]
      locus_name <- names(known_het)[i]
      
      K9me3_at_locus <- length(unique(queryHits(findOverlaps(K9me3_regions, locus))))
      K4me3_at_locus <- length(unique(queryHits(findOverlaps(K4me3_regions, locus))))
      
      nucs_at_locus <- length(unique(queryHits(findOverlaps(nuc_ranges, locus))))
      K9me3_nucs_at_locus <- length(unique(queryHits(findOverlaps(nuc_ranges[K9me3_nucs], locus))))
      
      locus_K9_rate <- ifelse(nucs_at_locus > 0, K9me3_nucs_at_locus / nucs_at_locus * 100, 0)
      genome_K9_rate <- sum(K9me3_status[nuclear]) / sum(nuclear) * 100
      
      cat(sprintf("%s (%s:%d-%d):\n", locus_name, 
                  as.character(seqnames(locus)), start(locus), end(locus)))
      cat(sprintf("  Nucleosomes in region: %d\n", nucs_at_locus))
      cat(sprintf("  K9me3 nucleosomes: %d (%.1f%%)\n", K9me3_nucs_at_locus, locus_K9_rate))
      cat(sprintf("  Genome average K9me3: %.1f%%\n", genome_K9_rate))
      cat(sprintf("  K4me3 regions overlapping: %d\n", K4me3_at_locus))
      
      if(locus_K9_rate > genome_K9_rate * 1.2) {
        cat(sprintf("  ✓ K9me3 enriched at %s (%.1fx above average)\n\n", locus_name, locus_K9_rate / genome_K9_rate))
      } else {
        cat(sprintf("  ~ K9me3 at genome average at %s\n\n", locus_name))
      }
    }
  }
  
  # ========== SUMMARY TABLE ==========
  cat("\n============================================================\n")
  cat("SPATIAL VALIDATION SUMMARY\n")
  cat("============================================================\n\n")
  
  cat("| Test | Type | Metric | Result |\n")
  cat("|------|------|--------|--------|\n")
  
  if(exists("depletion_fold")) {
    cat(sprintf("| K4 at K9 overlap | Rule 2 emergent | K4me3 depletion | %.1fx |\n", depletion_fold))
  }
  
  if(exists("fold_K4_subtelo")) {
    cat(sprintf("| K4 at subtelomeres | Rule 2 emergent | Fold enrichment | %.2fx |\n", fold_K4_subtelo))
  }
  
  if(exists("internal_rate") && exists("subtelo_rate")) {
    cat(sprintf("| K4 expression proxy | Emergent | Internal vs subtelo | %.1f%% vs %.1f%% |\n", internal_rate, subtelo_rate))
  }
  
  if(exists("fold_K9_subtelo")) {
    cat(sprintf("| K9 at subtelomeres | Targeting | Fold enrichment | %.2fx |\n", fold_K9_subtelo))
  }
  
  if(exists("cor_test") && !is.na(cor_test$estimate)) {
    cat(sprintf("| K9 telomeric gradient | Targeting | Spearman rho | %.3f |\n", cor_test$estimate))
  }
  
  cat("\n=== Spatial validation complete ===\n")
  
}

unlink(paste0(outputDir,"scLayerSetBoth.object.Rdata") )  # large object that is probably no longer needed.



