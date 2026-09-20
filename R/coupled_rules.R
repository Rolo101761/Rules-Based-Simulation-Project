# ============================================================
# Coupled 4-mark rule set (K4/K9/K27 writers + erasers, Rules 1-3) — binding-factor definitions
# ============================================================
# VERBATIM copy of lines 114-416 of R/model5_tf_active_dna_seq.R (the sacCer3 model, taken 2026-09-20),
# so the dm6 model uses exactly the rule wiring validated on sacCer3 (Rule 1/2 antagonism at initial
# targeting AND at me2->me3 promotion, Rule 3A K27 self-recruitment).
#
# Why a copy rather than a shared source(): the Phase 1 skopt campaign re-reads the sacCer3 script from
# disk on every trial, so that file must not be edited until the campaign finishes. TODO (after Phase 1):
# make model5_tf_active_dna_seq.R source() this file instead and delete the duplicated block.
#
# Requires `nucleosomeWidth` (147) to be defined before sourcing. Defines all bf_* objects and the
# bfList_H3K4/H3K9/H3K27_methylate / _demethyl lists. Does NOT define abundances (those depend on the run).
#
# NOTE: the K4/K9/K27 tf_*_active layers are seeded by the caller (sequence motifs); nothing here is
# genome-specific.

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
                                            profile.layers=c("tf_K4_active", "sampled_meUp", "H3K9me3"),
                                            patternLength = 6,
                                            profile.marks = c(1, 1, 0),
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
                                                     profile.layers = c("nucleosome", "H3K4me_promotion", "H3K9me3", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3"),
                                                     profile.marks = c(1,1,0,0,0,1,0),
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
                                                profile.layers=c("tf_K27_active", "sampled_K27_meUp", "H3K4me3"),
                                                patternLength = 6,
                                                profile.marks = c(1, 1, 0),
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
                                                         profile.layers = c("nucleosome", "H3K27me_promotion", "H3K4me3", "H3K27me1", "H3K27me2", "H3K27me3"),
                                                         profile.marks = c(1,1,0,0,1,0),
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

