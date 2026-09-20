# ============================================================
# dm6 sequence-anchored inputs: motif scan + empirical k-mer scan
# ============================================================
# Phase 2 step 1. Characterises the three sequence-anchored inputs of the coupled model on dm6:
#   K4  : TATA box                      (TATAWAWR, Drosophila core-promoter TATA)
#   K9  : satellite repeats             (AAGAG, AATAT, AACAC, dodeca; plus RepeatMasker satellites)
#   K27 : PRE motifs                    (GAF = GAGAG, Pho = GCCAT)
# and runs the same empirical 8-mer scan used for the yeast K9 motif (enrichment in
# heterochromatin vs euchromatin), then checks whether it independently rediscovers the
# annotated fly satellites (UCSC dm6 rmsk table = ground truth).
#
# Run from the repo root:  Rscript R/dm6_motif_scan.R
# Needs output/dm6_motifs/rmsk.txt.gz (UCSC goldenPath/dm6/database/rmsk.txt.gz).
# Outputs (all under output/dm6_motifs/): motif_summary.csv, kmer_top100.csv,
# kmer_rediscovery_families.csv, kmer_known_units.csv, TATA_TSS_profile.csv

suppressMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
})

outdir <- "output/dm6_motifs/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
genome <- BSgenome.Dmelanogaster.UCSC.dm6

chrs <- setdiff(seqnames(genome), "chrM")
main_chrs <- c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")
genome_bp <- sum(as.numeric(seqlengths(genome)[chrs]))
cat(sprintf("dm6 nuclear genome: %.1f Mb across %d sequences\n", genome_bp / 1e6, length(chrs)))

# ---------- helpers ----------

expand_iupac <- function(p) {
  opts <- lapply(strsplit(p, "")[[1]], function(ch) strsplit(IUPAC_CODE_MAP[[ch]], "")[[1]])
  unique(do.call(paste0, expand.grid(opts, stringsAsFactors = FALSE)))
}

rc_chr <- function(x) as.character(reverseComplement(DNAStringSet(x)))

# All exact occurrences of an IUPAC pattern on both strands -> strand-annotated GRanges
scan_motif <- function(pattern) {
  fwd <- expand_iupac(pattern)
  rev <- unique(rc_chr(fwd))
  pdF <- PDict(DNAStringSet(fwd))
  pdR <- PDict(DNAStringSet(rev))
  out <- list()
  for (chr in chrs) {
    s <- genome[[chr]]
    for (dir in c("+", "-")) {
      m <- unlist(matchPDict(if (dir == "+") pdF else pdR, s))
      if (length(m) > 0) out[[length(out) + 1]] <- GRanges(chr, IRanges(start(m), end(m)), strand = dir)
    }
  }
  with_genome_si(do.call(c, out))
}

# Tandem arrays of >= n copies of a repeat unit (both strands) -> merged array GRanges
scan_tandem <- function(unit, n = 3, mismatch = 0) {
  arr <- paste(rep(unit, n), collapse = "")
  out <- list()
  for (chr in chrs) {
    s <- genome[[chr]]
    for (pat in unique(c(arr, rc_chr(arr)))) {
      m <- matchPattern(DNAString(pat), s, max.mismatch = mismatch)
      if (length(m) > 0) out[[length(out) + 1]] <- GRanges(chr, IRanges(start(m), end(m)))
    }
  }
  if (length(out) == 0) return(GRanges())
  with_genome_si(reduce(do.call(c, out), min.gapwidth = 1))
}

pct <- function(bp) 100 * bp / genome_bp

# Give a GRanges the full dm6 seqinfo (adds unused seqlevels first, so the setter accepts it)
with_genome_si <- function(gr) {
  seqlevels(gr) <- seqlevels(genome)
  seqinfo(gr) <- seqinfo(genome)
  gr
}

# ---------- RepeatMasker ground truth ----------

rmsk <- read.delim(gzfile(paste0(outdir, "rmsk.txt.gz")), header = FALSE, stringsAsFactors = FALSE)
rmsk <- GRanges(rmsk$V6, IRanges(rmsk$V7 + 1, rmsk$V8), strand = "*",
                repName = rmsk$V11, repClass = rmsk$V12)
rmsk <- rmsk[as.character(seqnames(rmsk)) %in% chrs]
rmsk <- with_genome_si(rmsk)
sat <- rmsk[rmsk$repClass == "Satellite"]
simple <- rmsk[rmsk$repClass == "Simple_repeat"]
lowc <- rmsk[rmsk$repClass == "Low_complexity"]
cat(sprintf("rmsk: Satellite %d (%.2f Mb) | Simple_repeat %d (%.2f Mb) | Low_complexity %d (%.2f Mb)\n",
            length(sat), sum(width(reduce(sat))) / 1e6, length(simple), sum(width(reduce(simple))) / 1e6,
            length(lowc), sum(width(reduce(lowc))) / 1e6))

summary_rows <- list()
add_row <- function(input, motif, pattern, gr, note = "") {
  red <- reduce(gr, ignore.strand = TRUE)
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(
    input = input, motif = motif, pattern = pattern, n_hits = length(gr),
    hits_per_Mb = round(length(gr) / (genome_bp / 1e6), 1),
    hit_bp = sum(width(red)), pct_genome_by_hits = round(pct(sum(width(red))), 3),
    note = note, stringsAsFactors = FALSE)
}

# ============================================================
# 1. K4: TATA box, checked against annotated TSSs
# ============================================================
cat("\n=== K4: TATA box (TATAWAWR) ===\n")
tata <- scan_motif("TATAWAWR")
cat(sprintf("TATAWAWR hits: %d (%.0f per Mb)\n", length(tata), length(tata) / (genome_bp / 1e6)))

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
tx <- transcripts(txdb)
tx <- tx[as.character(seqnames(tx)) %in% chrs]
tss <- unique(resize(tx, 1, fix = "start"))
cat(sprintf("Unique TSSs: %d\n", length(tss)))

win <- function(from, to) {          # strand-aware window [TSS+from, TSS+to], positions relative to the TSS, 5'->3'
  plus <- as.character(strand(tss)) == "+"
  p <- start(tss)
  suppressWarnings(trim(GRanges(seqnames(tss), IRanges(ifelse(plus, p + from, p - to), ifelse(plus, p + to, p - from)),
                                strand = strand(tss), seqinfo = seqinfo(tss))))
}
rate_in <- function(w) {
  h <- countOverlaps(w, tata)                    # same-strand hits only (GRanges strand-aware); 8 bp hit vs 15 bp window
  sum(h) / (length(w) * median(width(w)))
}
profile <- data.frame(from = seq(-300, 100, by = 15))
profile$to <- profile$from + 14
profile$hits_per_tss_bp <- sapply(seq_len(nrow(profile)), function(i) rate_in(win(profile$from[i], profile$to[i])))
baseline <- mean(profile$hits_per_tss_bp[profile$to < -100 | profile$from > 50])
profile$fold_vs_baseline <- round(profile$hits_per_tss_bp / baseline, 2)
write.csv(profile, paste0(outdir, "TATA_TSS_profile.csv"), row.names = FALSE)
peak <- profile[which.max(profile$fold_vs_baseline), ]
cat(sprintf("TATA density peaks at TSS%+d..%+d with %.1fx the flanking baseline (expected ~ -35..-20)\n",
            peak$from, peak$to, peak$fold_vs_baseline))
core_w <- win(-35, -20)
has_tata <- countOverlaps(core_w, tata) > 0
cat(sprintf("TSSs with a TATA hit at -35..-20: %.1f%% (chance level ~ %.1f%%)\n",
            100 * mean(has_tata), 100 * (1 - exp(-baseline * 16))))
add_row("K4", "TATA box", "TATAWAWR", tata,
        sprintf("%.1fx enriched at TSS%+d..%+d; %.1f%% of TSSs carry one at -35..-20",
                peak$fold_vs_baseline, peak$from, peak$to, 100 * mean(has_tata)))

# ============================================================
# 2. K27: PRE motifs (GAF, Pho) — sites and clustering
# ============================================================
cat("\n=== K27: PRE motifs ===\n")
gaf <- scan_motif("GAGAG")
gaf9 <- scan_motif("GAGAGAGAG")
pho <- scan_motif("GCCAT")
pho_long <- scan_motif("GCCATHWT")
cat(sprintf("GAF GAGAG: %d | GAGAGAGAG: %d | Pho GCCAT: %d | Pho GCCATHWT: %d\n",
            length(gaf), length(gaf9), length(pho), length(pho_long)))

# PRE-like windows: 400 bp bins with >= 3 GAF sites AND >= 1 Pho site
bins <- unlist(tile(GRanges(chrs, IRanges(1, seqlengths(genome)[chrs]), seqinfo = seqinfo(genome)), width = 400))
n_gaf <- countOverlaps(bins, gaf); n_pho <- countOverlaps(bins, pho)
pre_like <- bins[n_gaf >= 3 & n_pho >= 1]
exp_bins <- mean(n_gaf >= 3) * mean(n_pho >= 1) * length(bins)
cat(sprintf("PRE-like 400bp bins (>=3 GAGAG and >=1 GCCAT): %d = %.2f%% of genome (independence expectation ~%.0f bins)\n",
            length(pre_like), pct(sum(width(pre_like))), exp_bins))
cat(sprintf("  ... of which on main arms: %s\n",
            paste(names(table(seqnames(pre_like)))[table(seqnames(pre_like)) > 0],
                  table(seqnames(pre_like))[table(seqnames(pre_like)) > 0], sep = "=", collapse = " ")))
add_row("K27", "GAF core", "GAGAG", gaf, "raw site: very frequent, needs clustering")
add_row("K27", "GAF 9-mer", "GAGAGAGAG", gaf9)
add_row("K27", "Pho core", "GCCAT", pho)
add_row("K27", "Pho ext", "GCCATHWT", pho_long)
add_row("K27", "PRE-like bin", ">=3 GAGAG & >=1 GCCAT / 400bp", pre_like,
        sprintf("%.1fx the independence expectation (%.0f bins)", length(pre_like) / exp_bins, exp_bins))

# ============================================================
# 3. K9: satellite repeat units
# ============================================================
cat("\n=== K9: satellite repeat units (tandem arrays, >=3 copies) ===\n")
sat_units <- list(AAGAG = "AAGAG", AATAT = "AATAT", AACAC = "AACAC")
sat_arrays <- list()
for (nm in names(sat_units)) {
  a <- scan_tandem(sat_units[[nm]], n = 3)
  sat_arrays[[nm]] <- a
  in_rmsk <- mean(countOverlaps(a, c(simple, sat, lowc)) > 0)
  cat(sprintf("(%s)3+ arrays: %d, %.2f Mb (%.2f%% of genome); %.0f%% overlap an rmsk repeat\n",
              nm, length(a), sum(width(a)) / 1e6, pct(sum(width(a))), 100 * in_rmsk))
  add_row("K9", paste0("(", nm, ")n"), paste0(sat_units[[nm]], " x>=3"), a,
          sprintf("%.0f%% of arrays overlap an rmsk repeat", 100 * in_rmsk))
}
# dodeca satellite (12-mer unit CCCGTACTCGGT), allow 1 mismatch per 24-mer
dodeca <- scan_tandem("CCCGTACTCGGT", n = 2, mismatch = 1)
cat(sprintf("dodeca (CCCGTACTCGGT x2, <=1 mm): %d arrays, %.3f Mb\n", length(dodeca), sum(width(dodeca)) / 1e6))
if (length(dodeca) > 0) cat(sprintf("  on: %s\n", paste(names(table(seqnames(dodeca)))[table(seqnames(dodeca)) > 0],
                                                        table(seqnames(dodeca))[table(seqnames(dodeca)) > 0], sep = "=", collapse = " ")))
add_row("K9", "dodeca", "CCCGTACTCGGT x2 (<=1 mm)", dodeca, "expected to be chr3 pericentromeric")

cat("\nRepeatMasker Satellite families (sequence not asserted; identity taken from rmsk names):\n")
fam <- do.call(rbind, lapply(split(sat, sat$repName), function(g) data.frame(
  family = g$repName[1], n = length(g), total_kb = round(sum(width(g)) / 1e3, 1),
  median_len = median(width(g)),
  top_chrs = paste(head(names(sort(table(as.character(seqnames(g))), decreasing = TRUE)), 3), collapse = ","))))
print(fam[order(-fam$total_kb), ], row.names = FALSE)
for (f in fam$family) add_row("K9", paste0("rmsk:", f), "RepeatMasker", sat[sat$repName == f],
                              sprintf("median unit %.0f bp", fam$median_len[fam$family == f]))
add_row("K9", "rmsk Satellite (all)", "RepeatMasker", sat)

# ============================================================
# 4. Empirical 8-mer scan (heterochromatin vs euchromatin) — as in the yeast K9 work
# ============================================================
cat("\n=== Empirical 8-mer scan: pericentromeric/heterochromatic vs euchromatic ===\n")
# Approximate eu/heterochromatin boundaries (from memory of published dm6 annotations; +-1 Mb) —
# only used to define foreground vs background for the enrichment ranking; real satellites are
# enriched >>10x so results are insensitive to a 1 Mb error.
het <- GRanges(c("chr2L", "chr2R", "chr3L", "chr3R", "chrX", "chr4", "chrY"),
               IRanges(c(22000001, 1, 23000001, 1, 22000001, 1, 1),
                       c(23513712, 6000000, 28110227, 4600000, 23542271, 1348131, 3667352)))
het <- c(het, GRanges(setdiff(chrs, main_chrs), IRanges(1, seqlengths(genome)[setdiff(chrs, main_chrs)])))
het <- with_genome_si(het)
het <- trim(het)
all_chr <- GRanges(chrs, IRanges(1, seqlengths(genome)[chrs]), seqinfo = seqinfo(genome))
eu <- setdiff(all_chr, het)
cat(sprintf("Foreground (het proxy): %.1f Mb | background (euchromatin proxy): %.1f Mb\n",
            sum(width(het)) / 1e6, sum(width(eu)) / 1e6))

K <- 8
all_k <- mkAllStrings(DNA_BASES, K)
rc_idx <- match(rc_chr(all_k), all_k)
count_kmers <- function(regions) {
  tot <- numeric(length(all_k))
  for (chr in unique(as.character(seqnames(regions)))) {
    r <- regions[seqnames(regions) == chr]
    tot <- tot + oligonucleotideFrequency(getSeq(genome, r), K, simplify.as = "collapsed")
  }
  tot
}
fg <- count_kmers(het); bg <- count_kmers(eu)
fgN <- sum(fg); bgN <- sum(bg)
keep <- seq_along(all_k) <= rc_idx                 # one representative per reverse-complement pair
fg_c <- ifelse(seq_along(all_k) == rc_idx, fg, fg + fg[rc_idx])
bg_c <- ifelse(seq_along(all_k) == rc_idx, bg, bg + bg[rc_idx])
kt <- data.frame(kmer = all_k, rc = all_k[rc_idx], fg = fg_c, bg = bg_c, stringsAsFactors = FALSE)[keep, ]
kt$log2_enrich <- log2(((kt$fg + 0.5) / fgN) / ((kt$bg + 0.5) / bgN))
kt <- kt[kt$fg >= 100, ]
kt <- kt[order(-kt$log2_enrich), ]
kt$rank <- seq_len(nrow(kt))
cat(sprintf("%d canonical 8-mers (of %d) pass fg count >= 100\n", nrow(kt), sum(keep)))

# Genomic positions of the top-N k-mers, then how much of each rmsk satellite/simple-repeat family they cover
TOPN <- 100
top <- head(kt, TOPN)
pdT <- PDict(DNAStringSet(c(top$kmer, top$rc[top$kmer != top$rc])))
hits <- list()
for (chr in chrs) {
  m <- unlist(matchPDict(pdT, genome[[chr]]))
  if (length(m) > 0) hits[[length(hits) + 1]] <- GRanges(chr, IRanges(start(m), end(m)))
}
hits <- with_genome_si(reduce(do.call(c, hits)))
cat(sprintf("Top-%d k-mers: %d loci (merged), %.2f Mb\n", TOPN, length(hits), sum(width(hits)) / 1e6))

ann <- c(sat, simple, lowc)
frac_in_ann <- sum(width(intersect(hits, reduce(ann)))) / sum(width(hits))
cat(sprintf("%.1f%% of top-%d k-mer bp fall inside rmsk Satellite/Simple_repeat/Low_complexity (genome-wide base rate: %.2f%%)\n",
            100 * frac_in_ann, TOPN, pct(sum(width(reduce(ann))))))

recall <- function(fam_gr) {
  g <- reduce(fam_gr)
  ov <- findOverlaps(g, hits)
  sum(width(pintersect(g[queryHits(ov)], hits[subjectHits(ov)]))) / sum(width(g))
}
fam_row <- function(cls, g) data.frame(class = cls, family = g$repName[1], total_kb = round(sum(width(reduce(g))) / 1e3, 1),
                                       pct_covered_by_top_kmers = round(100 * recall(g), 1))
simple_big <- split(simple, simple$repName)
simple_big <- simple_big[sapply(simple_big, function(g) sum(width(g)) >= 3000)]   # >= 3 kb total: skip thousands of tiny families
fam_tab <- rbind(do.call(rbind, lapply(split(sat, sat$repName), function(g) fam_row("Satellite", g))),
                 do.call(rbind, lapply(simple_big, function(g) fam_row("Simple_repeat", g))))
fam_tab <- fam_tab[order(-fam_tab$total_kb), ]
write.csv(fam_tab, paste0(outdir, "kmer_rediscovery_families.csv"), row.names = FALSE)
cat("\nRediscovery: % of each annotated repeat family covered by top-100 k-mer hits\n")
cat("(Satellite families, then the 12 biggest Simple_repeat families)\n")
print(rbind(fam_tab[fam_tab$class == "Satellite", ], head(fam_tab[fam_tab$class == "Simple_repeat", ], 12)), row.names = FALSE)

# Known units: is the best-ranked k-mer that is a substring of the (rotated) unit repeated?
known <- c(AAGAG = "AAGAG", AATAT = "AATAT", AACAC = "AACAC", AAGAC = "AAGAC", dodeca = "CCCGTACTCGGT")
kn <- do.call(rbind, lapply(names(known), function(nm) {
  u <- paste(rep(known[[nm]], ceiling(20 / nchar(known[[nm]]))), collapse = "")
  in_u <- function(km) grepl(km, u, fixed = TRUE) | grepl(rc_chr(km), u, fixed = TRUE) |
    grepl(km, rc_chr(u), fixed = TRUE)
  m <- kt[sapply(kt$kmer, in_u), ]
  data.frame(unit = nm, n_kmers_of_unit_in_top100 = sum(m$rank <= 100),
             best_rank = if (nrow(m)) min(m$rank) else NA_integer_,
             best_kmer = if (nrow(m)) m$kmer[which.min(m$rank)] else NA_character_,
             best_log2_enrich = if (nrow(m)) round(max(m$log2_enrich), 2) else NA_real_,
             n_kmers_of_unit_total = nrow(m))
}))
write.csv(kn, paste0(outdir, "kmer_known_units.csv"), row.names = FALSE)
cat("\nKnown fly satellite units vs the empirical scan:\n"); print(kn, row.names = FALSE)

write.csv(head(kt, 100), paste0(outdir, "kmer_top100.csv"), row.names = FALSE)
cat("\nTop 30 canonical 8-mers by log2 enrichment (het proxy vs euchromatin proxy):\n")
print(head(kt, 30), row.names = FALSE)

# ============================================================
# 5. Summary
# ============================================================
summ <- do.call(rbind, summary_rows)
write.csv(summ, paste0(outdir, "motif_summary.csv"), row.names = FALSE)
cat("\n=== MOTIF SUMMARY ===\n"); print(summ[, c("input", "motif", "n_hits", "pct_genome_by_hits", "note")], row.names = FALSE)
cat("\nDONE\n")
