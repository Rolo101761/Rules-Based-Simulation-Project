# TATAWAWR expands to 8 concrete 8-mers. Which are actually TSS-specific in dm6 (fold of hits at TSS -30..-16 vs flanks)?
# Also reports motif availability on the chr4 development slice. Run from repo root.
suppressMessages({library(Biostrings); library(GenomicRanges); library(BSgenome.Dmelanogaster.UCSC.dm6); library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)})
genome <- BSgenome.Dmelanogaster.UCSC.dm6
chrs <- setdiff(seqnames(genome), "chrM")
tx <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene); tx <- tx[as.character(seqnames(tx)) %in% chrs]
tss <- unique(resize(tx, 1, fix = "start")); plus <- as.character(strand(tss)) == "+"; p <- start(tss)
win <- function(from, to) suppressWarnings(trim(GRanges(seqnames(tss), IRanges(ifelse(plus, p + from, p - to), ifelse(plus, p + to, p - from)),
                                                          strand = strand(tss), seqinfo = seqinfo(tss))))
core <- win(-30, -16)
flanks <- c(win(-300, -101), win(51, 250))
g <- expand.grid(w1 = c("A","T"), w2 = c("A","T"), r = c("A","G"), stringsAsFactors = FALSE)
variants <- paste0("TATA", g$w1, "A", g$w2, g$r)
scan <- function(pat) {
  fw <- PDict(DNAStringSet(pat)); rv <- PDict(reverseComplement(DNAStringSet(pat)))
  out <- list()
  for (chr in chrs) for (d in c("+", "-")) {
    m <- unlist(matchPDict(if (d == "+") fw else rv, genome[[chr]]))
    if (length(m)) out[[length(out) + 1]] <- GRanges(chr, IRanges(start(m), end(m)), strand = d)
  }
  x <- do.call(c, out); seqlevels(x) <- seqlevels(genome); seqinfo(x) <- seqinfo(genome); x
}
res <- do.call(rbind, lapply(variants, function(v) {
  h <- scan(v)
  c_core <- sum(countOverlaps(core, h)) / (length(core) * 15)
  c_fl <- sum(countOverlaps(flanks, h)) / sum(width(flanks))
  c4 <- sum(as.character(seqnames(h)) == "chr4")
  data.frame(variant = v, genome_hits = length(h), hits_chr4 = c4, tss_core_fold = round(c_core / c_fl, 2),
             pct_TSS_with_hit = round(100 * mean(countOverlaps(core, h) > 0), 2))
}))
res <- res[order(-res$tss_core_fold), ]
print(res, row.names = FALSE)
cat("\nExpected per-TSS chance rate scales with genome hit density; fold is vs flanking TSS-proximal windows.\n")
