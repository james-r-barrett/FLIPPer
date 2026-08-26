#!/usr/bin/env Rscript
## Thin wrapper around DECIPHER::DetectRepeats() so FLIPPer.py can shell out to it the
## same way it used to shell out to xstream.jar. Reads a protein FASTA, runs DetectRepeats,
## and writes one row per detected tandem-repeat region to a CSV - Period/Copies are derived
## here (from the Left/Right per-copy position lists DetectRepeats returns) so the Python side
## never has to parse R list-columns, only a plain CSV.
##
## Optionally (when alignments_dir is non-empty) also extracts each hit's individual repeat
## units and aligns them with AlignSeqs(), writing one aligned FASTA per hit - this is what
## lets FLIPPer show repeat copies lined up against each other, the way XSTREAM's report did.
suppressMessages(library(DECIPHER))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: detect_repeats.R <fasta_path> <min_score> <max_period> <max_copies> <out_csv> <alignments_dir|\"\">")
}
fasta_path <- args[1]
min_score  <- as.numeric(args[2])
max_period <- as.numeric(args[3])
max_copies <- as.numeric(args[4])
out_csv    <- args[5]
alignments_dir <- args[6]

seqs <- readAAStringSet(fasta_path)

empty_result <- data.frame(ID = character(), Begin = integer(), End = integer(),
                            Period = numeric(), Copies = integer(), Score = numeric(),
                            RepeatIndex = integer(), UnitLefts = character(), UnitRights = character())

if (length(seqs) == 0) {
  write.csv(empty_result, out_csv, row.names = FALSE)
  quit(status = 0)
}

## DetectRepeats' processors argument does NOT parallelize the type="tandem" search (confirmed
## empirically: processors=NULL showed no speedup over processors=1 - it appears to only apply
## to the type="interspersed"/FindSynteny path). Real parallelism for large candidate sets is
## instead done at the process level, by run_detect_repeats() in FLIPPer_lib.py splitting the
## input across multiple concurrent invocations of this script - so this stays single-process.
hits <- DetectRepeats(seqs,
                       type = "tandem",
                       minScore = min_score,
                       maxPeriod = max_period,
                       maxCopies = max_copies,
                       processors = 1,
                       verbose = FALSE)

if (is.null(hits) || nrow(hits) == 0) {
  write.csv(empty_result, out_csv, row.names = FALSE)
  quit(status = 0)
}

## FASTA IDs here are always single tokens (no description) - FLIPPer.py only ever passes
## detect_repeats.R files it wrote itself using record.id, so this is a plain lookup, not a
## header-parsing step. sub("^>+", ...) strips any stray leading '>' characters so downstream
## Python-side ID matching against cleanly-read record.id values always lines up, regardless of
## how the upstream FASTA writer formatted its headers.
ids <- sub("^>+", "", names(seqs)[hits$Index])
periods <- mapply(function(l, r) mean(r - l + 1), hits$Left, hits$Right)
copies <- lengths(hits$Left)
## per-hit-row index within each ID, so a sequence with more than one detected repeat region
## gets a stable identifier tying its CSV row to its alignment file below
row_idx <- ave(seq_along(ids), ids, FUN = seq_along)
## semicolon-joined start/end of each individual repeat copy (1-based, inclusive) - kept
## alongside the summary Period/Copies columns so the Python side can draw a block per copy
## against the original sequence coordinates (not just the alignment-local ones)
unit_lefts <- sapply(hits$Left, paste, collapse = ";")
unit_rights <- sapply(hits$Right, paste, collapse = ";")

out <- data.frame(ID = ids,
                   Begin = hits$Begin,
                   End = hits$End,
                   Period = periods,
                   Copies = copies,
                   Score = hits$Score,
                   RepeatIndex = row_idx,
                   UnitLefts = unit_lefts,
                   UnitRights = unit_rights)

write.csv(out, out_csv, row.names = FALSE)

if (nzchar(alignments_dir)) {
  if (!dir.exists(alignments_dir)) dir.create(alignments_dir, recursive = TRUE)
  safe_ids <- gsub("[^A-Za-z0-9_.-]", "_", ids)
  for (i in seq_len(nrow(out))) {
    units <- extractAt(seqs[[hits$Index[i]]], IRanges(hits$Left[[i]], hits$Right[[i]]))
    if (length(units) < 2) next
    aligned <- tryCatch(AlignSeqs(units, verbose = FALSE), error = function(e) NULL)
    if (is.null(aligned)) next
    names(aligned) <- paste0("copy", seq_along(aligned))
    ## ConsensusSequence() is computed from this exact alignment, so what gets shown as the
    ## consensus is guaranteed to actually match the repeat structure in the alignment above it,
    ## rather than a separately-reimplemented (and potentially inconsistent) majority vote
    consensus <- tryCatch(ConsensusSequence(aligned, threshold = 0.3), error = function(e) NULL)
    if (!is.null(consensus)) {
      names(consensus) <- "consensus"
      aligned <- c(aligned, consensus)
    }
    fname <- file.path(alignments_dir, paste0(safe_ids[i], "__", row_idx[i], ".fasta"))
    writeXStringSet(aligned, fname)
  }
}
