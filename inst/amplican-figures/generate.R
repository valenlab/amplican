#!/usr/bin/env Rscript
# Regenerates the BLAST-style pairwise alignment SVGs embedded in the
# amplican walkthrough animation (inst/amplican_animation.html).
#
# Uses pwalign::pairwiseAlignment (the same Needleman-Wunsch aligner amplican
# uses internally) on the ID_1 amplicon. The "match" read is the perfect
# amplicon (no perfectly-clean read exists in the example data); the "deletion"
# read is the amplicon with positions 34..117 removed, mirroring read 4.
#
# Output: inst/amplican-figures/alignment_match.svg
#         inst/amplican-figures/alignment_deletion.svg
#
# Run from the package root:  Rscript inst/amplican-figures/generate.R

suppressMessages({
  library(pwalign)
  library(Biostrings)
})

out_dir <- dirname(normalizePath(".", mustWork = TRUE))
# locate the inst/amplican-figures directory relative to this script
this <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", this, value = TRUE)
if (length(script_arg)) {
  out_dir <- dirname(sub("--file=", "", script_arg))
} else {
  out_dir <- "inst/amplican-figures"
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cfg <- read.csv(system.file("extdata", "config.csv", package = "amplican"))
id1 <- cfg[cfg$ID == "ID_1", ]
amp <- toupper(id1$Amplicon)                 # 199 bp
mat <- nucleotideSubstitutionMatrix(match = 5, mismatch = -4, baseOnly = TRUE)

align <- function(read) {
  aln <- pairwiseAlignment(read, amp, type = "global",
                           substitutionMatrix = mat,
                           gapOpening = 25, gapExtension = 0)
  list(read = as.character(pattern(aln)), amp = as.character(subject(aln)),
       score = score(aln))
}

match_read   <- amp                                            # perfect match
del_read     <- paste0(substr(amp, 1, 33), substr(amp, 118, nchar(amp)))  # del 34..117

render <- function(al, title, file) {
  r <- strsplit(al$read, "")[[1]]
  a <- strsplit(al$amp,  "")[[1]]
  n <- length(r)
  cw <- 7        # char width
  fs <- 11       # font size
  x0 <- 56       # left padding for the "read/amplicon" labels
  W <- x0 + n * cw + 16
  H <- 104
  yR <- 44; yC <- 60; yA <- 80                  # read, consensus, amplicon baselines

  esc <- function(s) gsub("&", "&amp;", gsub("<", "&lt;", gsub(">", "&gt;", s)))

  lines <- c(
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" width="100%%" height="auto" font-family="ui-monospace,SFMono-Regular,Menlo,Consolas,monospace">', W, H),
    sprintf('<rect x="0" y="0" width="%.0f" height="%.0f" fill="#ffffff"/>', W, H),
    sprintf('<text x="10" y="22" font-size="13" font-weight="700" fill="#0f172a">%s</text>', esc(title)),
    sprintf('<text x="10" y="%g" font-size="10" fill="#64748b">read</text>', yR),
    sprintf('<text x="10" y="%g" font-size="10" fill="#64748b">amplicon</text>', yA),
    sprintf('<text x="%.0f" y="22" font-size="10" fill="#64748b" text-anchor="end">score %g</text>', W - 10, al$score)
  )

  # red highlight behind columns where the READ has a gap (a deletion in the read)
  i <- 1
  while (i <= n) {
    if (r[i] == "-") {
      j <- i
      while (j <= n && r[j] == "-") j <- j + 1
      x1 <- x0 + (i - 1) * cw - 1
      x2 <- x0 + (j - 1) * cw + 1
      lines <- c(lines, sprintf(
        '<rect x="%g" y="%g" width="%g" height="14" rx="2" fill="#fee2e2"/><text x="%g" y="%g" font-size="9" font-weight="700" fill="#ef4444" text-anchor="middle">deletion (%d bp)</text>',
        x1, yR - 10, x2 - x1, (x1 + x2) / 2, yR - 14, j - i))
      i <- j
    } else i <- i + 1
  }

  for (k in seq_len(n)) {
    x <- x0 + (k - 1) * cw + cw / 2
    rc <- r[k]; ac <- a[k]
    # consensus
    if (rc == "-" || ac == "-") {
      cc <- " "; ccol <- "#ffffff"
    } else if (rc == ac) {
      cc <- "|"; ccol <- "#cbd5e1"
    } else {
      cc <- "."; ccol <- "#ef4444"
    }
    rcol <- if (rc == "-") "#ef4444" else "#334155"
    acol <- if (ac == "-") "#ef4444" else "#334155"
    rch <- if (rc == "-") "-" else rc
    ach <- if (ac == "-") "-" else ac
    lines <- c(lines,
      sprintf('<text x="%g" y="%g" font-size="%g" fill="%s" text-anchor="middle">%s</text>', x, yR, fs, rcol, rch),
      sprintf('<text x="%g" y="%g" font-size="%g" fill="%s" text-anchor="middle">%s</text>', x, yC, fs - 2, ccol, cc),
      sprintf('<text x="%g" y="%g" font-size="%g" fill="%s" text-anchor="middle">%s</text>', x, yA, fs, acol, ach)
    )
  }
  lines <- c(lines, '</svg>')
  writeLines(lines, file)
  message("wrote ", file, " (", n, " aligned columns)")
}

render(align(match_read),   "Matching alignment (no edits)",
       file.path(out_dir, "alignment_match.svg"))
render(align(del_read),     "Deletion alignment (read vs amplicon)",
       file.path(out_dir, "alignment_deletion.svg"))
