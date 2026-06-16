# reproduce_bug.R --------------------------------------------------------------
#
# Minimal, self-contained demo of the findLQR() over-culling bug on fixture 01.
# Loads the in-development amplican source (so edits to findLQR are exercised),
# runs findLQR on the high-editing fixture, and shows that the CURRENT code culls
# ~62% of reads -- almost all of them edited -- dropping an apparent 90% editing
# rate to ~76%. Compare against the desired behaviour stored in expected.json.
#
# Usage:  Rscript reproduce_bug.R
# Expect: CURRENT prints a large cull (1866) and a FAIL vs desired; once findLQR
#         is fixed this should print a small/zero cull and PASS.

suppressMessages({library(data.table); library(jsonlite); library(amplican)})

AMP_PATH <- (function() {
  cand <- tryCatch({ a <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", a[grep("--file=", a)][1]); if (is.na(f)) stop(); dirname(normalizePath(f))
  }, error = function(e) getwd())
  for (i in 1:6) if (file.exists(file.path(cand, "DESCRIPTION"))) return(normalizePath(cand)) else cand <- dirname(cand)
  normalizePath(".", mustWork = TRUE)
})()
pkgload::load_all(AMP_PATH)

case_dir <- file.path(AMP_PATH, "findLQR_tuning", "fixtures", "01_high_edit_falsepos")
aln  <- as.data.frame(fread(file.path(case_dir, "events.csv")))
expd <- fromJSON(file.path(case_dir, "expected.json"))

# what does the in-development findLQR do?
flagged <- amplican::findLQR(aln)
onlyBR  <- unique(as.data.table(aln)[flagged, ], by = "read_id")
culled  <- if (nrow(onlyBR)) as.integer(sum(onlyBR$counts)) else 0L

# how many of the culled reads actually carried an indel (i.e. were edited)?
culled_ids <- onlyBR$read_id
all_ids    <- unique(as.data.table(aln)$read_id)
edited_ids <- unique(as.data.table(aln)[type %in% c("deletion", "insertion"), read_id])
nonwt_cull <- length(intersect(culled_ids, edited_ids))

total_reads <- length(all_ids)
cat(sprintf("Fixture:        01_high_edit_falsepos (nlgn4a @ 90%% editing)\n"))
cat(sprintf("findLQR culled: %d / %d unique reads  (%.1f%%)\n",
            length(culled_ids), total_reads, length(culled_ids) / total_reads * 100))
cat(sprintf("  ...of which %d carried an indel (edited) => the cull removes EDITED reads\n",
            nonwt_cull))
cat(sprintf("Low_Score equivalent (sum counts): %d  (pipeline: %d)\n",
            culled, expd$pipeline_low_score))

# apparent editing rate before vs after the filter
n_edited_all <- length(intersect(all_ids, edited_ids))
kept_ids <- setdiff(all_ids, culled_ids)
n_edited_kept <- length(intersect(kept_ids, edited_ids))
cat(sprintf("Apparent editing: %.1f%% (all reads) -> %.1f%% (after findLQR)\n",
            n_edited_all / total_reads * 100,
            if (length(kept_ids)) n_edited_kept / length(kept_ids) * 100 else NA))

desired_filter <- expd$desired$should_filter
cat(sprintf("\nDesired should_filter=%s  | current culled=%d  -> %s\n",
            desired_filter, culled,
            if (is.logical(desired_filter) && !desired_filter && culled == 0) "PASS" else "FAIL (bug present)"))
