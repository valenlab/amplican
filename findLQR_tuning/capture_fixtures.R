# capture_fixtures.R -----------------------------------------------------------
#
# Reconstructs the EXACT data.frame that amplican's findLQR() receives inside
# amplicanPipeline(), for a diverse set of cases, and writes each as a hermetic
# fixture under fixtures/<case>/events.csv plus an expected.json describing the
# intended behaviour.
#
# Faithfulness: amplicanPipeline() applies findEOP() then findPD() to the raw
# events BEFORE calling findLQR() (R/amplican.R lines ~289-319). raw_events.csv
# is written before those filters, so we re-apply findEOP+findPD here. The
# reconstructed Low_Score is asserted to equal the pipeline's reported Low_Score
# (from config_summary.csv) for every real case -- that equality is the
# acceptance criterion proving the fixture matches findLQR's true input.
#
# Usage:  Rscript capture_fixtures.R
# Edits to R/helpers_filters.R::findLQR are picked up via pkgload::load_all().

suppressMessages({library(data.table); library(jsonlite)})

AMP_PATH  <- (function() {
  cand <- tryCatch({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", a[grep("--file=", a)][1])
    if (is.na(f)) stop(); dirname(normalizePath(f))
  }, error = function(e) getwd())
  for (i in 1:6) if (file.exists(file.path(cand, "DESCRIPTION"))) return(normalizePath(cand)) else cand <- dirname(cand)
  normalizePath(".", mustWork = TRUE)
})()                          # this package's root (walks up to DESCRIPTION)
BENCH_ROOT <- "/media/ai/valenbackup/Projects/uib/CRISPR/efficiency/amplican_manuscript/analysis/indel_size/simulation/amplican_results"
FIX_DIR   <- file.path(AMP_PATH, "findLQR_tuning", "fixtures")

# load the in-development source so an edited findLQR is what we exercise
pkgload::load_all(AMP_PATH)

# -- reconstruction ------------------------------------------------------------
# Given a condition directory + ID, return the post-EOP/PD events for that ID,
# i.e. exactly the `aln_id` passed to findLQR() in amplicanPipeline().
# The full post-EOP/PD aln is memoized per condition (several cases share one).
.proc_cache <- new.env(parent = emptyenv())
processed_aln <- function(condition_dir) {
  if (!is.null(.proc_cache[[condition_dir]])) return(.proc_cache[[condition_dir]])
  aln  <- fread(file.path(condition_dir, "alignments", "raw_events.csv"))
  cfgT <- copy(fread(file.path(condition_dir, "config_summary.csv")))
  eop <- amplican::findEOP(aln, cfgT)            # mutates cfgT NA -> 1 / nchar
  aln <- aln[!eop, ]
  pd  <- amplican::findPD(aln, cfgT)
  onlyPD <- unique(aln[pd, ], by = c("seqnames", "read_id"))
  aln <- aln[!onlyPD, on = c("seqnames", "read_id")]
  .proc_cache[[condition_dir]] <- aln
  aln
}
reconstruct <- function(condition_dir, target_id) {
  as.data.frame(processed_aln(condition_dir)[seqnames == target_id])
}

# run current findLQR, report filtered read-counts (the Low_Score metric)
current_filtered <- function(aln_df) {
  l <- amplican::findLQR(aln_df)
  if (!any(l)) return(0L)
  onlyBR <- unique(as.data.table(aln_df)[l, ], by = "read_id")
  as.integer(sum(onlyBR$counts))
}

pipeline_low_score <- function(condition_dir, target_id) {
  cfgT <- fread(file.path(condition_dir, "config_summary.csv"))
  as.integer(cfgT[ID == target_id, Low_Score])
}

# `filtered` is precomputed by the caller so findLQR runs exactly once per case
write_fixture <- function(case_name, aln_df, expected, filtered) {
  odir <- file.path(FIX_DIR, case_name)
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)
  fwrite(as.data.table(aln_df), file.path(odir, "events.csv"))
  expected$current_filtered_counts <- filtered
  writeLines(toJSON(expected, auto_unbox = TRUE, pretty = TRUE),
             file.path(odir, "expected.json"))
  cat(sprintf("  %-26s rows=%6d  current_filtered=%5d  -> %s/%s\n",
              case_name, nrow(aln_df), filtered, "fixtures", case_name))
}

# -- real cases (from the synthetic benchmark) ---------------------------------
# 1freq = "No indels > 10bp"; 4freq = "Deletions > 10bp"
cond <- function(x) file.path(BENCH_ROOT, x)
real_cases <- list(
  list(case = "01_high_edit_falsepos",  cond = cond("270mut_30wt_1freq_150readlen"), id = "nlgn4a",
       truth = 90, note = "THE BUG: 90% editing, edited reads culled as 'junk'",
       desired_filter = FALSE, desired_reason = "edited reads are the majority; culling them is the failure"),
  list(case = "02_high_edit_escapee",   cond = cond("270mut_30wt_1freq_150readlen"), id = "cnsta",
       truth = 90, note = "borderline: k2 happens to beat k3 so nothing filtered",
       desired_filter = FALSE, desired_reason = "same population as case 01; fix must make ALL 90% samples behave like this"),
  list(case = "03_mid_edit_negative",   cond = cond("200mut_100wt_1freq_150readlen"), id = "nlgn4a",
       truth = 66.7, note = "66.7% editing: correctly not filtered",
       desired_filter = FALSE, desired_reason = "regression guard"),
  list(case = "04_low_edit_negative",   cond = cond("100mut_200wt_1freq_150readlen"), id = "nlgn4a",
       truth = 33.3, note = "33.3% editing: correctly not filtered",
       desired_filter = FALSE, desired_reason = "regression guard"),
  list(case = "05_pure_wt",             cond = cond("0mut_300wt_1freq_150readlen"),   id = "nlgn4a",
       truth = 0,   note = "WT only: nothing to filter",
       desired_filter = FALSE, desired_reason = "regression guard"),
  list(case = "07_large_del_high_edit", cond = cond("270mut_30wt_4freq_150readlen"), id = "nlgn4a",
       truth = 90, note = "bug also fires for large deletions -> editing-fraction driven",
       desired_filter = FALSE, desired_reason = "same root cause as case 01")
)

cat("Capturing real cases (reconstruct findEOP+findPD, verify vs pipeline):\n")
for (rc in real_cases) {
  aln_df <- reconstruct(rc$cond, rc$id)
  recon  <- current_filtered(aln_df)            # findLQR runs exactly once here
  pipe   <- pipeline_low_score(rc$cond, rc$id)
  status <- if (recon == pipe) "MATCH" else "*** MISMATCH ***"
  cat(sprintf("  %-26s reconstructed=%5d pipeline=%5d  %s\n", rc$case, recon, pipe, status))
  stopifnot(recon == pipe)   # refuse to ship a fixture that doesn't reproduce
  write_fixture(rc$case, aln_df, list(
    source = paste0("benchmark: ", basename(rc$cond), " / ID=", rc$id),
    editing_rate_truth = rc$truth,
    pipeline_low_score = pipe,
    description = rc$note,
    desired = list(should_filter = rc$desired_filter, reason = rc$desired_reason)
  ), filtered = recon)
}

# -- cases 06 & 08: synthetic, fully controlled (genuine junk must KEEP firing) -
# findLQR only fires when clara's k=3 silhouette beats k=2, i.e. the data needs
# THREE separable clusters. Real WT+edited reads are NOT well separated in the
# (score, events) space, so a hand-built 3-cluster design is the faithful way to
# probe the *legitimate* off-target/junk use case. We share mk_cluster() across
# the two cases and vary the editing majority:
#   case 06: WT majority  -> findLQR correctly culls the junk minority today
#   case 08: edited majority -> findLQR culls NOTHING today (false-negative on
#            junk under high editing) -- the complement of the over-culling bug
mk_cluster <- function(n, score_min, score_max, events_per_read, label) {
  epr <- rep_len(events_per_read, n)              # length n (scalar or vector ok)
  read_ids <- rep(paste0(label, "_", seq_len(n)), times = epr)
  nr <- length(read_ids)
  starts <- sample(40:200, nr, replace = TRUE)
  data.table(
    seqnames = "synthetic",
    start = starts,
    end   = starts + sample(0:5, nr, replace = TRUE),
    width = 1L,
    strand = "+",
    originally = "A",
    replacement = "C",
    type = "mismatch",
    read_id = read_ids,
    score = as.integer(round(runif(nr, score_min, score_max))),
    counts = 1L,
    readType = if (label == "junk") "offtarget" else "target"
  )
}

cat("\nCapturing synthetic case 06 (WT majority + genuine junk, must keep firing):\n")
set.seed(7)
case06 <- rbind(
  mk_cluster(2000, 690, 710, 1, "wt"),
  mk_cluster( 400, 620, 650, 2, "edited"),
  mk_cluster( 120, 280, 340, sample(3:4, 120, replace = TRUE), "junk")
)
f06 <- current_filtered(as.data.frame(case06))
write_fixture("06_genuine_junk", as.data.frame(case06), list(
  source = "hand-built: WT majority + edited minority + junk minority",
  editing_rate_truth = NA,
  description = "low-score/high-event junk minority; today findLQR correctly removes it",
  desired = list(should_filter = TRUE,
                 target_cluster = "junk",
                 reason = "minority junk must still be removed; the fix must not over-correct")
), filtered = f06)

# -- case 08: edited majority + junk (false-negative complement) ----------------
cat("\nCapturing synthetic case 08 (edited majority + junk):\n")
set.seed(11)
case08 <- rbind(
  mk_cluster(400,  690, 710, 1, "wt"),
  mk_cluster(2000, 620, 650, 2, "edited"),
  mk_cluster(500,  280, 340, sample(3:4, 500, replace = TRUE), "junk")
)
f08 <- current_filtered(as.data.frame(case08))
write_fixture("08_synthetic_boundary", as.data.frame(case08), list(
  source = "hand-built: edited majority + WT minority + junk minority",
  editing_rate_truth = NA,
  description = paste0("under high editing findLQR currently culls NOTHING,",
                       " i.e. it misses genuine junk (false-negative) -- the",
                       " complement of the over-culling bug. Desired: cull only junk."),
  desired = list(should_filter = TRUE,
                 target_cluster = "junk",
                 reason = "remove the small low-score junk cluster; the edited majority must survive")
), filtered = f08)

cat("\nDone. 8 fixtures under ", FIX_DIR, "\n", sep = "")
