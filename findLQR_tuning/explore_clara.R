# explore_clara.R --------------------------------------------------------------
#
# Workbench for tuning findLQR(). For each fixture it reconstructs the exact
# 2-D feature space findLQR sees -- (range01(score), range01(events/read_len)) --
# re-runs clara at k=2 and k=3, prints the silhouette decision, and draws the
# clusters + medoids + the cluster that the CURRENT rule would cull. This makes
# the failure mode visible: at high editing the edited reads form a separable
# cluster that the geometry-only rule deletes.
#
# Usage:
#   Rscript explore_clara.R [case_name]      # one case (opens/saves its plot)
#   Rscript explore_clara.R                   # all cases -> summary table + plots
# Plots are written to findLQR_tuning/plots/<case>.png

suppressMessages({library(data.table); library(cluster); library(ggplot2)})

AMP_PATH <- (function() {
  cand <- tryCatch({ a <- commandArgs(trailingOnly = FALSE)
    f <- sub("--file=", "", a[grep("--file=", a)][1]); if (is.na(f)) stop(); dirname(normalizePath(f))
  }, error = function(e) getwd())
  for (i in 1:6) if (file.exists(file.path(cand, "DESCRIPTION"))) return(normalizePath(cand)) else cand <- dirname(cand)
  normalizePath(".", mustWork = TRUE)
})()
FIX_DIR <- file.path(AMP_PATH, "findLQR_tuning", "fixtures")
PLOT_DIR <- file.path(AMP_PATH, "findLQR_tuning", "plots")
dir.create(PLOT_DIR, showWarnings = FALSE)

range01 <- function(x) { nx <- (x - min(x)) / diff(range(x)); nx[!is.finite(nx)] <- 0; nx }

# Mirror findLQR's per-read feature aggregation exactly.
feature_frame <- function(aln) {
  an <- as.data.table(aln)[, list(events = .N / max(end), score = max(score)),
                           by = c("read_id", "strand", "seqnames")]
  an <- an[, list(events = events / length(unique(strand)),
                  score = score / length(unique(strand))),
           by = c("read_id", "seqnames")]
  an[, `:=`(nscore = range01(score), nevents = range01(events))]
  an
}

# Run the k2-vs-k3 decision and report everything a tuner needs.
inspect_case <- function(aln, case_name) {
  an <- feature_frame(aln)
  if (nrow(an) < 1000) {
    out <- list(case = case_name, n = nrow(an), decision = "no filter (n<1000)")
    return(invisible(out))
  }
  ss <- min(1000, nrow(an))
  k2 <- clara(cbind(an$nscore, an$nevents), 2, samples = 500, sampsize = ss)
  k3 <- clara(cbind(an$nscore, an$nevents), 3, samples = 500, sampsize = ss)
  k2s <- mean(silhouette(k2)[, 3]); k3s <- mean(silhouette(k3)[, 3])
  fires <- k3s > k2s
  centers <- apply(k3$medoids, 1, function(z) sqrt((z[1] - 1) ^ 2 + z[2] ^ 2))
  cut_cl <- if (fires) which.max(centers) else NA
  an$cluster <- if (fires) k3$clustering else k2$clustering
  # how big is each cluster (fraction of reads)? -- motivates the size-guard idea
  sizes <- an[, .N, by = cluster][order(cluster)]
  sizes[, frac := N / sum(N)]
  if (fires) {
    culled_n <- sizes[cluster == cut_cl, N]
    culled_frac <- sizes[cluster == cut_cl, frac]
  } else { culled_n <- 0L; culled_frac <- 0 }

  cat(sprintf("\n=== %s === (n=%d reads)\n", case_name, nrow(an)))
  cat(sprintf("  silhouette: k2=%.4f  k3=%.4f  -> %s\n", k2s, k3s,
              if (fires) sprintf("FILTER (cull cluster %d)", cut_cl) else "no filter"))
  if (fires) {
    print(sizes)
    cat(sprintf("  -> culling %d reads (%.1f%% of all reads)\n", culled_n, culled_frac * 100))
  }

  # plot
  df <- as.data.frame(an)
  df$cluster <- factor(df$cluster)
  df$cut <- fires & (an$cluster == cut_cl)
  p <- ggplot(df, aes(nscore, nevents, colour = cluster)) +
    geom_point(alpha = 0.35, size = 1.2) +
    geom_point(data = as.data.frame(k3$medoids), aes(V1, V2), shape = 8, size = 5,
               stroke = 1.5, colour = "black", inherit.aes = FALSE) +
    scale_colour_grey(start = 0.75, end = 0.1) +
    labs(title = sprintf("%s  | k2=%.3f k3=%.3f -> %s",
                         case_name, k2s, k3s, if (fires) "FILTER" else "no filter"),
         x = "range01(alignment score)  <- high is good",
         y = "range01(events / read length)  <- high is 'busy'") +
    theme_bw() +
    annotate("text", x = 0.5, y = 0.97,
             label = if (fires) sprintf("culled cluster %d: %d reads (%.0f%%)",
                                        cut_cl, culled_n, culled_frac * 100) else "nothing culled",
             hjust = 0.5)
  ggsave(file.path(PLOT_DIR, paste0(case_name, ".png")), p, width = 6, height = 5, dpi = 120)

  invisible(list(case = case_name, n = nrow(an), k2s = k2s, k3s = k3s,
                 fires = fires, culled_frac = culled_frac))
}

cases <- list.dirs(FIX_DIR, recursive = FALSE, full.names = FALSE)
target <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(target)) cases <- cases[cases == target]
cat("findLQR workbench --", length(cases), "case(s)\n")
for (cs in cases) inspect_case(as.data.frame(fread(file.path(FIX_DIR, cs, "events.csv"))), cs)
cat("\nPlots written to ", PLOT_DIR, "\n", sep = "")
