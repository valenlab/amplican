# Run from the package root:
#   install.packages(c("pkgload", "testthat", "muttest"))
#   Rscript tools/mutation-pilot.R [report-directory] [target1,target2,...]
# muttest generates changes; this runner applies them only in memory, restoring
# both namespace and exported functions after each run. No source files change.
# Survivors need review. Test errors and setup errors do not count as kills.

for (package in c("pkgload", "testthat", "muttest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install the developer dependency: ", package)
  }
}
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else tempfile(
  "amplican-mutations-", tmpdir = dirname(tempdir()))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output)

suppressPackageStartupMessages(pkgload::load_all(
  ".", compile = FALSE, helpers = FALSE, export_all = FALSE, quiet = TRUE))
namespace <- asNamespace("amplican")
exports <- as.environment("package:amplican")
targets <- c(getEvents = "R/helpers_general.R",
             amplicanSummarize = "R/amplicanSummarize.R",
             amplicanNormalize = "R/amplicanNormalize.R",
             amplicanFilter = "R/amplicanFilter.R",
             findEOP = "R/helpers_filters.R", findPD = "R/helpers_filters.R",
             findLQR = "R/helpers_filters.R", goodBaseQuality = "R/helpers_filters.R",
             goodAvgQuality = "R/helpers_filters.R", alphabetQuality = "R/helpers_filters.R",
             amplicanOverlap = "R/amplicanSummarize.R",
             amplicanConsensus = "R/amplicanSummarize.R", getHits = "R/amplicanSummarize.R",
             is_hdr = "R/helpers_alignment.R", is_hdr_strict = "R/helpers_alignment.R",
             locate_pr_start = "R/helpers_alignment.R", get_seq = "R/helpers_general.R",
             get_left_primer = "R/helpers_general.R", get_right_primer = "R/helpers_general.R",
             upperGroups = "R/helpers_general.R", flipRanges = "R/helpers_general.R",
             amplicanMap = "R/helpers_general.R", getEventInfo = "R/helpers_general.R",
             extractEvents = "R/AlignmentsExperimentSet-class.R",
             makeAlignment = "R/helpers_alignment.R", amplicanAlign = "R/amplicanAlign.R",
             amplicanPipe = "R/amplican.R")
tests <- c(getEvents = "^getEvents$", amplicanSummarize = "^core_verbs$",
           amplicanNormalize = "^core_verbs$", amplicanFilter = "^core_verbs$",
           findEOP = "^(core_verbs|summary_helpers)$", findPD = "^core_verbs$",
           findLQR = "^(core_verbs|findLQR)$", goodBaseQuality = "^(core_verbs|filters_helpers)$",
           goodAvgQuality = "^(core_verbs|filters_helpers)$", alphabetQuality = "^(core_verbs|filters_helpers)$",
           amplicanOverlap = "^summary_helpers$", amplicanConsensus = "^summary_helpers$",
           getHits = "^(getHits|summary_helpers)$", is_hdr = "^hdr($|_)", is_hdr_strict = "^hdr($|_)")
upstream_helpers <- c("locate_pr_start", "get_seq", "get_left_primer", "get_right_primer",
                      "upperGroups", "flipRanges", "amplicanMap", "getEventInfo", "extractEvents")
tests[upstream_helpers] <- "^alignment_helpers$"
tests[c("get_seq", "flipRanges")] <- "^(alignment_helpers|summary_helpers|general_helpers|plot_helpers)$"
tests[c("get_left_primer", "get_right_primer")] <- "^plot_helpers$"
tests["amplicanMap"] <- "^general_helpers$"
tests["extractEvents"] <- "^amplican$"
tests[c("makeAlignment", "amplicanAlign", "amplicanPipe")] <- "^(amplican|AlignmentsExperimentSet)$"
regression_tests <- tests
if (file.exists("tests/testthat/test_upstream_helpers.R")) tests[upstream_helpers] <- "^upstream_helpers$"
if (file.exists("tests/testthat/test_upstream_assignment.R")) tests[c("makeAlignment", "amplicanAlign")] <- "^upstream_assignment$"
if (file.exists("tests/testthat/test_pipeline_invariants.R")) tests["amplicanPipe"] <- "^pipeline_invariants$"
# Fast cases run first; surviving changes also run the existing regression files.
if (file.exists("tests/testthat/test_hdr_decisions.R")) tests[c("is_hdr", "is_hdr_strict")] <- "^hdr_decisions$"
if (file.exists("tests/testthat/test_core_decisions.R")) {
  tests[c("amplicanNormalize", "amplicanFilter", "findEOP", "findPD", "findLQR",
          "goodBaseQuality", "goodAvgQuality", "alphabetQuality")] <- "^core_decisions$"
}
if (length(args) > 1L) {
  selected <- strsplit(args[2], ",", fixed = TRUE)[[1]]
  if (!all(selected %in% names(targets))) stop("Unknown mutation target: ", args[2])
  targets <- targets[selected]
}
bindings <- mget(unique(c(names(targets), "amplicanPipeline", "amplicanPipelineConservative")), namespace)
originals <- bindings[names(targets)]
if ("extractEvents" %in% names(targets)) {
  originals$extractEvents <- methods::selectMethod("extractEvents", "AlignmentsExperimentSet")
}
cat("Source:", normalizePath("."), "version",
    as.character(packageVersion("amplican")), "\nReports:", output, "\n")
writeLines(c(paste("Source:", normalizePath(".")), capture.output(sessionInfo())),
           file.path(output, "sessionInfo.txt"))

set_function <- function(name, fn) {
  for (env in list(namespace, exports)) {
    if (!exists(name, envir = env, inherits = FALSE)) next
    locked <- bindingIsLocked(name, env)
    if (locked) unlockBinding(name, env)
    assign(name, fn, envir = env)
    if (locked) lockBinding(name, env)
  }
}

run_case <- function(name, fn, filter = tests[[name]]) {
  probe <- new.env(parent = emptyenv())
  probe$calls <- 0L
  tracked <- function(...) {
    probe$calls <- probe$calls + 1L
    fn(...)
  }
  set_function(name, tracked)
  on.exit(set_function(name, bindings[[name]]))
  if (name == "amplicanPipe") {
    for (entry in c(amplicanPipeline = 0.01, amplicanPipelineConservative = 0.15)) {
      # Each closure needs its own function binding, rather than the loop's last value.
      pipeline <- local({
        frequency <- entry
        implementation <- fn(frequency)
        function(...) {
          probe$calls <- probe$calls + 1L
          implementation(...)
        }
      })
      set_function(if (entry == 0.01) "amplicanPipeline" else "amplicanPipelineConservative", pipeline)
    }
    on.exit({
      set_function("amplicanPipeline", bindings$amplicanPipeline)
      set_function("amplicanPipelineConservative", bindings$amplicanPipelineConservative)
    }, add = TRUE)
  }
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(testthat::test_dir(
    "tests/testthat", filter = filter, reporter = "silent",
    stop_on_failure = FALSE, package = "amplican", load_package = "none"),
    error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - start
  if (inherits(result, "error")) {
    return(data.frame(status = "setup_error", executions = probe$calls, assertions = 0L, failed = 0L,
                      errors = 1L, warnings = 0L, elapsed = elapsed,
                      detail = conditionMessage(result)))
  }
  counts <- colSums(as.data.frame(result)[, c("nb", "failed", "error", "warning")])
  status <- if (counts["error"] && counts["failed"]) "mixed_failure" else if (counts["error"]) "test_error" else if (counts["failed"]) {
    "killed"
  } else if (!counts["nb"]) "no_tests" else "survived"
  failures <- unlist(lapply(result, function(test) test$results), recursive = FALSE)
  failures <- Filter(function(x) inherits(x, if (counts["error"]) "expectation_error" else "expectation_failure"),
                     failures)
  data.frame(status = status, executions = probe$calls, assertions = unname(counts["nb"]),
             failed = unname(counts["failed"]), errors = unname(counts["error"]),
             warnings = unname(counts["warning"]), elapsed = elapsed,
             detail = if (length(failures)) conditionMessage(failures[[1]]) else "")
}

report <- list()
record <- function(name, kind, mutation, line, before, after, fn) {
  result <- run_case(name, fn)
  if (result$status == "survived" && regression_tests[[name]] != tests[[name]]) {
    regression <- run_case(name, fn, filter = regression_tests[[name]])
    regression$assertions <- regression$assertions + result$assertions
    regression$elapsed <- regression$elapsed + result$elapsed
    regression$warnings <- regression$warnings + result$warnings
    regression$executions <- regression$executions + result$executions
    result <- regression
  }
  if (kind == "baseline" && result$status == "survived") result$status <- "passed"
  if (kind != "baseline" && result$status == "survived" && result$executions == 0L) result$status <- "not_exercised"
  report[[length(report) + 1L]] <<- cbind(
    data.frame(target = name, kind = kind, mutation = mutation, line = line,
               original = before, replacement = after), result)
  write.csv(do.call(rbind, report), file.path(output, "mutations.csv"), row.names = FALSE)
  cat(sprintf("%-18s %-30s %-12s %.2fs\n", name, mutation, result$status, result$elapsed))
  result
}

for (name in names(targets)) {
  baseline <- record(name, "baseline", "unmodified source", NA_integer_, "", "",
                     originals[[name]])
  stopifnot(baseline$status == "passed", baseline$warnings == 0L)
  if (file.exists("tests/testthat/test_upstream_helpers.R")) stopifnot(baseline$executions > 0L)
}

seed <- function(name, label, replacements, kind = "known_fault") {
  if (!name %in% names(targets)) return(invisible(NULL))
  fn <- originals[[name]]
  code <- paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n")
  before <- code
  for (from in names(replacements)) {
    if (!grepl(from, code, fixed = TRUE)) stop("Fault seed needs updating: ", label)
    code <- sub(from, replacements[[from]], code, fixed = TRUE)
  }
  stopifnot(code != before)
  body(fn) <- parse(text = code)[[1]]
  record(name, kind, label, NA_integer_, before, code, fn)
}
seed("getEvents", "historical mismatch labels", c(
  'comparison[comparison != "+"]' = "comparison",
  'Biostrings::replaceAt(pattern, ins, value = "")' = "pattern",
  'Biostrings::replaceAt(subject, ins, value = "")' = "subject"))
seed("amplicanSummarize", "inverted frameshift predicate", setNames(
  deparse(quote(sum(net_width) %% 3 == 0)), deparse(quote(sum(net_width) %% 3 != 0))))
seed("amplicanNormalize", "cross-group normalization", c(
  'if (is.null(add))' = 'add <- setdiff(add, "Group"); if (is.null(add))'))
if ("amplicanNormalize" %in% names(targets) && grepl("aln_ctr_counts", paste(deparse(body(originals$amplicanNormalize)), collapse = " "))) {
  seed("amplicanNormalize", "double-count control strands", c(
    'unique(aln_ctr, by = unique(c(cols, "seqnames", "read_id")))' = 'aln_ctr'))
}
seed("amplicanFilter", "read ID without experiment", c(
  'on = .(seqnames, read_id)' = 'on = .(read_id)'))
if ("findEOP" %in% names(targets) && grepl("aln$start <= 0", paste(deparse(body(originals$findEOP)), collapse = " "), fixed = TRUE)) {
  seed("findEOP", "miss zero-start relative events", c('aln$start <= 0' = 'aln$start < 0'))
}
if ("findLQR" %in% names(targets) && grepl("sum(events)", paste(deparse(body(originals$findLQR)), collapse = " "), fixed = TRUE)) {
  seed("findLQR", "per-strand removal fraction", c(
    'sum(events)/length(unique(strand))' = 'events/length(unique(strand))',
    'sum(score)/length(unique(strand))' = 'score/length(unique(strand))'))
}
seed("amplicanConsensus", "prefer lower score", c(
  'aln_fwd$score[fi] >= aln_rve$score[ri]' = 'aln_fwd$score[fi] <= aln_rve$score[ri]'))
if ("amplicanConsensus" %in% names(targets) && grepl("only_b_fwd", paste(deparse(body(originals$amplicanConsensus)), collapse = " "))) {
  seed("amplicanConsensus", "asymmetric clipped-read sets", c(
    'only_b_fwd <- b_fwd[!b_rve, on = .(seqnames, read_id)]' = 'only_b_fwd <- b_fwd',
    'b_fwd <- only_b_fwd' = 'b_fwd <- only_b_fwd[!b_rve, on = .(seqnames, read_id)]'))
}
seed("is_hdr_strict", "accept partial donor", c(
  'hits$n == length(hdr_events)' = 'hits$n > 0'))
# After the overlap fix, replay its previous buffered-origin calculation.
if ("amplicanOverlap" %in% names(targets) && grepl("origins", paste(deparse(body(originals$amplicanOverlap)), collapse = " "))) {
  seed("amplicanOverlap", "buffered relative origin", c(
    'origins[i]' = 'IRanges::start(cutSites[[i]])[1]'))
}

if (file.exists("tests/testthat/test_upstream_helpers.R")) {
  seed("extractEvents", "wrong read multiplicity", c(
    'r_counts[[ID]][as.integer(tempGR$read_id)]' = 'rep(r_counts[[ID]][1], length(tempGR))'), "seeded_fault")
  seed("extractEvents", "swap mate HDR flags", c('fRT <- fwd_types[[ID]]' = 'fRT <- rve_types[[ID]]'), "seeded_fault")
  seed("extractEvents", "skip orientation conversion", c(
    'flipRanges(GenomicRanges::as.data.frame(finalGR, row.names = NULL), exp_data)' =
      'GenomicRanges::as.data.frame(finalGR, row.names = NULL)'), "seeded_fault")
  seed("makeAlignment", "require both primers in permissive mode", c(
    'unqT$forwardFound | unqT$reverseFound' = 'unqT$forwardFound & unqT$reverseFound'), "seeded_fault")
  seed("makeAlignment", "accept either primer in strict mode", c(
    'unqT$forwardFound & unqT$reverseFound' = 'unqT$forwardFound | unqT$reverseFound'), "seeded_fault")
  seed("makeAlignment", "omit primer trimming", c('start = IDunqT$fwdPrInReadPos' = 'start = 1L'), "seeded_fault")
  seed("makeAlignment", "omit reverse complementation", c(
    'rR <- Biostrings::reverseComplement(' = 'rR <- identity('), "seeded_fault")
  seed("makeAlignment", "lose duplicates across batches", c('Total = sum(Total)' = 'Total = max(Total)'), "seeded_fault")
  seed("makeAlignment", "collapse distinct reverse mates", c(
    '.(Total = .N), by = .(Forward, Reverse)' = '.(Total = .N, Reverse = Reverse[1]), by = .(Forward)'), "seeded_fault")
  seed("amplicanAlign", "lose configuration order", c('id_order <- match(ids, cfgT$ID)' = 'id_order <- seq_len(nrow(cfgT))'), "seeded_fault")
  seed("amplicanPipe", "skip off-target filtering", c('if (event_filter)' = 'if (FALSE)'), "seeded_fault")
  seed("amplicanPipe", "skip primer-dimer removal", c(
    'aln <- aln[!onlyPD, on = list(seqnames, read_id)]' = 'aln <- aln'), "seeded_fault")
  seed("amplicanPipe", "skip normalization", c(
    'aln <- amplicanNormalize(aln, cfgT, min_freq = min_freq, add = normalize)' = 'aln <- aln'), "seeded_fault")
  seed("amplicanPipe", "summarize both mates", c('aln$consensus & aln$overlaps' = 'aln$overlaps'), "seeded_fault")
}

mutators <- c(muttest::comparison_operators(), muttest::arithmetic_operators(),
              muttest::index_mutations(), muttest::logical_operators(), muttest::boolean_literals(),
              list(muttest::call_name("any", "all"), muttest::call_name("all", "any"),
                   muttest::call_name("sum", "max"), muttest::call_name("max", "sum")))
for (name in names(targets)) {
  if (name %in% c("makeAlignment", "amplicanAlign", "amplicanPipe")) next
  source <- parse(targets[[name]], keep.source = TRUE)
  definition <- which(vapply(source, function(x) {
    is.call(x) && if (name == "extractEvents") {
      identical(x[[1]], as.name("setMethod")) && identical(x[[2]], "extractEvents")
    } else identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
  }, logical(1)))
  span <- attr(source, "srcref")[[definition]]
  plan <- muttest::muttest_plan(mutators, source_files = targets[[name]])
  keep <- vapply(plan$mutation, function(m) {
    m$location$start$line >= span[1] && m$location$end$line <= span[3]
  }, logical(1))
  plan <- plan[keep, ]
  cat("Generated", nrow(plan), "mutations for", name, "\n")
  for (i in seq_len(nrow(plan))) {
    mutated <- parse(text = plan$mutated_code[[i]])
    fn <- eval(mutated[[definition]][[if (name == "extractEvents") 4L else 3L]], envir = namespace)
    stopifnot(!identical(c(deparse(formals(fn)), deparse(body(fn))),
                        c(deparse(formals(originals[[name]])), deparse(body(originals[[name]])))))
    line <- plan$mutation[[i]]$location$start$line
    record(name, "generated", paste(plan$mutator[[i]]$from, "to", plan$mutator[[i]]$to),
           line, plan$original_code[[i]][line], plan$mutated_code[[i]][line], fn)
  }
}

results <- do.call(rbind, report)
print(with(results, table(kind, status)))
cat("Results:", file.path(output, "mutations.csv"), "\n")
cat("Survivors were rechecked against existing regressions. Review survivors and errors separately.\n")
