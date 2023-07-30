library(testthat)
context("amplican main making analysis of example files")

# un-comment line below this comment and below test_that to make vignettes
# after any changes to rmarkdown templates
# changes paths so that system.file points to the actual development folder
# devtools::load_all()

config <- system.file("extdata", "config.csv", package = "amplican")
fastq_folder <- system.file("extdata", package = "amplican")
results_folder <- tempdir()
# results_folder <- system.file("extdata", "results", package = "amplican")
# vignettes_path <- system.file('vignettes', package = 'amplican')
#
# devtools::unload()
library(amplican)

test_that("amplican theme teplates are loaded properly", {
  expect_gte(
    nchar(system.file("rmarkdown", "templates",
                      "id_report", package = "amplican")), 1)
})

test_that("amplican runs through example files without any issues", {
  expect_warning(
    amplicanPipeline(config, fastq_folder, results_folder,
                     knit_reports = FALSE,
                     primer_mismatch = 0,
                     fastqfiles = 0,
                     continue = FALSE))
})

# # after successful amplican pipeline run, replace auto generated
# # system links in .Rmd reports with system.file type of links
# # also make vignettes out of .Rmd reports
# unlink(file.path(results_folder, "reports"), recursive = TRUE, force = FALSE)
#
# # run this when developing and willing to generate new vignettes
# # for visual quality test and for the users to see example reports
# # fix run param file
# runParam <- readLines(file.path(results_folder, "RunParameters.txt"))
# runParam <- gsub(config,
#                  "full/path/to/config/file/that/has/been/used.csv",
#                  runParam)
# cat(runParam,
#     file = file.path(results_folder, "RunParameters.txt"), sep = "\n")
# # fix paths to the summary config
# conf_file <- system.file('extdata',
#                          'results',
#                          'config_summary.csv',
#                          package = 'amplican')
# conf_df <- read.csv(conf_file,
#                     stringsAsFactors = FALSE)
# conf_df$Reverse_Reads_File <- basename(conf_df$Reverse_Reads_File)
# conf_df$Forward_Reads_File <- basename(conf_df$Forward_Reads_File)
# write.csv(conf_df, file = conf_file, row.names = FALSE)
#
# template_names <- c("id_report", "barcode_report", "group_report",
#                     "guide_report", "amplicon_report", "index_report")
# rmd_names <- c("id_report.Rmd", "barcode_report.Rmd", "group_report.Rmd",
#                "guide_report.Rmd", "amplicon_report.Rmd", "index.Rmd")
#
# for (i in seq_along(rmd_names)) {
#   rmd_file_path <- file.path(vignettes_path, paste0("example_", rmd_names[i]))
#   if (file.exists(rmd_file_path)) {
#     file.remove(rmd_file_path)
#   }
#   rmarkdown::draft(file = rmd_file_path,
#                    template = template_names[i],
#                    package = "amplican",
#                    edit = FALSE)
#
#   rmd_content <- readLines(rmd_file_path)
#   settings_end <- which(rmd_content == "---")[2]-1
#   fixed_rmd_content <- c(rmd_content[1:settings_end],
#                          "vignette: >",
#                          paste0("  %\\VignetteIndexEntry{example ",
#                                 sub("report_", "",
#                                     tools::file_path_sans_ext(rmd_names[i])),
#                                 " report}"),
#                          "  %\\VignetteEngine{knitr::rmarkdown}",
#                          "  %\\VignetteEncoding{UTF-8}",
#                          rmd_content[(settings_end+1):length(rmd_content)])
#
#   if (rmd_names[i] == "index.Rmd") {
#     fixed_rmd_content <-
#       c(fixed_rmd_content[1:13],
#         paste0('  links: \"',
#                '1. [Report by id](./example_id_report.html)\\n',
#                '2. [Report by barcode](./example_barcode_report.html)\\n',
#                '3. [Report by group](./example_group_report.html)\\n',
#                '4. [Report by guide](./example_guide_report.html)\\n',
#                '5. [Report by amplicon](./example_amplicon_report.html)\\n\"'),
#         fixed_rmd_content[15:length(fixed_rmd_content)])
#   }
#   cat(fixed_rmd_content, file = rmd_file_path, sep = "\n")
# }
#
# devtools::build_vignettes()
