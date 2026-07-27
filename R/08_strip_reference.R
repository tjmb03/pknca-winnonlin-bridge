#!/usr/bin/env Rscript
#' Remove reference-derived values from the committed platform records
#'
#' The records under reproducibility/ include a WinNonlin column and the
#' percent-difference columns computed from it. Those are numeric NCA
#' parameters from a source repository that declares no licence.
#'
#' This strips them, leaving only values this project computed itself. The
#' cross-platform figure (fig4) is unaffected -- it compares PKNCA against
#' PKNCA and NonCompart against NonCompart across machines. Figure 3 loses its
#' engines-vs-reference point and falls back to a note.
#'
#'   Rscript R/08_strip_reference.R

dirs <- list.dirs("reproducibility", recursive = FALSE)
if (length(dirs) == 0) stop("no reproducibility/ records found")

drop_cols <- c("WinNonlin", "pct_diff_PKNCA", "pct_diff_NonCompart",
               "max_abs_pct", "max_observed_pct")

for (d in dirs) {
  for (f in list.files(d, pattern = "\\.csv$", full.names = TRUE)) {
    x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    hit <- intersect(names(x), drop_cols)
    if (length(hit) == 0) next
    x <- x[, setdiff(names(x), hit), drop = FALSE]
    write.csv(x, f, row.names = FALSE)
    message("stripped ", paste(hit, collapse = ", "), " from ", f)
  }
}
message("\nDone. Re-run: Rscript R/07_figures.R")
