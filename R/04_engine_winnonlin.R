#' Phoenix WinNonlin arm
#'
#' Two ways to get WinNonlin numbers into this comparison:
#'
#'   (A) Published reference. A Phoenix WinNonlin 6.3/7.0 output for
#'       datasets::Theoph is public. Fetched at runtime. Requires no license and
#'       lets anyone reproduce the three-way comparison.
#'
#'   (B) Your own Phoenix run. Export "Final Parameters Pivoted" from a Phoenix
#'       NCA object to CSV, drop it in phoenix/reference/, and this module reads
#'       it with the same crosswalk. See phoenix/PHOENIX_SETUP.md for the exact
#'       GUI settings that correspond to NCA_RULES.
#'
#' Note on automation: Phoenix does have a documented command-line mode, but it
#' drives the NLME engine, not the NCA object. Scripted NCA runs go through the
#' separately-licensed AutoPilot Toolkit, or through Phoenix's built-in R Script
#' object inside a workflow. There is no supported `Phoenix.exe --run-nca`
#' equivalent. Plan the architecture around export/import, not shelling out.

source(file.path("R", "00_config.R"))

# ---- (A) Published reference ------------------------------------------------

fetch_winnonlin_reference <- function(which = "theoph_log",
                                      cache_dir = PATHS$phoenix_rf,
                                      force = FALSE) {

  url <- WINNONLIN_REFERENCE[[which]]
  if (is.null(url)) stop("Unknown reference: ", which,
                         ". Available: ",
                         paste(names(WINNONLIN_REFERENCE), collapse = ", "))

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(cache_dir, paste0("winnonlin_", which, ".csv"))

  if (!file.exists(dest) || force) {
    message("Fetching WinNonlin reference: ", which)
    ok <- tryCatch({
      utils::download.file(url, dest, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) {
      warning("Could not fetch reference (", conditionMessage(e),
              "). The three-way comparison will fall back to two-way.")
      FALSE
    })
    if (!ok) return(NULL)
  } else {
    message("Using cached WinNonlin reference: ", basename(dest))
  }

  read_winnonlin_csv(dest)
}

# ---- (B) Your own Phoenix export --------------------------------------------

read_winnonlin_csv <- function(path, subject_col = "Subject") {

  if (!file.exists(path)) {
    warning("No Phoenix output at ", path,
            " -- skipping the WinNonlin arm.")
    return(NULL)
  }

  w <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)

  if (!subject_col %in% names(w)) {
    # Phoenix sort keys are named after the source column. Try common variants
    # before giving up, since this is the most likely import friction point.
    candidates <- intersect(c("Subject", "ID", "USUBJID", "SUBJID"), names(w))
    if (length(candidates) == 0)
      stop("Cannot find a subject column in ", basename(path),
           ". Columns present: ", paste(names(w), collapse = ", "))
    subject_col <- candidates[1]
    message("Using '", subject_col, "' as the subject key")
  }

  param_cols <- intersect(names(w), PARAM_MAP$winnonlin)
  if (length(param_cols) == 0)
    stop("No recognised WinNonlin parameter columns in ", basename(path),
         ". Expected some of: ",
         paste(head(PARAM_MAP$winnonlin, 5), collapse = ", "), " ...")

  message("Read ", length(param_cols), " parameters for ",
          nrow(w), " subjects from ", basename(path))

  out <- do.call(rbind, lapply(param_cols, function(p) {
    data.frame(
      Subject  = w[[subject_col]],
      PPTESTCD = PARAM_MAP$PPTESTCD[match(p, PARAM_MAP$winnonlin)],
      value    = suppressWarnings(as.numeric(w[[p]])),
      engine   = "WinNonlin",
      stringsAsFactors = FALSE
    )
  }))

  out$Subject <- suppressWarnings(as.numeric(as.character(out$Subject)))
  rownames(out) <- NULL
  out
}

# ---- Export a Phoenix-ready input file --------------------------------------

export_for_phoenix <- function(phoenix_df,
                               path = file.path(PATHS$phoenix_in,
                                                "theoph_for_phoenix.csv")) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  # Phoenix is fussy about two things on import: it will not accept column
  # names with spaces or symbols in the mapping dialog without complaint, and
  # it treats an empty cell differently from a zero. Write clean and explicit.
  df <- phoenix_df
  names(df) <- gsub("[^A-Za-z0-9_]", "_", names(df))

  utils::write.csv(df, path, row.names = FALSE, na = "")
  message("Phoenix input written: ", path)
  message("Import this in Phoenix, then follow phoenix/PHOENIX_SETUP.md")
  invisible(path)
}
