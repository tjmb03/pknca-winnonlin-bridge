#!/usr/bin/env Rscript
#' Run the full comparison
#'
#'   Rscript run_all.R              # Theophylline validation track
#'   Rscript run_all.R --cdisc      # also run the CDISC SDTM -> PP track
#'   Rscript run_all.R --no-fetch   # skip the network fetch of the reference
#'   Rscript run_all.R --figures    # force-rebuild the committed README figures
#'
#' Outputs land in outputs/.

args      <- commandArgs(trailingOnly = TRUE)
do_cdisc  <- "--cdisc"    %in% args
no_fetch  <- "--no-fetch" %in% args

source(file.path("R", "00_config.R"))
source(file.path("R", "01_prepare_theoph.R"))
source(file.path("R", "02_engine_pknca.R"))
source(file.path("R", "03_engine_noncompart.R"))
source(file.path("R", "04_engine_winnonlin.R"))
source(file.path("R", "05_compare.R"))

# Fail fast if the two R engines are configured inconsistently.
assert_rules_consistent()

# ---------------------------------------------------------------- Track A ----
message("\n=== Track A: Theophylline engine validation ===\n")

theoph <- prepare_theoph()
validate_theoph(theoph)

write.csv(theoph$conc, file.path(PATHS$derived, "theoph_conc.csv"),
          row.names = FALSE)
write.csv(theoph$dose, file.path(PATHS$derived, "theoph_dose.csv"),
          row.names = FALSE)
export_for_phoenix(theoph$phoenix)

message("Running PKNCA ...")
res_pknca <- run_pknca(theoph$conc, theoph$dose)

message("Running NonCompart ...")
res_ncmp <- run_noncompart(theoph$conc, dose_mg = THEOPH_DOSE_MG)

# WinNonlin arm: prefer a local export, fall back to the published reference.
local_exports <- list.files(PATHS$phoenix_rf, pattern = "\\.csv$",
                            full.names = TRUE)
local_exports <- setdiff(local_exports,
                         list.files(PATHS$phoenix_rf,
                                    pattern = "^winnonlin_theoph_(linear|log)\\.csv$",
                                    full.names = TRUE))

res_wnl <- NULL
if (length(local_exports) > 0) {
  message("Reading local Phoenix export: ", basename(local_exports[1]))
  res_wnl <- read_winnonlin_csv(local_exports[1])
} else if (!no_fetch) {
  which_ref <- if (identical(NCA_RULES$noncompart_down, "Log"))
    "theoph_log" else "theoph_linear"
  res_wnl <- fetch_winnonlin_reference(which_ref)
}

if (is.null(res_wnl))
  message("No WinNonlin arm available; running PKNCA vs NonCompart only.")

# ---- Compare ----------------------------------------------------------------
cmp <- compare_engines(res_pknca, res_ncmp, res_wnl,
                       reference = if (is.null(res_wnl)) "NonCompart"
                                   else "WinNonlin")
agree <- summarise_agreement(cmp)
report_verdict(agree)

write.csv(cmp,   file.path(PATHS$outputs, "comparison_per_subject.csv"),
          row.names = FALSE)
write.csv(agree, file.path(PATHS$outputs, "comparison_summary.csv"),
          row.names = FALSE)

# ---- Precision attribution ---------------------------------------------------
# Separate "how closely do the engines agree with each other" from "how finely
# was the reference published". Conflating the two overstates engine error by
# about six orders of magnitude on this dataset.
ee <- engine_vs_engine(res_pknca, res_ncmp)
cat("\nENGINE vs ENGINE (no reference file involved)\n")
cat(sprintf("  PKNCA vs NonCompart : max |diff| %.3e %%  over %d comparisons\n",
            ee$max_pct_diff, ee$n))
cat(sprintf("  bit-identical       : %d/%d\n", ee$n_bit_identical, ee$n))
cat(sprintf("  double-precision floor: %.1e %%\n", ee$double_eps_pct))

ref_csv <- if (length(local_exports) > 0) local_exports[1] else
  file.path(PATHS$phoenix_rf,
            paste0("winnonlin_", if (identical(NCA_RULES$noncompart_down, "Log"))
                   "theoph_log" else "theoph_linear", ".csv"))

if (!is.null(res_wnl) && file.exists(ref_csv)) {
  prec <- attribute_precision(cmp, ref_csv, engine = "PKNCA")
  if (!is.null(prec)) {
    write.csv(prec, file.path(PATHS$outputs, "precision_attribution.csv"),
              row.names = FALSE)
    n_lim <- sum(prec$precision_limited, na.rm = TRUE)
    cat(sprintf(
      "\nPRECISION ATTRIBUTION vs reference\n  %d/%d parameters are precision-limited\n",
      n_lim, nrow(prec)))
    cat("  (observed difference <= what the reference's stored decimals can resolve)\n")
    worst <- prec[which.max(prec$max_observed_pct), ]
    cat(sprintf("  largest observed: %s  %.3e %%  vs rounding bound %.3e %% (%s dp)\n",
                worst$PPTESTCD, worst$max_observed_pct,
                worst$max_rounding_pct, worst$ref_decimals))
  }
}

ref <- attr(cmp, "reference")
for (e in setdiff(attr(cmp, "engines"), ref)) {
  plot_agreement(cmp, engine = e, reference = ref,
                 path = file.path(PATHS$outputs,
                                  paste0("agreement_", e, "_vs_", ref, ".png")))
}

# ---------------------------------------------------------------- Track B ----
if (do_cdisc) {
  message("\n=== Track B: CDISC SDTM -> NCA -> PP ===\n")
  ok <- tryCatch({
    source(file.path("R", "06_cdisc_track.R"))
    prep    <- prepare_cdisc_pk()
    analyte <- names(sort(table(prep$conc$PCTESTCD), decreasing = TRUE))[1]
    message("Analyte: ", analyte)

    nca_cdisc <- run_cdisc_nca(prep, analyte = analyte)
    pp        <- to_pp_domain(nca_cdisc)

    write.csv(pp, file.path(PATHS$outputs, "pp_domain.csv"), row.names = FALSE)
    message("PP-style dataset written: ", nrow(pp), " records, ",
            length(unique(pp$USUBJID)), " subjects")
    TRUE
  }, error = function(e) {
    message("CDISC track skipped: ", conditionMessage(e))
    FALSE
  })
}

# README figures are committed to the repository, and PNG output is not
# byte-reproducible (encoder metadata differs run to run). Regenerating them on
# every pipeline run would dirty the working tree every time and block the
# push pre-flight check. So: only build them when they are missing, or when
# explicitly asked with --figures.
figs <- file.path("figures",
                  c("fig1_agreement.png", "fig2_precision.png",
                    "fig3_scales.png", "fig4_platform.png"))
if (("--figures" %in% args) || !all(file.exists(figs))) {
  if (file.exists(file.path("reproducibility", "linux_x86_64",
                            "comparison_per_subject.csv"))) {
    tryCatch(source(file.path("R", "07_figures.R")),
             error = function(e) message("Figure generation skipped: ",
                                         conditionMessage(e)))
  }
}

capture_provenance()
message("\nDone. See outputs/\n")
