library(testthat)

# Tests assume they are run from the repo root:
#   testthat::test_dir("tests/testthat")

root <- normalizePath(file.path("..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R"))) root <- normalizePath(".")
withr_wd <- function(expr) { old <- setwd(root); on.exit(setwd(old)); expr }

withr_wd({
  source(file.path("R", "00_config.R"))
  source(file.path("R", "01_prepare_theoph.R"))
  source(file.path("R", "03_engine_noncompart.R"))
})

# ---- Configuration integrity ------------------------------------------------

test_that("AUC method crosswalk is consistent between engines", {
  withr_wd(expect_true(assert_rules_consistent()))
})

test_that("crosswalk has no duplicate or missing entries", {
  expect_equal(anyDuplicated(PARAM_MAP$PPTESTCD), 0)
  expect_equal(anyDuplicated(PARAM_MAP$pknca), 0)
  expect_equal(anyDuplicated(PARAM_MAP$winnonlin), 0)
  expect_false(any(is.na(PARAM_MAP$PPTESTCD)))
  expect_false(any(PARAM_MAP$label == ""))
})

test_that("crosswalk rows are all the same length", {
  lens <- vapply(PARAM_MAP, length, integer(1))
  expect_equal(length(unique(lens)), 1L)
})

test_that("tolerance is a sane fraction, not a percentage", {
  expect_gt(NCA_RULES$tolerance, 0)
  expect_lt(NCA_RULES$tolerance, 0.05)
})

# ---- Data preparation -------------------------------------------------------

test_that("Theophylline prepares to the expected shape", {
  d <- withr_wd(prepare_theoph())
  expect_equal(length(unique(d$conc$Subject)), 12L)
  expect_equal(nrow(d$conc), 132L)
  expect_equal(nrow(d$dose), 12L)
  expect_true(all(d$dose$Dose == THEOPH_DOSE_MG))
})

test_that("Subject is coerced through character, not factor level", {
  d <- withr_wd(prepare_theoph())
  # Theoph's Subject factor is ordered by peak concentration, so level index
  # and printed label differ. Subject 1 must have Cmax 10.50 -- if the coercion
  # were done via as.numeric(factor) this would silently be a different subject.
  s1 <- d$conc[d$conc$Subject == 1, ]
  expect_equal(max(s1$conc), 10.50, tolerance = 1e-8)
  expect_equal(s1$Time[which.max(s1$conc)], 1.12, tolerance = 1e-8)
})

test_that("validation catches malformed input", {
  d <- withr_wd(prepare_theoph())
  d$conc <- d$conc[d$conc$Subject <= 6, ]
  expect_error(validate_theoph(d), "12 subjects")
})

# ---- Known reference values -------------------------------------------------
# Anchored on the published Phoenix WinNonlin 6.3/7.0 output for Theoph.
# These are the numbers the whole repo exists to reproduce.

WNL_THEOPH_LOG <- data.frame(
  Subject = 1:12,
  CMAX    = c(10.50, 8.33, 8.20, 8.60, 11.40, 6.44,
              7.09, 7.56, 9.03, 10.21, 8.00, 9.75),
  TMAX    = c(1.12, 1.92, 1.02, 1.07, 1.00, 1.15,
              3.48, 2.02, 0.63, 3.55, 0.98, 3.52),
  LAMZHL  = c(14.3043776, 6.6593416, 6.7660874, 6.9812467,
              8.0022640, 7.8949979, 7.8466683, 8.5100379,
              8.4059988, 9.2469158, 7.2612365, 6.2865082),
  AUCLST  = c(147.2347485, 88.7312755, 95.8781978, 102.6336232,
              118.1793538, 71.6970150, 87.9692274, 86.8065635,
              83.9374360, 135.5760701, 77.8934723, 115.2202082)
)

test_that("observed Cmax and Tmax match the WinNonlin reference exactly", {
  d <- withr_wd(prepare_theoph())
  for (s in WNL_THEOPH_LOG$Subject) {
    sub <- d$conc[d$conc$Subject == s, ]
    i   <- which.max(sub$conc)
    expect_equal(sub$conc[i], WNL_THEOPH_LOG$CMAX[WNL_THEOPH_LOG$Subject == s],
                 tolerance = 1e-8,
                 info = paste("Cmax, subject", s))
    expect_equal(sub$Time[i], WNL_THEOPH_LOG$TMAX[WNL_THEOPH_LOG$Subject == s],
                 tolerance = 1e-8,
                 info = paste("Tmax, subject", s))
  }
})

# ---- Engine runs (skipped when packages absent) ------------------------------

test_that("PKNCA reproduces the WinNonlin reference within tolerance", {
  skip_if_not_installed("PKNCA")
  withr_wd({
    source(file.path("R", "02_engine_pknca.R"))
    d   <- prepare_theoph()
    res <- run_pknca(d$conc, d$dose)
    w   <- widen_engine(res)

    for (s in WNL_THEOPH_LOG$Subject) {
      row <- w[w$Subject == s, ]
      expect_equal(row$CMAX,
                   WNL_THEOPH_LOG$CMAX[WNL_THEOPH_LOG$Subject == s],
                   tolerance = NCA_RULES$tolerance,
                   info = paste("PKNCA Cmax, subject", s))
      expect_equal(row$LAMZHL,
                   WNL_THEOPH_LOG$LAMZHL[WNL_THEOPH_LOG$Subject == s],
                   tolerance = NCA_RULES$tolerance,
                   info = paste("PKNCA half-life, subject", s))
    }
  })
})

test_that("NonCompart reproduces the WinNonlin reference within tolerance", {
  skip_if_not_installed("NonCompart")
  withr_wd({
    d   <- prepare_theoph()
    res <- run_noncompart(d$conc, dose_mg = THEOPH_DOSE_MG)
    w   <- stats::reshape(res[, c("Subject", "PPTESTCD", "value")],
                          idvar = "Subject", timevar = "PPTESTCD",
                          direction = "wide")
    names(w) <- sub("^value\\.", "", names(w))

    for (s in WNL_THEOPH_LOG$Subject) {
      row <- w[w$Subject == s, ]
      expect_equal(row$AUCLST,
                   WNL_THEOPH_LOG$AUCLST[WNL_THEOPH_LOG$Subject == s],
                   tolerance = NCA_RULES$tolerance,
                   info = paste("NonCompart AUClast, subject", s))
    }
  })
})

# ---- Comparison machinery ---------------------------------------------------

test_that("compare_engines flags a deliberately injected discrepancy", {
  withr_wd(source(file.path("R", "05_compare.R")))

  a <- data.frame(Subject = 1:3, PPTESTCD = "CMAX",
                  value = c(10, 20, 30), engine = "A",
                  stringsAsFactors = FALSE)
  b <- a; b$engine <- "B"
  b$value[2] <- 20 * 1.05          # 5% off, well outside tolerance

  cmp <- compare_engines(a, b, reference = "A")
  expect_equal(sum(!cmp$pass_B), 1L)
  expect_equal(cmp$pct_diff_B[cmp$Subject == 2], 5, tolerance = 1e-8)
})

test_that("compare_engines survives a zero reference value", {
  withr_wd(source(file.path("R", "05_compare.R")))
  a <- data.frame(Subject = 1:2, PPTESTCD = "TLAG", value = c(0, 0),
                  engine = "A", stringsAsFactors = FALSE)
  b <- a; b$engine <- "B"
  expect_silent(cmp <- compare_engines(a, b, reference = "A"))
  expect_true(all(cmp$pass_B))
})

test_that("compare_engines refuses a single engine", {
  withr_wd(source(file.path("R", "05_compare.R")))
  a <- data.frame(Subject = 1, PPTESTCD = "CMAX", value = 1,
                  engine = "A", stringsAsFactors = FALSE)
  expect_error(compare_engines(a), "at least two")
})

# ---- Regression tests -------------------------------------------------------
# Both of these encode bugs found during the first real run of this suite.
# They looked like passes until the numbers were inspected.

test_that("REGRESSION: BLQ rule does not delete genuine time-zero measurements", {
  skip_if_not_installed("PKNCA")
  withr_wd({
    source(file.path("R", "02_engine_pknca.R"))
    d   <- prepare_theoph()
    res <- run_pknca(d$conc, d$dose)
    w   <- widen_engine(res)

    # PKNCA classifies conc == 0 as BLQ. Nine of twelve Theoph subjects have a
    # genuine pre-dose zero. With blq_first = "drop" those records vanish, the
    # interval start (0) then precedes the first measurement (0.25 h), and every
    # AUC-derived parameter comes back NA for those subjects.
    expect_equal(sum(is.na(w$AUCLST)), 0L,
                 info = "AUClast missing -- check NCA_RULES$blq_first")
    expect_equal(sum(is.na(w$AUCIFO)), 0L)
    expect_equal(sum(is.na(w$CLFO)),   0L)
    expect_equal(nrow(w), 12L)
  })
})

test_that("REGRESSION: PKNCA options are passed as a plain override list", {
  skip_if_not_installed("PKNCA")
  withr_wd({
    source(file.path("R", "02_engine_pknca.R"))
    d <- prepare_theoph()

    # PKNCA.options(...) mutates global state as a side effect. If run_pknca()
    # ever goes back to building its option list that way, the global default
    # will drift between engine runs and the comparison stops being reproducible.
    before <- PKNCA::PKNCA.options("auc.method")
    invisible(run_pknca(d$conc, d$dose))
    after  <- PKNCA::PKNCA.options("auc.method")
    expect_identical(before, after,
                     info = "run_pknca() leaked global PKNCA options")
  })
})

test_that("REGRESSION: a missing value is never counted as agreement", {
  withr_wd(source(file.path("R", "05_compare.R")))

  a <- data.frame(Subject = 1:3, PPTESTCD = "AUCLST",
                  value = c(100, 200, 300), engine = "A",
                  stringsAsFactors = FALSE)
  b <- a; b$engine <- "B"
  b$value[2] <- NA_real_          # engine B failed to produce a value

  cmp <- compare_engines(a, b, reference = "A")

  expect_false(cmp$pass_B[cmp$Subject == 2])
  expect_true(cmp$missing_B[cmp$Subject == 2])
  expect_equal(sum(cmp$missing_B), 1L)

  agree <- summarise_agreement(cmp)
  expect_equal(agree$n_missing, 1L)
  expect_equal(agree$n_compared, 2L)
  expect_equal(agree$n_subjects, 3L)

  # Wrap in capture.output: report_verdict prints by design, and an alarming
  # "NOT CLEAN" banner in the middle of a passing test run is confusing.
  invisible(capture.output(verdict <- report_verdict(agree)))
  expect_false(verdict$clean)
})

test_that("REGRESSION: full agreement with no gaps reports CLEAN", {
  withr_wd(source(file.path("R", "05_compare.R")))
  a <- data.frame(Subject = 1:3, PPTESTCD = "AUCLST",
                  value = c(100, 200, 300), engine = "A",
                  stringsAsFactors = FALSE)
  b <- a; b$engine <- "B"
  cmp <- compare_engines(a, b, reference = "A")
  invisible(capture.output(verdict <- report_verdict(summarise_agreement(cmp))))
  expect_true(verdict$clean)
})

# ---- Precision attribution --------------------------------------------------
# The headline "worst difference" against the reference measures the reference
# file's stored precision, not engine error. These tests pin that distinction
# so the README claim cannot silently become wrong.

test_that("engines agree with each other far more closely than with the reference", {
  skip_if_not_installed("PKNCA")
  skip_if_not_installed("NonCompart")
  withr_wd({
    source(file.path("R", "02_engine_pknca.R"))
    source(file.path("R", "05_compare.R"))
    d <- prepare_theoph()
    p <- run_pknca(d$conc, d$dose)
    n <- run_noncompart(d$conc, dose_mg = THEOPH_DOSE_MG)

    ee <- engine_vs_engine(p, n)

    # Engine-to-engine disagreement is floating-point noise, not method
    # difference: orders of magnitude below the reference's resolution.
    expect_lt(ee$max_pct_diff, 1e-10)
    expect_gt(ee$n_bit_identical, 50)
    expect_equal(ee$n, 192L)
  })
})

test_that("the reference's rounding bound exceeds every observed difference", {
  skip_if_not_installed("PKNCA")
  withr_wd({
    source(file.path("R", "05_compare.R"))
    ref <- file.path("phoenix", "reference", "winnonlin_theoph_log.csv")
    skip_if_not(file.exists(ref), "reference CSV not fetched")

    cmp  <- read.csv(file.path("outputs", "comparison_per_subject.csv"))
    skip_if_not("pct_diff_PKNCA" %in% names(cmp), "run run_all.R first")
    attr(cmp, "reference") <- "WinNonlin"

    prec <- attribute_precision(cmp, ref, engine = "PKNCA")
    expect_false(is.null(prec))

    # Every parameter must be precision-limited. If one is not, the engines
    # genuinely disagree with WinNonlin beyond what rounding can explain, and
    # that is a finding rather than a rounding artifact.
    expect_true(all(prec$precision_limited),
                info = paste("not precision-limited:",
                             paste(prec$PPTESTCD[!prec$precision_limited],
                                   collapse = ", ")))
  })
})

test_that("lambda-z is stored with fewer significant figures than other parameters", {
  withr_wd({
    ref <- file.path("phoenix", "reference", "winnonlin_theoph_log.csv")
    skip_if_not(file.exists(ref), "reference CSV not fetched")

    raw <- readLines(ref)
    hdr <- strsplit(raw[1], ",")[[1]]

    sigfig <- function(x) {
      d <- gsub("[.]", "", trimws(x))
      d <- sub("^0+", "", d)
      nchar(sub("0+$", "", d))
    }
    min_sf <- function(col) {
      j <- which(hdr == col)
      v <- vapply(raw[-1], function(l) strsplit(l, ",")[[1]][j], character(1),
                  USE.NAMES = FALSE)
      min(sigfig(v))
    }

    # The reference is written to a fixed ~11-character field. lambda-z is the
    # only parameter below 0.1, so its leading "0.0" costs two positions and it
    # ends up with markedly fewer significant figures than everything else.
    # That, not estimation accuracy, is why lambda-z shows the largest relative
    # deviation from the reference.
    #
    # NB: compare significant figures, not decimal places. Decimal counts are
    # misleading here because trailing zeros are stripped -- AUMCINF_obs prints
    # as "1313.951" while actually carrying 10 significant figures. An earlier
    # version of this test used decimal places and failed for exactly that
    # reason.
    expect_lt(min_sf("Lambda_z"), min_sf("AUClast"))
    expect_lt(min_sf("Lambda_z"), min_sf("HL_Lambda_z"))
    expect_lt(min_sf("Lambda_z"), min_sf("Vz_F_obs"))
    expect_lte(min_sf("Lambda_z"), 8L)
  })
})
