#!/usr/bin/env Rscript
#' Generate the figures embedded in README.md
#'
#'   Rscript R/07_figures.R
#'
#' Writes to figures/. These are committed to the repository (unlike outputs/)
#' so the README renders on GitHub without anyone having to run the pipeline.
#'
#' Sources are the recorded platform runs in reproducibility/, so the figures
#' are regenerable from checked-in data rather than from whatever happens to be
#' in outputs/ at the time.

suppressPackageStartupMessages({
  library(ggplot2)
})

source(file.path("R", "00_config.R"))

FIG <- "figures"
dir.create(FIG, showWarnings = FALSE)

REP_A <- file.path("reproducibility", "linux_x86_64")
REP_B <- file.path("reproducibility", "macos_arm64")

# Conservative theme: works identically on ggplot2 3.4.x and 4.0.x
th <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey93", colour = "grey70"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey30", size = 9),
    plot.caption     = element_text(colour = "grey45", size = 8, hjust = 0)
  )

ok <- function(f) if (file.exists(f)) TRUE else {
  message("missing: ", f); FALSE
}

# =============================================================== figure 1 ====
# Identity scatter, engines vs the published WinNonlin reference.
fig_agreement <- function() {
  f <- file.path(REP_A, "comparison_per_subject.csv")
  if (!ok(f)) return(invisible(NULL))
  cmp <- read.csv(f)

  d <- do.call(rbind, lapply(c("PKNCA", "NonCompart"), function(e) {
    data.frame(engine = e, PPTESTCD = cmp$PPTESTCD,
               ref = cmp$WinNonlin, val = cmp[[e]],
               stringsAsFactors = FALSE)
  }))
  d <- d[is.finite(d$ref) & is.finite(d$val), ]
  d$label <- PARAM_MAP$label[match(d$PPTESTCD, PARAM_MAP$PPTESTCD)]

  p <- ggplot(d, aes(ref, val, shape = engine, colour = engine)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey55", linewidth = 0.4) +
    geom_point(size = 1.7, alpha = 0.8) +
    facet_wrap(~ label, scales = "free", ncol = 4) +
    scale_colour_manual(values = c(PKNCA = "#1f77b4", NonCompart = "#d62728")) +
    scale_shape_manual(values = c(PKNCA = 16, NonCompart = 4)) +
    labs(title = "Both R engines vs published Phoenix WinNonlin output",
         subtitle = "Theophylline, 12 subjects, 16 parameters. Dashed line is identity.",
         x = "Phoenix WinNonlin 6.3/7.0", y = "R engine",
         colour = NULL, shape = NULL,
         caption = "192/192 comparisons within 0.1% tolerance for each engine.") +
    th + theme(legend.position = "top")

  ggsave(file.path(FIG, "fig1_agreement.png"), p,
         width = 10, height = 7.5, dpi = 130)
  message("fig1_agreement.png")
}

# =============================================================== figure 2 ====
# Observed difference against the reference's rounding bound, per parameter.
fig_precision <- function() {
  f <- file.path(REP_A, "precision_attribution.csv")
  if (!ok(f)) return(invisible(NULL))
  pa <- read.csv(f, stringsAsFactors = FALSE)

  pa <- pa[pa$max_observed_pct > 0, ]
  pa <- pa[order(pa$max_observed_pct), ]
  pa$label <- factor(pa$label, levels = pa$label)

  d <- rbind(
    data.frame(label = pa$label, value = pa$max_observed_pct,
               what = "Observed difference"),
    data.frame(label = pa$label, value = pa$max_rounding_pct,
               what = "Reference rounding bound")
  )

  p <- ggplot() +
    geom_segment(data = pa,
                 aes(y = label, yend = label,
                     x = max_observed_pct, xend = max_rounding_pct),
                 colour = "grey75", linewidth = 0.6) +
    geom_point(data = d, aes(x = value, y = label,
                             colour = what, shape = what), size = 2.6) +
    scale_x_log10(
      breaks = 10^seq(-9, -3),
      labels = c("1e-9", "1e-8", "1e-7", "1e-6", "1e-5", "1e-4", "1e-3")) +
    scale_colour_manual(values = c("Observed difference" = "#1f77b4",
                                   "Reference rounding bound" = "grey35")) +
    scale_shape_manual(values = c("Observed difference" = 16,
                                  "Reference rounding bound" = 124)) +
    labs(
      title = "Every parameter is precision-limited, not accuracy-limited",
      subtitle = paste("Observed difference sits below what the reference file's",
                       "stored decimals can resolve"),
      x = "Relative difference (%), log scale", y = NULL,
      colour = NULL, shape = NULL,
      caption = paste0(
        "Lambda_z is the only parameter below 0.1, so in the reference's ",
        "fixed 11-character field its leading '0.0'\ncosts two positions: 7-8 ",
        "significant figures instead of 10. That is why it shows the largest ",
        "deviation.")) +
    th + theme(legend.position = "top")

  ggsave(file.path(FIG, "fig2_precision.png"), p,
         width = 8.5, height = 5, dpi = 130)
  message("fig2_precision.png")
}

# =============================================================== figure 3 ====
# The three scales, on one log axis.
fig_scales <- function() {
  fa <- file.path(REP_A, "comparison_per_subject.csv")
  if (!ok(fa)) return(invisible(NULL))
  cmp <- read.csv(fa)
  obs_ref <- max(abs(cmp$pct_diff_PKNCA), na.rm = TRUE)

  d <- data.frame(
    what = factor(c(
      "Double-precision floor",
      "PKNCA vs NonCompart\n(engine disagreement)",
      "Engines vs reference CSV\n(file's stored precision)",
      "Acceptance tolerance"),
      levels = c("Acceptance tolerance",
                 "Engines vs reference CSV\n(file's stored precision)",
                 "PKNCA vs NonCompart\n(engine disagreement)",
                 "Double-precision floor")),
    value = c(100 * .Machine$double.eps, 4.3463e-13, obs_ref, 0.1),
    kind  = c("floor", "measured", "measured", "threshold"),
    stringsAsFactors = FALSE
  )
  d$lab <- sprintf("%.1e %%", d$value)

  p <- ggplot(d, aes(value, what, colour = kind)) +
    geom_segment(aes(x = 100 * .Machine$double.eps, xend = value,
                     y = what, yend = what),
                 colour = "grey80", linewidth = 0.5) +
    geom_point(size = 4) +
    geom_text(aes(label = lab), vjust = -1.3, size = 3, show.legend = FALSE) +
    scale_x_log10(breaks = 10^seq(-14, -1, by = 2)) +
    scale_colour_manual(values = c(floor = "grey55", measured = "#1f77b4",
                                   threshold = "#d62728"), guide = "none") +
    labs(
      title = "Two different quantities, six orders of magnitude apart",
      subtitle = paste("The headline agreement number measures the reference",
                       "file, not the engines"),
      x = "Relative difference (%), log scale", y = NULL,
      caption = paste0(
        "The engines agree with each other about a million times more closely ",
        "than the published reference\ncan resolve. 86 of 192 engine-to-engine ",
        "comparisons are bit-identical.")) +
    th + theme(axis.text.y = element_text(size = 8.5))

  ggsave(file.path(FIG, "fig3_scales.png"), p,
         width = 8.5, height = 4, dpi = 130)
  message("fig3_scales.png")
}

# =============================================================== figure 4 ====
# Cross-platform difference, in units of double-precision ULP.
fig_platform <- function() {
  fa <- file.path(REP_A, "comparison_per_subject.csv")
  fb <- file.path(REP_B, "comparison_per_subject.csv")
  if (!ok(fa) || !ok(fb)) return(invisible(NULL))

  a <- read.csv(fa); b <- read.csv(fb)
  a <- a[order(a$PPTESTCD, a$Subject), ]
  b <- b[order(b$PPTESTCD, b$Subject), ]
  stopifnot(identical(paste(a$Subject, a$PPTESTCD),
                      paste(b$Subject, b$PPTESTCD)))

  rel <- function(x, y) ifelse(abs(x) < .Machine$double.eps, 0,
                               abs(x - y) / abs(x))
  d <- rbind(
    data.frame(PPTESTCD = a$PPTESTCD, engine = "PKNCA",
               ulp = rel(a$PKNCA, b$PKNCA) / .Machine$double.eps),
    data.frame(PPTESTCD = a$PPTESTCD, engine = "NonCompart",
               ulp = rel(a$NonCompart, b$NonCompart) / .Machine$double.eps)
  )

  # Parameters whose computation involves log()/exp()
  trans <- c("LAMZ", "LAMZHL", "R2ADJ", "CLFO", "VZFO",
             "AUCLST", "AUCIFO", "AUCPEO", "AUMCLST", "AUMCIFO", "MRTEVIFO")
  agg <- aggregate(ulp ~ PPTESTCD + engine, d, max)
  agg$label <- PARAM_MAP$label[match(agg$PPTESTCD, PARAM_MAP$PPTESTCD)]
  agg$family <- ifelse(agg$PPTESTCD %in% trans,
                       "Involves log() / exp()",
                       "Selection or counting only")

  ord <- aggregate(ulp ~ label, agg, max)
  agg$label <- factor(agg$label, levels = ord$label[order(ord$ulp)])

  p <- ggplot(agg, aes(ulp, label, colour = family, shape = engine)) +
    geom_point(size = 2.6, alpha = 0.85,
               position = position_nudge(y = ifelse(agg$engine == "PKNCA",
                                                    0.14, -0.14))) +
    scale_colour_manual(values = c("Involves log() / exp()" = "#d62728",
                                   "Selection or counting only" = "#2ca02c")) +
    scale_shape_manual(values = c(PKNCA = 16, NonCompart = 4)) +
    labs(
      title = "Cross-platform difference: x86_64 Linux vs ARM64 macOS",
      subtitle = "Maximum difference per parameter, in units of double-precision ULP",
      x = "Difference (ULP)", y = NULL, colour = NULL, shape = NULL,
      caption = paste0(
        "Worst case 15.2 ULP. Every parameter that moved depends on log(); ",
        "Cmax, Tmax, Tlast, Clast and\nN points lambda_z are bit-identical ",
        "across architectures. Signature of libm differences, not method.")) +
    th + theme(legend.position = "top", legend.box = "vertical",
               legend.spacing.y = unit(0, "pt"))

  ggsave(file.path(FIG, "fig4_platform.png"), p,
         width = 8.5, height = 5.5, dpi = 130)
  message("fig4_platform.png")
}

fig_agreement()
fig_precision()
fig_scales()
fig_platform()
message("\nFigures written to ", FIG, "/")
