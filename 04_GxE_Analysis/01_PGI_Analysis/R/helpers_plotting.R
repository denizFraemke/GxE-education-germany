# helpers_plotting.R
# ggplot scales, theme, and plotting functions.

# ---- Plotting helpers ------------------------------------------------

# Colour + linetype scales (used in every plot)
.scale_col  <- scale_colour_manual(values = c(East = COL_EAST, West = COL_WEST))
.scale_fill <- scale_fill_manual(  values = c(East = COL_EAST, West = COL_WEST))
.scale_lt   <- scale_linetype_manual(values = c(female = LT_FEMALE, male = LT_MALE))

# Base theme
.theme <- theme_minimal(base_size = 12)

# Standard subtitle for any plot that draws non-parametric kernel-regression
# slopes as a dotted overlay alongside a GAM/LM fit.
.subtitle_with_dotted <- function(extra = NULL) {
  base <- "Solid: model fit (GAM/LM) with 95% CI ribbon. Dotted: non-parametric kernel-regression PGI slope (bandwidth = 8 y)."
  if (is.null(extra) || !nzchar(extra)) base else paste(extra, "|", base)
}

# Compute slope-based y-limits with a fractional pad. The lower bound is
# always at most -0.05 so the y = 0 reference line has visible breathing
# room and the dashed x-axis is comfortably above the plot floor.
# CI ribbons extending beyond are clipped visually by coord_cartesian.
.slope_ylim <- function(slope_vals, pad_frac = 0.30, neg_min = -0.05) {
  rng <- range(slope_vals, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(c(NA_real_, NA_real_))
  pad <- pad_frac * diff(rng)
  c(min(rng[1] - pad, neg_min),
    max(rng[2] + pad, 0))
}

# Pre-compute scale info for the histogram overlay + the matching count
# secondary axis. Histograms anchor at y = 0 (the x-axis line) and rise
# upward to fill hist_frac of the y > 0 visible range; the sec axis is
# the inverse of that scaling so right-edge sample-size labels run from
# 0 (at y = 0) up to max_count (at y = hist_frac * y_lim[2]).
.hist_scale_info <- function(hist_dat, y_lim, bins = 50, hist_frac = 1.0) {
  if (is.null(hist_dat) || !nrow(hist_dat) ||
      !all(is.finite(y_lim)) || diff(y_lim) == 0) return(NULL)
  if (y_lim[2] <= 0) return(NULL)   # need a positive vertical region
  groups <- if ("east_west" %in% names(hist_dat))
    unique(na.omit(as.character(hist_dat$east_west))) else NA_character_
  max_count <- 0
  for (g in groups) {
    d_sub <- if (!is.na(g))
      hist_dat[!is.na(hist_dat$east_west) &
                 as.character(hist_dat$east_west) == g, , drop = FALSE]
      else hist_dat
    if (!nrow(d_sub)) next
    h <- tryCatch(hist(d_sub$birth_year, breaks = bins, plot = FALSE),
                  error = function(e) NULL)
    if (is.null(h)) next
    max_count <- max(max_count, max(h$counts, na.rm = TRUE))
  }
  if (max_count == 0) return(NULL)
  list(
    scale_factor = hist_frac * y_lim[2] / max_count,
    y_baseline   = 0,                     # x-axis line
    max_count    = max_count,
    bins         = bins
  )
}

# Histogram of birth_year by region, drawn as pre-binned rectangles with
# explicit ymin = y_baseline so empty bins do not produce a phantom solid
# block.
.hist_layers_east_west <- function(hist_dat, scale_info, bins = 50,
                                   alpha = 0.10) {
  if (is.null(hist_dat) || is.null(scale_info)) return(list())
  cols <- c(West = COL_WEST, East = COL_EAST)
  Filter(Negate(is.null), lapply(c("West", "East"), function(reg) {
    d_sub <- hist_dat[!is.na(hist_dat$east_west) &
                       as.character(hist_dat$east_west) == reg, , drop = FALSE]
    if (!nrow(d_sub)) return(NULL)
    h <- tryCatch(hist(d_sub$birth_year, breaks = bins, plot = FALSE),
                  error = function(e) NULL)
    if (is.null(h)) return(NULL)
    binned <- data.frame(
      x_lo  = h$breaks[-length(h$breaks)],
      x_hi  = h$breaks[-1],
      count = h$counts,
      stringsAsFactors = FALSE
    )
    binned <- binned[binned$count > 0, , drop = FALSE]
    if (!nrow(binned)) return(NULL)
    binned$y_lo <- scale_info$y_baseline
    binned$y_hi <- scale_info$y_baseline + binned$count * scale_info$scale_factor
    geom_rect(
      data = binned,
      aes(xmin = x_lo, xmax = x_hi, ymin = y_lo, ymax = y_hi),
      fill = cols[reg], colour = NA, alpha = alpha,
      inherit.aes = FALSE
    )
  }))
}

# Sample-size companion (right) y-axis when a histogram is overlaid.
.scale_y_with_count <- function(scale_info) {
  if (is.null(scale_info))
    return(scale_y_continuous(expand = expansion(mult = c(0, 0.05))))
  brk <- pretty(c(0, scale_info$max_count), n = 4)
  scale_y_continuous(
    expand   = expansion(mult = c(0, 0.05)),
    sec.axis = sec_axis(
      trans  = ~ (. - scale_info$y_baseline) / scale_info$scale_factor,
      name   = NULL,
      breaks = brk
    )
  )
}

# Trim a curves data.frame to the per-cell birth_year range observed in
# hist_dat. Drops the long thin tails where there are very few observations
# (and where CIs balloon as a consequence). Default 2nd-98th percentile
# keeps the curve over 96 percent of the data and clips only the extreme
# boundary noise.
.trim_to_data_range <- function(curves, hist_dat,
                                group_vars = "east_west",
                                qlo = 0.02, qhi = 0.98) {
  if (is.null(hist_dat) || !nrow(hist_dat)) return(curves)
  if (!"birth_year" %in% names(curves) || !"birth_year" %in% names(hist_dat))
    return(curves)
  use_groups <- length(group_vars) > 0 &&
                all(group_vars %in% names(curves)) &&
                all(group_vars %in% names(hist_dat))
  if (!use_groups) {
    qs <- quantile(hist_dat$birth_year, c(qlo, qhi), na.rm = TRUE)
    return(curves[curves$birth_year >= qs[1] & curves$birth_year <= qs[2], ,
                  drop = FALSE])
  }
  cv_grp <- do.call(paste, c(curves[, group_vars, drop = FALSE], list(sep = "|")))
  hd_grp <- do.call(paste, c(hist_dat[, group_vars, drop = FALSE], list(sep = "|")))
  keep <- rep(FALSE, nrow(curves))
  for (g in unique(cv_grp)) {
    qs <- quantile(hist_dat$birth_year[hd_grp == g], c(qlo, qhi), na.rm = TRUE)
    if (any(is.na(qs))) next
    keep[cv_grp == g & curves$birth_year >= qs[1] & curves$birth_year <= qs[2]] <- TRUE
  }
  curves[keep, , drop = FALSE]
}

# Combined GAM/LM curve + non-parametric overlay (east_west only)
plot_combined_east_west <- function(gam_curves, np_slopes, title, ylab,
                                    hist_dat = NULL, hist_alpha = 0.10,
                                    hist_bins = 50) {
  gam_curves <- .trim_to_data_range(gam_curves, hist_dat, group_vars = "east_west")
  np_slopes  <- .trim_to_data_range(np_slopes,  hist_dat, group_vars = "east_west")
  y_lim <- .slope_ylim(c(gam_curves$slope, np_slopes$slope))
  scale_info <- .hist_scale_info(hist_dat, y_lim, bins = hist_bins)
  ggplot(gam_curves, aes(x = birth_year, y = slope, colour = east_west)) +
    .hist_layers_east_west(hist_dat, scale_info, bins = hist_bins, alpha = hist_alpha) +
    .scale_y_with_count(scale_info) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = east_west),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1) +
    geom_line(data = np_slopes,
              aes(x = birth_year, y = slope, colour = east_west),
              linewidth = 0.8, alpha = 0.75, linetype = "dotted",
              inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    .scale_col + .scale_fill +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Birth year", y = ylab, title = title,
         subtitle = .subtitle_with_dotted(),
         colour = "Region", fill = "Region") +
    .theme
}

# Combined GAM/LM curve only (east_west x gender) — np overlay dropped because
# the four (region x gender) cells are noisy and read like region averages
# when overlaid with the GAM ribbons.
plot_combined_east_west_gender <- function(gam_curves, np_slopes, title, ylab,
                                           hist_dat = NULL, hist_alpha = 0.10,
                                           hist_bins = 50) {
  gam_curves <- .trim_to_data_range(gam_curves, hist_dat,
                                    group_vars = c("east_west", "gender"))
  gam_curves$group <- paste(gam_curves$east_west, gam_curves$gender)
  y_lim <- .slope_ylim(gam_curves$slope)
  scale_info <- .hist_scale_info(hist_dat, y_lim, bins = hist_bins)

  ggplot(gam_curves, aes(x = birth_year, y = slope,
                          colour = east_west, linetype = gender, group = group)) +
    .hist_layers_east_west(hist_dat, scale_info, bins = hist_bins, alpha = hist_alpha) +
    .scale_y_with_count(scale_info) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi,
                    fill = east_west, group = group),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    .scale_col + .scale_fill + .scale_lt +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Birth year", y = ylab, title = title,
         subtitle = "Solid: model fit (GAM) with 95% CI ribbon, four Region x Gender cells.",
         colour = "Region", fill = "Region",
         linetype = "Gender") +
    .theme
}

# Gender-stratified slope plot, filtered to a single region
plot_gender_slopes_one_region <- function(gam_curves, np_slopes, title, ylab, region) {
  gc <- gam_curves[gam_curves$east_west == region, , drop = FALSE]
  ns <- np_slopes [np_slopes$east_west  == region, , drop = FALSE]
  if (nrow(gc) == 0) return(NULL)
  region_col <- if (region == "East") COL_EAST else COL_WEST

  ggplot(gc, aes(x = birth_year, y = slope,
                 linetype = gender, group = gender)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = 0.15, fill = region_col, colour = NA) +
    geom_line(linewidth = 1, colour = region_col) +
    geom_line(data = ns,
              aes(x = birth_year, y = slope, linetype = gender, group = gender),
              linewidth = 0.8, alpha = 0.75, colour = region_col,
              inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    .scale_lt +
    labs(x = "Birth year", y = ylab,
         title = sprintf("%s — %s", title, region),
         linetype = "Gender") +
    .theme
}

# Difference-smooth plot (with simultaneous CI). Two faceting modes:
#   gender_facet = TRUE  -> diffs is East - West, faceted by gender
#   region_facet = TRUE  -> diffs is female - male, faceted by east_west
plot_diff_smooth <- function(diffs, title,
                             gender_facet = FALSE,
                             region_facet = FALSE,
                             hist_dat = NULL, hist_alpha = 0.10,
                             hist_bins = 50) {
  trim_groups <- if (region_facet && "east_west" %in% names(diffs))
    "east_west" else if (gender_facet && "gender" %in% names(diffs))
    "gender" else character(0)
  diffs <- .trim_to_data_range(diffs, hist_dat, group_vars = trim_groups)
  y_lim <- .slope_ylim(diffs$diff)
  scale_info <- .hist_scale_info(hist_dat, y_lim, bins = hist_bins)
  ylab_expr <- if (region_facet)
    expression(Delta * beta ~ "(female - male)") else
    expression(Delta * beta ~ "(East - West)")
  p <- ggplot(diffs, aes(x = birth_year, y = diff)) +
    .hist_layers_east_west(hist_dat, scale_info, bins = hist_bins, alpha = hist_alpha) +
    .scale_y_with_count(scale_info) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.20, fill = "grey60") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 1, colour = "black") +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Birth year",
         y = ylab_expr,
         title = title,
         subtitle = "Solid: difference smooth with 95% simultaneous CI ribbon.") +
    .theme
  if (region_facet && "east_west" %in% names(diffs)) {
    p <- p + facet_wrap(~ east_west)
  } else if (gender_facet && "gender" %in% names(diffs)) {
    p <- p + facet_wrap(~ gender)
  }
  p
}

# Birth-year coverage histogram, by region
plot_cohort_coverage <- function(dat, title = "Birth-year coverage by region") {
  ggplot(dat, aes(x = birth_year, fill = east_west)) +
    geom_histogram(bins = 60, position = "identity", alpha = 0.5) +
    .scale_fill +
    labs(x = "Birth year", y = "Count", title = title, fill = "Region") +
    .theme
}

plot_attrition_flow <- function(attrition_df, title,
                                y_min = NULL, y_max = NULL) {
  plot_df <- attrition_df |>
    mutate(
      step_f = factor(step, levels = unique(step)),
      label = ifelse(n_dropped > 0, paste0("-", n_dropped), "")
    )

  fixed_y <- !is.null(y_min) || !is.null(y_max)
  p <- ggplot(plot_df, aes(x = step_order, y = n_after, group = 1)) +
    geom_line(linewidth = 0.8, colour = "grey35") +
    geom_point(size = 2.4, colour = "black") +
    geom_text(aes(label = label), nudge_y = 0.03 * max(plot_df$n_before, na.rm = TRUE),
              size = 3, colour = "grey20") +
    scale_x_continuous(
      breaks = plot_df$step_order,
      labels = plot_df$step
    ) +
    labs(x = "Exclusion step", y = "Remaining N", title = title) +
    facet_grid(pipeline ~ cohort, scales = if (fixed_y) "free_x" else "free_y") +
    .theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  if (fixed_y) p <- p + coord_cartesian(
    ylim = c(if (is.null(y_min)) NA_real_ else y_min,
             if (is.null(y_max)) NA_real_ else y_max)
  )
  p
}

plot_within_cohort_slopes <- function(curves, title, ylab,
                                      facet_by = "cohort",
                                      stratify_gender = FALSE,
                                      region_filter   = NULL) {
  if (is.null(curves) || nrow(curves) == 0) return(NULL)

  if (!is.null(region_filter)) {
    curves <- curves[!is.na(curves$east_west) & curves$east_west == region_filter, ,
                     drop = FALSE]
    if (nrow(curves) == 0) return(NULL)
  }

  has_east   <- "east_west" %in% names(curves) &&
                any(!is.na(curves$east_west)) && is.null(region_filter)
  has_gender <- stratify_gender && "gender" %in% names(curves) &&
                length(unique(na.omit(curves$gender))) > 1

  aes_base <- if (has_east && has_gender) {
    aes(x = birth_year, y = slope,
        colour = east_west, linetype = gender,
        group = interaction(east_west, gender))
  } else if (has_east) {
    aes(x = birth_year, y = slope, colour = east_west, group = east_west)
  } else if (has_gender) {
    aes(x = birth_year, y = slope, linetype = gender, group = gender)
  } else {
    aes(x = birth_year, y = slope, group = 1)
  }

  p <- ggplot(curves, aes_base) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40")
  if (has_east) {
    p <- p +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = east_west),
                  alpha = 0.15, colour = NA, show.legend = FALSE)
  } else {
    p <- p +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                  alpha = 0.20, fill = "grey60", colour = NA)
  }
  p <- p +
    geom_line(linewidth = 1) +
    labs(x = "Birth year", y = ylab, title = title,
         colour = "Region", fill = "Region", linetype = "Gender") +
    .theme

  if (has_east)  p <- p + .scale_col + .scale_fill
  if (has_gender) p <- p + .scale_lt

  p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free")
}
