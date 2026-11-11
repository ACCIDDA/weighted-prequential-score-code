# produce_all_figures.R
#
# Outputs
#   figure_1_intro.png        — Fig 1: Michigan intro forecast comparison
#   figure_2_et_weights.pdf   — Fig 2: ET weight distributions
#   figure_3_examples.png     — Fig 3: Alabama/Oklahoma illustrative examples
#   figure_4_rcrit_cdf.png    — Fig 4: Estimated CDFs of R_crit
#   figure_5_tradeoff.png     — Fig 5: Improvement/penalty tradeoff
#   table_1_rcrit_summary.tex — Table 1: R_crit summary (LaTeX booktabs)
#
# ==============================================================================

# ==============================================================================
# 0.  Setup
# ==============================================================================


library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)
library(purrr)
library(here)
library(MMWRweek)
library(ggrepel)


# get_script_dir <- function() {
#   args        <- commandArgs(trailingOnly = FALSE)
#   file_arg    <- "--file="
#   script_path <- sub(file_arg, "", args[grep(file_arg, args)])
#   if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
# }
# 
# script_dir   <- get_script_dir()
# project_root <- if (basename(script_dir) == "figure_scripts_fixed") dirname(script_dir) else script_dir
# scripts_dir  <- file.path(project_root, "figure_scripts_fixed")

# Main directory
# Works whether this directory is the repo root or sits inside a larger project
dir <- if (dir.exists(here("R"))) here() else here("figure_scripts_fixed")

add_season <- function(d) {
  mw <- MMWRweek(as.Date(d))
  # MMWR week 40 cutoff for season year
  start_year <- ifelse(mw$MMWRweek >= 40, mw$MMWRyear, mw$MMWRyear - 1)
  paste0(start_year, "/", start_year + 1)
}

source(file.path(dir, "R", "load_modern_data.R"))
source(file.path(dir, "R", "rank_analysis.R"))

# Create output folder to save figures
output_dir <- file.path(dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1.  Constants
# ==============================================================================

# Seasons and horizons used in Figs 3, 4, 5
SEASONS           <- c("2023/2024", "2024/2025", "2025/2026")
HORIZONS          <- c(2L, 3L)

# Event week definition
PEAK_WINDOW       <- 3L           # peak week event = peak ± 3 weeks (7 total)
GROWTH_QUANTILE   <- 0.50         # growth week event: threshold at median of positive increase

# Eligibility filters (identical across all figures)
COVERAGE          <- 0.9          # min fraction of season weeks submitted
THRESHOLD         <- 1.0          # min fraction of event weeks submitted
MIN_MODELS        <- 5L           # min eligible models per location
MIN_EVENTS_PEAK   <- 7L           # min event weeks for peak cost function
MIN_EVENTS_GROWTH <- 4L           # min event weeks for growth (after scored-date restriction)

# Fig 5 weight values
W_VALUES_FIG5 <- c(0.3, 0.5, 0.9)

# Season colors (Figs 4, 5)
SEASON_COLOURS <- c("2023/2024" = "#2166AC",
                    "2024/2025" = "#D6604D",
                    "2025/2026" = "#4DAF4A")

# ==============================================================================
# 2.  Shared helpers
# ==============================================================================

# ---- 2a. Nature theme (base_size = 14; used directly by Fig 3)
# Fig 5 adds explicit +2pt overrides inside make_tradeoff_panel to stay at 16.
figure_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.background   = element_rect(fill = "white", colour = NA),
      plot.background    = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.line.x        = element_line(colour = "black", linewidth = 0.4),
      axis.line.y        = element_line(colour = "black", linewidth = 0.4),
      axis.text          = element_text(size = 14, colour = "grey10"),
      axis.title         = element_text(size = 15, colour = "black", face = "plain"),
      axis.title.x       = element_blank(),
      axis.ticks         = element_line(colour = "black"),
      axis.ticks.length  = grid::unit(2.5, "pt"),
      legend.position    = "none",
      plot.title         = element_text(size = 14, face = "plain", colour = "black",
                                        margin = margin(b = 2))
    )
}

# ---- 2b. Forecast series extraction (Figs 1 and 3) -------------------------
# For a given location, season, horizon, and model, extracts three data frames
# used to draw one forecast panel:
#   observed — reported hospitalizations each week (the black line)
#   forecast — the model's median (q=0.5) prediction each week (the blue line)
#   interval — the model's 5th–95th percentile predictions each week (the blue ribbon)

get_series <- function(forecasts, truth, location, season, horizon, model) {
  # Extract the forecasts for given model, location, season, horizon
  fc <- forecasts |>
    dplyr::filter(
      location == .env$location, season == .env$season,
      horizon  == .env$horizon,  model  == .env$model
    )
  # Splits into observed, forecast, and interval dataframes
  observed <- truth |>
    dplyr::filter(location == .env$location, tgt_season == .env$season) |>
    dplyr::select(target_end_date, observed)

  forecast <- fc |>
    dplyr::filter(quantile_level == 0.5) |>
    dplyr::distinct(target_end_date, predicted)

  interval <- fc |>
    dplyr::filter(quantile_level %in% c(0.05, 0.95)) |>
    dplyr::distinct(target_end_date, quantile_level, predicted) |>
    tidyr::pivot_wider(names_from = quantile_level, values_from = predicted,
                       names_prefix = "q") |>
    dplyr::rename(lower_90 = q0.05, upper_90 = q0.95)
  # Returns list of the splits dataframes
  list(observed = observed, forecast = forecast, interval = interval)
}

# ---- 2c. Peak cost function event dates -------------------------------------
# Returns the 7 event dates centered on the seasonal peak for a given location
# and season: the peak week (highest observed hospitalizations) plus the 3 weeks
# before and after it. Returns NULL if the location has no truth data or if any
# of the 7 dates fall outside the observed record (e.g. early-season locations).

get_peak_event_dates <- function(truth, location, season, peak_window) {
  loc_tbl <- truth |>
    dplyr::filter(location == .env$location, tgt_season == .env$season) |>
    dplyr::arrange(target_end_date)
  # Guard against incorrect locations, seasons, that don't exist in the dataset
  if (nrow(loc_tbl) == 0L) return(NULL)
  
  dates    <- loc_tbl$target_end_date
  peak     <- dates[order(-loc_tbl$observed, loc_tbl$target_end_date)[1L]]
  # Dates are in days so gives peak +- *week* window by multiplying by 7
  window   <- peak + seq(-peak_window, peak_window) * 7L
  
  # Skip location is the window does not lie in the season dates
  if (!all(window %in% dates)) return(NULL)
  
  window
}

# truth_edge <- data.frame(
#   location        = "01",
#   tgt_season      = "2023/2024",
#   target_end_date = as.Date("2023-11-04") + seq(0, 8) * 7,
#   observed        = c(1, 5, 10, 100, 200, 30, 20, 10, 5)
# )
# get_peak_event_dates(truth_edge, "01", "2023/2024", 3)


# ---- 2d. Growth cost function event dates ------------------------------------
# Returns the high-growth event weeks for a given location and season.
# Computes week-over-week increases in hospitalizations, then flags weeks where
# the increase exceeds the median positive increase (or growth_quantile percentile).
# Returns NULL if there are fewer than 3 positive-growth weeks (not enough to
# define a meaningful threshold) or if no weeks exceed the threshold.
get_growth_event_dates <- function(truth, location, season,
                                   growth_quantile = GROWTH_QUANTILE) {
  loc_tbl <- truth |>
    dplyr::filter(location == .env$location, tgt_season == .env$season) |>
    dplyr::arrange(target_end_date)
  if (nrow(loc_tbl) < 2L) return(NULL)

  dates      <- loc_tbl$target_end_date
  deltas     <- diff(loc_tbl$observed)
  pos_deltas <- deltas[deltas > 0]
  if (length(pos_deltas) < 3L) return(NULL)

  thresh      <- quantile(pos_deltas, probs = growth_quantile, names = FALSE)
  event_dates <- dates[-1L][deltas > thresh]
  if (length(event_dates) == 0L) return(NULL)
  event_dates
}

# test function

# truth_growth <- data.frame(
#    location        = "01",
#    tgt_season      = "2023/2024",
#    target_end_date = as.Date("2023-11-04") + seq(0, 9) * 7,
#    observed        = c(10, 15, 25, 60, 40, 80, 150, 100, 50, 20)
# )
# 
# delta_test_raw <- diff(c(10, 15, 25, 60, 40, 80, 150, 100, 50, 20))
# delta_test <- delta_test_raw[delta_test_raw > 0]
# delta_test <- delta_test[delta_test > median(delta_test)]
# get_growth_event_dates(truth_growth, "01", "2023/2024")[delta_test %in% delta_test_raw]
#  "2023-12-09" "2023-12-16"
# 
# get_growth_event_dates(truth_growth, "01", "2023/2024")


# ---- 2e. Joint event weeks across both cost functions -------------------------
# Table 1 and Figs 4 and 5 require a model to have submitted for every event week
# under BOTH event definitions. A location is therefore analysed only when both
# definitions yield a valid event set, and the eligible pool is the intersection
# of the two. Returns NULL if either definition fails, in which case the location
# is dropped from every analysis.
get_both_event_dates <- function(truth, location, season, scored_dates) {
  peak <- get_peak_event_dates(truth, location, season, PEAK_WINDOW)
  if (is.null(peak) || length(peak) < MIN_EVENTS_PEAK) return(NULL)

  growth <- get_growth_event_dates(truth, location, season)
  if (is.null(growth)) return(NULL)
  growth <- growth[growth %in% scored_dates]
  if (length(growth) < MIN_EVENTS_GROWTH) return(NULL)

  list(peak = peak, growth = growth)
}

# Models meeting the submission criteria under BOTH event definitions.
get_joint_eligible_models <- function(engine, events, coverage, threshold) {
  scores_obs  <- engine$scores_obs
  n_all_dates <- dplyr::n_distinct(scores_obs$target_end_date)

  eligible_for <- function(event_dates) {
    scores_obs |>
      dplyr::mutate(E = dplyr::if_else(target_end_date %in% event_dates, 1L, 0L)) |>
      dplyr::group_by(model) |>
      dplyr::summarise(
        n_dates_present   = dplyr::n_distinct(target_end_date),
        n_flagged_present = dplyr::n_distinct(target_end_date[E == 1L]),
        .groups = "drop"
      ) |>
      dplyr::filter(
        n_dates_present   >= floor(coverage  * n_all_dates),
        n_flagged_present >= floor(threshold * length(event_dates))
      ) |>
      dplyr::pull(model)
  }

  intersect(eligible_for(events$peak), eligible_for(events$growth))
}

# ==============================================================================
# 3.  Data loading 
# ==============================================================================

# All quantile forecasts for the analysis seasons
modern_quantile <- load_modern_quantile(
  seasons  = SEASONS,
  horizons = 1:4
)

# Truth table one row per location and arget_end_date with observed value
truth_tbl <- modern_quantile |>
  dplyr::filter(!is.na(observed)) |>
  dplyr::distinct(location, target_end_date, observed) |>
  dplyr::mutate(tgt_season = add_season(target_end_date)) |>
  dplyr::filter(tgt_season %in% SEASONS) |>
  dplyr::arrange(location, target_end_date)

# Locations
locations <- sort(setdiff(unique(modern_quantile$location), "US"))

# ==============================================================================
# 4.  Fig 1 — Michigan intro: UMass-trends_ensemble vs UMass-flusion
# ==============================================================================

# Plotting code...

FIG1_LOC     <- "26"
FIG1_SEASON  <- "2023/2024"
FIG1_HORIZON <- 2L
FIG1_MODEL_A <- "UMass-trends_ensemble"
FIG1_MODEL_B <- "UMass-flusion"

fig1_series_a <- get_series(modern_quantile, truth_tbl, FIG1_LOC, FIG1_SEASON, FIG1_HORIZON, FIG1_MODEL_A)
fig1_series_b <- get_series(modern_quantile, truth_tbl, FIG1_LOC, FIG1_SEASON, FIG1_HORIZON, FIG1_MODEL_B)

fig1_x_lim <- range(fig1_series_a$observed$target_end_date, na.rm = TRUE)
fig1_y_lim <- range(c(
  fig1_series_a$observed$observed,
  fig1_series_a$forecast$predicted, fig1_series_b$forecast$predicted,
  fig1_series_a$interval$lower_90,  fig1_series_a$interval$upper_90,
  fig1_series_b$interval$lower_90,  fig1_series_b$interval$upper_90
), na.rm = TRUE)

make_intro_panel <- function(series, model_name, x_lim, y_lim) {
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = series$interval,
      ggplot2::aes(x = target_end_date, ymin = lower_90, ymax = upper_90),
      fill = "#2166AC", alpha = 0.15, inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = series$observed,
      ggplot2::aes(x = target_end_date, y = observed),
      colour = "black", linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = series$observed,
      ggplot2::aes(x = target_end_date, y = observed),
      colour = "black", size = 1.5
    ) +
    ggplot2::geom_line(
      data = series$forecast,
      ggplot2::aes(x = target_end_date, y = predicted),
      colour = "#2166AC", linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = series$forecast,
      ggplot2::aes(x = target_end_date, y = predicted),
      colour = "#2166AC", size = 1.5
    ) +
    ggplot2::coord_cartesian(ylim = y_lim) +
    ggplot2::labs(title = model_name, x = NULL, y = "Incident Hospitalizations") +
    ggplot2::theme_minimal(base_size = 22) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_blank(),
      axis.line.x        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.line.y        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.text          = ggplot2::element_text(size = 22, colour = "grey10"),
      axis.title.y       = ggplot2::element_text(size = 24, colour = "black", face = "plain"),
      axis.title.x       = ggplot2::element_blank(),
      axis.ticks         = ggplot2::element_line(colour = "black"),
      legend.position    = "none",
      plot.title         = ggplot2::element_text(size = 22, face = "plain", colour = "black")
    ) +
    ggplot2::scale_x_date(
      limits      = x_lim,
      date_breaks = "2 months",
      date_labels = "%b\n%Y"
    )
}

fig1_panel_a <- make_intro_panel(fig1_series_a, FIG1_MODEL_A, fig1_x_lim, fig1_y_lim)
fig1_panel_b <- make_intro_panel(fig1_series_b, FIG1_MODEL_B, fig1_x_lim, fig1_y_lim)

fig1 <- fig1_panel_b + fig1_panel_a +
  patchwork::plot_layout(axis_titles = "collect")

ggplot2::ggsave(
  file.path(output_dir, "figure_1_intro.png"),
  plot = fig1, width = 14, height = 5.5, dpi = 300
)

# ==============================================================================
# 5.  Fig 2 — ET weight distributions
# ==============================================================================

# 2022/2023 truth, read from the truth file directly. It cannot come from
# modern_quantile: that is forecasts joined to truth, and the hub holds no
# quantile forecasts before 2023, so the 2022/2023 season has no rows there.
us_obs_2223 <- load_hub_truth() |>
  dplyr::mutate(tgt_season = add_season(target_end_date)) |>
  dplyr::filter(location == "US", tgt_season == "2022/2023", !is.na(observed)) |>
  dplyr::distinct(target_end_date, observed) |>
  dplyr::arrange(target_end_date) |>
  dplyr::mutate(t = dplyr::row_number())

df_fig2 <- us_obs_2223
Tn       <- nrow(df_fig2)

# Plotting code...

peak_t   <- which.max(df_fig2$observed)[1L]
E_idx    <- seq(max(1L, peak_t - 3L), min(Tn, peak_t + 3L))
m        <- length(E_idx)
alpha_w  <- 0.50

bad_event_weights <- c(0.04, 0.08, 0.11, 0.14, 0.05, 0.05, 0.03)

df_fig2 <- df_fig2 |>
  dplyr::mutate(
    in_E   = t %in% E_idx,
    w_unif = 1 / Tn,
    w_et   = dplyr::if_else(in_E, alpha_w / m, (1 - alpha_w) / (Tn - m)),
    w_bad  = (1 - alpha_w) / (Tn - m)
  )
df_fig2$w_bad[E_idx] <- bad_event_weights

keep_idx      <- seq(max(1L, peak_t - 10L), min(Tn, peak_t + 10L))
df_plot_fig2  <- df_fig2 |> dplyr::slice(keep_idx)
df_plot_fig2  <- df_plot_fig2 |> dplyr::arrange(target_end_date) |> dplyr::mutate(plot_x = dplyr::row_number())

n_right_keep  <- 5L
right_after_E <- which(df_plot_fig2$t > max(E_idx))
df_plot_fig2  <- df_plot_fig2 |>
  dplyr::slice(c(which(df_plot_fig2$t <= max(E_idx)), head(right_after_E, n_right_keep))) |>
  dplyr::mutate(plot_x = dplyr::row_number())

y_top  <- 0.18
bg_lo  <- 0.065
bg_hi  <- 0.17
rng    <- range(df_fig2$observed, na.rm = TRUE)
df_plot_fig2 <- df_plot_fig2 |>
  dplyr::mutate(obs_bg = bg_lo + (bg_hi - bg_lo) * (observed - rng[1]) / (rng[2] - rng[1]))

e_x        <- df_plot_fig2 |> dplyr::filter(in_E) |> dplyr::pull(plot_x)
event_band <- tibble::tibble(xmin = min(e_x) - 0.5, xmax = max(e_x) + 0.5, ymin = 0, ymax = y_top)

make_weights_panel <- function(data, weight_col, title_txt, y_lab = "Weight") {
  bars_outside <- data |> dplyr::filter(!in_E) |> dplyr::transmute(plot_x, w = .data[[weight_col]])
  bars_E       <- data |> dplyr::filter(in_E)  |> dplyr::transmute(plot_x, w = .data[[weight_col]])

  ggplot2::ggplot(data, ggplot2::aes(x = plot_x)) +
    ggplot2::geom_rect(
      data = event_band,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scales::alpha("#F04E4E", 0.10), color = NA
    ) +
    ggplot2::geom_col(
      data = bars_outside |> dplyr::filter(w > 0),
      ggplot2::aes(x = plot_x, y = w), inherit.aes = FALSE,
      width = 0.72, fill = scales::alpha("#8DB1CC", 0.08), color = scales::alpha("#7CA3C2", 0.38),
      linewidth = 1.0
    ) +
    ggplot2::geom_col(
      data = bars_E |> dplyr::filter(w > 0),
      ggplot2::aes(x = plot_x, y = w), inherit.aes = FALSE,
      width = 0.92, fill = scales::alpha("#5F89AA", 0.08), color = NA
    ) +
    ggplot2::geom_col(
      data = bars_E |> dplyr::filter(w > 0),
      ggplot2::aes(x = plot_x, y = w), inherit.aes = FALSE,
      width = 0.78, fill = scales::alpha("#8DB1CC", 0.48), color = scales::alpha("#3E5C76", 0.88),
      linewidth = 1.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = obs_bg),
                        color = scales::alpha("#6A6A6A", 0.60), linewidth = 1.1,
                        linetype = "22", lineend = "round") +
    ggplot2::geom_point(
      data = data |> dplyr::filter(in_E),
      ggplot2::aes(x = plot_x, y = obs_bg), inherit.aes = FALSE,
      color = scales::alpha("#B12A1C", 0.80), size = 2.8
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_top), breaks = c(0, 0.05, 0.10, 0.15, 0.20),
      labels = scales::label_number(accuracy = 0.01),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(min(data$plot_x) - 0.5, max(data$plot_x) + 0.5),
      breaks = NULL, expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(x = NULL, y = y_lab, title = title_txt) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 19) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_blank(),
      axis.line.x        = ggplot2::element_blank(),
      axis.line.y        = ggplot2::element_blank(),
      axis.text          = ggplot2::element_text(size = 19, colour = "grey10"),
      axis.title         = ggplot2::element_text(size = 21, colour = "black", face = "plain"),
      axis.title.x       = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.ticks         = ggplot2::element_line(colour = "black"),
      axis.ticks.length  = grid::unit(2.5, "pt"),
      legend.position    = "none",
      plot.title         = ggplot2::element_text(size = 18, face = "plain", colour = "black",
                                                  hjust = 0.5, margin = ggplot2::margin(b = 2)),
      plot.margin        = ggplot2::margin(10, 10, 10, 10)
    )
}

p_arrows <- ggplot2::ggplot() +
  ggplot2::annotate("segment", x = 0.50, xend = 0.18, y = 0.30, yend = 0.30,
                     linewidth = 1.5, color = "#6C7A89",
                     arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed")) +
  ggplot2::annotate("segment", x = 0.50, xend = 0.82, y = 0.30, yend = 0.30,
                     linewidth = 1.5, color = "#6C7A89",
                     arrow = grid::arrow(length = grid::unit(0.18, "inches"), type = "closed")) +
  ggplot2::annotate("text", x = 0.50, y = 0.68, label = "Start From Uniform Weights",
                     fontface = "plain", size = 5.6, color = "grey10") +
  ggplot2::annotate("text", x = 0.18, y = 0.55, label = "Feasible Solution",
                     fontface = "plain", size = 5.2, color = "grey10") +
  ggplot2::annotate("text", x = 0.82, y = 0.55, label = "ET Solution",
                     fontface = "plain", size = 5.2, color = "grey10") +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 22, -6, 22),
                  plot.background = ggplot2::element_rect(fill = "white", color = NA))

# 18pt, not 15: Fig 2 is saved 12in wide against Fig 3's 10in, so a 15pt tag
# would render ~12.5pt once both are scaled to the same column width
p_left  <- make_weights_panel(df_plot_fig2, "w_bad",  "Alternative Weights", y_lab = "Weight") +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
    plot.tag.position = c(0, 0.99))
p_mid   <- make_weights_panel(df_plot_fig2, "w_unif", "Uniform Weights",  y_lab = NULL) +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                  axis.text.y  = ggplot2::element_blank(),
                  axis.ticks.y = ggplot2::element_blank(),
                  plot.tag          = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
                  plot.tag.position = c(0, 0.99))
p_right <- make_weights_panel(df_plot_fig2, "w_et",   "ET Weights",       y_lab = NULL) +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                  axis.text.y  = ggplot2::element_blank(),
                  axis.ticks.y = ggplot2::element_blank(),
                  plot.tag          = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
                  plot.tag.position = c(0, 0.99))

fig2 <- (p_left + patchwork::plot_spacer() + p_mid + patchwork::plot_spacer() + p_right +
     patchwork::plot_layout(ncol = 5, widths = c(1, 0.07, 1, 0.07, 1)))

ggplot2::ggsave(
  file.path(output_dir, "figure_2_et_weights.pdf"),
  plot = fig2, width = 12, height = 6, units = "in", device = grDevices::cairo_pdf
)

# ==============================================================================
# 6.  Alabama 2023/2024 and Oklahoma 2025/2026 examples in Figure 3
# ==============================================================================

# Plotting code...

CASE1 <- list(
  location        = "01",
  season          = "2023/2024",
  horizon         = 2L,
  model_a         = "MIGHTE-Nsemble",
  model_b         = "UMass-flusion",
  usual_rank_a    = 1L,
  weighted_rank_a = 3L,
  usual_rank_b    = 2L,
  weighted_rank_b = 1L,
  W               = 0.5,
  cost_type       = "growth",
  growth_quantile = 0.50
)

CASE2 <- list(
  location        = "40",
  season          = "2025/2026",
  horizon         = 2L,
  model_a         = "OHT_JHU-nbxd",
  model_b         = "FluSight-trained_mean",
  usual_rank_a    = NA_integer_,
  weighted_rank_a = NA_integer_,
  usual_rank_b    = NA_integer_,
  weighted_rank_b = NA_integer_,
  W               = 0.5,
  cost_type       = "peak",
  peak_window     = 3L
)

event_dates_1 <- get_growth_event_dates(truth_tbl, CASE1$location, CASE1$season, CASE1$growth_quantile)
event_dates_2 <- get_peak_event_dates(truth_tbl, CASE2$location, CASE2$season, CASE2$peak_window)

series_1a <- get_series(modern_quantile, truth_tbl, CASE1$location, CASE1$season, CASE1$horizon, CASE1$model_a)
series_1b <- get_series(modern_quantile, truth_tbl, CASE1$location, CASE1$season, CASE1$horizon, CASE1$model_b)
series_2a <- get_series(modern_quantile, truth_tbl, CASE2$location, CASE2$season, CASE2$horizon, CASE2$model_a)
series_2b <- get_series(modern_quantile, truth_tbl, CASE2$location, CASE2$season, CASE2$horizon, CASE2$model_b)

xlim1 <- range(series_1a$observed$target_end_date, na.rm = TRUE)
ylim1 <- range(c(series_1a$observed$observed,
                  series_1a$forecast$predicted, series_1b$forecast$predicted,
                  series_1a$interval$lower_90,   series_1a$interval$upper_90,
                  series_1b$interval$lower_90,   series_1b$interval$upper_90), na.rm = TRUE)

xlim2 <- range(series_2a$observed$target_end_date, na.rm = TRUE)
ylim2 <- range(c(series_2a$observed$observed,
                  series_2a$forecast$predicted, series_2b$forecast$predicted,
                  series_2a$interval$lower_90,   series_2a$interval$upper_90,
                  series_2b$interval$lower_90,   series_2b$interval$upper_90), na.rm = TRUE)

event_df_1 <- tibble::tibble(target_end_date = as.Date(event_dates_1))
event_df_2 <- tibble::tibble(target_end_date = as.Date(event_dates_2))

make_forecast_panel <- function(series, model_name, x_limits, y_limits, event_df,
                                 show_y_label = TRUE) {
  y_label <- if (show_y_label) "Incident Hospitalizations" else ""

  outside_event_df <- series$observed |>
    dplyr::inner_join(series$interval, by = "target_end_date") |>
    dplyr::filter(
      target_end_date %in% event_df$target_end_date,
      observed < lower_90 | observed > upper_90
    ) |>
    dplyr::select(target_end_date, observed)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data        = event_df,
      ggplot2::aes(xmin = target_end_date - 3.5, xmax = target_end_date + 3.5,
                    ymin = -Inf, ymax = Inf),
      fill        = "#F04E4E", alpha = 0.18, inherit.aes = FALSE
    ) +
    ggplot2::geom_ribbon(
      data        = series$interval,
      ggplot2::aes(x = target_end_date, ymin = lower_90, ymax = upper_90),
      fill        = "#2166AC", alpha = 0.15, inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data      = series$observed,
      ggplot2::aes(x = target_end_date, y = observed),
      colour    = "black", linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data   = series$observed |>
        dplyr::filter(!target_end_date %in% outside_event_df$target_end_date),
      ggplot2::aes(x = target_end_date, y = observed),
      colour = "black", size = 1.5
    ) +
    ggplot2::geom_line(
      data      = series$forecast,
      ggplot2::aes(x = target_end_date, y = predicted),
      colour    = "#2166AC", linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data   = series$forecast,
      ggplot2::aes(x = target_end_date, y = predicted),
      colour = "#2166AC", size = 1.5
    ) +
    ggplot2::scale_x_date(
      limits      = x_limits,
      date_breaks = "2 months",
      date_labels = "%b\n%Y"
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(title = model_name, x = "Date", y = y_label) +
    figure_theme() +
    ggplot2::geom_point(
      data        = outside_event_df,
      ggplot2::aes(x = target_end_date, y = observed),
      shape       = 4, colour = "#B12A1C", size = 1.8, stroke = 0.9,
      inherit.aes = FALSE
    )

  if (!show_y_label) p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  p
}

# the event-region key sits in the upper plot of the panel
plot_1a_left <- make_forecast_panel(series_1a, CASE1$model_a, xlim1, ylim1, event_df_1, show_y_label = TRUE) +
  # single key for the whole figure: event-week shading, forecast, observed
  ggplot2::annotate("rect",
    xmin = xlim1[2] - 60L, xmax = xlim1[2] - 39L,
    ymin = ylim1[2] * 0.93, ymax = ylim1[2] * 0.97,
    fill = "#F04E4E", alpha = 0.4) +
  ggplot2::annotate("text",
    x = xlim1[2] - 36L, y = ylim1[2] * 0.95,
    label = "Event Weeks", hjust = 0, size = 3.5, colour = "grey20") +
  ggplot2::annotate("segment",
    x = xlim1[2] - 60L, xend = xlim1[2] - 39L,
    y = ylim1[2] * 0.86, yend = ylim1[2] * 0.86,
    colour = "#2166AC", linewidth = 0.6) +
  ggplot2::annotate("text",
    x = xlim1[2] - 36L, y = ylim1[2] * 0.86,
    label = "2-Week Forecast", hjust = 0, size = 3.5, colour = "grey20") +
  ggplot2::annotate("segment",
    x = xlim1[2] - 60L, xend = xlim1[2] - 39L,
    y = ylim1[2] * 0.79, yend = ylim1[2] * 0.79,
    colour = "black", linewidth = 0.6) +
  ggplot2::annotate("text",
    x = xlim1[2] - 36L, y = ylim1[2] * 0.79,
    label = "Observed", hjust = 0, size = 3.5, colour = "grey20") +
  ggplot2::coord_cartesian(ylim = ylim1, clip = "off")
plot_1a_right <- make_forecast_panel(series_1b, CASE1$model_b, xlim1, ylim1, event_df_1, show_y_label = TRUE)

plot_2a_left  <- make_forecast_panel(series_2a, CASE2$model_a, xlim2, ylim2, event_df_2, show_y_label = TRUE)
plot_2a_right <- make_forecast_panel(series_2b, CASE2$model_b, xlim2, ylim2, event_df_2, show_y_label = TRUE)

subfig3a <- ((plot_1a_left + ggplot2::theme(axis.text.x = ggplot2::element_blank())) /
             plot_1a_right) +
  patchwork::plot_annotation(
    title = expression(bold("A") ~~~~~~~~~ "Alabama, 2023/2024 season"),
    theme = ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 15, face = "plain", colour = "black",
                                          margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(2, 8, 8, 8)
    )
  )

subfig3b <- ((plot_2a_left + ggplot2::theme(axis.text.x = ggplot2::element_blank())) /
             plot_2a_right) +
  patchwork::plot_annotation(
    title = expression(bold("B") ~~~~~~~~~ "Oklahoma, 2025/2026 season"),
    theme = ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 15, face = "plain", colour = "black",
                                          margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(2, 8, 8, 8)
    )
  )

fig3 <- (
  patchwork::wrap_elements(full = subfig3a) |
  patchwork::plot_spacer() |
  patchwork::wrap_elements(full = subfig3b)
) +
  patchwork::plot_layout(widths = c(1, 0.03, 1)) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))
  )

ggplot2::ggsave(
  file.path(output_dir, "figure_3_examples.png"),
  plot = fig3, width = 10, height = 8, dpi = 300
)

# ==============================================================================
# 7.  Fig 4 — Estimated CDFs of R_crit  +  Table 1
#
# Source: application/figure_rcrit_cdf.R
# R_crit = min_{b in B} (n_T - n_E)/n_E * (N_b - N_a)/(E_a - E_b)
# where a = PS winner, B = challengers with lower mean event-week WIS.
# R_crit = Inf when PS winner is unbeatable on event weeks.
# ==============================================================================

# ---- R_crit computation -----------------------------------------------------
compute_rcrit <- function(engine, event_dates, n_T, coverage, threshold, min_models,
                          allowed_models) {
  # engine returns list with (a) forecast data for a given location/horizon/season (b) a dataframe of the scores for these forecasts
  scores_obs <- engine$scores_obs
  # Number of event dates
  n_E        <- length(event_dates)
  
  # Obtain weekly scores for event dates
  scores_e <- scores_obs |>
    dplyr::mutate(E = dplyr::if_else(target_end_date %in% event_dates, 1L, 0L))

  # Total number of dates
  n_all_dates <- dplyr::n_distinct(scores_e$target_end_date)

  # This part filters the models to those satisfying the submission criteria
  eligible_models <- scores_e |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      n_dates_present   = dplyr::n_distinct(target_end_date), # number of total dates submitted
      n_flagged_present = dplyr::n_distinct(target_end_date[E == 1L]), # number of event dates submitted
      .groups = "drop"
    ) |>
    dplyr::filter(
      n_dates_present   >= floor(coverage  * n_all_dates), # require at least 0.90 of all dates
      n_flagged_present >= floor(threshold * n_E) # require ALL event week dates
    ) |>
    dplyr::pull(model)
  # Restrict to models also eligible under the other cost function
  eligible_models <- intersect(eligible_models, allowed_models)
  # If no models satisfy the filter then return nothing
  if (length(eligible_models) < min_models) return(NULL)

  # Compute event and non-event week WIS for each model
  model_stats <- scores_e |>
    dplyr::filter(model %in% eligible_models) |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      E_m       = mean(wis[E == 1L], na.rm = TRUE), # mean WIS on the event weeks,
      N_m       = mean(wis[E == 0L], na.rm = TRUE), # mean WIS on the non event weeks
      usual_wis = mean(wis,          na.rm = TRUE), # total WIS
      .groups   = "drop"
    )

  # a_idx = index of the top PS model
  a_idx       <- which.min(model_stats$usual_wis)
  E_a         <- model_stats$E_m[a_idx] # mean event week WIS for top PS model
  N_a         <- model_stats$N_m[a_idx] # mean non event week WIS for top PS model
  # return model scores for all models except top PS model; filter only to those that do better on the event weeks (if E_m > E_a then the model cannot do better in improvement)
  challengers <- model_stats[-a_idx, , drop = FALSE] |> dplyr::filter(E_m < E_a) 
  # Rcrit is infinite if no challengers
  if (nrow(challengers) == 0L) {
    return(list(r_crit = Inf, n_models = length(eligible_models), n_E = n_E, n_T = n_T))
  }
  # R_ab = Value at which top PS model is not the top under the WPS
  R_ab   <- ((n_T - n_E) / n_E) * (challengers$N_m - N_a) / (E_a - challengers$E_m)
  # Return min R_ab value needed to observe difference, 
  r_crit <- min(R_ab)

  list(r_crit = r_crit, n_models = length(eligible_models), n_E = n_E, n_T = n_T)
}

# ---- Main loop: compute R_crit per (cost_fn, season, horizon, location) -----
rcrit_rows <- list()

for (cost_fn in c("peak", "growth")) {
  for (seas in SEASONS) {
    for (horiz in HORIZONS) {
      message(sprintf("    %s | %s | h=%d", cost_fn, seas, horiz))

      n_skip_events <- 0L; n_skip_data <- 0L
      n_skip_engine <- 0L; n_skip_models <- 0L; n_computed <- 0L

      for (loc in locations) {

        # Gate 1: forecast data
        loc_data <- modern_quantile |> dplyr::filter(location == loc, horizon == horiz)
        if (nrow(loc_data) == 0L) { n_skip_data <- n_skip_data + 1L; next }

        # Gate 2: engine (needed before the event dates, to know the scored dates)
        engine <- tryCatch(
          suppressWarnings(build_rank_table_fast_engine(
            modern_all     = loc_data,
            location_value = loc,
            horizon_value  = horiz,
            season_value   = seas
          )),
          error = function(e) NULL
        )
        if (is.null(engine) || nrow(engine$scores_obs) == 0L) {
          n_skip_engine <- n_skip_engine + 1L; next
        }

        # Gate 3: both event definitions must be valid at this location
        scored_dates <- unique(as.Date(engine$scores_obs$target_end_date))
        events <- get_both_event_dates(truth_tbl, loc, seas, scored_dates)
        if (is.null(events)) { n_skip_events <- n_skip_events + 1L; next }
        event_dates <- events[[cost_fn]]

        # Models must satisfy the submission criteria under BOTH definitions
        allowed_models <- get_joint_eligible_models(engine, events, COVERAGE, THRESHOLD)

        # n_T: restrict 2025/2026 to post-Nov 22 2025
        scores_filtered <- engine$scores_obs
        if (seas == "2025/2026") {
          scores_filtered <- scores_filtered |>
            dplyr::filter(as.Date(target_end_date) >= as.Date("2025-11-22"))
        }
        n_T <- dplyr::n_distinct(scores_filtered$target_end_date)

        result <- compute_rcrit(
          engine = engine, event_dates = event_dates, n_T = n_T,
          coverage = COVERAGE, threshold = THRESHOLD, min_models = MIN_MODELS,
          allowed_models = allowed_models
        )

        if (is.null(result)) { n_skip_models <- n_skip_models + 1L; next }

        rcrit_rows[[length(rcrit_rows) + 1L]] <- tibble::tibble(
          cost_fn  = cost_fn, season = seas, horizon = horiz, location = loc,
          r_crit   = result$r_crit, n_models = result$n_models,
          n_E      = result$n_E,    n_T = result$n_T,
          is_inf   = !is.finite(result$r_crit)
        )
        n_computed <- n_computed + 1L
      }
      message(sprintf("      Computed: %d | Skip events: %d | Skip engine: %d | Skip models: %d",
                       n_computed, n_skip_events, n_skip_engine, n_skip_models))
    }
  }
}

rcrit_df <- dplyr::bind_rows(rcrit_rows)
message(sprintf("  Total R_crit: %d (finite: %d, Inf: %d)",
                 nrow(rcrit_df), sum(!rcrit_df$is_inf), sum(rcrit_df$is_inf)))

# ---- Manual ECDF (Inf in denominator, plateau below 1) ----------------------
ecdf_df <- rcrit_df |>
  dplyr::group_by(cost_fn, horizon, season) |>
  dplyr::group_modify(function(.x, .y) {
    total_n     <- nrow(.x)
    finite_vals <- sort(.x$r_crit[is.finite(.x$r_crit)])
    n_f         <- length(finite_vals)
    if (n_f == 0L) return(tibble::tibble(x = numeric(0), y = numeric(0)))
    tibble::tibble(
      x = c(1, finite_vals, 15.1),
      y = c(0, seq_len(n_f) / total_n, n_f / total_n)
    )
  }) |>
  dplyr::ungroup()

# ---- In-plot % R_crit = Inf annotations -------------------------------------
inf_label_df <- rcrit_df |>
  dplyr::group_by(cost_fn, horizon, season) |>
  dplyr::summarise(pct_inf = round(100 * mean(is_inf), 0), .groups = "drop") |>
  dplyr::arrange(cost_fn, horizon, season) |>
  dplyr::group_by(cost_fn, horizon) |>
  dplyr::mutate(y_pos = 0.97 - (dplyr::row_number() - 1L) * 0.08) |>
  dplyr::ungroup() |>
  dplyr::mutate(label = sprintf("%d%% (r_crit = ∞)", pct_inf))

# ---- Panel helper ------------------------------------------------------------
make_cdf_panel <- function(cost, horiz, show_legend = FALSE, show_y_label = TRUE) {
  df_plot   <- ecdf_df      |> dplyr::filter(cost_fn == cost, horizon == horiz)
  df_labels <- inf_label_df |> dplyr::filter(cost_fn == cost, horizon == horiz)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y, colour = season)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::geom_text(
      data        = df_labels,
      ggplot2::aes(x = 14.8, y = y_pos, label = label, colour = season),
      inherit.aes = FALSE, hjust = 1, vjust = 1,
      size        = 10 / ggplot2::.pt, fontface = "plain", show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = SEASON_COLOURS, name = "Season") +
    ggplot2::scale_x_continuous(breaks = c(1, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    ggplot2::coord_cartesian(xlim = c(1, 15), ylim = c(0, 1), expand = FALSE, clip = "on") +
    ggplot2::labs(
      subtitle = sprintf("Horizon = %d", horiz),
      x        = expression(italic(R)),
      y        = expression(Pr(r[crit] <= italic(r)))
    ) +
    figure_theme() +
    ggplot2::theme(
      plot.subtitle   = ggplot2::element_text(size = 11, colour = "black", hjust = 0,
                                               margin = ggplot2::margin(b = 2)),
      axis.text       = ggplot2::element_text(size = 12, colour = "grey10"),
      axis.title      = ggplot2::element_text(size = 13, colour = "black", face = "plain"),
      legend.position = if (show_legend) "right" else "none",
      legend.text     = ggplot2::element_text(size = 11, colour = "black"),
      legend.title    = ggplot2::element_text(size = 12, colour = "black", face = "plain"),
      plot.margin     = ggplot2::margin(2, 8, 8, 8)
    )

  if (!show_y_label) p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  p
}

panel_peak_h2   <- make_cdf_panel("peak",   2L, show_legend = FALSE, show_y_label = TRUE)
panel_peak_h3   <- make_cdf_panel("peak",   3L, show_legend = TRUE,  show_y_label = FALSE)
panel_growth_h2 <- make_cdf_panel("growth", 2L, show_legend = FALSE, show_y_label = TRUE)
panel_growth_h3 <- make_cdf_panel("growth", 3L, show_legend = TRUE,  show_y_label = FALSE)

subfig4a <- (panel_peak_h2 | panel_peak_h3) +
  patchwork::plot_annotation(
    title = expression(bold("A") ~~~~~~~~~~~ "Seasonal Peak Event Week"),
    theme = ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 13, face = "plain", colour = "black",
                                          margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(2, 8, 8, 8)
    )
  )

subfig4b <- (panel_growth_h2 | panel_growth_h3) +
  patchwork::plot_annotation(
    title = expression(bold("B") ~~~~~~~~~~~ "High-Growth Event Week"),
    theme = ggplot2::theme(
      plot.title  = ggplot2::element_text(size = 13, face = "plain", colour = "black",
                                          margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(2, 8, 8, 8)
    )
  )

fig4 <- (
  patchwork::wrap_elements(full = subfig4a) /
  patchwork::plot_spacer() /
  patchwork::wrap_elements(full = subfig4b)
) +
  patchwork::plot_layout(heights = c(1, 0.03, 1)) +
  patchwork::plot_annotation(
    caption = expression(italic(r)),
    theme   = ggplot2::theme(
      plot.margin  = ggplot2::margin(8, 8, 8, 8),
      plot.caption = ggplot2::element_text(size = 13, hjust = 0.5, colour = "black",
                                            margin = ggplot2::margin(t = 2))
    )
  )

ggplot2::ggsave(
  file.path(output_dir, "figure_4_rcrit_cdf.png"),
  plot = fig4, width = 8, height = 8, dpi = 300
)
message("Saved: figure_4_rcrit_cdf.png")

# ---- Table 1: R_crit summary (LaTeX) ----------------------------------------
# P(R_crit <= r) uses all eligible locations in the denominator (Inf counts as
# R_crit > r), consistent with the manual ECDF used for Figure 4.
summary_rcrit <- rcrit_df |>
  dplyr::group_by(cost_fn, season, horizon) |>
  dplyr::summarise(
    n_locs  = dplyr::n(),
    pct_inf = round(100 * mean(is_inf), 1),
    p_le_2  = round(mean(r_crit <= 2), 3),
    p_le_4  = round(mean(r_crit <= 4), 3),
    .groups = "drop"
  ) |>
  dplyr::arrange(cost_fn, season, horizon)

# Insert season and cost-function multirow labels manually
cost_labels <- c(peak = "Seasonal peak", growth = "High growth")

tex_body <- character(0)
cf_groups <- split(seq_len(nrow(summary_rcrit)), summary_rcrit$cost_fn)
cf_order  <- c("peak", "growth")

for (cf in cf_order) {
  idx      <- cf_groups[[cf]]
  seas_grp <- split(idx, summary_rcrit$season[idx])
  seas_ord <- SEASONS[SEASONS %in% names(seas_grp)]
  n_cf     <- length(idx)

  tex_body <- c(tex_body,
    sprintf("\\multirow{%d}{*}{%s}", n_cf, cost_labels[cf]))

  first_seas <- TRUE
  for (seas in seas_ord) {
    sidx <- seas_grp[[seas]]
    n_s  <- length(sidx)
    if (first_seas) {
      tex_body <- c(tex_body,
        sprintf("  & \\multirow{%d}{*}{%s}", n_s, seas))
      first_seas <- FALSE
    } else {
      tex_body <- c(tex_body, "  \\cmidrule{2-7}",
        sprintf("  & \\multirow{%d}{*}{%s}", n_s, seas))
    }
    for (i in sidx) {
      r <- summary_rcrit[i, ]
      tex_body <- c(tex_body,
        sprintf("  &                            & %d & %d & %.1f & %.3f & %.3f \\\\",
                r$horizon, r$n_locs, r$pct_inf, r$p_le_2, r$p_le_4))
    }
  }
  if (cf != tail(cf_order, 1)) tex_body <- c(tex_body, "\\midrule")
}

table_tex <- c(
  "\\begin{table}[ht]",
  "\\centering",
  "\\begin{tabular}{llccccc}",
  "\\toprule",
  "Cost function & Season & $h$ & $n$ & \\% $R_{\\text{crit}} = \\infty$ & $\\Pr(R_{\\text{crit}} \\leq 2)$ & $\\Pr(R_{\\text{crit}} \\leq 4)$ \\\\",
  "\\midrule",
  tex_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\caption{Summary of $R_{\\text{crit}}$ across locations by cost function, season, and horizon~$h$.",
  "$n$ is the number of eligible locations.",
  "The column \\% $R_{\\text{crit}} = \\infty$ gives the percentage of eligible locations where the",
  "PS top model is unbeatable under any degree of event-week upweighting.",
  "$\\Pr(R_{\\text{crit}} \\leq r)$ is computed over all eligible locations, with $R_{\\text{crit}} = \\infty$",
  "treated as $R_{\\text{crit}} > r$.}",
  "\\label{tab:rcrit_summary}",
  "\\end{table}"
)

writeLines(table_tex, file.path(output_dir, "table_1_rcrit_summary.tex"))
message("Saved: table_1_rcrit_summary.tex")

# ==============================================================================
# 8.  Fig 5 — Improvement and penalty tradeoff across W values
#
# Source: application/figure_event_tradeoff.R
# Peak cost function only. For each (season, horizon, location) where the top
# PS and WPS models differ, plots event-week WIS improvement vs. non-event-week
# WIS penalty, for w in {0.3, 0.5, 0.9}.
# ==============================================================================

message("\n=== Fig 5: Improvement/penalty tradeoff ===")

# ---- Build per-location WIS means and top-model IDs -------------------------
message("  Computing rank tables across W values...")

wis_model_rows <- list()
top_model_rows <- list()

for (seas in SEASONS) {
  for (horiz in HORIZONS) {
    for (loc in locations) {

      # Gate 1: forecast data
      loc_data <- modern_quantile |> dplyr::filter(location == loc, horizon == horiz)
      if (nrow(loc_data) == 0L) next

      # Gate 2: engine (built ONCE per location, before the event dates so the
      # scored dates are known)
      engine <- tryCatch(
        suppressWarnings(build_rank_table_fast_engine(
          modern_all     = loc_data,
          location_value = loc,
          horizon_value  = horiz,
          season_value   = seas
        )),
        error = function(e) NULL
      )
      if (is.null(engine) || nrow(engine$scores_obs) == 0L) next

      # Gate 3: both event definitions must be valid; peak drives this figure
      scored_dates <- unique(as.Date(engine$scores_obs$target_end_date))
      events <- get_both_event_dates(truth_tbl, loc, seas, scored_dates)
      if (is.null(events)) next
      event_dates <- events$peak

      # Models must satisfy the submission criteria under BOTH definitions
      allowed_models <- get_joint_eligible_models(engine, events, COVERAGE, THRESHOLD)
      if (length(allowed_models) < MIN_MODELS) next

      # Per-model WIS means (w independent, computed once)
      event_dates_d <- as.Date(event_dates)
      scores_e <- engine$scores_obs |>
        dplyr::mutate(
          target_end_date = as.Date(target_end_date),
          E               = dplyr::if_else(target_end_date %in% event_dates_d, 1L, 0L)
        )

      per_model_wis <- scores_e |>
        dplyr::filter(model %in% allowed_models) |>
        dplyr::group_by(model) |>
        dplyr::summarise(
          mean_event_wis    = mean(wis[E == 1L], na.rm = TRUE),
          mean_nonevent_wis = mean(wis[E == 0L], na.rm = TRUE),
          mean_all_wis      = mean(wis,          na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::mutate(season = seas, horizon = horiz, location = loc)

      wis_model_rows[[length(wis_model_rows) + 1L]] <- per_model_wis

      # Gate 4 + top-model identification: once per w value
      for (w in W_VALUES_FIG5) {
        rank_tbl_w <- tryCatch(
          suppressWarnings(make_rank_table_fast(
            target_end_dates = event_dates,
            W                = w,
            threshold        = THRESHOLD,
            coverage         = COVERAGE,
            engine           = engine
          )),
          error = function(e) NULL
        )
        if (is.null(rank_tbl_w)) next
        # restrict to the joint pool and re-rank within it
        rank_tbl_w <- rank_tbl_w |>
          dplyr::filter(model %in% allowed_models) |>
          dplyr::mutate(
            weighted_rank = dplyr::min_rank(weighted_wis),
            usual_rank    = dplyr::min_rank(usual_wis)
          )
        if (nrow(rank_tbl_w) < MIN_MODELS) next

        ps1_model  <- rank_tbl_w$model[rank_tbl_w$usual_rank    == min(rank_tbl_w$usual_rank)][1L]
        wps1_model <- rank_tbl_w$model[rank_tbl_w$weighted_rank == min(rank_tbl_w$weighted_rank)][1L]

        if (ps1_model == wps1_model) next

        top_model_rows[[length(top_model_rows) + 1L]] <- tibble::tibble(
          season = seas, horizon = horiz, location = loc,
          W = w, ps1_model = ps1_model, wps1_model = wps1_model
        )
      }
    }
  }
}

top_all <- dplyr::bind_rows(top_model_rows)
wis_all <- dplyr::bind_rows(wis_model_rows)

for (w in W_VALUES_FIG5) {
  message(sprintf("  w = %.1f: %d changed-top locations", w, sum(top_all$W == w)))
}

# ---- Build plot data per W value --------------------------------------------
build_tradeoff_df <- function(top_model_rows, wis_model_rows, w_val) {
  tops <- top_model_rows |>
    dplyr::filter(W == w_val) |>
    dplyr::select(season, horizon, location, ps1_model, wps1_model)
  if (nrow(tops) == 0L) return(tibble::tibble())

  wis <- wis_model_rows |>
    dplyr::select(season, horizon, location, model, mean_event_wis, mean_nonevent_wis, mean_all_wis)

  tops |>
    dplyr::left_join(
      wis |> dplyr::rename(mean_event_wis_PS_top    = mean_event_wis,
                             mean_nonevent_wis_PS_top = mean_nonevent_wis,
                             mean_all_wis_PS_top      = mean_all_wis),
      by = c("season", "horizon", "location", "ps1_model" = "model")
    ) |>
    dplyr::left_join(
      wis |> dplyr::select(season, horizon, location, model, mean_event_wis, mean_nonevent_wis) |>
        dplyr::rename(mean_event_wis_WPS_top    = mean_event_wis,
                       mean_nonevent_wis_WPS_top = mean_nonevent_wis),
      by = c("season", "horizon", "location", "wps1_model" = "model")
    ) |>
    dplyr::mutate(
      B = mean_all_wis_PS_top,
      x = (mean_nonevent_wis_WPS_top - mean_nonevent_wis_PS_top) / B,
      y = (mean_event_wis_PS_top     - mean_event_wis_WPS_top)   / B
    ) |>
    dplyr::filter(!is.na(x), !is.na(y))
}

plot_df_03 <- build_tradeoff_df(top_all, wis_all, 0.3)
plot_df_05 <- build_tradeoff_df(top_all, wis_all, 0.5)
plot_df_09 <- build_tradeoff_df(top_all, wis_all, 0.9)

# Shared axis range
all_vals <- c(plot_df_03$x, plot_df_03$y, plot_df_05$x, plot_df_05$y, plot_df_09$x, plot_df_09$y)
ax_pad   <- diff(range(all_vals, na.rm = TRUE)) * 0.05
ax_min   <- min(all_vals, na.rm = TRUE) - ax_pad
ax_max   <- max(all_vals, na.rm = TRUE) + ax_pad

# ---- Panel builder ----------------------------------------------------------
make_tradeoff_panel <- function(df, ax_min, ax_max, w_val, show_y_label = TRUE) {
  y_label     <- if (show_y_label) "Improvement" else ""
  bound_slope <- (1 - w_val) / w_val

  ribbon_x  <- seq(ax_min, ax_max, length.out = 300)
  ribbon_df <- data.frame(
    x    = ribbon_x,
    ymin = pmax(bound_slope * ribbon_x, ax_min),
    ymax = ax_max
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_ribbon(
      data = ribbon_df,
      ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
      fill = "grey85", alpha = 0.6, inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(slope = bound_slope, intercept = 0,
                          linetype = "solid", colour = "grey20", linewidth = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dotted", colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_point(colour = "#B12A1C", alpha = 0.65, size = 2.2) +
    ggplot2::scale_x_continuous(limits = c(ax_min, ax_max)) +
    ggplot2::scale_y_continuous(limits = c(ax_min, ax_max)) +
    ggplot2::labs(x = "Penalty", y = y_label) +
    figure_theme() +
    ggplot2::theme(
      axis.text    = ggplot2::element_text(size = 20, colour = "grey10"),
      axis.title.y = ggplot2::element_text(size = 21, colour = "black", face = "plain"),
      axis.title.x = ggplot2::element_text(size = 21, colour = "black", face = "plain")
    )

  if (!show_y_label) p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  p
}

n_swaps_03 <- nrow(plot_df_03)
n_swaps_05 <- nrow(plot_df_05)
n_swaps_09 <- nrow(plot_df_09)

panel_w03 <- make_tradeoff_panel(plot_df_03, ax_min, ax_max, 0.3, show_y_label = TRUE) +
  ggplot2::labs(tag = "A", subtitle = sprintf("w = 0.3   (n = %d)", n_swaps_03)) +
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 19, face = "bold", colour = "black"),
    plot.tag.position = c(0, 0.99),
    plot.subtitle     = ggplot2::element_text(size = 19, face = "plain", colour = "black",
                                               hjust = 0.5, margin = ggplot2::margin(b = 2)),
    plot.margin       = ggplot2::margin(2, 8, 8, 8)
  )

panel_w05 <- make_tradeoff_panel(plot_df_05, ax_min, ax_max, 0.5, show_y_label = FALSE) +
  ggplot2::labs(tag = "B", subtitle = sprintf("w = 0.5   (n = %d)", n_swaps_05)) +
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 19, face = "bold", colour = "black"),
    plot.tag.position = c(0, 0.99),
    plot.subtitle     = ggplot2::element_text(size = 19, face = "plain", colour = "black",
                                               hjust = 0.5, margin = ggplot2::margin(b = 2)),
    plot.margin       = ggplot2::margin(2, 8, 8, 8)
  )

panel_w09 <- make_tradeoff_panel(plot_df_09, ax_min, ax_max, 0.9, show_y_label = FALSE) +
  ggplot2::labs(tag = "C", subtitle = sprintf("w = 0.9   (n = %d)", n_swaps_09)) +
  ggplot2::theme(
    plot.tag          = ggplot2::element_text(size = 19, face = "bold", colour = "black"),
    plot.tag.position = c(0, 0.99),
    plot.subtitle     = ggplot2::element_text(size = 19, face = "plain", colour = "black",
                                               hjust = 0.5, margin = ggplot2::margin(b = 2)),
    plot.margin       = ggplot2::margin(2, 8, 8, 8)
  )

fig5 <- (panel_w03 | panel_w05 | panel_w09) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
  )

ggplot2::ggsave(
  file.path(output_dir, "figure_5_tradeoff.png"),
  plot = fig5, width = 14, height = 5, dpi = 300
)
message("Saved: figure_5_tradeoff.png")

# ==============================================================================
# 9.  Fig 6 — Weekly contribution to the PS and WPS (US, 2025/2026, h = 2)
#
#   A  observed trajectory               both event regions shaded
#   B  weekly share of the PS            both event regions shaded
#   C  weekly share of the WPS, peak     peak region shaded
#   D  weekly share of the WPS, growth   growth region shaded
#
# The WPS is one number per model, not a time series, so each week is shown as
# its contribution to the total — w_t * WIS_t for the WPS (w_t = W/|E| on event
# weeks, (1-W)/(T-|E|) elsewhere) and WIS_t/T for the PS — divided by that
# model's own total. Both panels then read as "% of this model's score", sum to
# 100% across weeks, and can share a y scale.
#
# Models are restricted to those eligible under BOTH cost functions and
# re-ranked within that pool. The PS is then identical across the two cost
# functions, which is why A and B each collapse to a single panel.
# ==============================================================================

FIG6_SEASON     <- "2025/2026"
FIG6_LOCATION   <- "US"
FIG6_HORIZON    <- 2L
FIG6_R          <- 2       # event week counts R times a non-event week
# W depends on |E| and T, so it is derived per cost function rather than fixed:
#   R = W/(1-W) * (T-|E|)/|E|   =>   W = R|E| / ((T-|E|) + R|E|)
fig6_w_from_R <- function(n_E, n_T) (FIG6_R * n_E) / ((n_T - n_E) + FIG6_R * n_E)
FIG6_TOP_N      <- 5L
FIG6_DROP_WEEKS <- 4L      # season opens with no forecasts, only whitespace

# y limit and colours are pinned rather than derived, so that changing the
# weighting changes only which models are shown and the curves themselves --
# not the axis range or which colour each model gets.
FIG6_YMAX <- 0.40

FIG6_PALETTE <- c(
  "CEPH-Rtrend_fluH"         = "#0072B2",
  "CMU-TimeSeries"           = "#D55E00",
  "MIGHTE-Nsemble"           = "#009E73",
  "MOBS-GLEAM_RL_FLUH"       = "#CC79A7",
  "NAU-vulPES"               = "#E69F00",
  "OHT_JHU-nbxd"             = "#7570B3",
  "UGA_flucast-Scenariocast" = "#A6761D",
  "UMass-flusion"            = "#56B4E9"
)

# Growth weeks sit inside the peak window, so the growth band is drawn slightly
# wider to leave a visible blue -> purple -> red gradient at the overlap. This
# is presentation only; the event weeks themselves are unchanged.
FIG6_PEAK_PAD   <- 3.5     # days either side of an event week
FIG6_GROWTH_PAD <- 5.5
# Same red and blue used elsewhere in this script: #F04E4E is the event-week
# shading in Figs 2 and 3, #2166AC the forecast blue in Figs 1, 3 and 5.
FIG6_SHADE_PEAK   <- list(fill = "#F04E4E", alpha = 0.18)
FIG6_SHADE_GROWTH <- list(fill = "#2166AC", alpha = 0.18)

# Hand-picked qualitative colours: the default hue palette places near-identical
# greens side by side, and black is reserved for the observed trajectory.
FIG6_COLS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#7570B3", "#A6761D", "#56B4E9", "#E7298A", "#1B9E77")

fig6_engine <- build_rank_table_fast_engine(
  modern_all     = modern_quantile |>
    dplyr::filter(location == FIG6_LOCATION, horizon == FIG6_HORIZON),
  location_value = FIG6_LOCATION,
  horizon_value  = FIG6_HORIZON,
  season_value   = FIG6_SEASON
)

fig6_scored    <- unique(as.Date(fig6_engine$scores_obs$target_end_date))
fig6_ev_peak   <- get_peak_event_dates(truth_tbl, FIG6_LOCATION, FIG6_SEASON, PEAK_WINDOW)
fig6_ev_growth <- get_growth_event_dates(truth_tbl, FIG6_LOCATION, FIG6_SEASON)
fig6_ev_growth <- fig6_ev_growth[fig6_ev_growth %in% fig6_scored]

fig6_n_T <- dplyr::n_distinct(fig6_engine$scores_obs$target_end_date)

# Models eligible under both cost functions
fig6_models <- intersect(
  make_rank_table_fast(target_end_dates = fig6_ev_peak,
                       W = fig6_w_from_R(length(fig6_ev_peak), fig6_n_T),
                       threshold = THRESHOLD, coverage = COVERAGE,
                       engine = fig6_engine)$model,
  make_rank_table_fast(target_end_dates = fig6_ev_growth,
                       W = fig6_w_from_R(length(fig6_ev_growth), fig6_n_T),
                       threshold = THRESHOLD, coverage = COVERAGE,
                       engine = fig6_engine)$model
)

# Weekly shares plus ranks recomputed within the common model pool
fig6_prepare <- function(event_dates) {
  sc <- fig6_engine$scores_obs |>
    dplyr::mutate(target_end_date = as.Date(target_end_date),
                  is_event = target_end_date %in% event_dates) |>
    dplyr::filter(model %in% fig6_models)

  n_T <- dplyr::n_distinct(sc$target_end_date)
  n_E <- length(event_dates)
  FIG6_W <- fig6_w_from_R(n_E, n_T)

  ranks <- sc |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      usual_wis    = mean(wis, na.rm = TRUE),
      weighted_wis = FIG6_W * mean(wis[is_event], na.rm = TRUE) +
                     (1 - FIG6_W) * mean(wis[!is_event], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::transmute(model,
                     usual_rank    = dplyr::min_rank(usual_wis),
                     weighted_rank = dplyr::min_rank(weighted_wis))

  sc |>
    dplyr::mutate(
      w_t         = dplyr::if_else(is_event, FIG6_W / n_E, (1 - FIG6_W) / (n_T - n_E)),
      contrib_wps = w_t * wis,
      contrib_ps  = wis / n_T
    ) |>
    dplyr::group_by(model) |>
    dplyr::mutate(
      pct_wps = contrib_wps / sum(contrib_wps, na.rm = TRUE),
      pct_ps  = contrib_ps  / sum(contrib_ps,  na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(ranks, by = "model") |>
    dplyr::mutate(in_top_wps = weighted_rank <= FIG6_TOP_N,
                  in_top_ps  = usual_rank    <= FIG6_TOP_N)
}

fig6_peak   <- fig6_prepare(fig6_ev_peak)
fig6_growth <- fig6_prepare(fig6_ev_growth)

fig6_obs <- truth_tbl |>
  dplyr::filter(location == FIG6_LOCATION, tgt_season == FIG6_SEASON) |>
  dplyr::arrange(target_end_date)

# x starts where every plotted model has a forecast, otherwise a model with an
# earlier first submission draws a lone line stretching left of all the others
fig6_shown <- unique(c(fig6_peak$model[fig6_peak$in_top_ps | fig6_peak$in_top_wps],
                       fig6_growth$model[fig6_growth$in_top_wps]))
fig6_x_lo <- max(
  min(fig6_obs$target_end_date) + FIG6_DROP_WEEKS * 7L,
  dplyr::bind_rows(fig6_peak, fig6_growth) |>
    dplyr::filter(model %in% fig6_shown) |>
    dplyr::group_by(model) |>
    dplyr::summarise(first = min(target_end_date), .groups = "drop") |>
    dplyr::pull(first) |> max()
)
fig6_xlim <- c(fig6_x_lo, max(fig6_obs$target_end_date))

fig6_ylim <- c(0, FIG6_YMAX)

fig6_palette <- FIG6_PALETTE[sort(fig6_shown)]
stopifnot(!anyNA(fig6_palette))   # every shown model must have a pinned colour

# One rectangle per contiguous run of event weeks. With a pad wider than half
# the 7-day spacing, per-week rectangles would overlap and the alpha would
# compound in the overlaps, rendering as vertical stripes.
fig6_band <- function(event_dates, pad) {
  ev  <- sort(event_dates)
  run <- cumsum(c(1, as.numeric(diff(ev)) > 7))
  data.frame(
    xmin = as.Date(tapply(ev, run, min), origin = "1970-01-01") - pad,
    xmax = as.Date(tapply(ev, run, max), origin = "1970-01-01") + pad
  )
}

fig6_shade <- function(p, which) {
  if (which %in% c("peak", "both")) {
    p <- p + ggplot2::geom_rect(
      data = fig6_band(fig6_ev_peak, FIG6_PEAK_PAD),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = FIG6_SHADE_PEAK$fill, alpha = FIG6_SHADE_PEAK$alpha)
  }
  if (which %in% c("growth", "both")) {
    p <- p + ggplot2::geom_rect(
      data = fig6_band(fig6_ev_growth, FIG6_GROWTH_PAD),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = FIG6_SHADE_GROWTH$fill, alpha = FIG6_SHADE_GROWTH$alpha)
  }
  p
}

# figure_theme() supplies the shared look; these are the Fig 6 additions:
# a bold panel letter, an in-panel legend, and x tick labels on the bottom row.
fig6_theme <- function(bottom) {
  figure_theme() +
    ggplot2::theme(
      plot.tag             = ggplot2::element_text(size = 15, face = "bold", colour = "black"),
      plot.tag.position    = c(0, 0.99),
      axis.text.x          = if (bottom) ggplot2::element_text(size = 14, colour = "grey10")
                             else ggplot2::element_blank(),
      plot.margin          = ggplot2::margin(10, 8, 10, 8),
      legend.position      = c(0.995, 0.97),
      legend.justification = c("right", "top"),
      legend.background    = ggplot2::element_rect(fill = scales::alpha("white", 0.75),
                                                   colour = NA),
      legend.key           = ggplot2::element_blank(),
      legend.key.height    = grid::unit(15, "pt"),
      legend.margin        = ggplot2::margin(3, 5, 3, 5),
      legend.title         = ggplot2::element_text(size = 13, colour = "black"),
      legend.text          = ggplot2::element_text(size = 12, colour = "black")
    )
}

fig6_x_scale <- ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")

fig6_score_panel <- function(scores, flag, y_col, rank_col, tag, y_lab,
                             legend_title, shading, bottom) {
  d   <- scores |> dplyr::filter(.data[[flag]])
  key <- d |> dplyr::distinct(model, rank = .data[[rank_col]]) |> dplyr::arrange(rank)

  fig6_shade(ggplot2::ggplot(), shading) +
    ggplot2::geom_line(data = d,
      ggplot2::aes(target_end_date, .data[[y_col]], colour = model), linewidth = 0.7) +
    ggplot2::geom_point(data = d,
      ggplot2::aes(target_end_date, .data[[y_col]], colour = model), size = 1.3) +
    ggplot2::scale_colour_manual(
      values = fig6_palette, breaks = key$model,
      labels = stats::setNames(paste0(key$rank, ". ", key$model), key$model)) +
    fig6_x_scale +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::coord_cartesian(xlim = fig6_xlim, ylim = fig6_ylim) +
    ggplot2::labs(tag = tag, x = NULL, y = y_lab, colour = legend_title) +
    fig6_theme(bottom)
}

# Panel A maps the two regions to a fill scale purely so the legend can name
# them; the keys are drawn more opaque than the bands, which at these alphas
# would be near-white at key size.
fig6_bands_legend <- rbind(
  transform(fig6_band(fig6_ev_peak,   FIG6_PEAK_PAD),   region = "Seasonal Peak"),
  transform(fig6_band(fig6_ev_growth, FIG6_GROWTH_PAD), region = "High Growth")
)
fig6_bands_legend$region <- factor(fig6_bands_legend$region,
                                   levels = c("Seasonal Peak", "High Growth"))

fig6_panel_a <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = fig6_bands_legend,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                 fill = region, alpha = region),
    inherit.aes = FALSE) +
  ggplot2::scale_fill_manual(
    values = c("Seasonal Peak" = FIG6_SHADE_PEAK$fill,
               "High Growth"   = FIG6_SHADE_GROWTH$fill), name = NULL) +
  ggplot2::scale_alpha_manual(
    values = c("Seasonal Peak" = FIG6_SHADE_PEAK$alpha,
               "High Growth"   = FIG6_SHADE_GROWTH$alpha), guide = "none") +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(
    alpha = c(FIG6_SHADE_PEAK$alpha, FIG6_SHADE_GROWTH$alpha) * 3))) +
  ggplot2::geom_line(data = fig6_obs,
    ggplot2::aes(target_end_date, observed), linewidth = 0.7) +
  ggplot2::geom_point(data = fig6_obs,
    ggplot2::aes(target_end_date, observed), size = 1.4) +
  fig6_x_scale +
  ggplot2::coord_cartesian(xlim = fig6_xlim) +
  ggplot2::labs(tag = "A", x = NULL, y = "Incident Hospitalizations") +
  fig6_theme(bottom = FALSE)

fig6_panel_b <- fig6_score_panel(fig6_peak, "in_top_ps", "pct_ps", "usual_rank",
                                 "B", "% PS", "Top 5 by PS",
                                 shading = "both", bottom = FALSE)

fig6_panel_c <- fig6_score_panel(fig6_peak, "in_top_wps", "pct_wps", "weighted_rank",
                                 "C", "% WPS", "Top 5 by Seasonal Peak WPS",
                                 shading = "peak", bottom = TRUE)

fig6_panel_d <- fig6_score_panel(fig6_growth, "in_top_wps", "pct_wps", "weighted_rank",
                                 "D", "% WPS", "Top 5 by High Growth WPS",
                                 shading = "growth", bottom = TRUE)

fig6 <- (fig6_panel_a | fig6_panel_b) / (fig6_panel_c | fig6_panel_d)

ggplot2::ggsave(
  file.path(output_dir, "figure_6_weekly_score_panels.png"),
  plot = fig6, width = 10, height = 6, dpi = 300
)
message("Saved: figure_6_weekly_score_panels.png")

# ==============================================================================
# 10. Fig 7 — Rank slopegraph: PS vs WPS (US, 2025/2026, h = 2)
#
# One line per model, from its standardised PS score to its standardised WPS
# score. Ranks alone hide how far apart models are; standardising within each
# column keeps the spacing, so a line's slope shows the rank change while its
# height shows how crowded that model is. Lower is better on both sides.
#
# Weighting is set by the emphasis ratio R rather than W: R is the weight of an
# event week relative to a non-event week, so it means the same thing in both
# panels, whereas a fixed W would imply different per-week emphasis because the
# two event definitions have different |E|. Inverting R = W/(1-W) * (T-|E|)/|E|
# gives the W passed to the scorer.
#
# Models are restricted to those eligible under BOTH event definitions and
# re-ranked within that pool, so the PS column is identical across panels.
# ==============================================================================

FIG7_SEASON      <- "2025/2026"
FIG7_LOCATION    <- "US"
FIG7_HORIZON     <- 2L
FIG7_R           <- 2      # event week carries R times a non-event week
FIG7_LABEL_SHIFT <- 3L     # label a model only if its rank moves this far

# Hand-picked qualitative colours, dealt out by vertical position below. The
# default hue palette spaces hues evenly, which hands near-identical colours to
# models that happen to sit next to each other.
FIG7_COLS <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
               "#7570B3", "#1B9E77", "#E7298A", "#A6761D", "#56B4E9",
               "#66A61E", "#B2182B", "#4575B4", "#762A83")

fig7_w_from_R <- function(R, n_E, n_T) (R * n_E) / ((n_T - n_E) + R * n_E)

# ---- scores, standardised within each panel ---------------------------------
fig7_engine <- build_rank_table_fast_engine(
  modern_all     = modern_quantile |>
    dplyr::filter(location == FIG7_LOCATION, horizon == FIG7_HORIZON),
  location_value = FIG7_LOCATION,
  horizon_value  = FIG7_HORIZON,
  season_value   = FIG7_SEASON
)

fig7_scored <- unique(as.Date(fig7_engine$scores_obs$target_end_date))
fig7_events <- get_both_event_dates(truth_tbl, FIG7_LOCATION, FIG7_SEASON, fig7_scored)
fig7_models <- get_joint_eligible_models(fig7_engine, fig7_events, COVERAGE, THRESHOLD)
fig7_n_T    <- dplyr::n_distinct(fig7_engine$scores_obs$target_end_date)

fig7_scores <- dplyr::bind_rows(lapply(c("peak", "growth"), function(cost_fn) {
  event_dates <- fig7_events[[cost_fn]]
  W <- fig7_w_from_R(FIG7_R, length(event_dates), fig7_n_T)

  make_rank_table_fast(
    target_end_dates = event_dates, W = W,
    threshold = THRESHOLD, coverage = COVERAGE, engine = fig7_engine
  ) |>
    dplyr::filter(model %in% fig7_models) |>
    dplyr::mutate(
      cost_fn  = cost_fn,
      PS_rank  = dplyr::min_rank(usual_wis),
      WPS_rank = dplyr::min_rank(weighted_wis),
      PS_z     = (usual_wis    - mean(usual_wis))    / stats::sd(usual_wis),
      WPS_z    = (weighted_wis - mean(weighted_wis)) / stats::sd(weighted_wis)
    )
}))

# ---- long form: one row per model per column --------------------------------
fig7_long <- fig7_scores |>
  dplyr::mutate(labelled = abs(WPS_rank - PS_rank) >= FIG7_LABEL_SHIFT |
                           PS_rank == 1L | WPS_rank == 1L) |>
  tidyr::pivot_longer(c(PS_z, WPS_z), names_to = "score", values_to = "z") |>
  dplyr::mutate(
    score = factor(dplyr::if_else(score == "PS_z", "PS", "WPS"), levels = c("PS", "WPS")),
    rank  = dplyr::if_else(score == "PS", PS_rank, WPS_rank),
    label = paste0(model, " (#", rank, ")")
  )

fig7_movers  <- fig7_long |> dplyr::filter(labelled)
fig7_stayers <- fig7_long |> dplyr::filter(!labelled)

# colours assigned by vertical position, so neighbours are maximally separated
fig7_order <- fig7_movers |>
  dplyr::filter(score == "PS") |>
  dplyr::group_by(model) |>
  dplyr::summarise(pos = mean(z), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(pos)) |>
  dplyr::pull(model)

fig7_palette <- stats::setNames(rep(FIG7_COLS, length.out = length(fig7_order)),
                                fig7_order)

# Labels sit in the gutters outside the plotting region; xlim bounds each side
# to its own margin so the two columns cannot collide.
fig7_label_layer <- function(side, nudge, label_xlim) {
  ggrepel::geom_text_repel(
    data        = fig7_movers |> dplyr::filter(score == side),
    ggplot2::aes(label = label, colour = model),
    hjust       = if (side == "PS") 1 else 0,
    nudge_x     = nudge, direction = "y", xlim = label_xlim,
    segment.size = 0.25, segment.alpha = 0.6,
    size = 4.6, box.padding = 0.28, point.padding = 0.15,
    force = 4, force_pull = 0.3, min.segment.length = 0,
    show.legend = FALSE, max.overlaps = Inf, seed = 1
  )
}

fig7 <- ggplot2::ggplot(fig7_long, ggplot2::aes(x = score, y = z, group = model)) +
  ggplot2::geom_line(data = fig7_stayers, colour = "grey80", linewidth = 0.4) +
  ggplot2::geom_point(data = fig7_stayers, colour = "grey75", size = 1.2) +
  ggplot2::geom_line(data = fig7_movers,
    ggplot2::aes(colour = model), linewidth = 0.7, show.legend = FALSE) +
  ggplot2::geom_point(data = fig7_movers,
    ggplot2::aes(colour = model), size = 2, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = fig7_palette) +
  fig7_label_layer("PS",  -0.12, c(NA, 0.93)) +
  fig7_label_layer("WPS",  0.12, c(2.07, NA)) +
  ggplot2::facet_wrap(
    ~ cost_fn, nrow = 2, axes = "all_x",
    labeller = ggplot2::labeller(
      cost_fn = c(peak = "Seasonal Peak", growth = "High Growth"))) +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 1.0)) +
  ggplot2::labs(
    title = sprintf("%s, %s Season", FIG7_LOCATION, FIG7_SEASON),
    x = NULL, y = "Standardized Score (Lower is Better)"
  ) +
  figure_theme() +
  ggplot2::theme(
    axis.text     = ggplot2::element_text(size = 16, colour = "grey10"),
    axis.title    = ggplot2::element_text(size = 17, colour = "black", face = "plain"),
    strip.text    = ggplot2::element_text(size = 16, colour = "black", face = "plain"),
    plot.title    = ggplot2::element_text(size = 16, face = "plain", colour = "black",
                                          margin = ggplot2::margin(b = 2)),
    panel.spacing = grid::unit(1.2, "lines"),
    plot.margin   = ggplot2::margin(10, 20, 10, 20)
  ) +
  ggplot2::coord_cartesian(clip = "off")   # let labels render in the gutters

ggplot2::ggsave(
  file.path(output_dir, "figure_7_rank_slopegraph.png"),
  plot = fig7, width = 11, height = 10, dpi = 300
)
message("Saved: figure_7_rank_slopegraph.png")

# ---- how many models actually moved -----------------------------------------
fig7_summary <- fig7_scores |>
  dplyr::group_by(cost_fn) |>
  dplyr::summarise(
    n_models      = dplyr::n(),
    n_moved       = sum(abs(WPS_rank - PS_rank) >= FIG7_LABEL_SHIFT),
    median_shift  = median(abs(WPS_rank - PS_rank)),
    max_shift     = max(abs(WPS_rank - PS_rank)),
    .groups = "drop"
  )

message(sprintf("\nFig 7 — models whose rank moved >= %d places (%s %s, h = %d, R = %g):",
                FIG7_LABEL_SHIFT, FIG7_LOCATION, FIG7_SEASON, FIG7_HORIZON, FIG7_R))
for (i in seq_len(nrow(fig7_summary))) {
  r <- fig7_summary[i, ]
  message(sprintf("    %-8s %2d of %d models   (median |shift| = %g, max = %g)",
                  r$cost_fn, r$n_moved, r$n_models, r$median_shift, r$max_shift))
}

message("\n=== All figures and table saved to final/outputs/ ===")
