suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(tibble)
})

nature_theme_plots <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.border       = ggplot2::element_blank(),
      axis.line.x        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.line.y        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.text          = ggplot2::element_text(size = 12, colour = "grey10"),
      axis.title         = ggplot2::element_text(size = 13, colour = "black", face = "plain"),
      axis.title.x       = ggplot2::element_blank(),
      axis.ticks         = ggplot2::element_line(colour = "black"),
      axis.ticks.length  = grid::unit(2.5, "pt"),
      legend.position    = "none",
      plot.title         = ggplot2::element_text(size = 12, face = "plain", colour = "black",
                                                  margin = ggplot2::margin(b = 2))
    )
}

build_plot_series <- function(case_df,
                              location,
                              season,
                              horizon,
                              model_name) {
  plot_df <- preprocess_rank_table_input_fast(
    modern_all = case_df,
    location_value = location,
    horizon_value = horizon,
    season_value = season
  )
  
  obs_df <- plot_df %>%
    filter(!is.na(.data$observed)) %>%
    distinct(.data$target_end_date, .keep_all = TRUE) %>%
    transmute(
      target_end_date = as.Date(.data$target_end_date),
      observed = as.numeric(.data$observed)
    ) %>%
    arrange(.data$target_end_date)
  
  model_df <- plot_df %>%
    filter(.data$model == model_name, .data$quantile_level == 0.5) %>%
    distinct(.data$target_end_date, .keep_all = TRUE) %>%
    transmute(
      target_end_date = as.Date(.data$target_end_date),
      predicted = as.numeric(.data$predicted)
    ) %>%
    arrange(.data$target_end_date)
  
  interval_df <- plot_df %>%
    filter(.data$model == model_name, .data$quantile_level %in% c(0.05, 0.95)) %>%
    transmute(
      target_end_date = as.Date(.data$target_end_date),
      quantile_level = as.numeric(.data$quantile_level),
      predicted = as.numeric(.data$predicted)
    ) %>%
    distinct(.data$target_end_date, .data$quantile_level, .keep_all = TRUE) %>%
    tidyr::pivot_wider(
      names_from = .data$quantile_level,
      values_from = .data$predicted,
      names_prefix = "q_"
    ) %>%
    transmute(
      target_end_date = .data$target_end_date,
      lower_90 = .data$q_0.05,
      upper_90 = .data$q_0.95
    ) %>%
    arrange(.data$target_end_date)
  
  list(
    observed = obs_df,
    forecast = model_df,
    interval = interval_df
  )
}

make_model_observed_panel <- function(observed_df,
                                      forecast_df,
                                      interval_df,
                                      model_name,
                                      location,
                                      season,
                                      window_start,
                                      window_end,
                                      x_limits,
                                      y_limits) {
  ggplot() +
    # annotate(
    #   "rect",
    #   xmin = as.Date(window_start),
    #   xmax = as.Date(window_end),
    #   ymin = -Inf,
    #   ymax = Inf,
    #   alpha = 0.12,
    #   fill = "#FDB863"
    # ) +
    geom_line(
      data = observed_df,
      aes(x = .data$target_end_date, y = .data$observed, color = "Observed"),
      linewidth = 0.9
    ) +
    geom_point(
      data = observed_df,
      aes(x = .data$target_end_date, y = .data$observed, color = "Observed"),
      size = 1.4
    ) +
    geom_ribbon(
      data = interval_df,
      aes(
        x = .data$target_end_date,
        ymin = .data$lower_90,
        ymax = .data$upper_90
      ),
      inherit.aes = FALSE,
      alpha = 0.20,
      fill = "#2C7FB8"
    ) +
    geom_line(
      data = forecast_df,
      aes(x = .data$target_end_date, y = .data$predicted, color = "Forecast"),
      linewidth = 0.9
    ) +
    geom_point(
      data = forecast_df,
      aes(x = .data$target_end_date, y = .data$predicted, color = "Forecast"),
      size = 1.4
    ) +
    scale_color_manual(
      values = c("Observed" = "#1F1F1F", "Forecast" = "#2C7FB8")
    ) +
    scale_x_date(
      limits = x_limits,
      date_breaks = "4 weeks",
      date_labels = "%b\n%Y"
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = model_name,
     # subtitle = paste(season, "Season"),
      x = "",
      y = "Incident Hospitalizations",
      color = NULL
    ) +
    nature_theme_plots()
}

create_example_plot_from_row <- function(example_row,
                                         case_data,
                                         save_path = NULL) {
  if (!is.data.frame(example_row) || nrow(example_row) != 1) {
    stop("example_row must be a one-row data frame or tibble.", call. = FALSE)
  }
  
  location <- example_row$location[[1]]
  season <- example_row$season[[1]]
  horizon <- example_row$horizon[[1]]
  window_start <- as.Date(example_row$window_start[[1]])
  window_end <- as.Date(example_row$window_end[[1]])
  model_a <- example_row$model_a[[1]]
  model_b <- example_row$model_b[[1]]
  example_id <- if ("example_id" %in% names(example_row)) example_row$example_id[[1]] else NA_integer_
  
  case_key <- paste(location, season, horizon, sep = ".")
  case_df <- case_data[[case_key]]
  if (is.null(case_df) || nrow(case_df) == 0) {
    stop("No case data found for this example row.", call. = FALSE)
  }
  
  series_a <- build_plot_series(case_df, location, season, horizon, model_a)
  series_b <- build_plot_series(case_df, location, season, horizon, model_b)
  
  if (nrow(series_a$observed) == 0 || nrow(series_a$forecast) == 0 || nrow(series_b$forecast) == 0) {
    stop("Insufficient observed or forecast data to build the plot.", call. = FALSE)
  }
  
  x_limits <- range(series_a$observed$target_end_date, na.rm = TRUE)
  y_limits <- range(
    c(
      series_a$observed$observed,
      series_a$forecast$predicted,
      series_b$forecast$predicted,
      series_a$interval$lower_90,
      series_a$interval$upper_90,
      series_b$interval$lower_90,
      series_b$interval$upper_90
    ),
    na.rm = TRUE
  )
  
  p_a <- make_model_observed_panel(
    observed_df = series_a$observed,
    forecast_df = series_a$forecast,
    interval_df = series_a$interval,
    model_name = model_a,
    location = location,
    season = season,
    window_start = window_start,
    window_end = window_end,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  p_b <- make_model_observed_panel(
    observed_df = series_b$observed,
    forecast_df = series_b$forecast,
    interval_df = series_b$interval,
    model_name = model_b,
    location = location,
    season = season,
    window_start = window_start,
    window_end = window_end,
    x_limits = x_limits,
    y_limits = y_limits
  )
  
  #title_text <- if (is.na(example_id)) {
  #  sprintf("%s vs %s", model_a, model_b)
  #} else {
  #  sprintf("Example %s: %s vs %s", example_id, model_a, model_b)
  #}
  
  title_text <- ""
  p <- p_a + p_b + patchwork::plot_layout(axis_titles = "collect")
  #+
  #  patchwork::plot_annotation(
  #    title = title_text,
  #    subtitle = sprintf(
  #      "%s | %s | Growth window: %s to %s",
  #      location,
  #      season,
  #      window_start,
  #      window_end
  #    )
   # )
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = 14,
      height = 5.5,
      dpi = 300
    )
  }
  
  p
}

save_top_example_plots <- function(example_pairs,
                                   case_data,
                                   output_dir) {
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  purrr::pwalk(
    example_pairs,
    function(example_id, ...) {
      example_row <- tibble::tibble(example_id = example_id, ...)
      create_example_plot_from_row(
        example_row = example_row,
        case_data = case_data,
        save_path = file.path(plot_dir, sprintf("example_%02d.png", example_id))
      )
      NULL
    }
  )
}
