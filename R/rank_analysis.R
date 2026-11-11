library(dplyr)
library(scoringutils)
library(MMWRweek)

weighted_score <- function(df, W, score_col = "wis", E_col = "E") {
  E_vals <- df[[E_col]]
  S_vals <- df[[score_col]]
  
  n_E  <- sum(E_vals == 1, na.rm = TRUE)
  n_Ec <- sum(E_vals == 0, na.rm = TRUE)
  
  if (n_E == 0 || n_Ec == 0) {
    return(mean(S_vals, na.rm = TRUE))
  }
  
  mean_E  <- mean(S_vals[E_vals == 1], na.rm = TRUE)
  mean_Ec <- mean(S_vals[E_vals == 0], na.rm = TRUE)
  
  W * mean_E + (1 - W) * mean_Ec
}

preprocess_rank_table_input <- function(modern_all,
                                        location_value,
                                        horizon_value,
                                        season_value) {
  modern_sub <- modern_all %>%
    filter(location == location_value, horizon == horizon_value) %>%
    filter(!is.na(predicted), !is.na(quantile_level), !is.na(forecast_date), !is.na(target_end_date)) %>%
    mutate(season = add_season(target_end_date)) %>%
    filter(season == season_value)
  
  if (season_value %in% c("2023/2024", "2024/2025")) {
    modern_sub <- modern_sub %>%
      filter(target == "wk inc flu hosp")
  }
  
  modern_sub %>%
    mutate(
      target = paste0("horizon_", horizon, "_rank_table_target")
    ) %>%
    arrange(
      model,
      forecast_date,
      target_end_date,
      quantile_level,
      desc(!is.na(observed))
    ) %>%
    distinct(
      model,
      forecast_date,
      location,
      target,
      target_end_date,
      horizon,
      quantile_level,
      .keep_all = TRUE
    )
}

keep_symmetric_quantile_rows <- function(modern_sub) {
  modern_sub %>%
    mutate(
      quantile_level_round = round(.data$quantile_level, 6)
    ) %>%
    group_by(
      .data$model,
      .data$forecast_date,
      .data$location,
      .data$target,
      .data$target_end_date,
      .data$horizon
    ) %>%
    group_modify(function(.x, .y) {
      q <- sort(unique(.x$quantile_level_round))
      valid_q <- q[q == 0.5]
      lower_q <- q[q < 0.5]
      paired_lower <- lower_q[round(1 - lower_q, 6) %in% q]
      valid_q <- sort(unique(c(valid_q, paired_lower, round(1 - paired_lower, 6))))
      .x %>%
        filter(.data$quantile_level_round %in% valid_q)
    }) %>%
    ungroup() %>%
    dplyr::select(-.data$quantile_level_round)
}

preprocess_rank_table_input_fast <- function(modern_all,
                                             location_value,
                                             horizon_value,
                                             season_value) {
  preprocess_rank_table_input(
    modern_all = modern_all,
    location_value = location_value,
    horizon_value = horizon_value,
    season_value = season_value
  ) %>%
    filter(.data$quantile_level >= 0, .data$quantile_level <= 1) %>%
    keep_symmetric_quantile_rows() %>%
    arrange(
      .data$model,
      .data$forecast_date,
      .data$location,
      .data$target,
      .data$target_end_date,
      .data$quantile_level
    ) %>%
    distinct(
      .data$model,
      .data$forecast_date,
      .data$location,
      .data$target,
      .data$target_end_date,
      .data$horizon,
      .data$quantile_level,
      .keep_all = TRUE
    )
}

build_rank_table_fast_engine <- function(modern_all,
                                         location_value,
                                         horizon_value,
                                         season_value) {
  modern_sub <- preprocess_rank_table_input_fast(
    modern_all = modern_all,
    location_value = location_value,
    horizon_value = horizon_value,
    season_value = season_value
  )
  
  if (nrow(modern_sub) == 0) {
    return(list(
      modern_sub = modern_sub,
      scores_obs = tibble::tibble()
    ))
  }
  
  fc <- as_forecast_quantile(
    modern_sub,
    model_col      = "model",
    date_col       = "forecast_date",
    target_col     = "target",
    quantile_col   = "quantile_level",
    prediction_col = "predicted",
    truth_col      = "observed",
    location_col   = "location"
  )
  
  scores <- score(
    fc,
    metrics = get_metrics(fc, select = c("wis", "overprediction", "underprediction", "dispersion")),
    na.rm = TRUE
  )
  
  needed_cols <- c("wis", "overprediction", "underprediction", "dispersion")
  if (!all(needed_cols %in% names(scores))) {
    return(list(
      modern_sub = modern_sub,
      scores_obs = tibble::tibble()
    ))
  }
  
  scores_obs <- scores %>%
    left_join(
      modern_sub %>%
        distinct(
          .data$model, .data$forecast_date, .data$location, .data$target, .data$target_end_date,
          .data$horizon, .data$output_type, .data$season, .data$observed
        ),
      by = c(
        "model", "forecast_date", "location", "target", "target_end_date",
        "horizon", "output_type", "season"
      )
    )
  
  list(
    modern_sub = modern_sub,
    scores_obs = scores_obs
  )
}

make_rank_table_fast <- function(modern_all = NULL,
                                 location_value = NULL,
                                 horizon_value = NULL,
                                 season_value = NULL,
                                 target_end_dates,
                                 W = 0.7,
                                 threshold = 0.8,
                                 coverage = 0.8,
                                 engine = NULL) {
  if (is.null(engine)) {
    engine <- build_rank_table_fast_engine(
      modern_all = modern_all,
      location_value = location_value,
      horizon_value = horizon_value,
      season_value = season_value
    )
  }
  
  scores_obs <- engine$scores_obs
  target_end_dates <- unique(as.Date(target_end_dates))
  n_flagged_dates <- length(target_end_dates)
  
  if (nrow(scores_obs) == 0) {
    return(tibble::tibble(
      model = character(),
      weighted_wis = numeric(),
      usual_wis = numeric(),
      weighted_rank = integer(),
      usual_rank = integer(),
      rank_change = integer()
    ))
  }
  
  scores_obs <- scores_obs %>%
    mutate(E = if_else(.data$target_end_date %in% target_end_dates, 1L, 0L))
  
  n_all_dates <- n_distinct(scores_obs$target_end_date)
  
  eligible_models <- scores_obs %>%
    group_by(.data$model) %>%
    summarise(
      n_dates_present = n_distinct(.data$target_end_date),
      overall_coverage = if (n_all_dates == 0) 1 else n_dates_present / n_all_dates,
      n_flagged_present = n_distinct(.data$target_end_date[.data$E == 1]),
      event_coverage = if (n_flagged_dates == 0) 1 else n_flagged_present / n_flagged_dates,
      .groups = "drop"
    ) %>%
    filter(
      .data$n_dates_present   >= floor(coverage  * n_all_dates),
      .data$n_flagged_present >= floor(threshold * n_flagged_dates)
    ) %>%
    pull(.data$model)
  
  scores_obs %>%
    filter(.data$model %in% eligible_models) %>%
    group_by(.data$model) %>%
    summarise(
      weighted_wis = weighted_score(pick(everything()), W = W),
      usual_wis    = mean(.data$wis, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_rank = min_rank(.data$weighted_wis),
      usual_rank    = min_rank(.data$usual_wis),
      rank_change   = .data$usual_rank - .data$weighted_rank
    ) %>%
    arrange(.data$weighted_rank, .data$weighted_wis)
}

make_rank_table <- function(modern_all,
                            location_value,
                            horizon_value,
                            season_value,
                            target_end_dates,
                            W = 0.7,
                            threshold = 0.8,
                            coverage = 0.8) {
  
  target_end_dates <- unique(as.Date(target_end_dates))
  n_flagged_dates <- length(target_end_dates)
  
  modern_sub <- preprocess_rank_table_input(
    modern_all = modern_all,
    location_value = location_value,
    horizon_value = horizon_value,
    season_value = season_value
  )
  
  fc <- as_forecast_quantile(
    modern_sub,
    model_col      = "model",
    date_col       = "forecast_date",
    target_col     = "target",
    quantile_col   = "quantile_level",
    prediction_col = "predicted",
    truth_col      = "observed",
    location_col   = "location"
  )
  
  scores <- score(
    fc,
    metrics = get_metrics(fc, select = c("wis", "overprediction", "underprediction", "dispersion"))
  )
  
  scores_obs <- scores %>%
    left_join(
      modern_sub %>%
        distinct(
          model, forecast_date, location, target, target_end_date,
          horizon, output_type, season, observed
        ),
      by = c(
        "model", "forecast_date", "location", "target", "target_end_date",
        "horizon", "output_type", "season"
      )
    ) %>%
    mutate(E = if_else(target_end_date %in% target_end_dates, 1L, 0L))
  
  n_all_dates <- n_distinct(scores_obs$target_end_date)
  
  eligible_models <- scores_obs %>%
    group_by(model) %>%
    summarise(
      n_dates_present = n_distinct(target_end_date),
      overall_coverage = if (n_all_dates == 0) 1 else n_dates_present / n_all_dates,
      n_flagged_present = n_distinct(target_end_date[E == 1]),
      event_coverage = if (n_flagged_dates == 0) 1 else n_flagged_present / n_flagged_dates,
      .groups = "drop"
    ) %>%
    filter(overall_coverage >= coverage, event_coverage >= threshold) %>%
    pull(model)
  
  scores_obs %>%
    filter(model %in% eligible_models) %>%
    group_by(model) %>%
    summarise(
      weighted_wis = weighted_score(pick(everything()), W = W),
      usual_wis    = mean(wis, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_rank = min_rank(weighted_wis),
      usual_rank    = min_rank(usual_wis),
      rank_change   = usual_rank - weighted_rank
    ) %>%
    arrange(weighted_rank, weighted_wis)
}
