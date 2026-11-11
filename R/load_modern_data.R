# load_modern_data.R
#
# Loads FluSight quantile forecasts and truth directly from the FluSight S3
# hub using hubData. No local clones required.

suppressPackageStartupMessages({
  library(dplyr)
  library(hubData)
  library(arrow)
  library(MMWRweek)
})

HUB_BUCKET <- "cdcepi-flusight-forecast-hub"

load_modern_quantile <- function(
    hub_bucket_name = HUB_BUCKET,
    horizons        = 1:4,
    seasons         = c("2023/2024", "2024/2025", "2025/2026")
) {
  # The hub is live: truth is revised weekly and late submissions appear
  # retroactively, so an unpinned pull gives different numbers run to run.
  # First run writes this cache; later runs reuse it. Delete the file to refresh.
  cache_key <- paste0(paste(gsub("/", "-", seasons), collapse = "_"),
                      "_h", paste(horizons, collapse = ""))
  cache_dir <- if (dir.exists(here::here("R"))) {
    here::here()
  } else here::here("figure_scripts_fixed")
  cache <- file.path(cache_dir, paste0("modern_quantile_", cache_key, ".rds"))
  if (file.exists(cache)) {
    message("Using cached hub pull: ", cache)
    return(readRDS(cache))
  }

  message("Connecting to FluSight hub (S3: ", hub_bucket_name, ") ...")
  bucket <- arrow::s3_bucket(hub_bucket_name)
  hub    <- hubData::connect_hub(bucket, file_format = "parquet", skip_checks = TRUE)

  message("Loading quantile forecasts ...")
  forecasts_raw <- hub |>
    dplyr::filter(
      output_type == "quantile",
      target      == "wk inc flu hosp",
      horizon     %in% horizons
    ) |>
    hubData::collect_hub()

  forecasts <- forecasts_raw |>
    dplyr::rename(
      model          = model_id,
      forecast_date  = reference_date,
      quantile_level = output_type_id,
      predicted      = value
    ) |>
    dplyr::mutate(
      forecast_date   = as.Date(forecast_date),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(quantile_level),
      output_type     = "quantile",
      season          = add_season(target_end_date)
    ) |>
    dplyr::filter(
      !is.na(quantile_level),
      quantile_level >= 0,
      quantile_level <= 1,
      !is.na(target_end_date),
      season %in% seasons
    )

  # ---- Truth ------------------------------------------------------------------
  message("Loading truth data from S3 ...")
  truth_path <- bucket$path("target-data/target-hospital-admissions.csv")
  truth_raw  <- arrow::read_csv_arrow(truth_path, as_data_frame = TRUE)

  truth <- truth_raw |>
    dplyr::select(location, date, value) |>
    dplyr::rename(target_end_date = date, observed = value) |>
    dplyr::mutate(
      target_end_date = as.Date(target_end_date),
      location        = trimws(location)
    ) |>
    dplyr::filter(!is.na(observed))

  # ---- Join -------------------------------------------------------------------
  message("Joining forecasts with truth ...")
  out <- forecasts |>
    dplyr::left_join(truth, by = c("location", "target_end_date")) |>
    dplyr::mutate(location = trimws(location)) |>
    dplyr::select(
      model, forecast_date, location, target, target_end_date,
      horizon, output_type, season, quantile_level, predicted, observed
    )

  message(sprintf("Loaded %d rows | %d models | %d locations | seasons: %s",
                  nrow(out),
                  dplyr::n_distinct(out$model),
                  dplyr::n_distinct(out$location),
                  paste(seasons, collapse = ", ")))
  saveRDS(out, cache)
  out
}

# Full historical truth — used by Fig 2 for the 2022/2023 US epidemic curve
load_hub_truth <- function(hub_bucket_name = HUB_BUCKET) {
  bucket     <- arrow::s3_bucket(hub_bucket_name)
  truth_path <- bucket$path("target-data/target-hospital-admissions.csv")
  arrow::read_csv_arrow(truth_path, as_data_frame = TRUE) |>
    dplyr::select(location, date, value) |>
    dplyr::rename(target_end_date = date, observed = value) |>
    dplyr::mutate(
      target_end_date = as.Date(target_end_date),
      location        = trimws(location)
    ) |>
    dplyr::filter(!is.na(observed))
}
