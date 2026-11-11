library(dplyr)
library(scoringutils)
library(MMWRweek)

add_season <- function(d) {
  mw <- MMWRweek(as.Date(d))
  start_year <- ifelse(mw$MMWRweek >= 40, mw$MMWRyear, mw$MMWRyear - 1)
  paste0(start_year, "/", start_year + 1)
}

# Helper function to obtain MMWR format date from filename
`%||%` <- function(a, b) if (!is.null(a)) a else b

parse_filename <- function(file, horizon = 1) {
  file_base <- basename(file)
  
  # 1) Extract EW
  ew <- as.integer(stringr::str_match(file_base, "EW(\\d{2})")[,2])
  
  # If EW is missing, return a fully NA row
  if (is.na(ew)) {
    return(tibble::tibble(
      file             = file,
      EW               = NA_integer_,
      EW_year          = NA_integer_,
      date_string      = NA_character_,
      submission_date  = NA_Date_,
      target_MMWRyear  = NA_integer_,
      target_MMWRweek  = NA_integer_
    ))
  }
  
  # 2) Extract explicit EW-year if present (e.g., EW52-2020-...)
  ew_year_raw <- as.integer(stringr::str_match(file_base, "EW\\d{2}-(\\d{4})-")[,2])
  
  # 3) Extract submission date if present
  date_string <- stringr::str_extract(file_base, "\\d{4}-\\d{2}-\\d{2}")
  submission_date <- suppressWarnings(lubridate::ymd(date_string))
  
  # 4) Infer EW_year
  if (!is.na(ew_year_raw)) {
    ew_year <- ew_year_raw
  } else if (!is.na(submission_date)) {
    issue_mmw <- MMWRweek::MMWRweek(submission_date)
    ew_year <- if (ew > issue_mmw$MMWRweek) issue_mmw$MMWRyear - 1 else issue_mmw$MMWRyear
  } else {
    # No EW-year and no submission date -> cannot determine
    ew_year <- NA_integer_
  }
  
  # 5) Compute target MMWR week/year
  target_week <- ew + horizon
  target_year <- ew_year
  
  if (!is.na(target_week) && target_week > 52) {
    target_week <- target_week - 52
    target_year <- target_year + 1
  }
  
  tibble::tibble(
    file             = file,
    EW               = ew,
    EW_year          = ew_year,
    date_string      = date_string,
    submission_date  = submission_date,
    target_MMWRyear  = target_year,
    target_MMWRweek  = target_week
  )
}



alphas <- c(
  0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4,
  0.5,
  0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99
)
