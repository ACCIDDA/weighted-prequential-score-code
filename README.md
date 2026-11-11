# Weighted prequential scoring — figures

Code to reproduce every figure and the summary table for the weighted
prequential score (WPS) analysis of FluSight influenza hospitalisation
forecasts.

## Running

```r
source("produce_all_figures.R")
```

One script produces everything. Outputs are written to `outputs/`:

| File | Content |
|---|---|
| `figure_1_intro.png` | Michigan 2023/2024, two models against the observed curve |
| `figure_2_et_weights.pdf` | Alternative / uniform / ET weight distributions |
| `figure_3_examples.png` | Alabama and Oklahoma illustrative rank changes |
| `figure_4_rcrit_cdf.png` | Estimated CDFs of the critical emphasis ratio *R*<sub>crit</sub> |
| `table_1_rcrit_summary.tex` | *R*<sub>crit</sub> summary, LaTeX booktabs |
| `figure_5_tradeoff.png` | Improvement / penalty tradeoff across *W* |
| `figure_6_weekly_score_panels.png` | Weekly contribution to the PS and WPS |
| `figure_7_rank_slopegraph.png` | PS vs WPS rank slopegraph |

## Data and reproducibility

Data comes from the CDC FluSight hub
(`cdcepi-flusight-forecast-hub`), which is read in with `hubData`.

To make figures reproducible, cached data is uploaded as well. 

## Layout

```
produce_all_figures.R       all figures and the tables
R/load_modern_data.R        loads data
R/rank_analysis.R           helper file to produce tables
outputs/                    stores figures and table
```

## Dependencies

```
arrow  dplyr  ggplot2  ggrepel  grid  here  hubData  MMWRweek
patchwork  purrr  scales  scoringutils  tibble  tidyr
```

`hubData` is from the hubverse:

```r
install.packages("hubData", repos = c("https://hubverse-org.r-universe.dev",
                                      "https://cloud.r-project.org"))
```


## Settings

Set near the top of `produce_all_figures.R`:

| Constant | Value | Meaning |
|---|---|---|
| `SEASONS` | 2023/24–2025/26 | seasons analysed |
| `HORIZONS` | 2, 3 | forecast horizons, in weeks |
| `PEAK_WINDOW` | 3 | peak event weeks = peak ± 3 |
| `GROWTH_QUANTILE` | 0.50 | growth event weeks exceed the median positive increase |
| `COVERAGE` | 0.9 | minimum share of scorable weeks a model must submit |
| `THRESHOLD` | 1.0 | a model must submit for **every** event week |
| `MIN_MODELS` | 5 | locations with fewer eligible models are dropped |
