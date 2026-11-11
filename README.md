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

Forecasts and truth come from the CDC FluSight hub on S3
(`cdcepi-flusight-forecast-hub`), read with `hubData`.

**The hub is live.** Truth is revised weekly — between two consecutive
vintages we observed 78 changed observations for the 2025/2026 season alone —
and late or corrected submissions appear retroactively. An unpinned pull
therefore gives different numbers run to run, and the effects compound: a
one-admission revision can move the seasonal peak by a week, which shifts the
event window, which changes which models meet the all-event-weeks eligibility
rule, which changes *n*, the rankings and *R*<sub>crit</sub>.

To make the figures reproducible, `load_modern_quantile()` caches its joined
pull to `modern_quantile_<seasons>_h<horizons>.rds` and reuses it on every
later run. **That file is committed on purpose.** Delete it to refresh against
current hub data, accepting that the numbers will move.

## Layout

```
produce_all_figures.R       all figures and the table, in order
R/load_modern_data.R        hub pull, truth join, cache
R/rank_analysis.R           WIS scoring engine and rank tables
outputs/                    generated figures and table
```

`R/rank_common.R` and `R/rank_example_plots.R` are not sourced by
`produce_all_figures.R` and are retained only for reference.

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

A `renv.lock` is included. Note it was generated from a system library rather
than an active renv project, so treat it as a record of versions rather than a
guaranteed-restorable environment.

## Key parameters

Set near the top of `produce_all_figures.R`:

| Constant | Value | Meaning |
|---|---|---|
| `SEASONS` | 2023/24–2025/26 | seasons analysed |
| `HORIZONS` | 2, 3 | forecast horizons, in weeks |
| `PEAK_WINDOW` | 3 | peak event weeks = peak ± 3 |
| `GROWTH_QUANTILE` | 0.50 | growth event weeks exceed the median positive weekly rise |
| `COVERAGE` | 0.9 | minimum share of scorable weeks a model must submit |
| `THRESHOLD` | 1.0 | a model must submit for **every** event week |
| `MIN_MODELS` | 5 | locations with fewer eligible models are dropped |

Table 1, Fig 4 and Fig 5 additionally require a model to be eligible under
*both* event definitions, so the eligible location set is identical across the
peak and growth rows.
