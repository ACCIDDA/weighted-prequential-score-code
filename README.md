# Weighted Prequential Score (WPS) Manuscript Figures

Code to reproduce figures for *Weighting the weeks that matter: entropic tilting for epidemic forecast evaluation*

## Producing Figures

```r
source("produce_all_figures.R")
```

The script above produces everything. Figures are written to `outputs/` directory:

| File | Content |
|---|---|
| `figure_1_intro.png` | Michigan 2023/2024 example |
| `figure_2_et_weights.pdf` | ET weights |
| `figure_3_examples.png` | Alabama and Oklahoma examples |
| `figure_4_rcrit_cdf.png` | Estimated CDFs of the critical emphasis ratio *r*<sub>crit</sub> |
| `table_1_rcrit_summary.tex` | *r*<sub>crit</sub> summary, LaTeX booktabs |
| `figure_5_tradeoff.png` | Improvement / penalty tradeoff  |
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

Important settings at the top of `produce_all_figures.R`:

| Constant | Value | Meaning |
|---|---|---|
| `SEASONS` | 2023/24–2025/26 | seasons analysed |
| `HORIZONS` | 2, 3 | forecast horizons, in weeks |
| `PEAK_WINDOW` | 3 | peak event weeks = peak ± 3 |
| `GROWTH_QUANTILE` | 0.50 | growth event weeks exceed the median positive increase |
| `COVERAGE` | 0.9 | minimum share of scorable weeks a model must submit |
| `THRESHOLD` | 1.0 | a model must submit for **every** event week |
| `MIN_MODELS` | 5 | locations with fewer eligible models are dropped |
