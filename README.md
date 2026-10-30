# Predicting ESG controversies - thesis code

R code for an MSc thesis (RSM) testing whether ESG scores, pillar scores and individual ESG indicators predict future ESG controversies.

## Scripts (run in this order)

| Script | What it does |
|---|---|
| `data_preparation.R` | LSEG Excel export → model-ready datasets (sample, indicator selection, imputation) |
| `appendix_indicator.R` | Appendix tables as LaTeX: the ESG indicators used (Table B.1) and how 65 raw indicators become the 27 used in the models (Appendix C) |
| `hypothesis_testing.R` | Hypothesis tests and robustness checks (logistic regression, elastic net, random forest, gradient boosting) |

`hypothesis_testing.R` is run one section at a time in RStudio (SETUP, then H1 to H5b, then the robustness checks R1 to R5), not sourced top to bottom. Comments inside the scripts still use their earlier file names (01_data_preparation_revised.R, 02/03_appendix_..., 04_hypothesis_tests_revised.R).

## Data

Not included. The input is an LSEG Workspace export (`Data_thesis.xlsx`, 20 March 2026), which is licensed and cannot be shared. Generated files (`.rds`, `.csv`, `.tex`) are excluded too; some contain firm-level LSEG data.

## Running

Each script sets its working directory with `setwd("~/Desktop/Thesis_ESG_data/Data")`; change this to the folder holding the data.

Requirements: R 4.4.2 with readxl, dplyr, tidyr, purrr, tibble (data preparation) and glmnet, ranger, xgboost, sandwich, lmtest (hypothesis tests); dplyr is also needed for R5, and fastshap only for the optional random-forest SHAP.
