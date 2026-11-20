# results/

## Folder content and structure

-   `models/` contains processed regression estimates from the primary and sensitivity analyses.

-   `tables/` contains descriptive summaries, sample sizes, model specifications, and regression results.

-   `figures/` contains descriptive HAI plots and adjusted regression figures.

Regression results are organized into overall and subtype-specific analyses.

The primary analyses use multiple imputation to account for left-censored HAI titers. Sensitivity analyses use either fixed substitution for censored titers or data restricted to the FluVac cohort.

Tables are generally saved as RDS files containing formatted flextable objects. Figures are saved as image files.

The manuscript and supplementary materials identify the relevant tables and figures.
