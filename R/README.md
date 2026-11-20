# R/

This directory contains the R functions used for data preparation, statistical analyses, and generation of tables and figures.

-   01-process-data.R applies eligibility and exclusion criteria to the processed study data and performs multiple imputation of left-censored HAI titers.

-   02-model-fits.R defines model covariates, fits overall and subtype-specific regression models, calculates participant-clustered robust variance estimates, and processes pooled estimates.

-   03-tables.R generates descriptive summaries, sample-size tables, model specifications, and regression results.

-   04-figures.R generates descriptive HAI antibody-response plots and adjusted regression figures.

The analysis workflow is defined in `_targets.R` in the repository root. Functions in this directory are sourced automatically by the `{targets}` pipeline.

For details on the statistical methods, including eligibility criteria, multiple imputation, model specifications, and sensitivity analyses, see the manuscript and supplementary materials.
