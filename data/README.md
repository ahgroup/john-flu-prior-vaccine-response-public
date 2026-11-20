# /data

This directory contains the datasets and reference files used in the analyses.

## Folder structure and contents

-   `raw/` contains study data used as input to the analysis pipeline.

-   `final/` contains regression-eligible analysis dataset generated after applying eligibility and exclusion criteria.

-   `reference/` contains reference files describing influenza vaccine formulations and homologous HAI assay strains.

The input data combine observations from the FluVac cohort and URIV clinical trial.

Each eligible vaccination pair represents vaccinations received by the same participant in consecutive influenza seasons, with the HAI response measured following the second vaccination. A vaccination pair may contribute observations for multiple influenza subtypes.

The analysis dataset is generated using func_reduce_data(), which applies vaccine-sequence and censoring exclusions and prepares the regression variables.

Censored pre- and post-vaccination HAI titers are subsequently imputed using func_impute_data(). The imputed datasets are managed by the {targets} pipeline.

Further details on eligibility criteria, censoring, and multiple imputation are provided in the manuscript and supplementary materials.
