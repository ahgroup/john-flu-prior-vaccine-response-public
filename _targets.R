###
# _targets.R
#
# Pipeline ====
#
# 00 external files
# Load raw data and reference files.
#
# 01 process data
# Process FluVac and URIV data, construct consecutive-season vaccination pairs, apply eligibility criteria, and impute censored titers.
#
# 02 fit models
# Primary analysis: overall models pooling available influenza subtypes, separate subtype-specific models.
# Sensitivity analyses: fixed substitution of censored titers at 1:5; FluVac-only models.
#
# 03 tables
# Generate descriptive, sample-size, censoring, clustering, covariate-specification, and regression-result tables.
#
# 04. figures
# Generate descriptive HAI plots and adjusted regression figures.
#
# Reproducibility ====
#
# Restore the project's R package environment first using: renv::restore()
# 
# Run the pipeline from the project root using: targets::tar_make()
#
# Target dependencies determine execution order and which results must be rebuilt after changes to inputs or functions.
###

# Setup ====

library(targets)
library(tarchetypes)

# number of cores to use
n_local_cores <- parallelly::availableCores(omit=2, constraints="connections")

# define local workers, store log for each worker
controller_local <- crew::crew_controller_local(
  name="local",
  workers=n_local_cores,
  options_local=crew::crew_options_local(log_directory="logs")
)

# set target options
tar_option_set(
  format="rds",
  controller=controller_local,
  resources=targets::tar_resources(
    crew=tar_resources_crew(controller="local")
  ),
  packages=c("dplyr")
)

# source functions
tar_source(files=list.files("R", pattern="\\.R$", full.names=TRUE, 
                              recursive=FALSE))

# Targets ====

list(
  
  ###
  ## 00 external files ====
  ###
  
  tar_target(
    file_raw_data,
    here::here("data/raw/dat_process.rds"),
    format = "file"
  ),
  
  tar_target(
    file_homologous_assays,
    here::here("data/reference/homologous_assay_list.xlsx"),
    format = "file"
  ),
  
  tar_target(
    file_vaccine_strains,
    here::here("data/reference/vaccine_list.xlsx"),
    format = "file"
  ),
  
  tar_target(
    file_vaccine_types,
    here::here("data/reference/vaccine_type_table.xlsx"),
    format = "file"
  ),
  
  ###
  ## 01 process data ====
  ###
  
  tar_target(
    reduce_data,
    func_reduce_data(
      datapath=file_raw_data,
      filepath=here::here("data/final/dat_final")
    )
  ),
  
  tar_target(
    impute_data,
    func_impute_data(
      data=reduce_data, 
      n_imputations=40, 
      iterations=100, 
      seed=333
    )
  ),
  
  tar_target(
    impute_data_fluvac,
    func_impute_data(
      data=reduce_data %>% filter(study == "FluVac"), 
      n_imputations=40, 
      iterations=100, 
      seed=333
    )
  ),
  
  ###
  ## 02 fit models ====
  ###
  
  ### primary ====
  
  #### overall ====
  
  tar_target(
    model_terms_primary_overall,
    func_model_terms_overall()
  ),
  
  tar_target(
    fit_models_primary_overall,
    func_fit_models_imputed_overall(
      imputed_data=impute_data,
      model_terms=model_terms_primary_overall
    )
  ),
  
  tar_target(
    pool_models_primary_overall,
    func_pool_models_overall(
      models=fit_models_primary_overall
    )
  ),
  
  tar_target(
    estimates_primary_overall,
    func_pool_estimates_overall(
      models=pool_models_primary_overall,
      filepath=here::here("results/models/primary/overall/estimates-primary-overall.rds")
    )
  ),
  
  #### subtype ====
  
  tar_target(
    model_terms_primary_subtype,
    func_model_terms_subtype()
  ),
  
  tar_target(
    fit_models_primary_subtype,
    func_fit_models_imputed_subtype(
      imputed_data=impute_data,
      model_terms=model_terms_primary_subtype
    )
  ),
  
  tar_target(
    pool_models_primary_subtype,
    func_pool_models_subtype(
      models=fit_models_primary_subtype
      )
  ),
  
  tar_target(
    estimates_primary_subtype,
    func_pool_estimates_subtype(
      models=pool_models_primary_subtype,
      filepath=here::here("results/models/primary/subtype/estimates-primary-subtype.rds")
    )
  ),
  
  ### sensitivity ====
  
  #### fixed lod ====
  
  ##### overall
  
  tar_target(
    fit_models_fixed_lod_overall,
    func_fit_models_overall(
      data=reduce_data,
      model_terms=model_terms_primary_overall
    )
  ),
  
  tar_target(
    estimates_fixed_lod_overall,
    func_process_estimates_fixed_lod(
      models=fit_models_fixed_lod_overall,
      by_subtype=FALSE,
      filepath=here::here("results/models/sensitivity/fixed-lod/overall/estimates-fixed-lod-overall.rds")
    )
  ),
  
  ##### subtype
  
  tar_target(
    fit_models_fixed_lod_subtype,
    func_fit_models_subtype(
      data=reduce_data,
      model_terms=model_terms_primary_subtype
    )
  ),
  
  tar_target(
    estimates_fixed_lod_subtype,
    func_process_estimates_fixed_lod(
      models=fit_models_fixed_lod_subtype,
      filepath=here::here("results/models/sensitivity/fixed-lod/subtype/estimates-fixed-lod-subtype.rds")
    )
  ),
  
  #### FluVac only ====
  
  ##### overall
  
  tar_target(
    model_terms_fluvac_overall,
    func_model_terms_overall_fluvac()
  ),
  
  tar_target(
    fit_models_fluvac_overall,
    func_fit_models_imputed_overall(
      imputed_data=impute_data_fluvac,
      model_terms=model_terms_fluvac_overall
    )
  ),
  
  tar_target(
    pool_models_fluvac_overall,
    func_pool_models_overall(
      models=fit_models_fluvac_overall
    )
  ),
  
  tar_target(
    estimates_fluvac_overall,
    func_pool_estimates_overall(
      models=pool_models_fluvac_overall,
      filepath=here::here("results/models/sensitivity/fluvac-only/overall/estimates-fluvac-overall.rds")
    )
  ),
  
  ##### subtype
  
  tar_target(
    model_terms_fluvac_subtype,
    func_model_terms_subtype_fluvac()
  ),
  
  tar_target(
    fit_models_fluvac_subtype,
    func_fit_models_imputed_subtype(
      imputed_data=impute_data_fluvac,
      model_terms=model_terms_fluvac_subtype
    )
  ),
  
  tar_target(
    pool_models_fluvac_subtype,
    func_pool_models_subtype(
      models=fit_models_fluvac_subtype
    )
  ),
  
  tar_target(
    estimates_fluvac_subtype,
    func_pool_estimates_subtype(
      models=pool_models_fluvac_subtype,
      filepath=here::here("results/models/sensitivity/fluvac-only/subtype/estimates-fluvac-subtype.rds")
    )
  ),
  
  ###
  ## 03 tables ====
  ###
  
  ### descriptive ====
  
  tar_target(
    tab_vaccine_descriptions,
    func_tab_vaccine_descriptions(
      vaccine_description_path=file_vaccine_types,
      filepath=here::here("results/tables/descriptive/vaccine-formulations.rds")
    )
  ),
  
  tar_target(
    tab_homologous_strains,
    func_tab_homologous_strains(
      homologous_assay_path=here::here("data/reference/homologous_assay_list.xlsx"),
      filepath=here::here("results/tables/descriptive/homologous-strains.rds")
    )
  ),
  
  tar_target(
    tab_demographics,
    func_tab_demographics(
      data=reduce_data,
      filepath=here::here("results/tables/descriptive/demographics.rds")
    )
  ),
  
  tar_target(
    tab_n_obs_all,
    func_tab_n_obs(
      data=reduce_data,
      stratify_by="all",
      filepath=here::here("results/tables/descriptive/eligible-observations/all/eligible-observations")
    )
  ),
  
  tar_target(
    tab_n_obs_study,
    func_tab_n_obs(
      data=reduce_data,
      stratify_by="study",
      filepath=here::here("results/tables/descriptive/eligible-observations/by-study/eligible-observations")
    )
  ),
  
  tar_target(
    tab_n_obs_site,
    func_tab_n_obs(
      data=reduce_data,
      stratify_by="site",
      filepath=here::here("results/tables/descriptive/eligible-observations/by-site/eligible-observations")
    )
  ),
  
  tar_target(
    tab_n_obs_subtype,
    func_tab_n_obs_subtype(
      data=reduce_data,
      filepath=here::here("results/tables/descriptive/eligible-observations/by-subtype/eligible-observations")
    )
  ),
  
  tar_target(
    tab_n_clusters_overall,
    func_tab_n_clusters(
      data=reduce_data,
      by_subtype=FALSE,
      filepath=here::here("results/tables/descriptive/cluster-counts/cluster-counts-overall.rds")
    )
  ),
  
  tar_target(
    tab_n_clusters_subtype,
    func_tab_n_clusters(
      data=reduce_data,
      by_subtype=TRUE,
      filepath=here::here("results/tables/descriptive/cluster-counts/cluster-counts-subtype.rds")
    )
  ),
  
  tar_target(
    tab_fc_small,
    func_tab_fc_small(
      data=reduce_data,
      filepath=here::here("results/tables/descriptive/small-fold-changes/small-fc-")
    )
  ),
  
  tar_target(
    tab_censored_counts,
    func_tab_censored_counts(
      data=reduce_data,
      filepath=here::here("results/tables/descriptive/censored-counts.rds")
    )
  ),
  
  tar_target(
    tab_strain_combs,
    func_tab_strain_combs(
      data=reduce_data,
      filepath=here::here("results/tables/descriptive/strain-combinations/combs-")
    )
  ),
  
  ### model spec ====
  
  tar_target(
    tab_model_terms_overall,
    func_tab_model_terms(
      model_terms=model_terms_primary_overall,
      model_type="overall",
      filepath=here::here("results/tables/model-specifications/model-terms-overall.rds")
    )
  ),
  
  tar_target(
    tab_model_terms_subtype,
    func_tab_model_terms(
      model_terms=model_terms_primary_subtype,
      model_type="subtype",
      filepath=here::here("results/tables/model-specifications/model-terms-subtype.rds")
    )
  ),
  
  ### regression ====
  
  #### primary 
  
  tar_target(
    tab_regression_primary_overall,
    func_tab_regression(
      estimates=estimates_primary_overall,
      filepath=here::here("results/tables/regression/primary/overall/regression-overall-")
    )
  ),
  
  tar_target(
    tab_regression_primary_subtype,
    func_tab_regression(
      estimates=estimates_primary_subtype,
      filepath=here::here("results/tables/regression/primary/subtype/regression-primary-")
    )
  ),
  
  #### fixed lod
  
  tar_target(
    tab_regression_fixed_lod_overall,
    func_tab_regression(
      estimates=estimates_fixed_lod_overall,
      filepath=here::here("results/tables/regression/sensitivity/fixed-lod/overall/regression-fixed-lod-overall-")
    )
  ),
  
  tar_target(
    tab_regression_fixed_lod_subtype,
    func_tab_regression(
      estimates=estimates_fixed_lod_subtype,
      filepath=here::here("results/tables/regression/sensitivity/fixed-lod/subtype/regression-fixed-lod-subtype-")
    )
  ),
  
  #### fluvac
  
  tar_target(
    tab_regression_fluvac_overall,
    func_tab_regression(
      estimates=estimates_fluvac_overall,
      filepath=here::here("results/tables/regression/sensitivity/fluvac-only/overall/regression-fluvac-overall-")
    )
  ),
  
  tar_target(
    tab_regression_fluvac_subtype,
    func_tab_regression(
      estimates=estimates_fluvac_subtype,
      filepath=here::here("results/tables/regression/sensitivity/fluvac-only/subtype/regression-fluvac-subtype-")
    )
  ),
  
  ###
  ## 04 figures ====
  ###
  
  ### descriptive ====
  
  tar_target(
    plot_current,
    func_plot_current(
      data=reduce_data,
      filepath=here::here("results/figures/descriptive/current-vaccine/plot-current.png"), 
      height=8, 
      width=10
    )
  ),
  
  tar_target(
    plot_prior_current,
    func_plot_prior_current(
      data=reduce_data,
      filepath=here::here("results/figures/descriptive/prior-current-vaccine/plot-prior-current-"),
      height=11,
      width=9
    )
  ),
  
  tar_target(
    plot_prior_current_study,
    func_plot_prior_current_by(
      data=reduce_data,
      facet="study",
      filepath=here::here("results/figures/descriptive/by-study/plot-prior-current-study-"),
      height=12,
      width=14
    )
  ),
  
  tar_target(
    plot_prior_current_site,
    func_plot_prior_current_by(
      data=reduce_data,
      facet="study_site",
      filepath=here::here("results/figures/descriptive/by-site/plot-prior-current-site-"),
      height=12,
      width=14
    )
  ),
  
  tar_target(
    plot_pre_post,
    func_plot_pre_post(
      data=reduce_data,
      filepath=here::here("results/figures/descriptive/pre-post/plot-pre-post-"),
      height=7,
      width=7
    )
  ),
  
  ### regression ====
  
  #### primary
  
  tar_target(
    plot_regression_primary_overall,
    func_plot_regression(
      estimates=estimates_primary_overall,
      filepath=here::here("results/figures/regression/primary/plot-regression-primary-overall.png"),
      height=7,
      width=3.5
    )
  ),
  
  tar_target(
    plot_regression_primary_subtype,
    func_plot_regression(
      estimates=estimates_primary_subtype,
      filepath=here::here("results/figures/regression/primary/plot-regression-primary-subtype.png"),
      height=7,
      width=8
    )
  ),
  
  #### fixed lod
  
  tar_target(
    plot_regression_fixed_lod_overall,
    func_plot_regression(
      estimates=estimates_fixed_lod_overall,
      filepath=here::here("results/figures/regression/sensitivity/fixed-lod/plot-regression-fixed-lod-overall.png"),
      height=7,
      width=3.5
    )
  ),
  
  tar_target(
    plot_regression_fixed_lod_subtype,
    func_plot_regression(
      estimates=estimates_fixed_lod_subtype,
      filepath=here::here("results/figures/regression/sensitivity/fixed-lod/plot-regression-fixed-lod-subtype.png"),
      height=7,
      width=8
    )
  ),
  
  #### fluvac
  
  tar_target(
    plot_regression_fluvac_overall,
    func_plot_regression(
      estimates=estimates_fluvac_overall,
      filepath=here::here("results/figures/regression/sensitivity/fluvac-only/plot-regression-fluvac-overall.png"),
      height=7,
      width=3.5
    )
  ),
  
  tar_target(
    plot_regression_fluvac_subtype,
    func_plot_regression(
      estimates=estimates_fluvac_subtype,
      filepath=here::here("results/figures/regression/sensitivity/fluvac-only/plot-regression-fluvac-subtype.png"),
      height=7,
      width=8
    )
  )
)
# END OF SCRIPT ====