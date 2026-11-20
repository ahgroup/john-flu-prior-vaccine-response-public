###
# 01-process-data.R
###

func_reduce_data <- function(datapath, filepath) {
  
  data <- readRDS(datapath)
  
  # define allowed previous vaccines per current vaccine (based on sample size)
  vaccine_rules <- list(
    "Fluzone-SD" = c("Fluzone-SD", "Fluzone-HD", "FluMist", "Flucelvax", 
                     "FluBlok", "Intradermal"),
    "Fluzone-HD" = c("Fluzone-SD", "Fluzone-HD", "Intradermal"),
    "Flucelvax" = c("Fluzone-SD", "Flucelvax", "FluBlok"),
    "FluMist" = c("Fluzone-SD", "FluMist"),
    "FluBlok" = c("Fluzone-SD", "Flucelvax", "FluBlok"),
    "FLUAD" = c("Fluzone-HD", "FLUAD")
  )
  
  allowed_combinations <- tibble::tibble(
    vaccine = rep(names(vaccine_rules), times = lengths(vaccine_rules)),
    prev_vaccine = unlist(vaccine_rules)
  )
  
  dat_combs <- data %>%
    semi_join(allowed_combinations, by = c("vaccine", "prev_vaccine")) %>%
    filter(!is.na(pretiter), !is.na(posttiter)) %>%
    mutate(
      prev_vaccine_abbr = factor(
        prev_vaccine, 
        levels = c("Fluzone-SD","Fluzone-HD","FluMist","Flucelvax","FluBlok",
                   "FLUAD","Intradermal"),
        labels = c("FZSD","FZHD","FM","FV","FB","FA","ID")
      ),
      vaccine = factor(
        vaccine, 
        levels = c("Fluzone-SD","Fluzone-HD","FluMist","Flucelvax","FluBlok",
                   "FLUAD")
      ),
      vaccine_sequence = paste0(prev_vaccine, ", ", vaccine)
    )
  
  # calculate censoring across pre/posttiter and subtype
  lod = 10
  censoring_threshold = .50
  
  dat_reduced <- dat_combs %>%
    group_by(vaccine, prev_vaccine, assay_subtype) %>%
    mutate(
      n_titer_measurements = sum(!is.na(pretiter)) + sum(!is.na(posttiter)),
      n_censored_measurements = sum(pretiter < lod, na.rm = TRUE) + sum(posttiter < lod, na.rm = TRUE),
      proportion_censored = n_censored_measurements / n_titer_measurements
    ) %>%
    # keep strata with no more than 50% censoring
    filter(proportion_censored <= censoring_threshold) %>%
    ungroup() %>%
    select(-n_titer_measurements, -n_censored_measurements, -proportion_censored) %>%
    # after censoring exclusion, drop vaccine-subtype strata which no longer have enough obs
    filter(
      !(vaccine == "FluMist" & assay_subtype == "H1N1"),
      !(vaccine == "Flucelvax" & assay_subtype == "H1N1"),
      !(vaccine == "FLUAD" & assay_subtype == "H1N1"),
      !(vaccine == "FLUAD" & assay_subtype == "B-Yamagata"), # HD no B-Yam
      !(vaccine == "FluBlok" & assay_subtype == "B-Yamagata") # only one group
    )
  
  dat_scale <- dat_reduced %>%
    group_by(vaccine) %>%
    mutate(
      age_sd_overall =
        (age - min(age, na.rm = TRUE)) /
        (max(age, na.rm = TRUE) - min(age, na.rm = TRUE)),
      log_pretiter_sd_overall =
        (log_pretiter - min(log_pretiter, na.rm = TRUE)) /
        (max(log_pretiter, na.rm = TRUE) -
           min(log_pretiter, na.rm = TRUE))
    ) %>%
    ungroup() %>%
    group_by(vaccine, assay_subtype) %>%
    mutate(
      age_sd_subtype =
        (age - min(age, na.rm = TRUE)) /
        (max(age, na.rm = TRUE) - min(age, na.rm = TRUE)),
      log_pretiter_sd_subtype =
        (log_pretiter - min(log_pretiter, na.rm = TRUE)) /
        (max(log_pretiter, na.rm = TRUE) -
           min(log_pretiter, na.rm = TRUE))
    ) %>%
    ungroup() %>%
    select(
      subject_id, person_year_id, prev_person_year_id, age, age_min, age_max,
      age_factor, age_sd_overall, age_sd_subtype, birth_year, gender, race, ethnicity, bmi, year, study, 
      study_site, total_times_participated, vaccine_sequence, vaccine, 
      prev_vaccine, prev_vaccine_abbr, prev_vaccine_2y, 
      h1n1_egg_equivalent, h1n1_strain_change, 
      h3n2_egg_equivalent, h3n2_strain_change,
      bvic_egg_equivalent, bvic_strain_change,
      byam_egg_equivalent, byam_strain_change,
      vaccine_strain, vaccine_strain_change, assay_subtype, assay_strain, 
      homologous_assay, pretiter, log_pretiter, log_pretiter_sd_overall, 
      log_pretiter_sd_subtype, posttiter, log_posttiter, fc, log_fc
    )
  
  dat_model <- dat_scale %>%
    filter(
      !is.na(subject_id),
      !is.na(prev_vaccine_abbr),
      !is.na(age),
      !is.na(log_pretiter),
      !is.na(log_posttiter)
    )
  
  # assign labels to variables
  labelled::var_label(dat_model) <- list(
    subject_id = "Subject ID",
    person_year_id = "Person-year ID",
    prev_person_year_id = "Previous person-year ID",
    year = "Year",
    study = "Study",
    study_site = "Study site",
    age = "Age",
    age_min = "Min Age",
    age_max = "Max Age",
    age_factor = "Age",
    age_sd_overall = "Age",
    age_sd_subtype = "Age",
    birth_year = "Birth year",
    gender = "Gender",
    race = "Race",
    prev_vaccine = "Previous vaccine",
    prev_vaccine_abbr = "Previous vaccine",
    prev_vaccine_2y = "3-year repeat vaccinator",
    vaccine = "Current vaccine",
    assay_subtype = "Assay subtype",
    assay_strain = "Assay strain",
    h1n1_egg_equivalent = "Id-level H1N1 standardized vaccine strain",
    h1n1_strain_change = "Id-level H1N1 vaccine strain change",
    h3n2_egg_equivalent = "H3N2 standardized vaccine strain",
    h3n2_strain_change = "H3N2 vaccine strain change",
    bvic_egg_equivalent = "B-Victoria standardized vaccine strain",
    bvic_strain_change = "B-Victoria vaccine strain change",
    byam_egg_equivalent = "B-Yamagata standardized vaccine strain",
    byam_strain_change = "B-Yamagata vaccine strain change",
    vaccine_strain = "Vaccine strain",
    vaccine_strain_change = "Vaccine strain change",
    log_pretiter = "Log pretiter",
    log_pretiter_sd_subtype = "Log pretiter",
    log_pretiter_sd_overall = "Log pretiter",
    fc = "Fold change",
    log_fc = "Log fold change",
    homologous_assay = "Assay strain is same as vaccine strain"
  )
  
  saveRDS(dat_model, paste0(filepath,".rds"))
  return(dat_model)
}

func_impute_data <- function(data, n_imputations, iterations, seed) {
  
  set.seed(seed)
  
  lod <- log2(10)
  
  dat_cens <- data %>%
    mutate(
      cens_pre = log_pretiter < lod,
      cens_post = log_posttiter < lod
    )
  
  # impute one complete dataset
  impute_one <- function(dat_cens) {
    
    dat_imputed <- dat_cens %>%
      group_split(vaccine, assay_subtype) %>%
      purrr::map_dfr(function(dat) {
        
        # identify censoring
        lod_pre <- dat$cens_pre
        lod_post <- dat$cens_post
        
        # get obs uncensoring at both timepoints
        dat_complete <- !lod_pre & !lod_post
        
        # estimate joint distn from complete cases for each current vaccine
        mu <- c(
          pre = mean(dat$log_pretiter[dat_complete], na.rm = TRUE),
          post = mean(dat$log_posttiter[dat_complete], na.rm = TRUE)
        )
        
        sigma <- cov(
          cbind(
            pre = dat$log_pretiter[dat_complete],
            post = dat$log_posttiter[dat_complete]
          ),
          use = "complete.obs"
        )
        
        # conditional distributions
        beta_pre_post <- 
          sigma["pre", "post"] / sigma["post", "post"]
        
        beta_post_pre <- 
          sigma["pre", "post"] / sigma["pre", "pre"]
        
        sd_pre_post <- sqrt(
          sigma["pre", "pre"] -
            sigma["pre", "post"]^2 / sigma["post", "post"]
        )
        
        sd_post_pre <- sqrt(
          sigma["post", "post"] -
            sigma["pre", "post"]^2 / sigma["pre", "pre"]
        )
        
        # start from observed values
        log_pre_assign <- dat$log_pretiter
        log_post_assign <- dat$log_posttiter
        
        # Gibbs sampling
        for (iter in seq_len(iterations)) {
          
          # impute censored pre-titers conditional on post-titer
          for (i in which(lod_pre)) {
            
            mu_pre_cond <- mu["pre"] +
              beta_pre_post *
              (log_post_assign[i] - mu["post"])
            
            log_pre_assign[i] <- truncnorm::rtruncnorm(
              n = 1,
              a = -Inf,
              b = lod,
              mean = mu_pre_cond,
              sd = sd_pre_post
            )
          }
          
          # impute censored post-titers conditional on pre-titer
          for (i in which(lod_post)) {
            
            mu_post_cond <- mu["post"] +
              beta_post_pre *
              (log_pre_assign[i] - mu["pre"])
            
            log_post_assign[i] <- truncnorm::rtruncnorm(
              n = 1,
              a = -Inf,
              b = lod,
              mean = mu_post_cond,
              sd = sd_post_pre
            )
          }
        }
        
        dat %>%
          mutate(
            log_pretiter = log_pre_assign,
            log_posttiter = log_post_assign,
            
            pretiter = 2^log_pre_assign,
            posttiter = 2^log_post_assign,
            
            fc = posttiter / pretiter,
            log_fc = log_posttiter - log_pretiter
          )
      })
    
    # recalculate scaled pretiter after imputation
    dat_imputed %>%
      group_by(vaccine) %>%
      mutate(
        log_pretiter_sd_overall =
          (log_pretiter - min(log_pretiter, na.rm = TRUE)) /
          (max(log_pretiter, na.rm = TRUE) -
             min(log_pretiter, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      group_by(vaccine, assay_subtype) %>%
      mutate(
        log_pretiter_sd_subtype =
          (log_pretiter - min(log_pretiter, na.rm = TRUE)) /
          (max(log_pretiter, na.rm = TRUE) -
             min(log_pretiter, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      select(-cens_pre, -cens_post)
  }
  
  # generate n complete datasets
  imputed_data <- vector("list", n_imputations)
  
  for (i in seq_len(n_imputations)) {
    imputed_data[[i]] <- impute_one(dat_cens)
  }
  
  names(imputed_data) <- paste0("imp_", seq_len(n_imputations))
  
  return(imputed_data)
}
# END OF SCRIPT ====