###
# 02-model-fits.R
###

# Model specifications ====
func_model_terms_overall <- function() {
  
  model_terms <- list(
    "Fluzone-SD" = c("study", "vaccine_strain", "vaccine_strain_change"),
    "Fluzone-HD" =c("study", "vaccine_strain", "vaccine_strain_change"),
    "FluMist" = c("vaccine_strain", "vaccine_strain_change"),
    "Flucelvax" = c("study", "vaccine_strain_change"),
    "FluBlok" = c("study", "vaccine_strain_change"), 
    "FLUAD" = character(0)
  )
  
  return(model_terms)
  
}

func_model_terms_subtype <- function() {
  
  model_terms <- list(
    "Fluzone-SD" = list(
      H1N1 = c("study", "vaccine_strain", "vaccine_strain_change"),
      H3N2 = c("study", "vaccine_strain", "vaccine_strain_change"),
      `B-Victoria` = c("study", "vaccine_strain", "vaccine_strain_change"),
      `B-Yamagata` = c("study", "vaccine_strain_change")
    ),
    "Fluzone-HD" = list(
      H1N1 = c("study", "vaccine_strain", "vaccine_strain_change"),
      H3N2 = c("study", "vaccine_strain", "vaccine_strain_change"),
      `B-Victoria` = c("study", "vaccine_strain", "vaccine_strain_change"),
      `B-Yamagata` = c("vaccine_strain")
    ),
    "FluMist" = list(
      # not enough complete H1N1 data,
      H3N2 = c("vaccine_strain_change"),
      `B-Victoria` = c("vaccine_strain"),
      `B-Yamagata` = character(0)
    ),
    "Flucelvax" = list(
      # not enough complete H1N1 data,
      H3N2 = c("study", "vaccine_strain_change"),
      `B-Victoria` = c("study", "vaccine_strain_change"),
      `B-Yamagata` = c("study")
    ),
    "FluBlok" = list(
      H1N1 = c("study"),
      H3N2 = c("study"),
      `B-Victoria` = c("study", "vaccine_strain_change"), 
      `B-Yamagata` = character(0) # doesnt run bc not enough B-Yamagata data
    ),
    "FLUAD" = list(
      # not enough complete H1N1 data,
      H3N2 = character(0),
      `B-Victoria` = character(0), # doesnt run bc no B-Yamagata data
      `B-Yamagata` = character(0)
    )
  )
  
  return(model_terms)
  
}

func_model_terms_overall_fluvac <- function() {
  
  purrr::map(
    func_model_terms_overall(),
    ~ setdiff(.x, "study")
  )
  
}

func_model_terms_subtype_fluvac <- function() {
  
  model_terms <- list(
    "Fluzone-SD" = list(
      H1N1 = c("vaccine_strain", "vaccine_strain_change"),
      H3N2 = c("vaccine_strain", "vaccine_strain_change"),
      `B-Victoria` = c("vaccine_strain", "vaccine_strain_change"),
      `B-Yamagata` = c("vaccine_strain_change")
    ),
    "Fluzone-HD" = list(
      H1N1 = c("vaccine_strain", "vaccine_strain_change"),
      H3N2 = c("vaccine_strain", "vaccine_strain_change"),
      `B-Victoria` = c("vaccine_strain", "vaccine_strain_change"),
      `B-Yamagata` = c("vaccine_strain")
    ),
    "FluMist" = list(
      # not enough complete H1N1 data,
      H3N2 = c("vaccine_strain_change"),
      `B-Victoria` = c("vaccine_strain"),
      `B-Yamagata` = character(0)
    ),
    "Flucelvax" = list(
      # not enough complete H1N1 data,
      H3N2 = c("vaccine_strain_change"),
      `B-Victoria` = character(0), #c("vaccine_strain_change") in model w uriv
      `B-Yamagata` = character(0)
    ),
    "FluBlok" = list(
      H1N1 = character(0),
      H3N2 = character(0),
      `B-Victoria` =  character(0), #c("vaccine_strain_change") in model w uriv
      `B-Yamagata` = character(0) # doesnt run bc not enough B-Yamagata data
    ),
    "FLUAD" = list(
      # not enough complete H1N1 data,
      H3N2 = character(0),
      `B-Victoria` = character(0), # doesnt run bc no B-Yamagata data
      `B-Yamagata` = character(0)
    )
  )
  return(model_terms)
}

# Model fitting ====
func_fit_models_overall <- function(data, model_terms) {
  
  library(tidymodels)
  
  lm_spec <- linear_reg(
    mode = "regression",
    engine = "lm"
  )
  
  vaccines <- unique(as.character(data$vaccine))
  
  models <- setNames(
    vector("list", length(vaccines)),
    vaccines
  )
  
  for (vac in vaccines) {
    
    dat_clean <- data %>%
      filter(vaccine == vac) %>%
      mutate(
        assay_subtype = droplevels(
          factor(assay_subtype)
        ),
        prev_vaccine_abbr = droplevels(
          factor(prev_vaccine_abbr)
        ),
        strain_year = as.numeric(
          stringr::str_extract(
            vaccine_strain, "\\d{4}"
          )
        ),
        participant_cluster = as.factor(subject_id)
      )
    
    # vaccine-level covariates
    terms <- model_terms[[vac]] %||%
      character(0)
    
    base_terms <- c("log_pretiter_sd_overall", "prev_vaccine_abbr",
      "age_sd_overall", "assay_subtype")
    
    rhs <- unique(c(base_terms, terms))
    
    formula <- stats::reformulate(
      rhs,
      response = "log_fc"
    )
    
    fit <- lm_spec %>% fit(formula, data = dat_clean)
    
    fit_result <- list(
      fit = fit,
      cluster = dat_clean$participant_cluster
    )
    
    models[[vac]] <- fit_result
  }
  
  return(models)
}

func_fit_models_subtype <- function(data, model_terms) {
  
  library(tidymodels)
  
  lm_spec <- linear_reg(mode="regression", engine="lm")
  
  combinations <- data %>%
    distinct(vaccine, assay_subtype)
  
  vaccines <- unique(as.character(data$vaccine))
  
  models <- setNames(
    lapply(vaccines, \(x) list()),
    vaccines
  )
  
  for (i in seq_len(nrow(combinations))) {
    
    vac <- as.character(combinations$vaccine[i])
    subtype <- as.character(combinations$assay_subtype[i])
      
      dat_clean <- data %>% 
        filter(
          vaccine == vac, 
          assay_subtype == subtype
        ) %>%
        mutate(
          strain_year = as.numeric(stringr::str_extract(vaccine_strain, "\\d{4}"))
        ) %>%
        mutate(
          vaccine_strain = forcats::fct_reorder(
            vaccine_strain, 
            strain_year, 
            .fun = min,
            na.rm = TRUE
          ),
          participant_cluster = as.factor(subject_id)
        )
      
      terms <- model_terms[[vac]][[subtype]] %||% character(0)
      
      base_terms <- c("log_pretiter_sd_subtype", "prev_vaccine_abbr",
                      "age_sd_subtype")
      
      rhs <- c(base_terms, terms)
      
      formula <- as.formula(
        paste(
          "log_fc ~",
          paste(rhs, collapse = " + ")
        )
      )
      
      fit <- lm_spec %>%
        fit(formula, data = dat_clean)
      
      fit_result <- list(
        fit = fit,
        cluster = dat_clean$participant_cluster
      )
      
      models[[vac]][[subtype]] <- fit_result
  }
  return(models)
}

func_fit_models_imputed_overall <- function(imputed_data, model_terms) {
  
  models <- purrr::map(
    imputed_data,
    ~ func_fit_models_overall(
      data = .x,
      model_terms = model_terms
    )
  )
  
  names(models) <- names(imputed_data)
  
  return(models)
}

func_fit_models_imputed_subtype <- function(imputed_data, model_terms) {
  
  models <- purrr::map(
    imputed_data,
    ~ func_fit_models_subtype(
      data = .x,
      model_terms = model_terms
    )
  )
  
  names(models) <- names(imputed_data)
  
  return(models)
}

# Pool/process fits ====
func_pool_models_overall <- function(models) {
  
  pooled <- list()
  
  vaccines <- names(models[[1]])
  expected_m <- length(models)
  
  for (vac in vaccines) {
    
    # get model from each imputation
    models_vac <- purrr::map(
      models,
      ~ .x[[vac]]
    )
    
    # n obs
    lm_fit <- parsnip::extract_fit_engine(models_vac[[1]]$fit)
    
    n_obs <- stats::model.frame(lm_fit) %>%
      count(
        prev_vaccine_abbr,
        name = "n"
      ) %>%
      mutate(
        term = paste0(
          "prev_vaccine_abbr",
          as.character(prev_vaccine_abbr)
        )
      ) %>%
      select(term, n)
    
    ref_vaccine <- if (vac == "FLUAD") "FZHD" else "FZSD"
    
    n_ref <- n_obs %>%
      filter(
        term == paste0("prev_vaccine_abbr", ref_vaccine)
      ) %>%
      pull(n)
    
    # get estimates and cluster robust SEs
    estimates <- purrr::map_dfr(
      models_vac,
      function(x) {
        
        lm_fit <- parsnip::extract_fit_engine(x$fit)
        
        cluster <- x$cluster
        
        V <- sandwich::vcovCL(
          lm_fit,
          cluster = cluster,
          type = "HC1"
        )
        
        coef_est <- stats::coef(lm_fit)
        
        tibble::tibble(
          term = names(coef_est),
          estimate = unname(coef_est),
          std.error = unname(sqrt(diag(V))[names(coef_est)]
          )
        )
      }
    )
    
    # Rubins rules
    pooled_model <- estimates %>%
      group_by(term) %>%
      summarise(
        m = n(),
        pooled_estimate = mean(estimate),
        within_var = mean(std.error^2),
        between_var = var(estimate),
        total_var = within_var + (1 + 1 / m) * between_var,
        pooled_se = sqrt(total_var),
        statistic = pooled_estimate / pooled_se,
        df = if (is.finite(between_var) && between_var > 0) {
          (m - 1) * (1 + within_var / ((1 + 1 / m) * between_var))^2
        } else {
          Inf
        },
        p.value = 2 * pt(abs(statistic), df = df, lower.tail = FALSE),
        conf.low = pooled_estimate - qt(0.975, df) * pooled_se,
        conf.high = pooled_estimate + qt(0.975, df) * pooled_se, 
        .groups = "drop"
      )
    
    pooled[[vac]] <- pooled_model %>%
      left_join(n_obs, by = "term") %>%
      mutate(n_ref = n_ref)
    
  }
  return(pooled)
}

func_pool_models_subtype <- function(models) {
  
  pooled <- list()
  
  vaccines <- names(models[[1]])
  
  for (vac in vaccines) {
    
    pooled[[vac]] <- list()
    
    subtypes <- names(models[[1]][[vac]])
    
    for (subtype in subtypes) {
      
      # get model from each imputation
      models_sub <- purrr::map(
        models,
        ~ .x[[vac]][[subtype]]
      )
      
      # number of observations in each prior-vaccine group
      lm_fit <- parsnip::extract_fit_engine(models_sub[[1]]$fit)
      
      n_obs <- stats::model.frame(lm_fit) %>%
        count(
          prev_vaccine_abbr,
          name = "n"
        ) %>%
        mutate(
          term = paste0(
            "prev_vaccine_abbr",
            as.character(prev_vaccine_abbr)
          )
        ) %>%
        select(term, n)
      
      # reference prior-season vaccine
      ref_vaccine <- if (vac == "FLUAD") "FZHD" else "FZSD"
      
      # number of observations in the reference group
      n_ref <- n_obs %>%
        filter(
          term == paste0("prev_vaccine_abbr", ref_vaccine)
        ) %>%
        pull(n)
      
      # get estimates and SEs
      estimates <- purrr::map_dfr(
        models_sub,
        function(x) {
          
          lm_fit <- parsnip::extract_fit_engine(x$fit)
          
          # check clujster ids align with fitted rows
          used_rows <- as.integer(
            rownames(stats::model.frame(lm_fit))
          )
          
          cluster <- x$cluster[used_rows]
          
          V <- sandwich::vcovCL(
            lm_fit,
            cluster = cluster,
            type = "HC1"
          )
          
          coef_est <- stats::coef(lm_fit)
          
          tibble::tibble(
            term = names(coef_est),
            estimate = unname(coef_est),
            std.error = unname(sqrt(diag(V))[names(coef_est)]
            )
          )
        }
      )
      
      # Rubins rules
      pooled_model <- estimates %>%
        group_by(term) %>%
        summarise(
          m = n(),
          pooled_estimate = mean(estimate),
          within_var = mean(std.error^2),
          between_var = if (m > 1) {
            var(estimate)
          } else {
            NA_real_
          },
          total_var = within_var + (1 + 1 / m) * between_var,
          pooled_se = sqrt(total_var),
          statistic = pooled_estimate / pooled_se,
          df = if (is.finite(between_var) && between_var > 0) {
            (m - 1) * (1 + within_var / ((1 + 1 / m) * between_var))^2
          } else {
            Inf
          },
          p.value = 2 * pt(abs(statistic), df = df, lower.tail = FALSE),
          conf.low = pooled_estimate - qt(0.975, df) * pooled_se,
          conf.high = pooled_estimate + qt(0.975, df) * pooled_se, 
          .groups = "drop"
        )
      
      pooled[[vac]][[subtype]] <- pooled_model %>%
        left_join(n_obs, by = "term") %>%
        mutate(n_ref = n_ref)
      
    }
  }
  return(pooled)
}

func_pool_estimates_overall <- function(models, filepath) {
  
  fits <- purrr::imap_dfr(models, function(fit_obj, vaccine_name) {
    
    fit_obj %>%
      mutate(
        estimate_exp = 2^pooled_estimate,
        conf.low_exp = 2^conf.low,
        conf.high_exp = 2^conf.high,
        vaccine = vaccine_name
      ) %>%
      select(
        vaccine, term, n, n_ref, estimate=pooled_estimate, se=pooled_se, df, 
        conf.low, conf.high, estimate_exp, conf.low_exp, conf.high_exp, p.value
      ) %>%
      # add clean term names
      mutate(
        term_clean = case_when(
          term == "(Intercept)" ~ "Intercept",
          term %in% c(
            "age_sd", "age_sd_overall", "age_sd_subtype"
          ) ~ "Age",
          
          term %in% c("log_pretiter_sd_overall", "log_pretiter_sd_subtype") ~ "Log pretiter",
          term == "studyURIV" ~ "Study: URIV",
          term == "vaccine_strain_change" ~ "Vaccine strain change",
          
          stringr::str_detect(term, "^prev_vaccine_abbr") ~
            recode(
              stringr::str_remove(term, "^prev_vaccine_abbr"),
              FZSD = "Previous vaccine: Fluzone-SD",
              FZHD = "Previous vaccine: Fluzone-HD",
              FM   = "Previous vaccine: FluMist",
              FV   = "Previous vaccine: Flucelvax",
              FB   = "Previous vaccine: FluBlok",
              FA   = "Previous vaccine: FLUAD",
              ID   = "Previous vaccine: Intradermal"
            ),
          
          stringr::str_detect(term, "^assay_subtype") ~
            recode(
              stringr::str_remove(term, "^assay_subtype"),
              H3N2 = "Assay subtype: H3N2",
              `B-Yamagata` = "Assay subtype: B-Yamagata",
              `B-Victoria` = "Assay subtype: B-Victoria",
            ),
          
          stringr::str_detect(term, "^vaccine_strain[BH]") ~
            paste0(
              "Strain: ",
              stringr::str_remove(term, "^vaccine_strain")
            ),
          
          TRUE ~ term
        )
      )
  })
  # add multiple comparisons correction to prev vaccine pvalues
  fits_correct <- fits %>%
    ungroup() %>%
    mutate(
      p_fdr = p.adjust(
        if_else(
          stringr::str_detect(term, "^prev_vaccine_abbr"),
          p.value,
          NA_real_
        ),
        method = "BH"
      )
    )
  
  fits_round <- fits_correct %>%
    mutate(
      across(
        where(is.numeric) &
          !any_of(c("p.value", "p_fdr")),
        ~ round(.x, digits = 2)
      ),
      across(
        c(p.value, p_fdr),
        ~ round(.x, digits = 3)
      )
    )
  
  saveRDS(fits_round, file = filepath)
  return(fits_round)
}

func_pool_estimates_subtype <- function(models, filepath) {
  
  fits <- purrr::imap_dfr(models, function(sublist, vaccine_name) {
    purrr::imap_dfr(sublist, function(fit_obj, subtype_name) {
      
      fit_obj %>%
        mutate(
          estimate_exp = 2^pooled_estimate,
          conf.low_exp = 2^conf.low,
          conf.high_exp = 2^conf.high,
          vaccine = vaccine_name,
          subtype = subtype_name
        ) %>%
        select(
          vaccine, subtype, term, n, n_ref, estimate=pooled_estimate, 
          se=pooled_se, df, conf.low, conf.high, estimate_exp, conf.low_exp, 
          conf.high_exp, p.value
        ) %>%
        # add clean term names
        mutate(
          term_clean = case_when(
            term == "(Intercept)" ~ "Intercept",
            term %in% c(
              "age_sd", "age_sd_overall", "age_sd_subtype"
            ) ~ "Age",
            
            term %in% c("log_pretiter_sd_overall", "log_pretiter_sd_subtype") ~ "Log pretiter",
            term == "studyURIV" ~ "Study: URIV",
            term == "vaccine_strain_change" ~ "Vaccine strain change",
            
            stringr::str_detect(term, "^prev_vaccine_abbr") ~
              recode(
                stringr::str_remove(term, "^prev_vaccine_abbr"),
                FZSD = "Previous vaccine: Fluzone-SD",
                FZHD = "Previous vaccine: Fluzone-HD",
                FM   = "Previous vaccine: FluMist",
                FV   = "Previous vaccine: Flucelvax",
                FB   = "Previous vaccine: FluBlok",
                FA   = "Previous vaccine: FLUAD",
                ID   = "Previous vaccine: Intradermal"
              ),
            
            stringr::str_detect(term, "^vaccine_strain[BH]") ~
              paste0(
                "Strain: ",
                stringr::str_remove(term, "^vaccine_strain")
              ),
            
            TRUE ~ term
          )
        )
    })
  })
  
  # add multiple comparisons correction to prev vaccine pvalues
  fits_correct <- fits %>%
    mutate(
      p_fdr = p.adjust(
        if_else(
          stringr::str_detect(term, "^prev_vaccine_abbr"),
          p.value,
          NA_real_
        ),
        method = "BH"
      )
    )
  
  fits_round <- fits_correct %>%
    mutate(
      across(
        where(is.numeric) &
          !any_of(c("p.value", "p_fdr")),
        ~ round(.x, digits = 2)
      ),
      across(
        c(p.value, p_fdr),
        ~ round(.x, digits = 3)
      )
    )
  
  saveRDS(fits_round, file = filepath)
  return(fits_round)
}

func_process_estimates_fixed_lod <- function(models, by_subtype = TRUE, filepath) {
  
  fits <- purrr::imap_dfr(models, function(sublist, vaccine_name) {
    
    if (by_subtype) {
      model_list <- sublist
    } else {
      model_list <- list(overall = sublist)
    }
    
    purrr::imap_dfr(model_list, function(fit_obj, subtype_name) {
      
      lm_fit <- parsnip::extract_fit_engine(fit_obj$fit)
      
      # number of observations per prior-vaccine group
      n_obs <- stats::model.frame(lm_fit) %>%
        count(
          prev_vaccine_abbr,
          name = "n"
        ) %>%
        mutate(
          term = paste0(
            "prev_vaccine_abbr",
            as.character(prev_vaccine_abbr)
          )
        ) %>%
        select(term, n)
      
      ref_vaccine <- if (vaccine_name == "FLUAD") "FZHD" else "FZSD"
      
      n_ref <- n_obs %>%
        filter(
          term == paste0("prev_vaccine_abbr", ref_vaccine)
        ) %>%
        pull(n)
      
      cluster <- fit_obj$cluster
      
      V <- sandwich::vcovCL(
        lm_fit,
        cluster = cluster,
        type = "HC1"
      )
      
      n_clusters <- n_distinct(cluster)
      
      estimates <- broom::tidy(lm_fit) %>%
        mutate(
          std.error = sqrt(diag(V)[term]),
          statistic = estimate / std.error,
          
          # cluster-based degrees of freedom
          df = n_clusters - 1,
          
          p.value = 2 * stats::pt(
            abs(statistic),
            df = df,
            lower.tail = FALSE
          ),
          
          conf.low = estimate -
            stats::qt(0.975, df) * std.error,
          
          conf.high = estimate +
            stats::qt(0.975, df) * std.error
        )
      
      estimates <- estimates %>%
        left_join(
          n_obs,
          by = "term"
        ) %>%
        mutate(n_ref = n_ref)
      
      estimates %>%
        mutate(
          estimate_exp = 2^estimate,
          conf.low_exp = 2^conf.low,
          conf.high_exp = 2^conf.high,
          vaccine = vaccine_name,
          subtype = if (by_subtype) subtype_name else NULL
        ) %>%
        mutate(
          term_clean = case_when(
            term == "(Intercept)" ~ "Intercept",
            term %in% c(
              "age_sd", "age_sd_overall", "age_sd_subtype"
            ) ~ "Age",
            
            term %in% c("log_pretiter_sd_overall", "log_pretiter_sd_subtype") ~ "Log pretiter",
            term == "studyURIV" ~ "Study: URIV",
            term == "vaccine_strain_change" ~ "Vaccine strain change",
            
            stringr::str_detect(term, "^prev_vaccine_abbr") ~
              recode(
                stringr::str_remove(term, "^prev_vaccine_abbr"),
                FZSD = "Previous vaccine: Fluzone-SD",
                FZHD = "Previous vaccine: Fluzone-HD",
                FM   = "Previous vaccine: FluMist",
                FV   = "Previous vaccine: Flucelvax",
                FB   = "Previous vaccine: FluBlok",
                FA   = "Previous vaccine: FLUAD",
                ID   = "Previous vaccine: Intradermal"
              ),
            
            stringr::str_detect(term, "^vaccine_strain[BH]") ~
              paste0(
                "Strain: ",
                stringr::str_remove(term, "^vaccine_strain")
              ),
            
            TRUE ~ term
          )
        )
    })
  })
  
  # add multiple comparisons correction to prev vaccine pvalues
  fits_correct <- fits %>%
    ungroup() %>%
    mutate(
      p_fdr = p.adjust(
        if_else(
          stringr::str_detect(term, "^prev_vaccine_abbr"),
          p.value,
          NA_real_
        ),
        method = "BH"
      )
    )
  
  fits_round <- fits_correct %>%
    mutate(
      across(
        where(is.numeric) &
          !any_of(c("p.value", "p_fdr")),
        ~ round(.x, digits = 2)
      ),
      across(
        c(p.value, p_fdr),
        ~ round(.x, digits = 3)
      )
    )
  
  saveRDS(fits_round, file = filepath)
  return(fits_round)
}

# END OF SCRIPT ====