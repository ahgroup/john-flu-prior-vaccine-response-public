###
# 03-tables.R
###

# Descriptive ====
func_tab_vaccine_descriptions <- function(vaccine_description_path, filepath) {
  
  dat_types <- readxl::read_xlsx(path = vaccine_description_path)
  tab <- flextable::flextable(dat_types) %>% 
    flextable::width(
      j = "Vaccine",
      width = 1.19
    ) %>%
    flextable::width(
      j = "Manufacturer",
      width = 1.06
    ) %>%
    flextable::width(
      j = "Formulation",
      width = 1.25
    ) %>%
    flextable::width(
      j = "Production",
      width = 1.0
    ) %>%
    flextable::width(
      j = "Age range",
      width = 0.81
    ) %>%
    flextable::width(
      j = "Administration",
      width = 1.19
    )
  
  saveRDS(tab, file = filepath)
  return(tab)
}

func_tab_homologous_strains <- function(homologous_assay_path, filepath) {
  
  homologous_assay_list <- readxl::read_xlsx(homologous_assay_path)
  
  dat_tab <- homologous_assay_list %>% 
    select(
      Season=season, H1N1=H1N1_full_name, H3N2=H3N2_full_name,
      `B-Yamagata`=Yamagata_full_name,`B-Victoria`=Bvictoria_full_name
    ) %>%
    mutate(
      Season = paste0(Season,"-",(as.numeric(Season)+1))
    )
  
  tab <- flextable::flextable(dat_tab) %>% flextable::autofit()
  
  saveRDS(tab, file = filepath)
  
}

func_tab_n_obs <- function(data, stratify_by = c("all", "study", "site"), filepath) {
  
  library(gtsummary)
  
  stratify_by <- match.arg(stratify_by)
  
  if (stratify_by == "study") {
    groups <- c("FluVac", "URIV")
    group_var <- "study"
  } else if (stratify_by == "site") {
    groups <- c("UGA", "PA", "FL", "CCF")
    group_var <- "study_site"
  } else {
    groups <- c("all")
    group_var <- "study_site"
  }
  
  tables <- list()
  
  for (group in groups) {
    
    dat_group <- data %>%
      distinct(person_year_id, .keep_all = TRUE)
    
    if (group != "all") {
      dat_group <- dat_group %>%
        filter(.data[[group_var]] == group)
    }
    
    tab <- dat_group %>% 
      droplevels() %>% # removes 0s column
      tbl_summary(
        by = vaccine,
        include = c(prev_vaccine),
        label = list(
          prev_vaccine = "Prior vaccine",
          vaccine = "Vaccine"
        ),
        statistic = list(all_categorical() ~ "{n}")
      ) %>%
      modify_header(
        label ~ "**Current vaccine**",
        all_stat_cols() ~ "**{level} (n={n})**"
      ) %>%
      modify_table_body(
        ~ .x %>%
          mutate(across(starts_with("stat_"),
                        ~ ifelse(. == "0", "-", .)))
      ) %>%
      as_flex_table() %>%
      flextable::fontsize(size = 10, part = "footer") %>%
      flextable::autofit()
    
    if (group != "all") {
      header_label <- if (group_var == "study") {
        "Study:"
      } else {
        "Study site:"
      }
      
      tab <- flextable::add_header_lines(
        tab,
        values = paste(header_label, group)
      )
    }
    
    tab <- flextable::bg(tab, bg = "white", part = "all")
    saveRDS(tab, file=paste0(filepath, "-", group, ".rds"))
    tables[[group]] <- tab
    
  }
  return(tables)
}

func_tab_n_obs_subtype <- function(data, filepath) {
  
  library(gtsummary)
  
  subtypes <- unique(data$assay_subtype)
  
  tables <- list()
  
  for (subtype in subtypes) {
    
    dat_subtype <- data %>%
      filter(assay_subtype == subtype) %>%
      # drop flublok and fluad for yamagata
      filter(
        !(subtype == "B-Yamagata" & vaccine %in% c("FluBlok", "FLUAD"))
      )
    
    tab <- dat_subtype %>%
      droplevels() %>%
      tbl_summary(
        by = vaccine,
        include = c(prev_vaccine),
        label = list(
          prev_vaccine = "Previous vaccine",
          vaccine = "Vaccine"
        ),
        statistic = list(all_categorical() ~ "{n}")
      ) %>%
      modify_header(
        label ~ "**Current vaccine**",
        all_stat_cols() ~ "**{level} (n={n})**"
      ) %>%
      modify_table_body(
        ~ .x %>%
          mutate(
            across(
              starts_with("stat_"),
              ~ ifelse(. == "0", "-", .)
            )
          )
      ) %>%
      as_flex_table() %>%
      flextable::fontsize(size = 10, part = "footer") %>%
      flextable::autofit()
    
    tab <- flextable::add_header_lines(
      tab,
      values = paste("Influenza subtype:", subtype)
    )
    
    tab <- flextable::bg(tab, bg = "white", part = "all")
    saveRDS(tab, file = paste0(filepath, "-", subtype, ".rds"))
    tables[[subtype]] <- tab
  }
  return(tables)
}

func_tab_n_clusters <- function(data, by_subtype = TRUE, filepath) {
  
  group_vars <- if (by_subtype) {
    c("vaccine", "assay_subtype")
  } else {
    "vaccine"
  }
  
  tab <- data %>%
    group_by(across(all_of(c(group_vars, "subject_id")))) %>%
    summarise(
      n_obs = n(),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_observations = sum(n_obs),
      n_participants = n(),
      n_repeated = sum(n_obs > 1),
      max_repeated = max(n_obs),
      .groups = "drop"
    )
    
  labelled::var_label(tab) <- list(
    n_observations = "Observations",
    n_participants = "Participants",
    n_repeated = "Participants with >1 observation",
    max_repeated = "Maximum observations per participant"
  )
  
  tab_flex <- tab %>% 
    flextable::flextable() %>%
    flextable::width(
      j = c("n_observations", "n_participants"),
      width = 0.95
    )
  
  saveRDS(tab_flex, filepath)
  return(tab_flex)
  
}

func_tab_fc_small <- function(data, filepath) {
  
  library(gtsummary)
  library(flextable)
  
  dat_small <- data %>% 
    filter(fc < 0.25)
  
  vaccines <- unique(dat_small$vaccine)
  subtypes <- c("H1N1", "H3N2", "B-Yamagata", "B-Victoria")
  
  tabs <- list()
  
  for (st in subtypes) {
    
    dat_st_all <- dat_small %>% filter(assay_subtype == st)
    
    if (nrow(dat_st_all) == 0) next
    
    sub_tabs <- list()
    sub_vaccines <- character(0)
    
    for (vax in vaccines) {
      
      dat_filter <- dat_st_all %>% filter(vaccine == vax)
      
      if (nrow(dat_filter) == 0) next
      
      keep_levels <- dat_filter %>%
        count(prev_vaccine) %>%
        filter(n > 0) %>%
        pull(prev_vaccine)
      
      dat_vax <- dat_filter %>%
        mutate(
          fc = factor(round(fc, 4)),
          prev_vaccine = factor(prev_vaccine, levels = keep_levels)
        )
      
      tab <- dat_vax %>%
        tbl_summary(
          by = prev_vaccine,
          include = fc,
          label = list(fc = paste0(st, " fold change")),
          statistic = list(all_categorical() ~ "{n}")
        ) %>%
        modify_header(label ~ "")
      
      sub_tabs[[vax]] <- tab
      sub_vaccines <- c(sub_vaccines, vax)
    }
    
    if (length(sub_tabs) == 0) next
    
    # merge all vaccine sub-tables side by side, spanning header = vaccine name
    merged_tab <- tbl_merge(
      tbls       = sub_tabs,
      tab_spanner = sub_vaccines
    ) %>%
      modify_table_body(
        ~ .x %>%
          mutate(across(starts_with("stat_"),
                        ~ ifelse(is.na(.) | . == "0", "-", .)))
      ) %>%
      as_flex_table() %>%
      flextable::autofit()
    
    saveRDS(merged_tab, paste0(filepath, st, ".rds"))
    tabs[[st]] <- merged_tab
  }
  return(tabs)
}

func_tab_demographics <- function(data, filepath) {
  
  library(gtsummary)
  
  # create new dataset of eligible observations by unique subject_id
  dat_dem <- data %>%
    distinct(subject_id, .keep_all = TRUE) %>%
    select(subject_id, study_site, age_min, gender, race) %>%
    mutate(
      race = case_when(
        race %in% c("Not Specified", "Other", "Unknown") ~ NA_character_,
        TRUE ~ race
      ),
      study_site = factor(study_site, levels = c("FL","PA","UGA","CCF","URIV"))
    )
  
  # create table 
  tab <- tbl_summary(
    data = dat_dem,
    by = study_site,
    include = c(age_min, gender, race),
    label = list(
      age_min = "Age at first enrollment",
      gender = "Gender",
      race = "Race"
    )
  ) %>%
    add_overall() %>%
    modify_table_body(
      ~ .x %>%
        mutate(across(starts_with("stat_"),
                      ~ ifelse(. == "0 (NA%)", "-", .)))
    ) %>%
    as_flex_table() %>%
    flextable::autofit()
  
  tab <- flextable::bg(tab, bg = "white", part = "all")
  saveRDS(tab, file=filepath)
  return(tab)
}

func_tab_censored_counts <- function(data, filepath) {
  
  library(gtsummary)
  
  dat_lod <- data %>% 
    group_by(vaccine, prev_vaccine, assay_subtype) %>%
    summarise(
      n_lod = sum(pretiter < 10, na.rm = TRUE) + sum(posttiter < 10, na.rm = TRUE),
      n_total = sum(!is.na(pretiter)) + sum(!is.na(posttiter)),
      percent_lod = round(100 * n_lod / n_total, 2),
      .groups = "drop"
    ) %>%
    arrange(vaccine, prev_vaccine, assay_subtype)
  
  tab <- dat_lod %>% 
    flextable::flextable() %>%
    flextable::set_header_labels(
      #study_site = "Study",
      vaccine = "Current vaccine",
      assay_subtype = "Subtype",
      timepoint = "Timepoint",
      n_lod = "LOD (n)",
      n_total = "Total (n)",
      percent_lod = "LOD (%)"
    ) %>%
    flextable::colformat_num(
      j = "percent_lod",
      digits = 1
    ) %>%
    flextable::autofit()
  
  saveRDS(tab, filepath)
  return(tab)
  
}

func_tab_strain_combs <- function(data, filepath) {
  
  # produces cross-tabulations for every current-vaccine/subtype combination
  vaccines <- unique(data$vaccine)
  subtypes <- unique(data$assay_subtype)
  results <- list()
  
  for (vac in vaccines) {
    results[[vac]] <- list()
    
    for (subtype in subtypes) {
      
      dat_clean <- data %>%
        filter(
          vaccine == vac,
          assay_subtype == subtype
        ) %>%
        droplevels()
      
      tab_prior_strain <- table(
        dat_clean$prev_vaccine,
        dat_clean$vaccine_strain,
        useNA = "ifany"
      )
      
      tab_prior_change <- table(
        dat_clean$prev_vaccine,
        dat_clean$vaccine_strain_change,
        useNA = "ifany"
      )
      
      tab_strain_change <- table(
        dat_clean$vaccine_strain,
        dat_clean$vaccine_strain_change,
        useNA = "ifany"
      )
      
      # store tables
      results[[vac]][[subtype]] <- list(
        prior_vs_strain = tab_prior_strain,
        prior_vs_change = tab_prior_change,
        strain_vs_change = tab_strain_change
      )
      
      saveRDS(
        tab_prior_strain,
        file = paste0(filepath, vac, "-", subtype, "-prior-vs-strain.rds")
      )
      
      saveRDS(
        tab_prior_change,
        file = paste0(filepath, vac, "-", subtype, "-prior-vs-change.rds")
      )
      
      saveRDS(
        tab_strain_change,
        file = paste0(filepath, vac, "-", subtype, "-strain-vs-change.rds")
      )
    }
  }
  return(results)
}

# Model specifications ====
func_tab_model_terms <- function(model_terms, model_type, filepath) {
  
  model_spec_table <- purrr::map_dfr(
    names(model_terms),
    function(vac) {
      
      if (model_type == "subtype") {
        
        purrr::map_dfr(
          names(model_terms[[vac]]),
          function(subtype) {
            
            terms <- model_terms[[vac]][[subtype]]
            
            tibble::tibble(
              `Current-season vaccine` = vac,
              `Influenza subtype` = subtype,
              Study = if ("study" %in% terms) "Yes" else "No",
              `Vaccine strain` =
                if ("vaccine_strain" %in% terms) "Yes" else "No",
              `Vaccine strain change` =
                if ("vaccine_strain_change" %in% terms) "Yes" else "No"
            )
          }
        )
        
      } else {
        
        terms <- model_terms[[vac]]
        
        tibble::tibble(
          `Current-season vaccine` = vac,
          Study = if ("study" %in% terms) "Yes" else "No",
          `Vaccine strain` =
            if ("vaccine_strain" %in% terms) "Yes" else "No",
          `Vaccine strain change` =
            if ("vaccine_strain_change" %in% terms) "Yes" else "No"
        )
      }
    }
  )
  
  # Apply subtype-specific exclusions only
  if (model_type == "subtype") {
    
    model_spec_table <- model_spec_table %>%
      dplyr::filter(
        !(
          `Current-season vaccine` %in%
            c("FluBlok", "FLUAD") &
            `Influenza subtype` == "B-Yamagata"
        )
      )
  }
  
  # Preserve current-season vaccine ordering
  model_spec_table <- model_spec_table %>%
    dplyr::mutate(
      `Current-season vaccine` = factor(
        `Current-season vaccine`,
        levels = c(
          "Fluzone-SD",
          "Fluzone-HD",
          "FluMist",
          "Flucelvax",
          "FluBlok",
          "FLUAD"
        )
      )
    ) %>%
    dplyr::arrange(`Current-season vaccine`)
  
  tab <- model_spec_table %>%
    flextable::flextable() %>%
    flextable::autofit()
  
  saveRDS(tab, filepath)
  
  return(tab)
}

# Regression estimates ====
func_tab_regression <- function(estimates, filepath) {

  by_subtype <- "subtype" %in% names(estimates)
  
  if (by_subtype) {
    
    fits_grouped <- estimates %>%
      group_by(vaccine, subtype) %>%
      group_split()
    
  } else {
    
    fits_grouped <- estimates %>%
      group_by(vaccine) %>%
      group_split()
    
  }
  
  tables <- purrr::map(fits_grouped, function(dat) {
      
      vaccine <- unique(dat$vaccine)
      
      if (by_subtype) {
        subtype <- unique(as.character(dat$subtype))
        table_name <- paste(vaccine, subtype, sep = "-")
        table_header <- paste(vaccine, "-", subtype)
      } else {
        subtype <- NULL
        table_name <- vaccine
        table_header <- vaccine
      }
      
      table <- dat %>%
        filter(!is.na(estimate_exp)) %>%
        select(
          term_clean, estimate_exp, conf.low_exp, conf.high_exp, p.value, p_fdr
        ) %>%
        mutate(
          `Exponentiated estimate (95% CI)` = paste0(
            estimate_exp," (",conf.low_exp, "-",conf.high_exp, ")"
          ),
          `P-value` = format.pval(
            p.value,
            digits = 3,
            eps = 0.001
          ),
          `BH-adjusted p` = p_fdr
        ) %>%
        select(
          Term = term_clean, `Exponentiated estimate (95% CI)`, `P-value`, 
          `BH-adjusted p`
        ) %>%
        flextable::flextable() %>%
        flextable::add_header_lines(values = table_header) %>%
        flextable::autofit() %>%
        flextable::width(
          j = "Term",
          width = 2.63
        ) %>%
        flextable::width(
          j = "Exponentiated estimate (95% CI)",
          width = 2.19
        ) %>%
        flextable::width(
          j = "P-value",
          width = 0.63
        ) %>%
        flextable::width(
          j = "BH-adjusted p",
          width = 1.06
        )
      
        list(
          name = table_name,
          vaccine = vaccine,
          subtype = subtype,
          table = table
        )
    })
  
  purrr::iwalk(
    tables,
    ~ saveRDS(
      .x$table,
      file = paste0(filepath, .x$name, ".rds")
    )
  )
  return(tables)
}
# END OF SCRIPT ====