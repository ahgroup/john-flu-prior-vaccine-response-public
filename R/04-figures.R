###
# 04-figures.R
###

# Descriptive ====
func_plot_current <- function(data, filepath, height, width) {
  
  library(ggplot2)
  set.seed(333)
  
  ybreak_list <- list(
    pretiter = c(5,10,20,40,80,160,320,640,1280,2560),
    posttiter = c(5,10,20,40,80,160,320,640,1280,2560,5120),
    fc = c(0.25,0.5,1,2,4,8,16,32,64,128,256,512)
  )
  
  ylimits <- list(
    pretiter = c(5,2560),
    posttiter = c(5,5120),
    fc = c(0.25, 512)
  )
  
  dat_long <- data %>%
    mutate(
      vaccine = factor(
        vaccine,
        levels = c("Fluzone-SD","Fluzone-HD","FluMist", "Flucelvax","FluBlok",
                   "FLUAD")
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c(pretiter, posttiter, fc),
      names_to = "outcome",
      values_to = "value"
    ) %>%
    mutate(
      outcome = factor(
        outcome,
        levels = c("pretiter","posttiter","fc"),
        labels = c("Pre-titer","Post-titer","Fold change")
      )
    )
  
  p <- ggplot(
    dat_long, 
    aes(
      x = vaccine, y = value, shape = study, color = age_factor)
  ) +
    geom_jitter(size = 1.5, alpha = 0.5, width = 0.25, height = 0) +
    geom_hline(
      data = dat_long %>% filter(outcome == "Fold change"),
      aes(yintercept = 1),
      linetype = "dashed", 
      color = "grey30"
    ) +
    geom_hline(
      data = dat_long %>% filter(outcome == "Post-titer"),
      aes(yintercept = 40),
      linetype = "dashed", 
      color = "grey30"
    ) +
    scale_color_manual(
      values = c(
        "<65"     = ggokabeito::palette_okabe_ito()[1],
        "≥65"     = ggokabeito::palette_okabe_ito()[2],
        "Unknown" = ggokabeito::palette_okabe_ito()[8]
      )
    ) +
    stat_summary(
      inherit.aes = FALSE,
      aes(x = vaccine, y = value, group = vaccine),
      fun.data = function(x) data.frame(
        ymin = quantile(x, 0.25, na.rm = TRUE),
        ymax = quantile(x, 0.75, na.rm = TRUE)
      ),
      geom = "errorbar", width = 0.2, color = "black"
    ) +
    stat_summary(
      inherit.aes = FALSE,
      aes(x = vaccine, y = value, group = vaccine),
      fun = median,
      geom = "point",
      color = "black", shape = 18, size = 3
    ) +
    ggh4x::facet_grid2(
      outcome ~ assay_subtype,
      scales = "free",
      switch = "y",
      drop = TRUE
    ) +
    ggh4x::facetted_pos_scales(
      y = list(
        outcome == "Pre-titer"   ~ scale_y_continuous(
          trans = "log2", 
          breaks = ybreak_list$pretiter,  
          limits = ylimits$pretiter
        ),
        outcome == "Post-titer"  ~ scale_y_continuous(
          trans = "log2", 
          breaks = ybreak_list$posttiter, 
          limits = ylimits$posttiter
        ),
        outcome == "Fold change" ~ scale_y_continuous(
          trans = "log2", 
          breaks = ybreak_list$fc,
          limits = ylimits$fc,
          labels = ybreak_list$fc
        )
      )
    ) +
    labs(
      x = "Current vaccine",
      y = NULL,
      color = "Age",
      shape = "Study"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    plot = p,
    filename = filepath,
    height = height,
    width = width
  )
}

func_plot_prior_current <- function(data, filepath, height, width) {
  
  library(ggplot2)
  set.seed(333)

  outcomes <- c(
    pretiter = "Pre-titer",
    posttiter = "Post-titer",
    fc = "Fold change"
  )
  
  ybreaks <- list(
    pretiter = c(5,10,20,40,80,160,320,640,1280,2560),
    posttiter = c(5,10,20,40,80,160,320,640,1280,2560,5120),
    fc = c(0.25,0.5,1,2,4,8,16,32,64,128,256,512)
  )
  
  ylimits <- list(
    pretiter = c(5,2560),
    posttiter = c(5,5120),
    fc = c(0.25,512)
  )
  
  plots <- list()
  
  for (outcome in names(outcomes)) {
    
    dat_plot <- data %>%
      filter(
        !is.na(.data[[outcome]]),
        .data[[outcome]] > 0
      ) %>%
      mutate(
        vaccine = factor(
          vaccine,
          levels = c(
            "Fluzone-SD", "Fluzone-HD", "FluMist", "Flucelvax", "FluBlok", 
            "FLUAD"
          )
        ),
        assay_subtype = factor(
          assay_subtype,
          levels = c("H1N1", "H3N2","B-Yamagata", "B-Victoria")
        ),
        prev_vaccine = factor(
          prev_vaccine,
          levels = levels(data$prev_vaccine)
        )
      )
    
    p <- ggplot(
      dat_plot,
      aes(
        x = prev_vaccine,
        y = .data[[outcome]],
        color = age_factor,
        shape = study
      )
    ) +
      geom_jitter(width = 0.25, height = 0, alpha = 0.6, size = 1.5) +
      stat_summary(
        aes(x = prev_vaccine, y = .data[[outcome]], group = prev_vaccine),
        inherit.aes = FALSE,
        fun.data = function(x) data.frame(
          ymin = quantile(x, 0.25),
          ymax = quantile(x, 0.75)
        ),
        geom = "errorbar",
        width = 0.2,
        color = "black"
      ) +
      stat_summary(
        aes(x = prev_vaccine, y = .data[[outcome]]),
        inherit.aes = FALSE,
        fun = median,
        geom = "point",
        color = "black",
        shape = 18,
        size = 3
      ) +
      ggh4x::facet_grid2(
        rows = vars(vaccine),
        cols = vars(assay_subtype),
        scales = "free_x",
        #independent = "x",
        switch = "y"
      ) +
      scale_color_manual(
        values = c(
          "<65" = ggokabeito::palette_okabe_ito()[1],
          "≥65" = ggokabeito::palette_okabe_ito()[2],
          "Unknown" = ggokabeito::palette_okabe_ito()[8]
        )
      ) +
      scale_shape_manual(
        values = c("FluVac" = 16, "URIV" = 17)
      ) +
      scale_y_continuous(
        trans = "log2",
        breaks = ybreaks[[outcome]],
        limits = ylimits[[outcome]],
        labels = ybreaks[[outcome]]
      ) +
      labs(
        x = "Prior-season vaccine",
        y = outcomes[[outcome]],
        color = "Age",
        shape = "Study"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
    
    if (outcome == "posttiter") {
      p <- p +
        geom_hline(
          yintercept = 40,
          linetype = "dashed",
          color = "grey40"
        )
    }
    
    if (outcome == "fc") {
      p <- p +
        geom_hline(
          yintercept = 1,
          linetype = "dashed",
          color = "grey40"
        )
    }
    
    plots[[outcome]] <- p
    
    ggsave(
      filename = paste0(filepath, outcome, ".png"),
      plot = p,
      height = height,
      width = width
    )
  }
}

func_plot_prior_current_by <- function(data, facet, filepath, height, width) {
  
  library(ggplot2)
  set.seed(333)
  
  # outcomes
  yvars <- c("pretiter", "posttiter", "fc")
  ylabels <- c(
    pretiter = "Pre-titer", 
    posttiter = "Post-titer", 
    fc = "Fold change"
  )
  
  ybreak_list <- list(
    pretiter = c(5,10,20,40,80,160,320,640,1280,2560),
    posttiter = c(5,10,20,40,80,160,320,640,1280,2560,5120),
    fc = c(0.25,0.5,1,2,4,8,16,32,64,128,256,512)
  )
  
  ylimits <- list(
    pretiter = c(5,2560),
    posttiter = c(5,5120),
    fc = c(0.25, 512)
  )
  
  plots <- list()
  
  for (y in yvars) {
    
    p <- ggplot(
      data,
      aes(x = prev_vaccine, y = .data[[y]])
    ) +
      geom_jitter(
        aes(shape = study,  color = age_factor),
        size = 1.5, alpha = 0.6, width = .25, height = 0
      ) +
      stat_summary(
        fun.data = function(x) data.frame(
          ymin = quantile(x, 0.25),
          ymax = quantile(x, 0.75)
        ),
        geom = "errorbar", 
        width = 0.2, 
        color = "black"
      ) +
      stat_summary(
        fun = median, 
        geom = "point", 
        color = "black",
        shape = 18,
        size = 3
      ) +
      scale_color_manual(
        values = c(
          "<65" = ggokabeito::palette_okabe_ito()[1],
          "≥65" = ggokabeito::palette_okabe_ito()[2],
          "Unknown" = ggokabeito::palette_okabe_ito()[8]
        )
      ) +
      scale_shape_manual(
        values = c(
          "FluVac" = 16,
          "URIV" = 17
        )
      ) +
      scale_y_continuous(
        trans = "log2",
        breaks = ybreak_list[[y]],
        limits = ylimits[[y]],
        labels = ybreak_list[[y]]
      ) +
      ggh4x::facet_nested(
        rows = vars(vaccine),
        cols = vars(assay_subtype, !!rlang::sym(facet)),
        scales = "free_x",
        space = "free_x",
        switch = "y",
        nest_line = TRUE,
        drop = TRUE
      ) +
      scale_x_discrete(drop = TRUE) +
      labs(
        x = "Prior-season vaccine",
        y = ylabels[[y]],
        color = "Age",
        shape = "Study"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom", 
        panel.grid.minor = element_blank()
      )
    
    if (y == "posttiter") {
      p <- p +
        geom_hline(
          yintercept = 40,
          linetype = "dashed",
          color = "grey30"
        )
    }
    
    if (y == "fc") {
      p <- p +
        geom_hline(
          yintercept = 1,
          linetype = "dashed",
          color = "grey30"
        )
    }
    
    plots[[y]] <- p
    
    ggsave(
      plot = p,
      filename = paste0(filepath, y, ".png"),
      height = height,
      width = width
    )
  }
}

func_plot_pre_post <- function(data, filepath, height, width) {
  
  library(ggplot2)
  library(patchwork)
  set.seed(333)
  
  subtypes <- unique(data$assay_subtype)
  plots <- list()
  
  for (sub in subtypes) {
    
    dat_sub <- data %>%
      filter(assay_subtype == sub)
    
    vaccines <- levels(droplevels(dat_sub$vaccine))
    
    plots[[sub]] <- list()
    
    for (vac in vaccines) {
      
      dat_filter <- data %>% filter(vaccine == vac, assay_subtype == sub)
      
      dat_long <- dat_filter %>%
        mutate(prev_vaccine = factor(prev_vaccine)) %>%
        tidyr::pivot_longer(
          cols = c(pretiter, posttiter),
          names_to = "titer_type",
          values_to = "titer"
        ) %>%
        mutate(
          titer_type = recode(
            titer_type, pretiter = "Pretiter", posttiter = "Posttiter"
          ),
          titer_type = factor(titer_type, levels = c("Pretiter", "Posttiter"))
        )
      
      # get manual jitter
      dodge_width <- 0.25
      jitter_width <- 0.1
      dat_jitter <- dat_long %>%
        mutate(
          prev_vac_numeric = as.numeric(prev_vaccine),
          dodge_offset = ifelse(
            titer_type == "Pretiter", -dodge_width/2, dodge_width/2
          ),
          x_jittered = prev_vac_numeric + dodge_offset + 
            runif(n(), -jitter_width/2, jitter_width/2)
        )
      
      # medians
      dat_median <- dat_long %>%
        group_by(prev_vaccine, titer_type) %>%
        summarise(med = median(titer, na.rm = TRUE), .groups = "drop") %>%
        mutate(
          prev_vac_numeric = as.numeric(prev_vaccine),
          dodge_offset = ifelse(
            titer_type == "Pretiter", -dodge_width/2, dodge_width/2
          ),
          x = prev_vac_numeric + dodge_offset
        )
      
      p <- ggplot(
        dat_jitter, 
        aes(x = x_jittered, y = titer, color = titer_type)
      ) +
        # lines connecting pre/post for each person
        geom_line(aes(group = person_year_id), color = "grey40", alpha = 0.1) +
        # points
        geom_point(size = 1, alpha = 0.5) +
        geom_point(
          data = dat_median,
          aes(x = x, y = med),
          size = 3, shape = 18, color = "black", alpha = 1
        ) +
        ggokabeito::scale_color_okabe_ito(order = c(3, 5))+
        scale_x_continuous(
          breaks = seq_along(levels(dat_jitter$prev_vaccine)),
          labels = levels(dat_jitter$prev_vaccine)
        ) +
        scale_y_continuous(
          trans = "log2",
          breaks = c(5,10,20,40,80,160,320,640,1280,2560,5120),
          limits = c(5,5120)
        ) +
        labs(
          title = paste0(vac),
          x = "Prior-season vaccine",
          y = "HAI titer",
          color = NULL
        ) +
        theme_bw() +
        theme(
          legend.position = "bottom", 
          text = element_text(size=12),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      
      plots[[sub]][[vac]] <- p
    }
    
    sub_plots <- plots[[sub]]
    
    combined <- patchwork::wrap_plots(sub_plots, ncol = 3) +
      plot_layout(guides = "collect") +
      plot_layout(axis_titles = "collect") +
      plot_annotation(title = sub) &
      theme(legend.position = "bottom")
    
    ggsave(
      plot = combined,
      filename = paste0(filepath, sub, ".png"),
      height = if (sub == "H1N1") height / 1.8 else height,
      width = width
    )
  }
}

# Regression estimates ====
func_plot_regression <- function(estimates, filepath, height, width) {
  
  library(ggplot2)
  library(patchwork)
  
  by_subtype <- "subtype" %in% names(estimates)
  
  dat_plot <- estimates %>%
    filter(
      stringr::str_starts(term, "prev_vaccine_abbr")
    ) %>%
    mutate(
      prev_vaccine = factor(
        stringr::str_remove(
          term, "^prev_vaccine_abbr"
        ),
        levels = c(
          "FZSD", "FZHD", "FM", "FV",
          "FB", "FA", "ID"
        ),
        labels = c(
          "Fluzone-SD", "Fluzone-HD",
          "FluMist", "Flucelvax",
          "FluBlok", "FLUAD",
          "Intradermal"
        )
      ),
      vaccine = factor(
        vaccine,
        levels = c(
          "Fluzone-SD", "Fluzone-HD",
          "FluMist", "Flucelvax",
          "FluBlok", "FLUAD"
        )
      )
    )
  if (by_subtype) {
    dat_plot <- dat_plot %>%
      mutate(
        subtype = factor(
          subtype,
          levels = c(
            "H1N1", "H3N2","B-Yamagata", "B-Victoria"
          )
        )
      )
  }
  
  dat_ref <- dat_plot %>%
    distinct(
      vaccine,
      across(
        any_of("subtype")
      ),
      n_ref
    ) %>%
    mutate(
      ref_label = paste0("Ref: n=", n_ref)
    )
  
  xlimits <- c(
    min(dat_plot$conf.low_exp, na.rm = TRUE),
    max(dat_plot$conf.high_exp, na.rm = TRUE)
  )  
  
  p <- ggplot(
    dat_plot,
    aes(
      x = estimate_exp, 
      y = prev_vaccine, 
      xmin = conf.low_exp, 
      xmax = conf.high_exp,
    )
  ) +
    geom_pointrange(
      aes(color=prev_vaccine)
      )+
    geom_text(
      aes(x=Inf, label = paste0("n=", .data$n)),
      nudge_y = -0.25,
      hjust = 1.05,
      size = 3.5,
      color = "grey35",
      show.legend = FALSE
    ) +
    geom_text(
      data = dat_ref,
      aes(x = Inf, y = Inf, label = ref_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.2,
      size = 2.8,
      color = "grey40"
    ) +
    scale_color_manual(
      values = c(
        "Fluzone-SD" = "#2166AC", # blue
        "Fluzone-HD" = "#CC79A7", # pink
        "FluMist" = "#8C564B", # brown
        "Flucelvax" = "#9467BD", # purple
        "FluBlok" = "#E6AB02", # yellow
        "FLUAD" = "#009E73", # green
        "Intradermal" = "#D55E00" # red-orange
      ),
      name = "Prior-season vaccine"
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    scale_y_discrete(
      position = "right",
      limits = rev,
      expand = expansion(add = c(0.5, 0.9))
    ) +
    scale_x_continuous(
      transform = "log2", 
      limits = xlimits,
      breaks = c(0.25, 0.5, 1, 2, 4),
      labels = c("0.25", "0.5", "1.0", "2.0", "4.0")
    ) +
    labs(
      x = "Adjusted fold-change ratio (95% CI)",
      y = "Prior-season vaccine"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  if (by_subtype) {
    
    p <- p +
      facet_grid(
        rows = vars(vaccine),
        cols = vars(subtype),
        scales = "free_y",
        space = "free_y",
        switch = "y"
      )
    
  } else {
    
    p <- p +
      facet_grid(
        rows = vars(vaccine),
        scales = "free_y",
        space = "free_y",
        switch = "y"
      )
    
  }
  
  ggsave(
    plot = p,
    filename = filepath,
    height = height,
    width = width
  )
}
# END OF SCRIPT ====