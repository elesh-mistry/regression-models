

# functions ---------------------------------------------------------------

#' @title Kaplan-Meier survival function
#' 
#' Executes a Kaplan-Meier survival model. Returns survival curves, risk table, censor plots, 
#' and log-rank test if the curves are stratified.
#' 
#' @param df dataframe to analyse.
#' @param period.start first period the observed unit appears
#' @param period.stop last period the observed unit appears
#' @param period.event a binary for whether an event occurs. Ideally 0/1, but also accepts 1/2
#' @param strata variables by which to stratify the curve
#' @param model.name.short Shorter name of the model. Used to create objects assigned to environment.
#' @param conf.int Show confidence intervals for the curves
#' @param risk.table Show on the chart the number of units still at risk
#' @param ncensor.plot Show on the chart the number of units censored
reg_km <- function(df, period.start, period.stop, period.event,
                   strata = NULL, model.name.short = "surv.km",
                   conf.int = TRUE, risk.table = FALSE, ncensor.plot = FALSE, 
                   ...) {

  # input validation
  test.binary <- unique(df[[period.event]])
  if (length(test.binary) != 2) {
    stop("period.event must be a binary outcome variable of pairings 0/1, 1/2, or FALSE/TRUE.")
  } else if ((!setequal(test.binary, c(0, 1)) &&
              !setequal(test.binary, c(1, 2)) &&
              !setequal(test.binary, c(TRUE, FALSE)))) {
    stop("y.var must be made up pairings of 1/0", 
         " (also acceptable are 1/2, and TRUE/FALSE)")
  }

  # count unique values in strata. If only 1, will cause errors in producing some plots
  if (is.null(strata)) {
    i <- 1
  } else {
    i <- 0
    for(s in strata) {
      j <- length(unique(df[[s]]))
      i <- i + j
    }
    if (i == 1) message("strata only has one unique value and will be ignored.")
  }
  
  # create formula string to feed into model function
  rhs <- if (i == 1) "1" else paste(strata, collapse = "+")
  formula_str <<- paste0(
    "survival::Surv(time = ", period.start,
    ", time2 = ", period.stop,
    ", event = ", period.event, ") ~ ",
    rhs
  )
  
  # run model
  model <- survival::survfit(as.formula(formula_str), data = df)
  assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
  res <- tidy(model)
  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  
  # plot model. Issues with survminer package handling no strata for ncensor.plot means need if statement
  if (i == 1) {
    
    plt <- survminer::ggsurvplot(fit = model, 
                                 data = df,
                                 size = 1,
                                 conf.int = conf.int, 
                                 risk.table = risk.table,
                                 ncensor.plot = ncensor.plot,
                                 palette = "#F8766D")
    
  } else {
    
    plt <- survminer::ggsurvplot(fit = model, 
                                 data = df,
                                 size = 1,
                                 conf.int = conf.int, 
                                 risk.table = risk.table,
                                 ncensor.plot = ncensor.plot)
    
  }
  
  print(plt)
  assign(paste0(model.name.short, ".plt"), plt, envir = parent.frame())
  
  # run log-rank test if multiple strata
  if (i != 1) {
    formula_diff_str <<- paste0(
      "survival::Surv(", period.stop, ", ", period.event, ") ~ ",
      rhs
    )
    compare <- survival::survdiff(as.formula(formula_diff_str), data = df)
    message("Results from log-rank test on survival curves by strata:")
    print(compare)
  }
  
  message(paste0("Table of results stored in dataframe '", model.name.short, ".res'. ", 
                 "Plot stored as '", model.name.short, ".plt'. ", 
                 "Model stored as '", model.name.short, ".model'."))
  
}



#' @title Cox proportional hazards survival function
#' 
#' Executes a Cox proportional hazards survival model. Returns survival curves, risk table, censor plots, 
#' and regression results.
#' 
#' @param df dataframe to analyse.
#' @param period.start first period the observed unit appears
#' @param period.stop last period the observed unit appears
#' @param period.event a binary for whether an event occurs. Ideally 0/1, but also accepts 1/2
#' @param covariates covariates for regression
#' @param strata covariates for regression by which to stratify the modelled survival curve. Curves are fitted at means for numeric and last for non-numeric
#' @param model.name.short Shorter name of the model. Used to create objects assigned to environment.
#' @param conf.int Show confidence intervals for the curves
#' @param risk.table Show on the chart the number of units still at risk
#' @param ncensor.plot Show on the chart the number of units censored
reg_cox <- function(df, period.start, period.stop, period.event, covariates,
                    strata = NULL, cluster.var = NULL, model.name.short = "surv.cox",
                    conf.int = TRUE, risk.table = FALSE, ncensor.plot = FALSE,
                    ...) {
  
  # input validation
  test.binary <- unique(df[[period.event]])
  if (length(test.binary) != 2) {
    stop("period.event must be a binary outcome variable of pairings 0/1, 1/2, or FALSE/TRUE.")
  } else if ((!setequal(test.binary, c(0, 1)) &&
              !setequal(test.binary, c(1, 2)) &&
              !setequal(test.binary, c(TRUE, FALSE)))) {
    stop("y.var must be made up pairings of 1/0", 
         " (also acceptable are 1/2, and TRUE/FALSE)")
  }
  
  if (!is.null(covariates) && !is.null(strata)) {
    overlap <- intersect(covariates, strata)
    if (length(overlap) > 0) {
      stop(paste0("Variables found in both 'covariates' and 'strata': ",
                  paste(overlap, collapse = ", ")))
    }
  }
  
  # count unique values in strata. If only 1, will cause errors in producing some plots
  if (is.null(strata)) {
    i <- 1
  } else {
    i <- 0
    for(s in strata) {
      j <- length(unique(df[[s]]))
      i <- i + j
    }
    if (i == 1) message("strata only has one unique value and will be ignored.")
  }
  
  # create formula string to feed into model function
  covar_term <- paste(covariates, collapse = "+")
  strata_term <- if (i == 1) "" else paste("+", strata, collapse = "+")
  cluster_term <- if (!is.null(cluster.var)) paste0("+ cluster(", cluster.var, ")") else ""
  
  formula_str <<- paste0(
    "survival::Surv(time = ", period.start,
    ", time2 = ", period.stop,
    ", event = ", period.event, ") ~ ",
    covar_term,
    strata_term,
    cluster_term
  )
  
  # run model
  model <- survival::coxph(as.formula(formula_str), data = df)
  assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
  
  res <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
  res <- res %>%
    mutate(sig.level = case_when(p.value < 0.01 ~ "1% level",
                                 p.value < 0.05 ~ "5% level",
                                 p.value < 0.1  ~ "10% level",
                                 TRUE ~ "Not significant")) %>%
    mutate_if(is.numeric, ~round(., 3))
  
  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  model_summary <- summary(model)
  
  # model diagnostics
  
  ## Concordance = if take a pair of patients, expect failure in higher risk patient first
  test.concordance <- paste0(round(model_summary$concordance[1], 4), 
                             " (s.e.: ", round(model_summary$concordance[2], 4), ")")
  
  # rough estimate of Uno's C, using end periods
  formula_uno_str <<- paste0(
    "survival::Surv(", period.stop, ", ", period.event, ") ~ ",
    covar_term,
    strata_term,
    cluster_term
  )
  model_uno <- survival::coxph(as.formula(formula_uno_str), data = df)
  concordance_uno <- survival::concordance(model_uno, timewt = "n/G2")
  test.concordance.uno <- paste0(round(concordance_uno$concordance, 4),
                                 " (s.e.: ", round(sqrt(concordance_uno$var), 4), ")")
  
  ## proportional hazards assumption
  test.ph <- cox.zph(model)
  test.ph.global <- round(test.ph$table["GLOBAL", "p"], 4)
  test.ph.fail <- test.ph$table[
    rownames(test.ph$table) != "GLOBAL" & test.ph$table[, "p"] <= 0.1, 
    "p", 
    drop = FALSE
  ]
  test.ph.fail <- if (nrow(test.ph.fail) == 0) {
    "No individual variable with p-value > 0.1"
  } else {
    paste0(rownames(test.ph.fail), " (", round(test.ph.fail[, "p"], 4), ")", 
           collapse = ", ")
  }
  
  # PH as graph
  if(nrow(df) < 10000) {
    
    message("Schoenfeld's test and test of influential observations (using dfbetas) printed in 'plots'.")
    print(ggcoxzph(test.ph))
    print(ggcoxdiagnostics(model, type = "dfbetas"))
    
  } else {
    
    message("Due to large number of rows, test of constatnt PH using Schoenfeld residuals is printed for each variable individually in 'plots'.",
            "Test of influential observations (using dfbetas) has not been printed.")
    
    print(plot(cox.zph(model)))
  }
  
  # fitting model
  
  if (i == 1) {
    
    surv_at_means <- survival::survfit(model)
    plt <- survminer::ggsurvplot(fit = surv_at_means, 
                                 data = df,
                                 size = 1,
                                 conf.int = conf.int, 
                                 risk.table = risk.table,
                                 ncensor.plot = ncensor.plot,
                                 color = "#F8766D")
    
  } else {
    
    df_plt <- df %>%
      mutate(across(where(is.factor), as.character))
    covar_numeric <- colnames(df_plt %>% select(all_of(covariates)) %>% select_if(is.numeric))
    covar_character <- colnames(df_plt %>% select(all_of(covariates)) %>% select_if(is.character))
    
    predict_at <- df_plt %>%
      group_by(across(all_of(strata))) %>%
      summarise(
        across(all_of(covar_numeric), mean),
        across(all_of(covar_character), last)
      ) %>%
      ungroup()
    
    strata_labels <- predict_at %>%
      unite("label", all_of(strata), sep = ", ") %>%
      pull(label)
    
    surv_at_predict <- survival::survfit(model, newdata = predict_at)
    
    plt <- survminer::ggsurvplot(fit = surv_at_predict, 
                                 data = df,
                                 size = 1,
                                 legend.labs = strata_labels,
                                 conf.int = conf.int, 
                                 risk.table = risk.table,
                                 ncensor.plot = ncensor.plot)
    
  }
  
  print(plt)
  
  # print diagnostics
  test.results <- data.frame(
    'Model' = model.name.short,
    'Concordance' = test.concordance, 
    'Concordance Uno' = test.concordance.uno,
    'PH global p' = test.ph.global,
    'PH low p' = test.ph.fail,
    stringsAsFactors = FALSE
  )
  assign(paste0(model.name.short, ".test"), test.results, envir = parent.frame())
  
  message(paste0("Table of exponentiate coefficients stored in dataframe '", model.name.short, ".res'. ", 
                 "Model stored as '", model.name.short, ".model'. ",
                 "Results from model fit tests, stored in dataframe '", paste0(model.name.short, ".test'")))
  t(test.results)
  
}


# examples ----------------------------------------------------------------

## set up ----

options(scipen=999)

packages <- c(
  # utility
  "zoo", 
  # plotting
  "cowplot", 
  # stats
  "tidymodels", "survival", "survminer",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}


patient_list <- readr::read_csv("patient_mortality.csv")

patient_list <- patient_list %>% 
  mutate(intervention = factor(intervention, levels = c(0,1)),
         sex = factor(sex, levels = c("F", "M")))

patient_list_aggregated <- patient_list %>% 
  group_by(patient, intervention) %>% 
  summarise(start = min(start),
            stop = max(stop),
            event = max(event),
            age = mean(age),
            sex = last(sex),
            n_chronic_cond = mean(n_chronic_cond),
            n_medication = mean(n_medication),
            n_ae_atts = sum(n_ae_atts),
            n_nel_admissions = sum(n_nel_admissions),
            geography = last(geography)
  ) %>% 
  ungroup() %>% 
  filter(stop > start)

## Kaplan-Meier examples ----
reg_km(df = patient_list_aggregated,
       period.start = "start",
       period.stop = "stop",
       period.event = "event",
       strata = NULL,
       model.name.short = "surv.km")

reg_km(df = patient_list_aggregated,
       period.start = "start",
       period.stop = "stop",
       period.event = "event",
       strata = "intervention",
       model.name.short = "surv.km.split",
       conf.int = TRUE,
       risk.table = TRUE,
       ncensor.plot = TRUE)


## Cox examples ----

# interesting to note the differences between the intervention and non-intervention group before reviewing the models:
patient_list_aggregated %>% 
  mutate(female = ifelse(sex == "F", 1, 0)) %>% 
  group_by(intervention) %>% 
  summarise(age = mean(age),
            female = mean(female),
            n_chronic_cond = mean(n_chronic_cond),
            n_medication = mean(n_medication),
            n_ae_atts = mean(n_ae_atts),
            n_nel_admissions = mean(n_nel_admissions)) %>% 
  ungroup()


reg_cox(df = patient_list_aggregated,
        period.start = "start",
        period.stop = "stop",
        period.event = "event",
        covariates = c("sex", "age", "n_chronic_cond", "n_medication", "n_ae_atts", "n_nel_admissions"),
        strata = "intervention",
        cluster.var = "geography",
        model.name.short = "surv.cox.agg",
        conf.int = TRUE,
        risk.table = TRUE,
        ncensor.plot = TRUE)


# toy example using the non-aggregated data. Note the censor plot wouldn't make sense here.
reg_cox(df = patient_list,
        period.start = "start",
        period.stop = "stop",
        period.event = "event",
        covariates = c("sex", "age", "n_chronic_cond", "n_medication", "n_ae_atts", "n_nel_admissions"),
        strata = "intervention",
        cluster.var = "geography",
        model.name.short = "surv.cox.full",
        conf.int = TRUE,
        risk.table = TRUE,
        ncensor.plot = FALSE)
