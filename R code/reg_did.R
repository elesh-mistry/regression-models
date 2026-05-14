
# function ----------------------------------------------------------------

#' @title Difference-in-differences regression function
#' 
#' Executes a simple difference-in-differences model using OLS regression.
#' Returns regression outputs
#' 
#' @param df dataframe to analyse.
#' @param y.var dependent variable.
#' @param treat.field Field containing treatment flag. Binary: 0 for untreated period, 1 for treated period.
#' @param time.field Field containing time periods. Ideally binary: 0 for untreated period, 1 for treated period.
#' @param x.vars optional covariates.
#' @param model.name Name of the model. Will be included in the title of diagnostics plots.
#' @param model.name.short Shorter name of the model. Used to create dataframes assigned to environment.
#' @param fixed.effects Apply fixed effects. Choice of 'none' (default), 'individual', 'time', 'twoways', or 'nested'
#' @param group.field If applying group-level fixed effects, specify which field represents the group
#' @param fe.method If applying group-level fixed effects, what method to use. See ?plm::plm for more. Default is "within" estimators.
#' @param cluster.se.var For clustered standard errors, which field to cluster by. If using fixed effects, it will only be based on group.field
reg_did <- function(df, y.var, treat.field, time.field, x.vars = NULL,
                    model.name = "Difference-in-differences model",
                    model.name.short = "did.reg",
                    fixed.effects = "none", group.field = NULL, fe.method = "within",
                    cluster.se.var = NULL, ...) {

  # input validation
  model.name.short <- gsub(" ", ".", model.name.short)
  
  test.binary <- unique(df[[treat.field]])
  if (length(test.binary) != 2) {
    warning("'treat.field' is not a binary outcome variable.", 
            "Treat the DID estimator with caution as only a pooled estimate is produced.", 
            "Model will continue with current specification.")
  } else if ((!setequal(test.binary, c(1, 0)) &&
              !setequal(test.binary, c(TRUE, FALSE)) &&
              !setequal(test.binary, c("Y", "N")) &&
              !setequal(test.binary, c("Yes", "No")))) {
    warning("treat.field should ideally be made up pairings of 1/0 or TRUE/FALSE", 
         " (also acceptable are 'Y'/'N' and 'Yes'/'No'). Model will continue with current specification.")
  }
  
  test.binary <- unique(df[[time.field]])
  if (length(test.binary) != 2 && !fixed.effects %in% c("time", "twoways", "nested")) {
    warning(paste("time.field is not a binary variable.",
                  "For data with multiple treatment periods, consider using fixed.effects = 'time', 'twoways' or 'nested'.",
                  "Alternatively, especially for staggered treatment timing, consider using doubly-robust DID.",
                  "Continuing with a pooled estimator."))
  } else if ((!setequal(test.binary, c(1, 0)) &&
              !setequal(test.binary, c(TRUE, FALSE)) &&
              !setequal(test.binary, c("Y", "N")) &&
              !setequal(test.binary, c("Yes", "No")))) {
    stop("time.field must be made up pairings of 1/0 or TRUE/FALSE", 
         " (also acceptable are 'Y'/'N' and 'Yes'/'No')")
  }
  
  fixed.effects <- tolower(fixed.effects)
  fixed.effects <- case_when(
    fixed.effects == "twfe" ~ "twoways",
    fixed.effects == "unit" ~ "individual",
    is.null(fixed.effects) ~ "none",
    TRUE ~fixed.effects)
  if(!fixed.effects %in% c("none", "individual", "time", "twoways", "nested")) {
    stop("'fixed.effects' must be one of 'none', 'individual', 'time', 'twoways', or 'nested'")
  }
  
  if(fixed.effects %in% c("individual", "twoways", "nested") && is.null(group.field)) {
    stop(paste0("With 'fixed.effects' = ", fixed.effects, "must specify 'group.field' for calculating fixed effects."))
  }
  
  fe.method <- tolower(fe.method)
  if(!fe.method %in% c("within", "random", "ht", "between", "pooling", "fd")) {
    stop("'fe.method' must be one of 'within', 'random', 'ht', 'between', 'pooling', 'fd'.")
  }
  
  # get variables
  interaction_regs <- grep("\\*", x.vars, value = TRUE)
  interaction_replace <- unique(unlist(strsplit(interaction_regs, "\\*")))
  allvar <- c(y.var, 
              treat.field,
              time.field,
              group.field,
              unlist(lapply(x.vars, function(x) {
                if (x %in% interaction_regs) {
                  return(interaction_replace)
                } else {
                  return(x)
                }
              }))
  )
  df.reg <- df %>% 
    select(all_of(allvar)) %>% 
    mutate(treat.indicator = ifelse(!!sym(treat.field) == 1
                                    | !!sym(treat.field) == TRUE
                                    | !!sym(treat.field) == "Y"
                                    | !!sym(treat.field) == "Yes", 
                                    1, 0),
           time.indicator = ifelse(!!sym(time.field) == 1
                                    | !!sym(time.field) == TRUE
                                    | !!sym(time.field) == "Y"
                                    | !!sym(time.field) == "Yes", 
                                    1, 0)
    )
  
  
  # run model
  if(fixed.effects == "none") {
    
    reg.vars <- c(x.vars, "time.indicator*treat.indicator")
    
    model <- lm(reformulate(reg.vars, y.var), 
                data = df.reg)
    assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
    
    if (!is.null(cluster.se.var)) {
      
      res <- tidy(lmtest::coeftest(model, 
                                   vcov = vcovCL(model, type = 'HC0', 
                                                 cluster = df.reg[[cluster.se.var]]))) %>% 
        select(-c('statistic')) %>%
        mutate(sig.level = case_when(p.value < 0.01 ~ "1% level",
                                     p.value < 0.05 ~ "5% level",
                                     p.value < 0.1  ~ "10% level",
                                     TRUE ~ "Not significant")
        ) %>%
        mutate_if(is.numeric, ~round(.,3)) %>% 
        mutate(term = ifelse(term == "time.indicator:treat.indicator", "DID estimate", term))
      
    } else {
      
      res <- tidy(lmtest::coeftest(model, 
                                   vcov = vcovHC(model, type = 'HC0'))) %>% 
        select(-c('statistic')) %>%
        mutate(sig.level = case_when(p.value < 0.01 ~ "1% level",
                                     p.value < 0.05 ~ "5% level",
                                     p.value < 0.1  ~ "10% level",
                                     TRUE ~ "Not significant")
        ) %>%
        mutate_if(is.numeric, ~round(.,3)) %>% 
        mutate(term = ifelse(term == "time.indicator:treat.indicator", "DID estimate", term))
      
    }

    # tests
    test.bg <- try(round(bgtest(model)$p.value, 3))
    test.pcd <- NA_real_
    
  } else {
    
    reg.vars <- c(x.vars, "time.indicator*treat.indicator")
    
    panel <- plm::pdata.frame(df.reg, index = c(group.field, "time.indicator"))
    model <- plm::plm(reformulate(reg.vars, y.var), 
                      data = panel,
                      effect = fixed.effects,
                      model = fe.method)
    
    if (!is.null(cluster.se.var)) {
      warning("'cluster.se.var' has been specified, but package 'plm' only clusters on 'group.field'.")
    }
      
    res <- tidy(lmtest::coeftest(model, 
                                 vcov = plm::vcovHC(model, type = 'HC0', 
                                                    cluster = "group"))) %>% 
      select(-c('statistic')) %>%
      mutate(sig.level = case_when(p.value < 0.01 ~ "1% level",
                                   p.value < 0.05 ~ "5% level",
                                   p.value < 0.1  ~ "10% level",
                                   TRUE ~ "Not significant")
      ) %>%
      mutate_if(is.numeric, ~round(.,3)) %>% 
      mutate(term = ifelse(term == "time.indicator1:treat.indicator", "DID estimate", term))
      
  
    # tests
    test.bg <- try(round(plm::pbgtest(model)$p.value, 3))
    test.pcd <- try(round(plm::pcdtest(model)$p.value, 3))
    
  }
  
  res.did <- res %>% filter(term == "DID estimate")
  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  message(paste0("Regression results stored in dataframe '", model.name.short, ".res'. ",
                 "DID estimate printed below. ",
                 "Model stored as '", model.name.short, ".model'."))
  print(res.did)
  
  # print diagnostic p values
  test.results <- data.frame(
    'Model' = model.name,
    'Breusch–Godfrey correlation' = test.bg,
    'Pesaran CD test' = test.pcd,
    stringsAsFactors = FALSE
  )
  assign(paste0(model.name.short, ".test"), test.results, envir = parent.frame())
  message("Results from model fit tests, stored in dataframe '", paste0(model.name.short, ".test"), "'. ")
  print(t(test.results))

}


# examples ----------------------------------------------------------------

## set up ----
options(scipen=999)

packages <- c(
  # utility
  "zoo", 
  # stats
  "stats", "tidymodels", "lmtest", "sandwich", "plm",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}

patient_list <- readr::read_csv("patient_treatment.csv") 

patient_list <- patient_list %>% 
  mutate(patient = as.character(patient),
         sex = factor(sex, levels = c("F", "M")),
         imd_decile = factor(imd_decile, levels = sort(unique(patient_list$imd_decile))),
         ethnicity = factor(ethnicity, levels = sort(unique(patient_list$ethnicity))),
         geography = factor(geography, levels = sort(unique(patient_list$geography)))
  )


## running regressions ----
reg_did(df = patient_list %>% filter(time %in% c(0,1)),
        y.var = "n_pc_cons",
        treat.field = "treated",
        time.field = "time",
        model.name = "Simple DID model",
        model.name.short = "did.simple",
        cluster.se.var = "geography")


reg_did(df = patient_list %>% filter(time %in% c(0,1)),
        y.var = "n_pc_cons",
        treat.field = "treated",
        time.field = "time",
        x.vars = c("age", "sex", "imd_decile", "ethnicity", "n_chronic_cond", 
                   "n_medication", "n_ae_atts", "geography"),
        model.name = "TWFE DID model with covariates",
        model.name.short = "did.twfe",
        fixed.effects = "twoways",
        group.field = "patient",
        fe.method = "within")

