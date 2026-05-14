
# function ----------------------------------------------------------------

#' @title Propensity score matching model function
#' 
#' Executes a propensity score matching model. Returns matched data, significant variables, std mean differences, 
#' and propensity score overlap tests.
#' 
#' @param df dataframe to analyse.
#' @param y.var binary dependent variable.
#' @param x.vars independent variables.
#' @param model.name Name of the model. Will be included in the title of diagnostics plots.
#' @param model.name.short Shorter name of the model. Used to create dataframes assigned to environment.
#' @param method most common are "nearest" for nearest neighbour, "optimal" for optimal pair matching, 
#               "full" for optimal full matching. See ?MatchIt::matchit for more options.
#' @param distance See ?MatchIt::distance for details. Default is logit. Choosing some options may prompt installation of extra packages.
#' @param ratio Number of controls to match to one treated. Must be a positive integer. Default = 1.
#' @param replace Sample controls with (TRUE) or without (FALSE, default) replacement
#' @param caliper Caliper to apply to match distances. Based on standard deviation units. 
reg_psm <- function(df, y.var, x.vars, 
                    model.name = "Propensity score matching",
                    model.name.short = "psm",
                    method = "nearest", distance = "logit",
                    ratio = 1, replace = FALSE, caliper = NULL,
                    ...) {
  
  # input validation
  model.name.short <- gsub(" ", ".", model.name.short)
  
  if(ratio != round(ratio,0) | ratio < 1) {
    stop("ratio must be an integer greater than or equal to 1.")
  }
  
  test.binary <- unique(df[[y.var]])
  if (length(test.binary) != 2) {
    stop("y.var must be a binary outcome variable.")
  } else if ((!setequal(test.binary, c(1, 0)) &&
              !setequal(test.binary, c(TRUE, FALSE)) &&
              !setequal(test.binary, c("Y", "N")) &&
              !setequal(test.binary, c("Yes", "No")))) {
    stop("y.var must be made up pairings of 1/0 or TRUE/FALSE", 
         " (also acceptable are 'Y'/'N' and 'Yes'/'No')")
  }
  
  # create the table for which the regression will be run
  total.rows <- nrow(df %>% select(all_of(c(y.var, x.vars))))
  complete.rows <- nrow(df %>% select(all_of(c(y.var, x.vars))) %>% filter(complete.cases(.)))
  message(paste0("Filtering for complete cases. Number of complete cases: ", 
                 complete.rows, " / ", total.rows, ", ", round(100*complete.rows/total.rows, 2), "%"))
  
  df.reg <- df %>% 
    filter(complete.cases(select(., all_of(c(y.var, x.vars))))) %>% 
    mutate(positive = factor(ifelse(!!sym(y.var) == 1
                                    | !!sym(y.var) == TRUE
                                    | !!sym(y.var) == "Y"
                                    | !!sym(y.var) == "Yes", 
                                    1, 0),
                             levels = c(0,1)))
  y.var <- "positive"
  
  # execute the model
  set.seed(147) # if tree-based method, BART, or other similar gets selected, or if 2+ controls have the same propensity score 
  model <- MatchIt::matchit(formula = reformulate(x.vars, y.var),
                            data = df.reg,
                            method = method,
                            distance = distance,
                            ratio = ratio,
                            replace = replace,
                            caliper = caliper)
  assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
  message(paste0("Model stored as '", model.name.short, ".model'."))
  
  # check how many did match, then extract the matches
  res_summary <- summary(model)
  message(paste0("Number of matches printed below. ",
                 "Matched data stored as '", model.name.short, ".matched'."))
  print(res_summary$nn)
  
  matched_data <- MatchIt::get_matches(object = model)
  assign(paste0(model.name.short, ".matched"), matched_data, envir = parent.frame())
  
  # review underlying model, depending on distance type
  if (distance %in% c("glm", "logit")) {
    
    matching_model <- summary(model$model)$coefficients %>% 
      as.data.frame() %>% 
      mutate(across(is.numeric), round(., 5),
             sig.level = case_when(
               `Pr(>|z|)` < 0.01 ~ "***",
               `Pr(>|z|)` < 0.05 ~ "**",
               `Pr(>|z|)` < 0.1  ~ "*",
               TRUE ~ ""
             ))
    assign(paste0(model.name.short, ".dist.model"), matching_model, envir = parent.frame())
    message(paste0("Coefficients for matching model (", distance, "). ",
                   "Results stored as '", model.name.short, ".dist.model'."))
    print(matching_model)
    
  } else if (distance == "gam") {
    
    summary <- summary(model$model)
    matching_model <- data.frame(
      Estimate = summary$p.coeff, 
      Std.Error = summary$se, 
      t.value = summary$p.t, 
      Pr = summary$p.pv
    ) %>% 
      mutate(across(is.numeric), round(., 5),
             sig.level = case_when(
               Pr < 0.01 ~ "***",
               Pr < 0.05 ~ "**",
               Pr < 0.1  ~ "*",
               TRUE ~ ""
             ))
    assign(paste0(model.name.short, ".dist.model"), matching_model, envir = parent.frame())
    message(paste0("Coefficients for matching model (parametric components of ", distance, "). ",
                   "Results stored as '", model.name.short, ".dist.model'."))
    print(matching_model)
    
  } else {
    
    message(paste0("No coefficients (or equivalent) to show drivers of matching model for distance option ", distance))
    
  }
  
  # diagnostics - SMD / covariate balance
  plt.smd <- cobalt::love.plot(model,
                               binary = "std",
                               threshold = 0.2) +
    labs(subtitle = paste("Model:", model.name)) +
    theme(
      plot.title = element_text(color = "black", face = "bold", hjust = 0.5, size = 10),
      plot.subtitle = element_text(color = "black", hjust = 0.5, size = 9),
    )
  
  # diagnostics - PS overlap
  if (distance %in% c("glm", "logit", "gam")) {
    plt.overlap <- cobalt::bal.plot(model,
                                    var.name = "distance",
                                    which = "both",
                                    type = "density") +
      labs(title = "Distributional Balance",
           subtitle = paste("Model:", model.name)) +
      theme(
        plot.title = element_text(color = "black", face = "bold", hjust = 0.5, size = 10),
        plot.subtitle = element_text(color = "black", hjust = 0.5, size = 9),
      )
  }
  
  # combine plots if needed
  if (distance %in% c("glm", "logit", "gam")) {
    plt.combined <- cowplot::plot_grid(plt.smd, plt.overlap, nrow = 2, rel_heights = c(1.5, 1))
    print(plt.combined)
    message("See Plots pane for standardised mean differences and propensity score overlap tests.")
  } else {
    print(plt.smd)
    message("See Plots pane for standardised mean differences test.")
  }
  
}


# examples ----------------------------------------------------------------

## set up ----
# note that MatchIt may prompt you to install additional packages for certain options

options(scipen=999)

packages <- c(
  # utility
  "zoo", 
  # plotting
  "cowplot", "cobalt", 
  # stats
  "tidymodels", "MatchIt",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}

# note this a toy example using the frailty flag in place of a treatment flag
patient_list <- readr::read_csv("frail_patients.csv") 

patient_list <- patient_list %>% 
  mutate(patient = as.character(patient),
         frailty_flag = factor(frailty_flag, levels = c(0,1)),
         sex = factor(sex, levels = c("F", "M")),
         imd_decile = factor(imd_decile, levels = sort(unique(patient_list$imd_decile))),
         ethnicity = factor(ethnicity, levels = sort(unique(patient_list$ethnicity))),
         geography = factor(geography, levels = sort(unique(patient_list$geography)))
  )


## running models ----
reg_psm(df = patient_list,
        y.var = "frailty_flag",
        x.vars = c("age", "sex", "imd_decile", "ethnicity", "n_chronic_cond", 
                   "n_medication", "n_ae_atts", "n_pc_cons", "geography"),
        model.name = "Propensity score matching model",
        model.name.short = "psm",
        method = "nearest",
        distance = "logit",
        ratio = 1,
        replace = FALSE,
        caliper = NULL)

reg_psm(df = patient_list,
        y.var = "frailty_flag",
        x.vars = c("age"),
        model.name = "PSM model reduced specification",
        model.name.short = "psm.reduced",
        method = "nearest",
        distance = "logit",
        ratio = 1,
        replace = FALSE,
        caliper = NULL)

reg_psm(df = patient_list,
        y.var = "frailty_flag",
        x.vars = c("age", "sex", "ethnicity", "n_chronic_cond", "geography"),
        model.name = "PSM model alternative methods",
        model.name.short = "psm.alt",
        method = "optimal",
        distance = "randomforest",
        ratio = 2)
