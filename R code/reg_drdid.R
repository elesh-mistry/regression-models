
# function ----------------------------------------------------------------

#' @title Doubly-robust difference-in-differences regression function
#' 
#' Executes a doubly-robust difference-in-differences regression for group-time averaged treatment effects.
#' For background, see https://bcallaway11.github.io/did/reference/att_gt.html and https://www.sciencedirect.com/science/article/abs/pii/S0304407620303948 
#' Returns group- and time-averaged treatment effects by period.
#' 
#' @param df dataframe to analyse.
#' @param y.var dependent variable.
#' @param x.vars covariates.
#' @param model.name Name of the model. Will be included in the title of diagnostics plots.
#' @param model.name.short Shorter name of the model. Used to create dataframes assigned to environment.
#' @param time.field Field containing time period. Should be consecutive integers.
#' @param id.field Field containing unit IDs. Ideally consecutive integers.
#' @param treat.period Field stating which period (in time.field) the first treatment took place. 0 if untreated.
#' @param cluster.var Field to cluster standard errors by. Default is id.field.
#' @param method which method to use. Default is "dr" for doubly robust. See ?did::att_gt for more.
#' @param unbalanced Allow for an unbalanced panel. Default is FALSE. Note processing time is much slower if TRUE.
#' @param control.group what to use as control group. Default is "nevertreated". See ?did::att_gt for more.
#' @param anticipation The number of time periods before participating in the treatment where units can anticipate participating in the treatment. Default = 0.
#' @param alpha Set alpha for confidence intervals. Default = 0.05.
#' @param n.iter number of iterations for producing bootstrap confidence intervals Default = 1,000, but reduce if working with very large dataframes.
#' @param nthreads Number of processors used for parallel calculations.

reg_drdid <- function(df, y.var, x.vars,
                      model.name = "Doubly-robust difference-in-differences model",
                      model.name.short = "drdid.reg",
                      time.field, id.field, treat.period, cluster.var = id.field,
                      method = "dr", unbalanced = FALSE,
                      control.group = "nevertreated", anticipation = 0,
                      alpha = 0.05,
                      n.iter = 1000, nthreads = parallel::detectCores() - 1) {
    
  # execute the model
  set.seed(147)
  model <- did::att_gt(yname = y.var,
                       tname = time.field,
                       idname = id.field,
                       gname = treat.period,
                       xformla = reformulate(x.vars),
                       data = df,
                       panel = TRUE,
                       allow_unbalanced_panel = unbalanced,
                       control_group = control.group,
                       est_method = method,
                       anticipation = anticipation,
                       clustervars = cluster.var,
                       bstrap = TRUE,
                       biters = n.iter,
                       pl = nthreads > 1,
                       cores = nthreads,
                       alp = alpha)
  assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
  
  # combined results table
  agg.att.group <- aggte(model, type = "group", na.rm = TRUE)
  res.group <- tidy(agg.att.group) %>% 
    select(type, group, estimate, std.error, conf.low, conf.high) %>%
    rename(time.group.name = group) %>% 
    mutate(time.group.name = as.character(time.group.name))
  
  agg.att.time <- aggte(model, type = "dynamic", na.rm = TRUE)
  res <- tidy(agg.att.time) %>% 
    select(type, event.time, estimate, std.error, conf.low, conf.high) %>%
    rename(time.group.name = event.time) %>% 
    mutate(time.group.name = as.character(time.group.name)) %>% 
    bind_rows(res.group) %>% 
    mutate_if(is.numeric, round, digits = 4) %>% 
    mutate(significant = ifelse(conf.low > 0 | conf.high < 0, TRUE, FALSE))

  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  print(res)
    
  # results plots
  plt.time <- ggdid(agg.att.time) +
    labs(x = "Periods pre/post", 
         y = "Average Treatment Effect", 
         title = "Time-averaged treatment effect",
         subtitle = paste("Model:", model.name),
         color = "Period") +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 12, hjust = 0),
          plot.subtitle = element_text(size = 11, hjust = 0),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
    )
  
  plt.group <- ggdid(agg.att.group) +
    labs(x = "Average Treatment Effect", 
         y = "Group", 
         title = "Group-averaged treatment effect", 
         subtitle = paste("Model:", model.name),
         color = "Period") +
    coord_flip() +
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          plot.title = element_text(size = 12, hjust = 0),
          plot.subtitle = element_text(size = 11, hjust = 0),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
    )
  
  res.plt <- cowplot::plot_grid(plt.time, plt.group,
                                nrow = 2)
  print(res.plt)
  assign(paste0(model.name.short, ".plt"), res.plt, envir = parent.frame())
  
  message(paste0("Regression results, stored in dataframe '", model.name.short, ".res'. ", 
                 "Plot stored as '", model.name.short, ".plt'. ", 
                 "Model stored as '", model.name.short, ".model'."))
  
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
  "tidymodels", "did",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}

mpdta <- mpdta # in-built data from package did

## running regressions ----
reg_drdid(df = mpdta,
          y.var = "lemp",
          x.vars = c("lpop"),
          model.name = "Base DRDID model",
          model.name.short = "drdid.reg",
          time.field = "year",
          id.field = "countyreal",
          treat.period = "first.treat")


reg_drdid(df = mpdta,
          y.var = "lemp",
          x.vars = c("lpop"),
          model.name = "Not-yet-treated control",
          model.name.short = "drdid.notyet",
          time.field = "year",
          id.field = "countyreal",
          treat.period = "first.treat",
          control.group = "notyettreated",
          anticipation = 1)
