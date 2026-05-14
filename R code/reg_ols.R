
# function ----------------------------------------------------------------

#' @title OLS regression function
#' 
#' Executes an OLS regression. Returns 'tidy' results table, fit tests, and normality tests.
#' Current version of model does not return Ljung-Box test for autocorrelated errors, but will be added as option.
#' 
#' @param df dataframe to analyse.
#' @param y.var dependent variable.
#' @param x.vars independent variables. Can include interaction terms.
#' @param model.name Name of the model. Will be included in the title of diagnostics plots.
#' @param model.name.short Shorter name of the model. Used to create objects assigned to environment.
#' @param cluster.se.var For clustered standard errors, which field to cluster by.
reg_ols <- function(df, y.var, x.vars, 
                    model.name = "OLS regression model",
                    model.name.short = "ols.reg",
                    cluster.se.var = NULL,
                    ...) {
  
  # helpers
  model.name.short <- gsub(" ", ".", model.name.short)
  
  f_wrap_strings <- function(vector_of_strings, width){
    as.character(sapply(vector_of_strings, FUN = function(x){
      paste(strwrap(x,width=width), collapse="\n")}
    ))
  }
  
  # get variables
  interaction_regs <- grep("\\*", x.vars, value = TRUE)
  interaction_replace <- unique(unlist(strsplit(interaction_regs, "\\*")))
  allvar <- c(y.var, 
              cluster.se.var,
              unlist(lapply(x.vars, function(x) {
                if (x %in% interaction_regs) {
                  return(interaction_replace)
                } else {
                  return(x)
                }
              }))
  )
  df.reg <- df %>% select(all_of(allvar))
  
  # run model
  model <- lm(reformulate(x.vars, y.var), data = df.reg)
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
      mutate_if(is.numeric, ~round(.,3))
    
  } else {
    
    res <- tidy(lmtest::coeftest(model, 
                                 vcov = vcovHC(model, type = 'HC0'))) %>% 
      select(-c('statistic')) %>%
      mutate(sig.level = case_when(p.value < 0.01 ~ "1% level",
                                   p.value < 0.05 ~ "5% level",
                                   p.value < 0.1  ~ "10% level",
                                   TRUE ~ "Not significant")
      ) %>%
      mutate_if(is.numeric, ~round(.,3))
    
  }
  
  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  message(paste0("Regression results, stored in dataframe '", model.name.short, ".res'. ", 
                 "Model stored as '", model.name.short, ".model'."))
  print(res)
  
  # vectors for testing
  gap <- nrow(df.reg) - length(model$residuals)
  actual <- df.reg[[y.var]] %>% {if (gap > 0) tail(.,-gap) else .}
  pred <- model$fitted.values
  resids <- model$residuals
  stdresids <- rstandard(model)
  cook <- cooks.distance(model)
  
  test.plt.df <- cbind(actual, pred, resids, stdresids, cook) %>%
    as.data.frame() %>% 
    rename_with(~c('Actual', 'Fitted', 'Residuals', 'Std. residuals', 'Cook')) %>% 
    mutate(id = row_number())
  
  # model fit
  ## R-squared. Higher R2 = better fit. Crude measure.
  test.r2 <- round(summary(model)$r.squared, 3)
  ## ANOVA. Which variables drive the model? If not significant, consider removing from model
  test.anova <- anova(model, test="Chisq") %>%
    mutate(`Pr(>F)` = round(`Pr(>F)`, 3))
  anova.insignificant <- rownames(test.anova[test.anova$'Pr(>Chi)' > 0.05, ])
  if (is_empty(anova.insignificant)) {
    anova.insignificant <- "N/A"
  }
  message("Results from ANOVA test.")
  print(test.anova)
  ## Root MSE. Lower = better fit. An option for comparing models.
  test.rmse <- round(ModelMetrics::rmse(actual, pred),3)
  ## Akaike Information Criteria. For comparing similar models. Lower = better.
  test.aic <- stats::AIC(model)
  ## Mincer-Zarnowitz test. More for in/out-of-sample, but useful check. How well does model fit actuals? Should sit on 45-deg line
  # null = intercept = 0 and slope = 1 (good fit)
  mz.fit <- lm(actual ~ pred)
  mz.b <- coef(mz.fit)
  mz.cov <- vcov(mz.fit)  
  mz.s <- (mz.b - c(0, 1)) %*% solve(mz.cov) %*% c(mz.b - c(0, 1))
  test.mz <- pchisq(mz.s, 2, lower.tail = FALSE)
  ## plot actuals vs fitted
  plt.title <- f_wrap_strings(paste0("Actual vs fitted activity for ", model.name, sep = ""), 100)
  plt.test.scatter.fit <- ggplot(test.plt.df, aes(x = Actual)) +
    geom_smooth(aes(y = Fitted, color = "Fitted"), method = "lm", se = TRUE, linewidth = 1.25) +
    geom_point(aes(y = Fitted), color = "black", shape = 1, size = 1.25, alpha = 0.25, stroke = 0, show.legend = FALSE) +
    geom_line(aes(y = Actual, color = "45-degree"), linetype = "dashed", linewidth = 1.1, show.legend = TRUE) +
    scale_color_manual(values = c("Fitted" = "black", "45-degree" = "red")) +
    labs(x = "Actual", y = "Fitted", color = "Model",
         title = plt.title,
         subtitle = paste0("Mincer-Zarnowitz test p-value: ", test.mz)) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, face = "italic"),
          axis.title = element_text(size = 10))
  
  # normality tests
  ## Jarque-Bera test (null = normally distributed)
  test.jb <- try(round(jarque.bera.test(stdresids)$p.value, 3))
  ## Q-Q plot
  plt.title <- f_wrap_strings(paste0("QQ-plot for ", model.name, sep = ""), 100)
  plt.test.qq <- ggplot(test.plt.df, aes(sample = stdresids)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = plt.title,
         subtitle = paste0("Jarque-Bera test p-value: ", test.jb)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, face = "italic"),
          axis.title = element_text(size = 10))
  
  ## autocorrelated errors
  # use Breusch-Godfrey test (null = no autocorrelation)
  test.bg <- try(round(bgtest(model)$p.value, 3))
  
  # residuals plots
  ## residuals vs fitted values. Want a "cloud" with no clear pattern
  plt.title <- f_wrap_strings(paste0("Residuals vs fitted values for ", model.name, sep = ""), 100)
  plt.test.scatter.resid <- ggplot(test.plt.df, aes(x = Fitted, y = Residuals)) +
    geom_point(color = "black", shape = 1, size = 2, alpha = 1, stroke = 0, show.legend = FALSE) +
    labs(title = plt.title) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))
  ## Cook's distances. Want low values and need to investigate outliers
  plt.title <- f_wrap_strings(paste0("Cook's distances for ", model.name, sep = ""), 100)
  rule_of_thumb = 4 / (nrow(df.reg) - length(x.vars) - 1)
  plt.test.scatter.cook <- ggplot(test.plt.df, aes(x = id, y = Cook)) +
    geom_point(color = "black", shape = 1, size = 2, alpha = 1, stroke = 0, show.legend = FALSE) +
    geom_hline(yintercept = rule_of_thumb, color = "red", linetype = "dashed") +
    geom_text(
      data = subset(test.plt.df, Cook > rule_of_thumb),
      aes(label = id),
      vjust = -0.5, size = 3, color = "red"
    ) +
    labs(x = "row id", y = "Cook's Distance",
         title = plt.title,
         subtitle = paste0("Red line = rule of thumb threshold 4 / (N - k - 1) = ", round(rule_of_thumb, 3))
         ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9, face = "italic"),
          axis.title = element_text(size = 10))
  
  # print a grid of diagnostic plots
  plt.test.all <- cowplot::plot_grid(plt.test.scatter.fit, plt.test.qq,
                                     plt.test.scatter.resid, plt.test.scatter.cook, 
                                     nrow = 2
  )
  print(plt.test.all)
  
  # print diagnostic p values
  test.results <- data.frame(
    'Model' = model.name,
    'R-squared' = test.r2, 
    'RMSE' = test.rmse,
    'AIC' = test.aic,
    'ANOVA insignificant' = paste(unlist(anova.insignificant), collapse = ", "),
    'Mincer-Zarnowitz fit' = test.mz,
    'Jacque-Bera normality' = test.jb, 
    'Breusch–Godfrey correlation' = test.bg,
    stringsAsFactors = FALSE
  )
  assign(paste0(model.name.short, ".test"), test.results, envir = parent.frame())
  message("Results from model fit tests, stored in dataframe '", paste0(model.name.short, ".test"), "'. ",
          "See Plots pane for QQ plot, residuals plots, and actual vs pred.")
  print(t(test.results))
  
}

# examples ----------------------------------------------------------------

## set up ----
options(scipen=999)

packages <- c(
  # utility
  "zoo", 
  # plotting
  "scales", "cowplot", 
  # stats
  "stats", "tidymodels", "lmtest", "sandwich", "tseries",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}

data(iris) 
iris <- iris %>%
  mutate(Species = as.factor(Species))

## running regressions ----
reg_ols(df = iris,
        y.var = "Sepal.Length",
        x.vars = c("Sepal.Width", "Petal.Length*Petal.Width", "Species"),
        model.name = "Model 1: All variables, interact petal variables",
        model.name.short = "mod1",
        cluster.se.var = "Species")

reg_ols(df = iris,
        y.var = "Sepal.Length",
        x.vars = "Species",
        model.name = "Model 2: under-specified model",
        model.name.short = "mod2")
