
# function ----------------------------------------------------------------

#' @title Logistic regression function
#' 
#' Executes a simple logistic regression. Returns 'tidy' results table, ROC curve, f-beta at different cut-offs,
#' and a confusion matrix at optimal f-beta.
#' 
#' @param df dataframe to analyse.
#' @param y.var binary dependent variable.
#' @param x.vars independent variables. Can include interaction terms.
#' @param model.name Name of the model. Will be included in the title of diagnostics plots.
#' @param model.name.short Shorter name of the model. Used to create dataframes assigned to environment.
#' @param filter.complete Filter out all NAs. Recommend keeping as TRUE, as the model may fail to predict values otherwise.
#' @param cluster.se.var For clustered standard errors, which field to cluster by.
#' @param std.error.corrected For non-clustered standard errors, re-estimate corrected confidence intervals.
#' @param f.beta Beta value for f-beta scores.
#' @param train.pc Percentage of data for training dataset (vs testing).

reg_logit <- function(df, y.var, x.vars, 
                      model.name = "Logit regression model",
                      model.name.short = "logit.reg",
                      filter.complete = TRUE,
                      cluster.se.var = NULL, 
                      std.error.corrected = FALSE,
                      f.beta = 1, train.pc = 1, ...) {
  
  # input validation
  model.name.short <- gsub(" ", ".", model.name.short)
  
  if (!is.null(train.pc) && (train.pc <= 0 | train.pc > 1)) stop("train.pc must be between 0 and 1 (inclusive).")
  
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
  interaction_regs <- grep("\\*", x.vars, value = TRUE)
  interaction_replace <- unique(unlist(strsplit(interaction_regs, "\\*")))
  allvar <- c(y.var, 
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
    mutate(positive = factor(ifelse(!!sym(y.var) == 1
                                    | !!sym(y.var) == TRUE
                                    | !!sym(y.var) == "Y"
                                    | !!sym(y.var) == "Yes", 
                                    1, 0),
                             levels = c(0,1)))
  y.var <- "positive"
  
  # check data completeness and actual positive rate
  total.rows <- nrow(df.reg)
  complete.rows <- nrow(df.reg %>% filter(complete.cases(.)))
  total.positive <- nrow(df.reg %>% 
                           {if(filter.complete) filter(., complete.cases(.)) else .} %>% 
                           filter(positive == 1)
                         )
  
  if(filter.complete) {
    message(paste0("Number of complete cases: ", complete.rows, " / ", total.rows, ", ", round(100*complete.rows/total.rows, 2), "%", "\n",
                   "Number of positives: ", total.positive, " / ", complete.rows, ", ", round(100*total.positive/complete.rows, 2), "%"))
  } else {
    message(paste0("Number of complete cases: ", complete.rows, " / ", total.rows, ", ", round(100*complete.rows/total.rows, 2), "%", "\n",
                   "Number of positives: ", total.positive, " / ", total.rows, ", ", round(100*total.positive/total.rows, 2), "%"))
  }

  # filter down to complete cases
  if(filter.complete) {
    df.reg <- df.reg %>% filter(complete.cases(.))
  }
  
  # train/test split
  do.split <- !is.null(train.pc) && train.pc < 1
  
  if (do.split) {
    set.seed(147)
    train.idx <- sample(seq_len(nrow(df.reg)), size = floor(train.pc * nrow(df.reg)))
    df.train <- df.reg[train.idx, ]
    df.test  <- df.reg[-train.idx, ]
    message(paste0("Train/test split applied: ", nrow(df.train), " train rows / ", nrow(df.test), " test rows"))
  } else {
    df.train <- df.reg
    df.test  <- NULL
    message("No train/test split applied. Using full dataset for training and diagnostics.")
  }
  
  # helper functions for ROC and f-beta
  ## build ROC plot
  build_roc_plot <- function(actual, pred, split.label) {
    pr  <- ROCR::prediction(pred, actual)
    prf <- ROCR::performance(pr, measure = "tpr", x.measure = "fpr")
    auc <<- ROCR::performance(pr, measure = "auc")@y.values[[1]]
    
    df.plt <- data.frame(
      `False positive rate` = unlist(prf@x.values),
      `True positive rate`  = unlist(prf@y.values),
      check.names = FALSE
    )
    
    ggplot(df.plt, aes(x = `False positive rate`, y = `True positive rate`)) +
      geom_line(color = "black", linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal() +
      labs(title = paste0("ROC curve, ", split.label, ". AUC: ", round(auc, 3)),
           subtitle = paste("Model:", model.name)) +
      theme(plot.title    = element_text(size = 10),
            plot.subtitle = element_text(size = 9))
  }
  
  ## calculate F-beta score
  f_beta_score <- function(y.true, y.prob, cutoff, beta) {
    y.hat     <- ifelse(y.prob >= cutoff, 1, 0)
    MLmetrics::FBeta_Score(y_true = y.true, 
                           y_pred = y.hat, 
                           positive = 1, 
                           beta = f.beta)
  }
  
  ## build F-beta plot
  build_fbeta_plot <- function(actual, pred, split.label, cutoff.train = NULL) {
    cutoffs <- seq(0, 1, length.out = 201)
    fb_df <- data.frame(
      cutoff = cutoffs,
      f.beta = sapply(cutoffs, function(x) f_beta_score(actual, pred, x, f.beta))
    )
    
    # if no cutoff supplied (i.e. train call), find the optimal and store it
    # if cutoff supplied (i.e. test call), use the train-derived cutoff instead
    if (is.null(cutoff.train)) {
      best.row <- fb_df[which.max(fb_df$f.beta), ]
      cutoff.train <<- best.row$cutoff
      active.cutoff <<- best.row$cutoff
    } else {
      active.cutoff <<- cutoff.train
    }
    
    active.fb <<- f_beta_score(actual, pred, active.cutoff, f.beta)
    cutoff.lbl <- sprintf("%.3f", active.cutoff)
    fb.lbl <- sprintf("%.3f", active.fb)
    
    ggplot(fb_df, aes(cutoff, f.beta)) +
      geom_line(color = "black", linewidth = 1) +
      geom_vline(xintercept = active.cutoff, linetype = "dashed", linewidth = 0.5, color = "blue") +
      geom_hline(yintercept = active.fb,     linetype = "dashed", linewidth = 0.5, color = "blue") +
      annotate("text", x = active.cutoff, y = 0,
               label = paste0("Cutoff = ", cutoff.lbl, "\nF-score = ", fb.lbl),
               vjust = 0, hjust = -0.1, color = "blue") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title    = paste0("F-", f.beta, " Score, ", split.label),
           subtitle = paste("Model:", model.name),
           x = "Probability Cutoff",
           y = paste0("F-", f.beta, " Score")) +
      theme_minimal() +
      theme(plot.title    = element_text(size = 10),
            plot.subtitle = element_text(size = 9))
  }
  
  
  ## confusion matrix at optimal cut-off
  confusion_optimal <- function(y.true, y.prob, cutoff, beta) {
    y.hat     <- ifelse(y.prob >= cutoff, 1, 0)
    TP        <- sum(y.hat == 1 & y.true == 1)
    TN        <- sum(y.hat == 0 & y.true == 0)
    FP        <- sum(y.hat == 1 & y.true == 0)
    FN        <- sum(y.hat == 0 & y.true == 1)
    cm <- data.frame(
      Actual = c("Positive", "Positive", "Negative", "Negative"),
      Prediction = c("Positive", "Negative", "Positive", "Negative"),
      N = c(TP, FN, FP, TN)
    ) %>% 
      as_tibble()
    return(cm)
  }
  
  ## confusion matrix plot
  confusion_plot <- function(cm, split.label) {
    # note: warning suppression due to package "ggimage" being needed for the plot to work in full
    # however, the plot in its current format is sufficient for this wrapper
    suppressWarnings(
      cvms::plot_confusion_matrix(cm, target_col = "Actual", prediction_col = "Prediction", counts_col = "N") +
        labs(title = paste0("Confusion matrix using optimal train F-beta cutoff, ", split.label),
             subtitle = paste("Model:", model.name)) +
        theme(plot.title    = element_text(size = 10),
              plot.subtitle = element_text(size = 9))
    )
  }
  
  # execute the regression on training data
  model <- glm(reformulate(x.vars, y.var), family = binomial(link = 'logit'), data = df.train)
  assign(paste0(model.name.short, ".model"), model, envir = parent.frame())
  
  # create 'tidy' results table
  res <- tidy(model)
  
  if (!gtools::invalid(cluster.se.var)) {
    res <- broom::tidy(lmtest::coeftest(model, vcov. = vcovCL(model, 
                                                              cluster = df.train[[cluster.se.var]], 
                                                              type = "HC0")),
                       conf.int = TRUE)
  } else if (std.error.corrected) {
    ci <- confint(model)
    ci <- as.data.frame(ci, row.names = FALSE) %>% 
      rename_with(~c("conf.low", "conf.high")) %>%
      mutate(ci.low = as.numeric(ci.low),
             ci.high = as.numeric(ci.high))
    res <- res %>% cbind(ci)
  } else {
    res <- res %>% 
      mutate(conf.low = estimate - 2*std.error,
             conf.high = estimate + 2*std.error)
  }
  
  res <- res %>%
    mutate(odds.ratio = exp(estimate),
           sig.level = case_when(p.value < 0.01 ~ "1% level",
                                 p.value < 0.05 ~ "5% level",
                                 p.value < 0.1  ~ "10% level",
                                 TRUE ~ "Not significant")) %>%
    mutate_if(is.numeric, ~round(., 3))
  res$odds.ratio.low  <- round(exp(res$conf.low), 3)
  res$odds.ratio.high <- round(exp(res$conf.high), 3)
  
  assign(paste0(model.name.short, ".res"), res, envir = parent.frame())
  message(paste0("Regression results, stored in dataframe '", model.name.short, ".res'. ", 
                 "Model stored as '", model.name.short, ".model'."))
  print(res)
  
  # model fit
  ## ANOVA. Which variables drive the model? If not significant, consider removing from model
  test.anova <- anova(model, test = "Chisq") %>% 
    mutate_if(is.numeric, ~round(., 3))
  anova.insignificant <- rownames(test.anova[test.anova$'Pr(>Chi)' > 0.05, ])
  if (is_empty(anova.insignificant)) {
    anova.insignificant <- "N/A"
  }
  message("Results from ANOVA test.")
  print(test.anova)
  
  ## Akaike Information Criteria. For comparing similar models. Lower = better.
  test.aic <- model$aic
  
  # generate predictions and plots
  pred.train  <- predict(model, newdata = df.train, type = "response")
  actual.train <- as.numeric(as.character(df.train[[y.var]]))
  
  plt.roc.train <- build_roc_plot(actual.train,  pred.train,  "training data")
  auc.train <- auc
  
  plt.fbeta.train <- build_fbeta_plot(actual.train, pred.train, "training data")
  fbeta.train <- active.fb
  
  cm.train <- confusion_optimal(actual.train,  pred.train, active.cutoff, best.fb)
  plt.cm.train <- confusion_plot(cm.train, "training data")
  
  # apply train/test split if needed
  if (do.split) {
    pred.test   <- predict(model, newdata = df.test, type = "response")
    actual.test <- as.numeric(as.character(df.test$positive))
    
    plt.roc.test <- build_roc_plot(actual.test,   pred.test,   "test data")
    auc.test <- auc
    
    plt.fbeta.test <- build_fbeta_plot(actual.test, pred.test, "Test", cutoff.train = cutoff.train)
    fbeta.test <- active.fb
    
    cm.test <- confusion_optimal(actual.test,  pred.test, active.cutoff, best.fb)
    plt.cm.test <- confusion_plot(cm.test, "test data")
    
    print(cowplot::plot_grid(plt.roc.train, plt.fbeta.train, plt.cm.train,
                    plt.roc.test, plt.fbeta.test, plt.cm.test,
                    nrow = 2, ncol = 3))
  } else {
    print(cowplot::plot_grid(plt.roc.train, plt.fbeta.train, plt.cm.train, nrow = 1))
  }
  
  # model diagnostics summary row
  test.results <- data.frame(
    'Model.Name' = model.name,
    'AIC' = round(test.aic, 3),
    'ANOVA.insignificant' = paste(unlist(anova.insignificant), collapse = ", "),
    'Train.AUC' = round(auc.train, 3),
    'Test.AUC' = ifelse(do.split, round(auc.test, 3), NA),
    'Beta' = f.beta,
    'Train.FBeta' = round(fbeta.train, 3),
    'Test.FBeta' = ifelse(do.split, round(fbeta.test, 3), NA),
    stringsAsFactors = FALSE
  )
  assign(paste0(model.name.short, ".test"), test.results, envir = parent.frame())
  message("Results from model fit tests, stored in dataframe '", paste0(model.name.short, ".test"), "'. ",
          "See Plots pane for ROC and F-beta results.")
  print(t(test.results))
}

# examples ----------------------------------------------------------------

## set up ----
options(scipen=999)

packages <- c(
  # utility
  "zoo", 
  # plotting
  "cowplot", "cvms", 
  # stats
  "tidymodels", "lmtest", "sandwich", "ROCR", "MLmetrics",
  # tidyverse
  "tidyverse"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
for (pack in packages) {library(pack, character.only = TRUE)}

patient_list <- readr::read_csv("frail_patients.csv") 

patient_list <- patient_list %>% 
  mutate(patient = as.character(patient),
         frailty_flag = factor(frailty_flag, levels = c(0,1)),
         sex = factor(sex, levels = c("F", "M")),
         imd_decile = factor(imd_decile, levels = sort(unique(patient_list$imd_decile))),
         ethnicity = factor(ethnicity, levels = sort(unique(patient_list$ethnicity))),
         geography = factor(geography, levels = sort(unique(patient_list$geography)))
         )

## running regressions ----
reg_logit(df = patient_list,
          y.var = "frailty_flag",
          x.vars = c("age", "sex", "imd_decile", "ethnicity", "n_chronic_cond", 
                     "n_medication", "n_ae_atts", "n_pc_cons", "geography"),
          model.name = "Frailty logit model",
          model.name.short = "logit.frail",
          cluster.se.var = "geography",
          filter.complete = TRUE,
          std.error.corrected = FALSE,
          f.beta = 2,
          train.pc = 0.8)

reg_logit(df = patient_list,
          y.var = "frailty_flag",
          x.vars = c("age"),
          model.name = "Under-specified model",
          model.name.short = "logit.underspec",
          cluster.se.var = NULL,
          std.error.corrected = FALSE,
          f.beta = 2,
          train.pc = 0.8)
