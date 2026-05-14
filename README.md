# Regression models
This project contains wrapper functions for regression models. It is designed to standardise outputs across projects. Formula construction, model fitting, assumption testing, and plot generation all happen in a single call. Outputs (model objects, tidy results, and plots) are assigned directly to the calling environment.

The sections below list out the functions available. The pre-requisite packages are packages that need to be called for the functions to work. **These packages are needed for all models: tidyverse, scales, cowplot, and gtools.**

Models included in this project:
 - [Ordinary Least Squares (OLS)](#ols)
 - [Logistic regression (logit)](#logistic-regression)
 - [XGBoost](#xgboost)
 - [Propensity score matching (PSM)](#propensity-score-matching)
 - [Difference-in-differences (DID)](#difference-in-differences)
 - [Doubly-robust difference-in-differences (DRDID)](doubly-robust-difference-in-differences)
 - [Hazard models (Kaplan-Meier and Cox proportional hazards)](#hazard-models)

Additional models to be added:
 - ARIMAX
 - Short-term time-series forecasting (e.g., Prophet)
 - Dynamic panel data approaches
 - Poisson

---

## OLS

See script `reg_ols`. This holds a wrapper function for a basic OLS regression, ideally for cross-sectional data. This is built around the function `stats::lm`.

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages(c("lmtest", "sandwich", "tseries"))
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for dependent variable |
| `x.vars` | character vector | — | Columns for independent variable(s) |
| `model.name` | string | `"OLS regression model"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"ols.reg"` | Prefix for output objects assigned to parent environment |
| `cluster.se.var` | string | `NULL` | Column name for clustering variable for robust standard errors |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`lm`) |
| `[prefix].res` | Tidy results table |
| `[prefix].test` | Vector of diagnostic test results |

The following model diagnostics are performed:
- R-squared
- Root mean squared error (RMSE)
- Akaike Information Criteria (AIC)
- ANOVA test of significance of independent variables
- Breusch–Godfrey correlated residuals test
- Plot fitted vs actual. Test quality of fit using a Mincer-Zarnowitz test (normally reserved for time-series forecasting, but used as general sensecheck)
- QQ plot test for normally distributed residuals. Accompanied by Jacques-Bera test p-value
- Plot of residuals vs fitted values
- Plot of Cook's distances

Note that all plots are amalgamated into a single plot in the console using cowplot.


### Performance and Limitations

#### Robust standard errors
There is no test for heteroskedastic errors. Instead, the wrapper automatically calculates robust standard errors using `sandwich::vcovHC` or `sandwich::vcovCL` (if clustered) with type set to HC0 (White's estimator).

#### Train/test split
No train/test split option is available in this model. It is meant purely for guidance on the whole dataset, but can be adapted to apply to split data.

### Future developments

- Apply test/train split
- A wrapper for time-series specific data. To include additional checks (e.g., Ljung-Box Q test)

</details>

---

## Logistic regression

See script `reg_logit`. This holds a wrapper function for a basic logit regression. The model uses AUC and F-beta as the main assessment criteria. This is built around the function `stats::glm`.

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages(c("lmtest", "sandwich", "cvms", "ROCR", "MLmetrics"))
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for binary dependent variable |
| `x.vars` | character vector | — | Columns for independent variable(s) |
| `model.name` | string | `"Logit regression model"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"logit.reg"` | Prefix for output objects assigned to parent environment |
| `filter.complete` | logical | `TRUE` | Filter out all NAs. Recommend keeping as TRUE, as the model may fail to predict values otherwise |
| `cluster.se.var` | string | `NULL` | Column name for clustering variable for robust standard errors |
| `std.error.corrected` | logical | `FALSE` | For non-clustered standard errors, re-estimate corrected confidence intervals |
| `f.beta` | numeric | `1` | Beta value for f-beta scores |
| `train.pc` | numeric | `1` | Percentage of data for training dataset (vs testing). Set to 1 if no split desired |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`glm`) |
| `[prefix].res` | Tidy results table |
| `[prefix].test` | Vector of diagnostic test results |

The following model diagnostics are performed, done for a train/test split unless stated otherwise:
- Akaike Information Criteria (AIC), training data only
- ANOVA test of significance of independent variables, training data only
- ROC curve and associated AUC
- F-beta scores at 201 probability cutoffs from 0 to 1 inclusive
- Confusion matrix at the optimal F-beta/probability cutoff calculated on the training data

Note that all plots are amalgamated into a single plot in the console using cowplot.


### Performance and Limitations

#### Model assessment criteria
The wrapper is focused on AUC and F-beta to assess the models. F-beta is in particular has been selected, as in healthcare settings we regularly assess cases where there is high class imbalance in the data with fewer positives than negatives. See the [Wiki for F-score](https://en.wikipedia.org/wiki/F-score) for more. This means the confusion matrix is presented at the optimal probability cutoff to maximise the F-beta score rather than the usual 0.5. This can be changed in the wrapper by changing the helper function `confusion_optimal` by setting `cutoff` to 0.5 and changing the title in the helper function `confusion_plot`.

#### Corrected standard errors
If setting `std.error.corrected = TRUE`, note that this will slow down the performance of the model. We recommend first running multiple models to decide the best specification based on the assessment criteria, then re-run the "best" model with corrected standard errors.

### Future developments

- Mixed-effects logit regression (to resolve slow performance issues)

</details>
 
---

## XGBoost

See script `reg_xgboost`. This holds a wrapper function for a binary XGBoost. This is a form of gradient boosted tree-based methods which often outperforms others in its class for classification exercises. This is built around the package `xgboost`. Full documentation can be found [here](https://xgboost.readthedocs.io/). The model uses AUC and F-beta as the main assessment criteria.
Currently, the wrapper only works for binary classification problems.

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages(c("xgboost", "cvms", "ROCR", "MLmetrics"))
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for dependent variable |
| `x.vars` | character vector | — | Columns for independent variable(s) |
| `model.name` | string | `"XGBoost model"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"xgb"` | Prefix for output objects assigned to parent environment |
| `filter.complete` | logical | `FALSE` | Filter out all NAs. XGBoost is robust to a small number of NAs, so default is FALSE |
| `nthreads ` | integer | `parallel::detectCores() - 1` | Number of processors used for parallel calculations |
| Various tuning parameters: `nrounds`, `max.depth`, `learning.rate`, `reg.lambda`, `min.split.loss`, `subsample`, `colsample.bytree`, `scale.pos.weight` | numeric | `nrounds = 100`, `max.depth = 6`, `learning.rate = 0.5`, `reg.lambda = 1`, `min.split.loss = 0`, `subsample = 1`, `colsample.bytree = 1`, `scale.pos.weight = 1` | For detailed understanding, [see online documentation]( https://xgboost.readthedocs.io/en/latest/parameter.html #parameters-for-tree-booster) |
| `f.beta` | numeric | `1` | Beta value for f-beta scores |
| `train.pc` | numeric | `1` | Percentage of data for training dataset (vs testing). Set to 1 if no split desired |
| `importance.max ` | integer | `20` | The maximum number of bars to show for feature importance graph. Selects top 'x' based on gain |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`xgb.Booster`) |
| `[prefix].importance` | Table of feature importance using metrics cover, frequency, and gain. The top ‘x’ are included in the plots pane |
| `[prefix].test` | Vector of diagnostic test results |

The following model diagnostics are performed, done for a train/test split unless stated otherwise:
- ROC curve and associated AUC
- F-beta scores at 201 probability cutoffs from 0 to 1 inclusive
- Confusion matrix at the optimal F-beta/probability cutoff calculated on the training data

Note that all plots are amalgamated into a single plot in the console using cowplot.

### Performance and Limitations

#### Model tuning
There are no hard and fast rules for setting tuning parameters, as each model will behave differently depending on the underlying data. The [package documentation]( https://xgboost.readthedocs.io/en/stable/tutorials/param_tuning.html) provides a short discussion on the matter. You can consider using grid search and cross-validation to find optimal parameters. But in general, getting the right features in `x.vars` will give you greater gains than tweaking hyperparameters, especially when introducing new data.

#### Model assessment criteria
The wrapper is focused on AUC and F-beta to assess the models. F-beta is in particular has been selected, as in healthcare settings we regularly assess cases where there is high class imbalance in the data with fewer positives than negatives. See the [Wiki for F-score](https://en.wikipedia.org/wiki/F-score) for more. This means the confusion matrix is presented at the optimal probability cutoff to maximise the F-beta score rather than the usual 0.5. This can be changed in the wrapper by changing the helper function `confusion_optimal` by setting `cutoff` to 0.5 and changing the title in the helper function `confusion_plot`.

#### SHAP scores
SHAP is often preferred to other feature importance metrics like “gain” because it is more consistent and can be calculated for each observation. It is possible to produce SHAP scores for XGBoost, but we have had trouble implementing packages `SHAPforxgboost` and `shapviz` into this wrapper. 

### Future developments

- Non-binary classification
- SHAP scores

</details>
 
---

## Propensity score matching

See script `reg_psm`. This holds a wrapper function for propensity score matching. This is built around the package `MatchIt`.

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages(c("cobalt", "MatchIt"))
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for dependent variable |
| `x.vars` | character vector | — | Columns for independent variable(s) |
| `model.name` | string | `"Propensity score matching"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"psm"` | Prefix for output objects assigned to parent environment |
| `method` | string | `"nearest"` | Method for finding matches.  Most common are "nearest" for nearest neighbour, "optimal" for optimal pair matching, "full" for optimal full matching. See ?MatchIt::matchit for more options. |
| `distance` | string | `"logit"` | How to calculate distances to matches. See ?MatchIt::matchit for more options. Choosing some options may prompt installation of extra packages. |
| `ratio` | integer | `1` | Number of controls to match to one treated |
| `replace` | logical | `FALSE` | Sample controls with (TRUE) or without (FALSE) replacement |
| `caliper` | numeric | `NULL` | Caliper to apply to match distances. Based on standard deviation units |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`matchit`) |
| `[prefix].dist.model` | Where a parametric model is selected in `distance` ("glm", "logit", or “gam”), the model coefficients and p-values are outputted to a dataframe |
| `[prefix].matched` | A dataframe with the matched units. This is the original dataframe passed through the wrapper, plus three columns: `id` for row number in the original dataframe, `subclass` for a group number, and `weights` which equals 1 for treated unit and 1/ratio for matched controls |

The following model diagnostics are performed:
- Standardised mean differences (SMD) / covariate balance. For each independent variable/factor, we plot the SMD for the unmatched and matched data. A target threshold of +/-0.2 is shown on the chart, though SMDs of up to 0.5 might be acceptable.
- Propensity score overlap. Only where `distance` is one of "glm", "logit", or “gam”. Presents density plots of propensity scores for treatment and control groups, before and after matching.

Note that all plots are amalgamated into a single plot in the console using cowplot.


### Performance and Limitations

#### Difference-in-differences calculations
The main output from this process is getting the matches. It does not include an option to calculate, for example, difference-in-differences. However, the output `[prefix].matched` noted above is designed to be used directly in these sorts of calculations. There are additional wrapper for this.

#### Distance options
`distance` includes machine learning approaches to find matching controls. Though this might improve the quality of matches, it has two limitations. First, these algorithms may run more slowly than parametric options like logit. Second, the MatchIt package does not provide indicators of the most influential covariates (e.g., SHAP), so you will not be able to discern what is driving the matches. In most cases, “logit” should be sufficient.

### Future developments

- Other clustering algorithms not centred around treatment variables
- Additional difference-in-differences calculations

</details>
 
---

## Difference-in-differences

See script `reg_did`. This holds a wrapper function for difference-in-differences regressions. This includes flexibility for group- and time-fixed effects. Where fixed effects are applied, this is built around the package `plm`. Else, it simply uses OLS via `stats::lm`. 

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages("lmtest", "sandwich", "plm")
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for dependent variable |
| `treat.field` | string | — | Column name containing binary treatment flag. |
| `time.field` | string | — | Column name containing time periods. Ideally binary: 0 for untreated period, 1 for treated period. |
| `x.vars` | character vector | `NULL` | Optional columns for covariates |
| `model.name` | string | `"Difference-in-differences model"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"did.reg"` | Prefix for output objects assigned to parent environment |
| `fixed.effects` | string | `"none"` | Apply fixed effects. Choice of 'none' (default), 'individual', 'time', 'twoways', or 'nested' |
| `group.field` | string | `NULL` | If applying group-level fixed effects, specify which field represents the group |
| `fe.method` | string | `"within"` | If applying group-level fixed effects, what method to use. See ?plm::plm for more |
| `cluster.se.var` | string | `NULL` | Column name for clustering variable for robust standard errors. If applying fixed-effects, package plm will automatically select group.field |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`lm` without fixed-effects or `plm` with fixed effects) |
| `[prefix].res` | A dataframe with estimate and confidence intervals. DID estimator is labelled and printed in the console |
| `[prefix].test` | Vector of diagnostic test results |

The following model diagnostics are performed:
- Breusch–Godfrey correlated residuals test
- Pesaran’s cross-sectional dependence (CD) test, if fixed effects are applied

### Performance and Limitations
##### Model diagnostics
There is no direct test of the parallel trends assumption. If you have multiple pre-treatment periods, we recommend plotting means to see changes over time.
Also be sure to consider other issues with treatment design:
- Selection biases (e.g., the patient self-selecting or clinical choice)
- Lack of adherence to treatment being investigated
- Unobserved confounding factors (e.g., other actions taken by either group of patients impacting their health)
- Ensure enough periods post-treatment have occurred, especially where seasonality is a factor (even though this wrapper allows for time-fixed effects)

#### Multiple treatment periods
We recommend only using this for simple 2x2 designs. For multiple time periods, especially with staggered treatment timing, we recommend using the doubly-robust approach in the next section.

### Future developments

- Other clustering algorithms not centred around treatment variables
- Additional difference-in-differences calculations (dynamic panels for long time-series)

</details>

---

## Doubly robust difference-in-differences

See script `reg_drdid`. This holds a wrapper function for doubly-robust difference-in-differences regression for group-time averaged treatment effects. This is preferred to the "vanilla" DID model if there is staggered treatment timing. This is built around the package `did`. For background, see [the package documentation](https://bcallaway11.github.io/did/reference/att_gt.html) and [the original paper](https://www.sciencedirect.com/science/article/abs/pii/S0304407620303948).

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages("did")
```

### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `y.var` | string | — | Column name for dependent variable |
| `x.vars` | character vector | — | Columns for independent variable(s) |
| `model.name` | string | `"Doubly-robust difference-in-differences model"` | Model name to be included in the title of plots and diagnostics tables |
| `model.name.short` | string | `"drdid.reg"` | Prefix for output objects assigned to parent environment |
| `time.field` | string | — | Field containing time period. Should be integers, ideally consecutive |
| `id.field` | string | — | Field containing unit IDs. Ideally consecutive integers |
| `treat.period` | string | — | Field stating which period (in time.field) the first treatment took place. 0 if untreated |
| `cluster.var ` | string | `”id.field”` | Field to cluster standard errors by. Default is `id.field` |
| `method` | string | `"dr"` | Which method to use. Default is "dr" for doubly robust. See ?did::att_gt for more |
| `unbalanced ` | logical | `FALSE` | Allow for an unbalanced panel (e.g.,, attrition). Setting to FALSE will filter for complete cases across all periods |
| `control.group` | string | `"nevertreated"` | What to use as control group: “nevertreated” or “notyettreated”. See ?did::att_gt for more |
| `anticipation` | integer | `0` | The number of time periods before participating in the treatment where units can anticipate participating in the treatment. If unsure or expected to conform to treatment, use default of 0 |
| `alpha` | numeric | `0.05` | Set alpha for confidence intervals |
| `n.iter` | integer | `1000` | Number of iterations for producing bootstrap confidence intervals Default = 1,000, but reduce if working with very large dataframes |
| `nthreads ` | integer | `parallel::detectCores() - 1` | Number of processors used for parallel calculations for bootstrap confidence intervals |

### Outputs

The wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`MP`) |
| `[prefix].res` | A dataframe with estimate and confidence intervals for average treatment effects for each time-period (in `time.field`) and group (in `first.treat`) |
| `[prefix].plt` | A plot of the results above |


### Performance and Limitations
##### Model diagnostics
No model diagnostics are performed, although pre-treatment parallel trends assumption test is effectively included within the results table. Simply ensure that negative time periods are not statistically significantly different from zero.
We recommend before running the algorithm to check whether inverse propensity score matching is suitable for your data. If not, the “doubly robust” nature of the model will instead only rely on regression model and control group being correctly specified.
Also be sure to consider other issues with treatment design:
- Selection biases (e.g., the patient self-selecting or clinical choice)
- Lack of adherence to treatment being investigated
- Unobserved confounding factors (e.g., other actions taken by either group of patients impacting their health)
- Ensure enough periods post-treatment have occurred, especially where seasonality is a factor (even though this is an implementation of two-way fixed effects)

#### Unbalanced panels
The drdid function allows for unbalanced panels. However, the calculation takes much longer. Consider whether the model would still be robust if filtering out observations not present in all periods.

### Future developments

- Attempt to unpick the underlying model in the package to add Breusch-Godfrey and Pesaran's CD tests

</details>
 
---

## Hazard models

See script `reg_hazard` and data patient_mortality.csv. This holds two wrapper functions for two common survival analysis models: Kaplan-Meier and Cox proportional hazards. These are built around the packages `survival` and `survminer`.

<details>
 
### Additional prerequisites

In addition to those mentioned above, install the following required packages from CRAN:

```r
install.packages(c("survival", "survminer"))
```

### Parameters

#### Shared parameters (`reg_km`, `reg_cox`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `df` | data.frame | — | Input dataset |
| `period.start` | string | — | Column name for entry/start time |
| `period.stop` | string | — | Column name for exit/stop time |
| `period.event` | string | — | Column name for event indicator (1 = event, 0 = censored) |
| `strata` | character vector | `NULL` | Grouping variables for stratified analysis. If `NULL`, fits overall curve |
| `model.name.short` | string | `"surv.km"` | Prefix for output objects assigned to parent environment |
| `conf.int` | logical | `TRUE` | Show confidence intervals on plot |
| `risk.table` | logical | `FALSE` | Show risk table below plot |
| `ncensor.plot` | logical | `FALSE` | Show censoring plot |

#### `reg_cox` only

| Parameter | Type | Default | Description |
|---|---|---|---|
| `covariates` | character vector | `NULL` | Adjustment covariates. Must not overlap with `strata` |
| `cluster.var` | string | `NULL` | Column name for clustering variable for robust standard errors |

### Outputs

Each wrapper assigns the following objects to the calling environment, prefixed by `model.name.short`:

| Object | Description |
|---|---|
| `[prefix].model` | Fitted model object (`survfit` or `coxph`) |
| `[prefix].res` | Tidy results table |
| `[prefix].plt` | `ggsurvplot` plot object |

For `reg_cox`, additional diagnostics are printed to the console:
- Concordance (Weighted by 'N' and Uno's C statistic weighted by 'N/G^2')
- Proportional hazards test results (`cox.zph`), including variables failing at p < 0.1

### Performance and Limitations

#### Counting process format
Both wrappers use the three-argument `Surv(start, stop, event)` counting process format, which supports delayed entry and time-varying covariates. Note that `survdiff()` does not support this format, which is why the log-rank test in `reg_km()` uses only the stop time. Entry times are not accounted for in that test.

#### Proportional hazards
`reg_cox()` tests the PH assumption via `cox.zph()` and flags violations at p < 0.1. This should be interpreted alongside the Schoenfeld residual plots. Variables flagged as violating PH may require a time interaction term or stratification in the model. Where there are more than 1,000 rows of data, Schoenfeld residual plots will be produced in base R instead of the survival package function to speed up performance.

#### Concordance
Uno's C statistic may be used in preference to the standard concordance index, as it is more robust when patients have short observed time windows and many pairs are therefore incomparable. It is implemented via `survival::concordance(..., timewt = "uno")`.

#### Clustering
When `cluster.var` is specified, robust (sandwich) standard errors are used. This accounts for within-cluster correlation (e.g. repeated patient episodes) but does not model the clustering structure explicitly.


### Future developments

- Warnings about deprecated ggplot functions. These occur when producing plots that are parts of the survival packages and cannot be addressed by the wrapper directly.
- Gompertz survival wrapper
- Automatic PH remediation suggestions (time interaction or stratification)

</details>
