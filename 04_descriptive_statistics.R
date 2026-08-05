# Cannabis Use and Everyday Functioning Among People Living With HIV
# Script: 04_descriptive_statistics.R
# Purpose: 04 Descriptive Statistics


# ============================================================
# ADL PAPER: OUTCOME-SELECTION AND MODEL-BUILDING SCRIPT
# ============================================================

## Purpose:
## This script helps decide which ADL / functional variable should be the
## primary outcome before committing to final cannabis + HIV disease-severity models.
##
## The script:
## 1. Loads the HIV+ dataset.
## 2. Recreates shared candidacy variables.
## 3. Creates ADL candidate outcome versions.
## 4. Checks distributions, zero/floor effects, and correlations.
## 5. Tests whether the ADL denominator may be biased.
## 6. Screens candidate ADL outcomes against demographic/clinical variables.
## 7. Runs candidate cannabis/disease-severity models after outcome checks.
##
## Important:
## This script creates new variables and does not intentionally overwrite raw variables.


# ============================================================
# 0. PACKAGES
# ============================================================

## Install packages manually if needed.
## I am not including automatic install.packages() inside the script because
## package version changes can make dissertation analyses harder to reproduce.

library(haven)
library(labelled)
library(lavaan)
library(lmtest)
library(sandwich)
library(psych)
library(dplyr)
library(tidyr)
library(broom)
library(emmeans)
library(ggplot2)
library(effectsize)
library(performance)
library(parameters)
library(car)


# ============================================================
# 1. LOAD DATA
# ============================================================

## Update this path if your .sav file is stored somewhere else.
## Do not name the object "data" because data() is also a base R function.

adl_raw <- read_sav("data/HIV+113.sav")

## Create a working copy. Raw imported data remain untouched.
adl_work <- adl_raw

## Optional view.
View(adl_raw)


# ============================================================
# 2. NON-OVERWRITE CHECK
# ============================================================

## This protects your raw dataset. If any planned derived variable already
## exists in the source dataset, the script stops instead of overwriting it.

planned_derived_vars <- c(
  "mst_ct_rev",
  "gpt_time_rev",
  "wcpe_rev",
  "bvtr_z",
  "hvtt_z",
  "hvdr_z",
  "bvdr_z",
  "mstct_z",
  "gpttim_z",
  "wcpe_z",
  "igtnt_z",
  "cd_z",
  "ss_z",
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z",
  "num_tests_available",
  "global_z",
  "cd4sqrt",
  "nadirsqrt",
  "log_vl2",
  "DiseaseSev",
  "DiseaseSev_c",
  "cd4sqrt_c",
  "nadirsqrt_c",
  "log_vl2_c",
  "du_mar4_12m_aBin",
  "du_mar4_12m_aBin_ord",
  "du_mar6_30d_aBin",
  "du_mar6_30d_aBin_ord",
  "du_mar6_30d_aBin30_lowincl",
  "du_mar6_30d_aBin30_lowincl_ord",
  "du_mar2_life_aBin",
  "du_mar2_life_aBin_ord",
  "du_mar2_life_aBin_q75",
  "du_mar2_life_aBin_q75_ord",
  "du_mar2_life_a_log1p",
  "phq_4_latin_num",
  "phq_5_race_num",
  "race_eth_binary",
  "race_eth_binary_covnum",
  "race_eth_binary_covfac",
  "sex",
  "sex_covnum",
  "sex_covfac",
  "phq_2_age_c",
  "phq_7_degree_covfac",
  "phq_7_degree_num",
  "phq_7_degree_c",
  "du_alc4_12m_a_log1p",
  "du_alc6_30d_a_log1p",
  "du_alc4_12m_a_log1p_c",
  "du_alc6_30d_a_log1p_c",
  "adl_total_any",
  "adl_total_any_num",
  "adl_total_positive",
  "adl_total_log1p",
  "adl_b_sum_now_log1p",
  "mac_misr_log1p"
)

existing_planned_vars <- intersect(planned_derived_vars, names(adl_work))

if (length(existing_planned_vars) > 0) {
  stop(
    paste0(
      "Non-overwrite rule triggered. These planned derived variable names ",
      "already exist in the source dataset: ",
      paste(existing_planned_vars, collapse = ", "),
      ". Rename the derived variables before continuing."
    )
  )
}


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

## Check whether expected variables exist in the dataset.

check_vars <- function(var_names, dataset = adl_work) {
  missing_vars <- setdiff(var_names, names(dataset))
  if (length(missing_vars) > 0) {
    message("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  invisible(missing_vars)
}


## Find variables by partial name.

find_vars <- function(pattern, dataset = adl_work) {
  grep(pattern, names(dataset), ignore.case = TRUE, value = TRUE)
}


## Center a variable without scaling it.

center_only <- function(x) {
  as.numeric(scale(x, center = TRUE, scale = FALSE))
}


## Run linear model with HC3 robust standard errors.

robust_lm <- function(formula, dataset = adl_work, vcov_type = "HC3") {
  model <- lm(formula, data = dataset)
  
  robust_results <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = vcov_type)
  )
  
  list(
    formula = formula,
    n = nobs(model),
    model = model,
    robust_results = robust_results
  )
}


## Extract robust linear model results into a clean table.

extract_robust_lm <- function(model_list, model_set_name) {
  output <- data.frame()
  
  for (i in seq_along(model_list)) {
    robust_i <- model_list[[i]]$robust_results
    
    temp <- data.frame(
      model_set = model_set_name,
      outcome = model_list[[i]]$outcome,
      n = model_list[[i]]$n,
      term = rownames(robust_i),
      estimate = robust_i[, "Estimate"],
      robust_se = robust_i[, "Std. Error"],
      t_value = robust_i[, "t value"],
      p_value = robust_i[, "Pr(>|t|)"],
      row.names = NULL
    )
    
    output <- rbind(output, temp)
  }
  
  output
}


## Run the same model structure across multiple outcomes.

run_outcome_models <- function(outcomes, rhs, model_set_name, dataset = adl_work) {
  
  model_list <- lapply(outcomes, function(y) {
    f <- as.formula(paste(y, "~", rhs))
    result <- robust_lm(f, dataset = dataset)
    result$outcome <- y
    result
  })
  
  list(
    models = model_list,
    table = extract_robust_lm(model_list, model_set_name)
  )
}


## Spearman correlation block.
## Useful for non-normal, ordinal, skewed, or zero-heavy variables.

spearman_block <- function(row_vars, col_vars, dataset = adl_work) {
  
  vars_needed <- unique(c(row_vars, col_vars))
  check_vars(vars_needed, dataset)
  
  dat <- dataset[, vars_needed]
  
  corr <- psych::corr.test(
    dat,
    method = "spearman",
    use = "pairwise.complete.obs",
    adjust = "none"
  )
  
  r <- corr$r[row_vars, col_vars, drop = FALSE]
  p <- corr$p[row_vars, col_vars, drop = FALSE]
  
  table <- as.data.frame(as.table(r))
  names(table) <- c("row_variable", "column_variable", "rho")
  table$p_value <- as.vector(p)
  table <- table %>% arrange(p_value)
  
  list(r = r, p = p, table = table)
}


# ============================================================
# 4. SHARED CANDIDACY VARIABLES
# ============================================================

# ------------------------------------------------------------
# 4a. Neurocognitive composites retained for background checks
# ------------------------------------------------------------

## These are not the primary ADL outcomes, but they are retained because
## your ADL paper may need to check whether functional outcomes relate to cognition.

adl_work$mst_ct_rev <- -1 * adl_work$mst_s_ct1
adl_work$gpt_time_rev <- -1 * adl_work$gpt_nd_ttime
adl_work$wcpe_rev <- -1 * adl_work$wc_pe_rs

adl_work$bvtr_z   <- as.numeric(scale(adl_work$bv_tr_rs))
adl_work$hvtt_z   <- as.numeric(scale(adl_work$hv_tt_rs))
adl_work$hvdr_z   <- as.numeric(scale(adl_work$hv_dr_rs))
adl_work$bvdr_z   <- as.numeric(scale(adl_work$bv_dr_rs))
adl_work$mstct_z  <- as.numeric(scale(adl_work$mst_ct_rev))
adl_work$gpttim_z <- as.numeric(scale(adl_work$gpt_time_rev))
adl_work$wcpe_z   <- as.numeric(scale(adl_work$wcpe_rev))
adl_work$igtnt_z  <- as.numeric(scale(adl_work$igt_nt_rs))
adl_work$cd_z     <- as.numeric(scale(adl_work$wai_cd_rs))
adl_work$ss_z     <- as.numeric(scale(adl_work$wai_ss_rs))

test_z_vars <- c(
  "bvtr_z", "hvtt_z", "hvdr_z", "bvdr_z",
  "mstct_z", "gpttim_z", "wcpe_z", "igtnt_z",
  "cd_z", "ss_z"
)

adl_work$learning_z <- rowMeans(adl_work[, c("bvtr_z", "hvtt_z")], na.rm = TRUE)
adl_work$memory_z <- rowMeans(adl_work[, c("hvdr_z", "bvdr_z")], na.rm = TRUE)
adl_work$motor_z <- rowMeans(adl_work[, c("mstct_z", "gpttim_z")], na.rm = TRUE)
adl_work$executive_z <- rowMeans(adl_work[, c("wcpe_z", "igtnt_z")], na.rm = TRUE)
adl_work$processing_speed_z <- rowMeans(adl_work[, c("cd_z", "ss_z")], na.rm = TRUE)

adl_work$num_tests_available <- rowSums(!is.na(adl_work[, test_z_vars]))

adl_work$global_z <- rowMeans(adl_work[, test_z_vars], na.rm = TRUE)
adl_work$global_z[adl_work$num_tests_available < 9] <- NA

neurocog_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "motor_z",
  "executive_z",
  "processing_speed_z"
)

summary(adl_work[, neurocog_vars])
colSums(is.na(adl_work[, neurocog_vars]))


# ------------------------------------------------------------
# 4b. HIV disease-severity factor
# ------------------------------------------------------------

## Higher DiseaseSev should represent greater HIV disease burden.
## Current CD4 and nadir CD4 are reverse-scored through negative square roots.
## Viral load is log-transformed.

adl_work$cd4sqrt <- -sqrt(adl_work$lab_thel)
adl_work$nadirsqrt <- -sqrt(adl_work$mhq_4_cd4_lowest)
adl_work$log_vl2 <- log1p(adl_work$lab_hiv)

DiseaseSev_model <- '
  DiseaseSev =~ nadirsqrt + cd4sqrt + log_vl2
  cd4sqrt ~~ 0.001*cd4sqrt
'

fit <- lavaan::cfa(
  DiseaseSev_model,
  data = adl_work,
  estimator = "MLR",
  missing = "fiml"
)

summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

adl_work$DiseaseSev <- as.numeric(lavaan::lavPredict(fit))
adl_work$DiseaseSev_c <- center_only(adl_work$DiseaseSev)

adl_work$cd4sqrt_c <- center_only(adl_work$cd4sqrt)
adl_work$nadirsqrt_c <- center_only(adl_work$nadirsqrt)
adl_work$log_vl2_c <- center_only(adl_work$log_vl2)

disease_marker_cor <- cor(
  adl_work[, c("cd4sqrt", "nadirsqrt", "log_vl2", "DiseaseSev")],
  use = "pairwise.complete.obs"
)

round(disease_marker_cor, 3)


# ------------------------------------------------------------
# 4c. Cannabis bins from candidacy syntax
# ------------------------------------------------------------

## Past-year cannabis is the primary cannabis exposure.
## None / low / high are based on zero use and the median among nonzero users.

nonzero_vals_12m <- adl_work$du_mar4_12m_a[adl_work$du_mar4_12m_a > 0]
median_12m <- median(nonzero_vals_12m, na.rm = TRUE)

adl_work$du_mar4_12m_aBin <- NA_real_
adl_work$du_mar4_12m_aBin[adl_work$du_mar4_12m_a == 0] <- 0
adl_work$du_mar4_12m_aBin[
  adl_work$du_mar4_12m_a > 0 &
    adl_work$du_mar4_12m_a <= median_12m
] <- 1
adl_work$du_mar4_12m_aBin[adl_work$du_mar4_12m_a > median_12m] <- 2

adl_work$du_mar4_12m_aBin_ord <- ordered(
  adl_work$du_mar4_12m_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low", "high")
)

contrasts(adl_work$du_mar4_12m_aBin_ord) <- contr.poly(3)

table(adl_work$du_mar4_12m_aBin_ord, useNA = "ifany")


## Past-30-day cannabis: median split sensitivity.

nonzero_vals_30d <- adl_work$du_mar6_30d_a[adl_work$du_mar6_30d_a > 0]
median_30d <- median(nonzero_vals_30d, na.rm = TRUE)

adl_work$du_mar6_30d_aBin <- NA_real_
adl_work$du_mar6_30d_aBin[adl_work$du_mar6_30d_a == 0] <- 0
adl_work$du_mar6_30d_aBin[
  adl_work$du_mar6_30d_a > 0 &
    adl_work$du_mar6_30d_a <= median_30d
] <- 1
adl_work$du_mar6_30d_aBin[adl_work$du_mar6_30d_a > median_30d] <- 2

adl_work$du_mar6_30d_aBin_ord <- ordered(
  adl_work$du_mar6_30d_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low_recent", "high_recent")
)

contrasts(adl_work$du_mar6_30d_aBin_ord) <- contr.poly(3)

table(adl_work$du_mar6_30d_aBin_ord, useNA = "ifany")


## Past-30-day cannabis: 30g threshold sensitivity.

adl_work$du_mar6_30d_aBin30_lowincl <- NA_real_
adl_work$du_mar6_30d_aBin30_lowincl[adl_work$du_mar6_30d_a == 0] <- 0
adl_work$du_mar6_30d_aBin30_lowincl[
  adl_work$du_mar6_30d_a > 0 &
    adl_work$du_mar6_30d_a <= 30
] <- 1
adl_work$du_mar6_30d_aBin30_lowincl[adl_work$du_mar6_30d_a > 30] <- 2

adl_work$du_mar6_30d_aBin30_lowincl_ord <- ordered(
  adl_work$du_mar6_30d_aBin30_lowincl,
  levels = c(0, 1, 2),
  labels = c("none", "lower_recent_0to30", "heavier_recent_over30")
)

contrasts(adl_work$du_mar6_30d_aBin30_lowincl_ord) <- contr.poly(3)

table(adl_work$du_mar6_30d_aBin30_lowincl_ord, useNA = "ifany")


## Lifetime cannabis: median split sensitivity.

nonzero_vals_life <- adl_work$du_mar2_life_a[adl_work$du_mar2_life_a > 0]
median_life <- median(nonzero_vals_life, na.rm = TRUE)

adl_work$du_mar2_life_aBin <- NA_real_
adl_work$du_mar2_life_aBin[adl_work$du_mar2_life_a == 0] <- 0
adl_work$du_mar2_life_aBin[
  adl_work$du_mar2_life_a > 0 &
    adl_work$du_mar2_life_a <= median_life
] <- 1
adl_work$du_mar2_life_aBin[adl_work$du_mar2_life_a > median_life] <- 2

adl_work$du_mar2_life_aBin_ord <- ordered(
  adl_work$du_mar2_life_aBin,
  levels = c(0, 1, 2),
  labels = c("none", "low_lifetime", "high_lifetime")
)

contrasts(adl_work$du_mar2_life_aBin_ord) <- contr.poly(3)

table(adl_work$du_mar2_life_aBin_ord, useNA = "ifany")


## Lifetime cannabis: q75 split sensitivity.

q75_life <- quantile(
  adl_work$du_mar2_life_a[adl_work$du_mar2_life_a > 0],
  probs = .75,
  na.rm = TRUE
)

adl_work$du_mar2_life_aBin_q75 <- NA_real_
adl_work$du_mar2_life_aBin_q75[adl_work$du_mar2_life_a == 0] <- 0
adl_work$du_mar2_life_aBin_q75[
  adl_work$du_mar2_life_a > 0 &
    adl_work$du_mar2_life_a < q75_life
] <- 1
adl_work$du_mar2_life_aBin_q75[adl_work$du_mar2_life_a >= q75_life] <- 2

adl_work$du_mar2_life_aBin_q75_ord <- ordered(
  adl_work$du_mar2_life_aBin_q75,
  levels = c(0, 1, 2),
  labels = c("none", "lower_moderate_lifetime", "high_lifetime")
)

contrasts(adl_work$du_mar2_life_aBin_q75_ord) <- contr.poly(3)

table(adl_work$du_mar2_life_aBin_q75_ord, useNA = "ifany")


## Lifetime cannabis continuous sensitivity.

adl_work$du_mar2_life_a_log1p <- log1p(adl_work$du_mar2_life_a)


## Cannabis bin decision table.

cannabis_decision_table_adl <- data.frame(
  exposure = c(
    "Past-year median split",
    "30-day median split",
    "30-day 30g threshold",
    "Lifetime median split",
    "Lifetime q75 split"
  ),
  raw_variable = c(
    "du_mar4_12m_a",
    "du_mar6_30d_a",
    "du_mar6_30d_a",
    "du_mar2_life_a",
    "du_mar2_life_a"
  ),
  bin_variable = c(
    "du_mar4_12m_aBin",
    "du_mar6_30d_aBin",
    "du_mar6_30d_aBin30_lowincl",
    "du_mar2_life_aBin",
    "du_mar2_life_aBin_q75"
  ),
  n_none = c(
    sum(adl_work$du_mar4_12m_aBin == 0, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin == 0, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin30_lowincl == 0, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin == 0, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin_q75 == 0, na.rm = TRUE)
  ),
  n_low_or_lower = c(
    sum(adl_work$du_mar4_12m_aBin == 1, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin == 1, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin30_lowincl == 1, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin == 1, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin_q75 == 1, na.rm = TRUE)
  ),
  n_high_or_heavier = c(
    sum(adl_work$du_mar4_12m_aBin == 2, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin == 2, na.rm = TRUE),
    sum(adl_work$du_mar6_30d_aBin30_lowincl == 2, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin == 2, na.rm = TRUE),
    sum(adl_work$du_mar2_life_aBin_q75 == 2, na.rm = TRUE)
  )
)

cannabis_decision_table_adl


# ------------------------------------------------------------
# 4d. Demographic and clinical covariates
# ------------------------------------------------------------

## First, inspect labels before recoding race/ethnicity and sex.

labelled::val_labels(adl_work$phq_4_latin)
labelled::val_labels(adl_work$phq_5_race)
labelled::val_labels(adl_work$phq_6_gender)
labelled::val_labels(adl_work$phq_7_degree)

table(adl_work$phq_4_latin, useNA = "ifany")
table(adl_work$phq_5_race, useNA = "ifany")
table(adl_work$phq_6_gender, useNA = "ifany")
table(adl_work$phq_7_degree, useNA = "ifany")


## Remove SPSS labels for numeric recoding.

adl_work$phq_4_latin_num <- as.numeric(labelled::remove_labels(adl_work$phq_4_latin))
adl_work$phq_5_race_num <- as.numeric(labelled::remove_labels(adl_work$phq_5_race))


## Race/ethnicity coding.
## IMPORTANT: This assumes:
## - phq_5_race_num == 5 means White
## - phq_4_latin_num == 0 means non-Latinx
##
## Confirm this with the label tables above before interpreting results.

adl_work$race_eth_binary <- factor(
  dplyr::case_when(
    adl_work$phq_5_race_num == 5 & adl_work$phq_4_latin_num == 0 ~ "White non-Hispanic",
    !is.na(adl_work$phq_5_race_num) | !is.na(adl_work$phq_4_latin_num) ~ "Minoritized",
    TRUE ~ NA_character_
  ),
  levels = c("Minoritized", "White non-Hispanic")
)

adl_work$race_eth_binary_covnum <- as.numeric(adl_work$race_eth_binary) - 1
adl_work$race_eth_binary_covfac <- as.factor(adl_work$race_eth_binary)

table(adl_work$race_eth_binary_covfac, useNA = "ifany")


## Safer sex/gender recoding.
## This assumes:
## - phq_6_gender == 1 means Male
## - phq_6_gender == 2 means Female
##
## Confirm this with the label table above before interpreting.

adl_work$sex <- factor(
  dplyr::case_when(
    adl_work$phq_6_gender == 1 ~ "Male",
    adl_work$phq_6_gender == 2 ~ "Female",
    TRUE ~ NA_character_
  ),
  levels = c("Male", "Female")
)

adl_work$sex_covnum <- as.numeric(adl_work$sex) - 1
adl_work$sex_covfac <- as.factor(adl_work$sex)

table(adl_work$sex_covfac, useNA = "ifany")


## Age and education covariates.

adl_work$phq_2_age_c <- center_only(adl_work$phq_2_age)

adl_work$phq_7_degree_covfac <- as.factor(adl_work$phq_7_degree)
adl_work$phq_7_degree_num <- as.numeric(adl_work$phq_7_degree)
adl_work$phq_7_degree_c <- center_only(adl_work$phq_7_degree_num)

summary(adl_work[, c("phq_2_age", "phq_2_age_c", "phq_7_degree", "phq_7_degree_c")])


## Alcohol covariates retained for sensitivity models.

adl_work$du_alc4_12m_a_log1p <- log1p(adl_work$du_alc4_12m_a)
adl_work$du_alc6_30d_a_log1p <- log1p(adl_work$du_alc6_30d_a)

adl_work$du_alc4_12m_a_log1p_c <- center_only(adl_work$du_alc4_12m_a_log1p)
adl_work$du_alc6_30d_a_log1p_c <- center_only(adl_work$du_alc6_30d_a_log1p)


# ============================================================
# 5. ADL AND FUNCTIONAL OUTCOME SETUP
# ============================================================

## Core ADL variables:
## adl_a_no        = number of scorable ADL items
## adl_b_sum_now   = current ADL difficulty
## adl_c_sum_best  = best-functioning ADL difficulty
## adl_total       = standardized ADL decline score, usually (B - C) / A

adl_core_vars <- c(
  "adl_a_no",
  "adl_b_sum_now",
  "adl_c_sum_best",
  "adl_total"
)

## Secondary functional variables.

functional_secondary_vars <- c(
  "ft_total",
  "mac_td",
  "mac_3_forget",
  "mac_6_worse",
  "mac_tmr4",
  "mac_misd",
  "mac_misr",
  "mmt_ts"
)

adl_all_vars <- c(adl_core_vars, functional_secondary_vars)

check_vars(adl_all_vars, adl_work)

summary(adl_work[, adl_core_vars])
summary(adl_work[, intersect(functional_secondary_vars, names(adl_work))])
colSums(is.na(adl_work[, intersect(adl_all_vars, names(adl_work))]))

