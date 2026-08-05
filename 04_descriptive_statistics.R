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



# ============================================================
# ADL PAPER MODEL SEQUENCE
# cannabis-only -> core adjusted -> demographic-expanded -> HIV-adjusted
# ============================================================

library(dplyr)
library(broom)
library(MASS)
library(lmtest)
library(sandwich)

# ============================================================
# 0. Confirm cannabis variable is ordered
# ============================================================

adl_work$du_mar4_12m_aBin_ord <- factor(
  adl_work$du_mar4_12m_aBin_ord,
  levels = c("none", "low", "high"),
  ordered = TRUE
)

contrasts(adl_work$du_mar4_12m_aBin_ord)

# Linear contrast = .L
# Quadratic contrast = .Q

# ============================================================
# 1. PRIMARY: ORDINAL LOGISTIC MODEL
# Outcome: declined ADL domains
#0 declined domains / 1 declined domain / 2+ declined domains
# ============================================================

# Make sure outcome is ordered correctly
adl_work$adl_declined_domain_ord <- factor(
  adl_work$adl_declined_domain_ord,
  levels = c(
    "0_no_decline",
    "1_domain_decline",
    "2plus_domain_decline"
  ),
  ordered = TRUE
)

table(adl_work$adl_declined_domain_ord, useNA = "ifany")

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_ord_0_cannabis <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# cannabis + age + BDI + sex
# ------------------------------------------------------------

model_ord_1_core <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# core + education + race/ethnicity
# ------------------------------------------------------------

model_ord_2_demo <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# demographic-expanded + nadir CD4
# ------------------------------------------------------------

model_ord_3_hiv <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Tidy ordinal model results
# ------------------------------------------------------------
# Note:
# broom::tidy() for MASS::polr() does not always return p.value.
# We calculate an approximate two-tailed p-value from the z statistic.

tidy_polr_results <- function(model, model_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      model = model_name,
      p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    ) %>%
    dplyr::select(
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
}

ordinal_results <- dplyr::bind_rows(
  tidy_polr_results(model_ord_0_cannabis, "0 Cannabis only"),
  tidy_polr_results(model_ord_1_core, "1 Core adjusted"),
  tidy_polr_results(model_ord_2_demo, "2 Demographic-expanded"),
  tidy_polr_results(model_ord_3_hiv, "3 HIV-adjusted")
)

ordinal_results

# Cannabis terms only
ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))
# ------------------------------------------------------------
# Model fit comparison
# ------------------------------------------------------------

AIC(
  model_ord_0_cannabis,
  model_ord_1_core,
  model_ord_2_demo,
  model_ord_3_hiv
)

BIC(
  model_ord_0_cannabis,
  model_ord_1_core,
  model_ord_2_demo,
  model_ord_3_hiv
)

# ============================================================
# Figure: Cannabis odds ratios across ordinal model steps
# ============================================================

library(dplyr)
library(ggplot2)

ordinal_cannabis_plot_df <- ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    )
  )

ggplot(
  ordinal_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    group = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08
  ) +
  labs(
    title = "Cannabis Effects Across Ordinal ADL Model Steps",
    subtitle = "Outcome: 0, 1, or 2+ declined ADL domains",
    x = "Model step",
    y = "Odds ratio",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# ============================================================
# Predicted probabilities from ordinal HIV-adjusted model
# ============================================================

library(emmeans)
library(dplyr)
library(ggplot2)

emm_ord_hiv <- emmeans(
  model_ord_3_hiv,
  ~ du_mar4_12m_aBin_ord | adl_declined_domain_ord,
  mode = "prob"
)

emm_ord_hiv_df <- as.data.frame(emm_ord_hiv)

emm_ord_hiv_plot_df <- emm_ord_hiv_df %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    ),
    adl_decline_category = factor(
      adl_declined_domain_ord,
      levels = c(
        "0_no_decline",
        "1_domain_decline",
        "2plus_domain_decline"
      ),
      labels = c(
        "No declined domains",
        "1 declined domain",
        "2+ declined domains"
      ),
      ordered = TRUE
    )
  )

ggplot(
  emm_ord_hiv_plot_df,
  aes(
    x = cannabis_group,
    y = prob,
    group = adl_decline_category,
    linetype = adl_decline_category
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.08
  ) +
  labs(
    title = "Predicted Probability of ADL Decline by Cannabis Exposure",
    subtitle = "HIV-adjusted ordinal model; outcome = 0, 1, or 2+ declined ADL domains",
    x = "Past-year cannabis exposure",
    y = "Predicted probability",
    linetype = "ADL decline category"
  ) +
  theme_classic(base_size = 13)

# ============================================================
# Common-sample ordinal models for fair AIC/BIC comparison
# ============================================================

adl_ord_common <- adl_work %>%
  dplyr::select(
    adl_declined_domain_ord,
    du_mar4_12m_aBin_ord,
    phq_2_age_c,
    bdi_total,
    sex_covfac,
    phq_7_degree_c,
    race_eth_binary_covfac,
    nadirsqrt_c
  ) %>%
  tidyr::drop_na()

nrow(adl_ord_common)

model_ord_0_cannabis_common <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord,
  data = adl_ord_common,
  Hess = TRUE
)

model_ord_1_core_common <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac,
  data = adl_ord_common,
  Hess = TRUE
)

model_ord_2_demo_common <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac +
    phq_7_degree_c + race_eth_binary_covfac,
  data = adl_ord_common,
  Hess = TRUE
)

model_ord_3_hiv_common <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac +
    phq_7_degree_c + race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_ord_common,
  Hess = TRUE
)

AIC(
  model_ord_0_cannabis_common,
  model_ord_1_core_common,
  model_ord_2_demo_common,
  model_ord_3_hiv_common
)

BIC(
  model_ord_0_cannabis_common,
  model_ord_1_core_common,
  model_ord_2_demo_common,
  model_ord_3_hiv_common
)

sapply(
  list(
    model_ord_0_cannabis_common,
    model_ord_1_core_common,
    model_ord_2_demo_common,
    model_ord_3_hiv_common
  ),
  nobs
)

# ============================================================
# Clean cannabis summary table: ordinal ADL models
# ============================================================

ordinal_cannabis_summary <- ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    OR_CI = paste0(
      round(estimate, 2),
      " [",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      "]"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    OR_CI,
    p_value
  )

ordinal_cannabis_summary

# ============================================================
# Common-sample model fit table
# ============================================================

ordinal_model_fit_common <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_ord_0_cannabis_common),
    nobs(model_ord_1_core_common),
    nobs(model_ord_2_demo_common),
    nobs(model_ord_3_hiv_common)
  ),
  df = c(
    attr(logLik(model_ord_0_cannabis_common), "df"),
    attr(logLik(model_ord_1_core_common), "df"),
    attr(logLik(model_ord_2_demo_common), "df"),
    attr(logLik(model_ord_3_hiv_common), "df")
  ),
  AIC = c(
    AIC(model_ord_0_cannabis_common),
    AIC(model_ord_1_core_common),
    AIC(model_ord_2_demo_common),
    AIC(model_ord_3_hiv_common)
  ),
  BIC = c(
    BIC(model_ord_0_cannabis_common),
    BIC(model_ord_1_core_common),
    BIC(model_ord_2_demo_common),
    BIC(model_ord_3_hiv_common)
  )
) %>%
  dplyr::mutate(
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  )

ordinal_model_fit_common

# ============================================================
# Proportional odds assumption check
# ============================================================

# install.packages("brant") # if needed
library(brant)

brant::brant(model_ord_1_core)
brant::brant(model_ord_2_demo)
brant::brant(model_ord_3_hiv)

# ============================================================
# Multicollinearity check
# ============================================================

performance::check_collinearity(model_ord_3_hiv)

# ============================================================
# Multicollinearity check: core adjusted ordinal model
# ============================================================

performance::check_collinearity(model_ord_1_core)

# ============================================================
# Multicollinearity check: demographic-expanded ordinal model
# ============================================================

performance::check_collinearity(model_ord_2_demo)

# ============================================================
# ANALYSIS 2: BINARY ANY-ADL-DECLINE MODEL
# Outcome: no decline vs any ADL decline
# ============================================================

library(dplyr)
library(broom)
library(lmtest)
library(sandwich)

# ------------------------------------------------------------
# 0. Confirm binary outcome coding
# ------------------------------------------------------------

adl_work$adl_total_any_explore <- factor(
  adl_work$adl_total_any_explore,
  levels = c("no_decline", "any_decline")
)

table(adl_work$adl_total_any_explore, useNA = "ifany")

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_any_0_cannabis <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# cannabis + age + BDI + sex
# ------------------------------------------------------------

model_any_1_core <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# core + education + race/ethnicity
# ------------------------------------------------------------

model_any_2_demo <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# demographic-expanded + nadir CD4
# ------------------------------------------------------------

model_any_3_hiv <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Clean results table
# ------------------------------------------------------------

any_decline_results <- dplyr::bind_rows(
  broom::tidy(model_any_0_cannabis, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "0 Cannabis only"),
  
  broom::tidy(model_any_1_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "1 Core adjusted"),
  
  broom::tidy(model_any_2_demo, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "2 Demographic-expanded"),
  
  broom::tidy(model_any_3_hiv, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "3 HIV-adjusted")
) %>%
  dplyr::select(
    model,
    term,
    estimate,
    conf.low,
    conf.high,
    statistic,
    p.value
  )

any_decline_results

# Cannabis terms only
any_decline_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

# ============================================================
# Cannabis summary table: binary any-ADL-decline model
# ============================================================

any_decline_cannabis_summary <- any_decline_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    OR_CI = paste0(
      round(estimate, 2),
      " [",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      "]"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    OR_CI,
    p_value
  )

any_decline_cannabis_summary

# ============================================================
# Common-sample binary any-decline models
# ============================================================

adl_any_common_full <- adl_work %>%
  dplyr::select(
    adl_total_any_explore,
    du_mar4_12m_aBin_ord,
    phq_2_age_c,
    bdi_total,
    sex_covfac,
    phq_7_degree_c,
    race_eth_binary_covfac,
    nadirsqrt_c
  ) %>%
  tidyr::drop_na()

nrow(adl_any_common_full)

model_any_0_cannabis_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord,
  data = adl_any_common_full,
  family = binomial(link = "logit")
)

model_any_1_core_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_any_common_full,
  family = binomial(link = "logit")
)

model_any_2_demo_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common_full,
  family = binomial(link = "logit")
)

model_any_3_hiv_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_any_common_full,
  family = binomial(link = "logit")
)

any_model_fit_common <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_any_0_cannabis_common),
    nobs(model_any_1_core_common),
    nobs(model_any_2_demo_common),
    nobs(model_any_3_hiv_common)
  ),
  df = c(
    attr(logLik(model_any_0_cannabis_common), "df"),
    attr(logLik(model_any_1_core_common), "df"),
    attr(logLik(model_any_2_demo_common), "df"),
    attr(logLik(model_any_3_hiv_common), "df")
  ),
  AIC = c(
    AIC(model_any_0_cannabis_common),
    AIC(model_any_1_core_common),
    AIC(model_any_2_demo_common),
    AIC(model_any_3_hiv_common)
  ),
  BIC = c(
    BIC(model_any_0_cannabis_common),
    BIC(model_any_1_core_common),
    BIC(model_any_2_demo_common),
    BIC(model_any_3_hiv_common)
  )
) %>%
  dplyr::mutate(
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  )

any_model_fit_common

# ============================================================
# Figure 1: Cannabis ORs across binary any-ADL-decline model steps
# ============================================================

library(dplyr)
library(ggplot2)

any_cannabis_plot_df <- any_decline_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    )
  )

ggplot(
  any_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    group = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08
  ) +
  labs(
    title = "Cannabis Effects Across Binary ADL Decline Model Steps",
    subtitle = "Outcome: any ADL decline vs no decline",
    x = "Model step",
    y = "Odds ratio",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# ============================================================
# Figure 2: Predicted probability of any ADL decline
# Core adjusted binary logistic model
# ============================================================

library(emmeans)
library(dplyr)
library(ggplot2)

emm_any_core <- emmeans(
  model_any_1_core,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_any_core_df <- as.data.frame(emm_any_core) %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

emm_any_core_df

ggplot(
  emm_any_core_df,
  aes(
    x = cannabis_group,
    y = prob,
    group = 1
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.08
  ) +
  labs(
    title = "Predicted Probability of Any ADL Decline by Cannabis Exposure",
    subtitle = "Core adjusted model: age, BDI, and sex",
    x = "Past-year cannabis exposure",
    y = "Predicted probability of any ADL decline"
  ) +
  theme_classic(base_size = 13)
# ============================================================
# Figure 3: Observed proportion with any ADL decline by cannabis group
# ============================================================

observed_any_df <- adl_work %>%
  dplyr::filter(
    !is.na(du_mar4_12m_aBin_ord),
    !is.na(adl_total_any_explore)
  ) %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_any_decline = sum(adl_total_any_explore == "any_decline"),
    prop_any_decline = n_any_decline / n,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

observed_any_df

ggplot(
  observed_any_df,
  aes(
    x = cannabis_group,
    y = prop_any_decline
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        n_any_decline,
        "/",
        n,
        "\n",
        round(prop_any_decline * 100, 1),
        "%"
      )
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  ylim(0, 1) +
  labs(
    title = "Observed Any ADL Decline by Cannabis Exposure",
    subtitle = "Descriptive proportions before covariate adjustment",
    x = "Past-year cannabis exposure",
    y = "Proportion with any ADL decline"
  ) +
  theme_classic(base_size = 13)

# ============================================================
# ANALYSIS 3: CURRENT ADL DIFFICULTY BINARY MODEL
# Outcome: no/minimal current difficulty vs any current difficulty
# ============================================================

summary(adl_work$adl_b_sum_now)
table(adl_work$adl_b_sum_now, useNA = "ifany")

# ------------------------------------------------------------
# Create binary current ADL difficulty outcome
# ------------------------------------------------------------

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_now_any_binary = dplyr::case_when(
      is.na(adl_b_sum_now) ~ NA_character_,
      adl_b_sum_now == 1 ~ "no_current_difficulty",
      adl_b_sum_now > 1 ~ "any_current_difficulty"
    ),
    adl_now_any_binary = factor(
      adl_now_any_binary,
      levels = c("no_current_difficulty", "any_current_difficulty")
    )
  )

table(adl_work$adl_now_any_binary, useNA = "ifany")

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_now_0_cannabis <- glm(
  adl_now_any_binary ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# ------------------------------------------------------------

model_now_1_core <- glm(
  adl_now_any_binary ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# ------------------------------------------------------------

model_now_2_demo <- glm(
  adl_now_any_binary ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# ------------------------------------------------------------

model_now_3_hiv <- glm(
  adl_now_any_binary ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work,
  family = binomial(link = "logit")
)

# ------------------------------------------------------------
# Clean results table
# ------------------------------------------------------------

current_adl_results <- dplyr::bind_rows(
  broom::tidy(model_now_0_cannabis, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "0 Cannabis only"),
  
  broom::tidy(model_now_1_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "1 Core adjusted"),
  
  broom::tidy(model_now_2_demo, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "2 Demographic-expanded"),
  
  broom::tidy(model_now_3_hiv, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "3 HIV-adjusted")
) %>%
  dplyr::select(
    model,
    term,
    estimate,
    conf.low,
    conf.high,
    statistic,
    p.value
  )

current_adl_results

# Cannabis terms only
current_adl_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

# ============================================================
# Cannabis summary table: current ADL difficulty model
# ============================================================

current_adl_cannabis_summary <- current_adl_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    OR_CI = paste0(
      round(estimate, 2),
      " [",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      "]"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    OR_CI,
    p_value
  )

current_adl_cannabis_summary

# ============================================================
# Figure: Observed current ADL difficulty by cannabis group
# ============================================================

observed_now_df <- adl_work %>%
  dplyr::filter(
    !is.na(du_mar4_12m_aBin_ord),
    !is.na(adl_now_any_binary)
  ) %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_any_current = sum(adl_now_any_binary == "any_current_difficulty"),
    prop_any_current = n_any_current / n,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

observed_now_df

ggplot(
  observed_now_df,
  aes(x = cannabis_group, y = prop_any_current)
) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        n_any_current,
        "/",
        n,
        "\n",
        round(prop_any_current * 100, 1),
        "%"
      )
    ),
    vjust = -0.4,
    size = 3.5
  ) +
  ylim(0, 1) +
  labs(
    title = "Observed Current ADL Difficulty by Cannabis Exposure",
    subtitle = "Descriptive proportions before covariate adjustment",
    x = "Past-year cannabis exposure",
    y = "Proportion with any current ADL difficulty"
  ) +
  theme_classic(base_size = 13)

# ============================================================
# Figure: Cannabis effects across current ADL difficulty models
# ============================================================

current_cannabis_plot_df <- current_adl_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    )
  )

ggplot(
  current_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    group = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08
  ) +
  labs(
    title = "Cannabis Effects Across Current ADL Difficulty Model Steps",
    subtitle = "Outcome: any current ADL difficulty vs no current difficulty",
    x = "Model step",
    y = "Odds ratio",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )


# ============================================================
# CREATE ORDINAL CURRENT ADL DIFFICULTY OUTCOME
# Outcome: no / low / higher current ADL difficulty burden
# ============================================================

library(dplyr)

# First inspect original distribution again
summary(adl_work$adl_b_sum_now)
table(adl_work$adl_b_sum_now, useNA = "ifany")

# Create ordinal current ADL burden variable
# Assumption:
#   adl_b_sum_now == 1 reflects no/minimal current difficulty
#   2-3 reflects low current difficulty
#   >=4 reflects higher current difficulty

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_now_ordinal_burden = dplyr::case_when(
      is.na(adl_b_sum_now) ~ NA_character_,
      adl_b_sum_now == 1 ~ "no_current_difficulty",
      adl_b_sum_now %in% c(2, 3) ~ "low_current_difficulty",
      adl_b_sum_now >= 4 ~ "higher_current_difficulty"
    ),
    adl_now_ordinal_burden = factor(
      adl_now_ordinal_burden,
      levels = c(
        "no_current_difficulty",
        "low_current_difficulty",
        "higher_current_difficulty"
      ),
      ordered = TRUE
    )
  )

table(adl_work$adl_now_ordinal_burden, useNA = "ifany")
prop.table(table(adl_work$adl_now_ordinal_burden))

# ============================================================
# Check score ranges and values by ordinal category
# ============================================================

adl_now_burden_ranges <- adl_work %>%
  dplyr::filter(
    !is.na(adl_b_sum_now),
    !is.na(adl_now_ordinal_burden)
  ) %>%
  dplyr::group_by(adl_now_ordinal_burden) %>%
  dplyr::summarise(
    n = dplyr::n(),
    min_score = min(adl_b_sum_now, na.rm = TRUE),
    max_score = max(adl_b_sum_now, na.rm = TRUE),
    mean_score = mean(adl_b_sum_now, na.rm = TRUE),
    median_score = median(adl_b_sum_now, na.rm = TRUE),
    .groups = "drop"
  )

adl_now_burden_ranges

# Exact observed values by category
adl_now_values_by_category <- adl_work %>%
  dplyr::filter(
    !is.na(adl_b_sum_now),
    !is.na(adl_now_ordinal_burden)
  ) %>%
  dplyr::count(adl_now_ordinal_burden, adl_b_sum_now)

adl_now_values_by_category

# ============================================================
# ANALYSIS 3B: ORDINAL CURRENT ADL DIFFICULTY BURDEN MODEL
# Outcome: no / low / higher current difficulty
# ============================================================

library(MASS)
library(broom)
library(dplyr)

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_now_ord_0_cannabis <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# cannabis + age + BDI + sex
# ------------------------------------------------------------

model_now_ord_1_core <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# core + education + race/ethnicity
# ------------------------------------------------------------

model_now_ord_2_demo <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# demographic-expanded + nadir CD4
# ------------------------------------------------------------

model_now_ord_3_hiv <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work,
  Hess = TRUE
)

# ============================================================
# Tidy ordinal current ADL burden model results
# ============================================================

tidy_polr_results <- function(model, model_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      model = model_name,
      p.value = 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    ) %>%
    dplyr::select(
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
}

current_ordinal_results <- dplyr::bind_rows(
  tidy_polr_results(model_now_ord_0_cannabis, "0 Cannabis only"),
  tidy_polr_results(model_now_ord_1_core, "1 Core adjusted"),
  tidy_polr_results(model_now_ord_2_demo, "2 Demographic-expanded"),
  tidy_polr_results(model_now_ord_3_hiv, "3 HIV-adjusted")
)

current_ordinal_results

# Cannabis terms only
current_ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

# ============================================================
# Cannabis summary table: ordinal current ADL burden model
# ============================================================

current_ordinal_cannabis_summary <- current_ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    OR_CI = paste0(
      round(estimate, 2),
      " [",
      round(conf.low, 2),
      ", ",
      round(conf.high, 2),
      "]"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    OR_CI,
    p_value
  )

current_ordinal_cannabis_summary

# ============================================================
# Common-sample model fit: ordinal current ADL burden
# ============================================================

adl_now_ord_common <- adl_work %>%
  dplyr::select(
    adl_now_ordinal_burden,
    du_mar4_12m_aBin_ord,
    phq_2_age_c,
    bdi_total,
    sex_covfac,
    phq_7_degree_c,
    race_eth_binary_covfac,
    nadirsqrt_c
  ) %>%
  tidyr::drop_na()

nrow(adl_now_ord_common)

model_now_ord_0_common <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord,
  data = adl_now_ord_common,
  Hess = TRUE
)

model_now_ord_1_common <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_now_ord_common,
  Hess = TRUE
)

model_now_ord_2_common <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_now_ord_common,
  Hess = TRUE
)

model_now_ord_3_common <- MASS::polr(
  adl_now_ordinal_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_now_ord_common,
  Hess = TRUE
)

current_ordinal_model_fit_common <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_now_ord_0_common),
    nobs(model_now_ord_1_common),
    nobs(model_now_ord_2_common),
    nobs(model_now_ord_3_common)
  ),
  df = c(
    attr(logLik(model_now_ord_0_common), "df"),
    attr(logLik(model_now_ord_1_common), "df"),
    attr(logLik(model_now_ord_2_common), "df"),
    attr(logLik(model_now_ord_3_common), "df")
  ),
  AIC = c(
    AIC(model_now_ord_0_common),
    AIC(model_now_ord_1_common),
    AIC(model_now_ord_2_common),
    AIC(model_now_ord_3_common)
  ),
  BIC = c(
    BIC(model_now_ord_0_common),
    BIC(model_now_ord_1_common),
    BIC(model_now_ord_2_common),
    BIC(model_now_ord_3_common)
  )
) %>%
  dplyr::mutate(
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  )

current_ordinal_model_fit_common

# ============================================================
# FIGURE 2: Cannabis ORs across ordinal current ADL burden steps
# ============================================================

library(dplyr)
library(ggplot2)

current_ordinal_cannabis_plot_df <- current_ordinal_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    )
  )

ggplot(
  current_ordinal_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    group = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.08
  ) +
  labs(
    title = "Cannabis Effects Across Current ADL Burden Model Steps",
    subtitle = "Outcome: no, low, or higher current ADL difficulty burden",
    x = "Model step",
    y = "Odds ratio",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# ============================================================
# OPTIONAL FIGURE 3: Observed current ADL burden proportions
# by cannabis group
# ============================================================

observed_now_ord_df <- adl_work %>%
  dplyr::filter(
    !is.na(du_mar4_12m_aBin_ord),
    !is.na(adl_now_ordinal_burden)
  ) %>%
  dplyr::count(
    du_mar4_12m_aBin_ord,
    adl_now_ordinal_burden
  ) %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::mutate(
    group_n = sum(n),
    prop = n / group_n
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    ),
    current_adl_burden_category = factor(
      adl_now_ordinal_burden,
      levels = c(
        "no_current_difficulty",
        "low_current_difficulty",
        "higher_current_difficulty"
      ),
      labels = c(
        "No current difficulty",
        "Low current difficulty",
        "Higher current difficulty"
      ),
      ordered = TRUE
    )
  )

observed_now_ord_df

ggplot(
  observed_now_ord_df,
  aes(
    x = cannabis_group,
    y = prop,
    fill = current_adl_burden_category
  )
) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = paste0(n, "/", group_n)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3
  ) +
  ylim(0, 1) +
  labs(
    title = "Observed Current ADL Burden by Cannabis Exposure",
    subtitle = "Descriptive proportions before covariate adjustment",
    x = "Past-year cannabis exposure",
    y = "Proportion",
    fill = "Current ADL burden"
  ) +
  theme_classic(base_size = 13)

# ============================================================
# ANALYSIS 4: RAW ADL_TOTAL LINEAR MODEL WITH HC3 ROBUST SEs
# Outcome: continuous adl_total
# ============================================================

summary(adl_work$adl_total)

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_raw_0_cannabis <- lm(
  adl_total ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# ------------------------------------------------------------

model_raw_1_core <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# ------------------------------------------------------------

model_raw_2_demo <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# ------------------------------------------------------------

model_raw_3_hiv <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work
)

# ============================================================
# Function for HC3 robust linear model results
# More stable version using broom::tidy()
# ============================================================

tidy_hc3 <- function(model, model_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  broom::tidy(hc3) %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(
      model,
      term,
      estimate,
      std.error,
      statistic,
      p.value
    ) %>%
    dplyr::rename(
      std_error_hc3 = std.error
    )
}

# ============================================================
# HC3 robust results for raw adl_total models
# ============================================================

raw_adl_hc3_results <- dplyr::bind_rows(
  tidy_hc3(model_raw_0_cannabis, "0 Cannabis only"),
  tidy_hc3(model_raw_1_core, "1 Core adjusted"),
  tidy_hc3(model_raw_2_demo, "2 Demographic-expanded"),
  tidy_hc3(model_raw_3_hiv, "3 HIV-adjusted")
)

raw_adl_hc3_results

# Cannabis terms only
raw_adl_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

# ============================================================
# Cannabis summary table: raw adl_total HC3 model
# ============================================================

raw_adl_cannabis_summary <- raw_adl_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    B_SE = paste0(
      round(estimate, 3),
      " (",
      round(std_error_hc3, 3),
      ")"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    B_SE,
    p_value
  )

raw_adl_cannabis_summary

# ============================================================
# FIGURE: Raw ADL total sensitivity model
# Cannabis linear/quadratic effects across model steps
# HC3 robust standard errors
# ============================================================

library(dplyr)
library(ggplot2)
library(stringr)

# ------------------------------------------------------------
# 1. Create clean plotting dataset from HC3 results
# ------------------------------------------------------------

raw_adl_cannabis_plot_df <- raw_adl_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast",
      TRUE ~ term
    ),
    
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    ),
    
    cannabis_term = factor(
      cannabis_term,
      levels = c(
        "Linear cannabis contrast",
        "Quadratic cannabis contrast"
      )
    ),
    
    ci_low = estimate - 1.96 * std_error_hc3,
    ci_high = estimate + 1.96 * std_error_hc3
  )

raw_adl_cannabis_plot_df

# ------------------------------------------------------------
# 2. Coefficient plot
# ------------------------------------------------------------

figure_raw_adl_hc3_coefficients <- ggplot(
  raw_adl_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    ymin = ci_low,
    ymax = ci_high,
    group = cannabis_term,
    shape = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.45),
    linewidth = 0.6
  ) +
  geom_line(
    position = position_dodge(width = 0.45),
    linewidth = 0.6
  ) +
  labs(
    title = "Raw ADL Decline Sensitivity Model",
    subtitle = "Cannabis effects across sequential adjustment models; HC3 robust 95% CIs",
    x = "Model step",
    y = "Regression coefficient for raw ADL decline score",
    shape = "Cannabis term",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

figure_raw_adl_hc3_coefficients

# ============================================================
# OPTIONAL FIGURE: Predicted raw ADL total by cannabis group
# Core adjusted model
# ============================================================

library(emmeans)

emm_raw_core <- emmeans(
  model_raw_1_core,
  ~ du_mar4_12m_aBin_ord
)

emm_raw_core_df <- as.data.frame(emm_raw_core) %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

emm_raw_core_df


figure_raw_adl_predicted_core <- ggplot(
  emm_raw_core_df,
  aes(
    x = cannabis_group,
    y = emmean,
    ymin = lower.CL,
    ymax = upper.CL,
    group = 1
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_errorbar(width = 0.10, linewidth = 0.6) +
  labs(
    title = "Predicted Raw ADL Decline Score by Cannabis Group",
    subtitle = "Core adjusted model: age, depressive symptoms, and sex",
    x = "Past-year cannabis exposure",
    y = "Predicted raw ADL decline score"
  ) +
  theme_classic(base_size = 12)

figure_raw_adl_predicted_core

# ============================================================
# ANALYSIS 5: CURRENT ADL BURDEN AS CONTINUOUS OUTCOME
# Outcome: adl_b_sum_now
# Model: linear regression with HC3 robust SEs
# Purpose: sensitivity/comparability model
# ============================================================

library(dplyr)
library(broom)
library(lmtest)
library(sandwich)
library(ggplot2)

# ------------------------------------------------------------
# 0. Inspect outcome distribution
# ------------------------------------------------------------

summary(adl_work$adl_b_sum_now)
table(adl_work$adl_b_sum_now, useNA = "ifany")

mean(adl_work$adl_b_sum_now, na.rm = TRUE)
var(adl_work$adl_b_sum_now, na.rm = TRUE)
var(adl_work$adl_b_sum_now, na.rm = TRUE) / mean(adl_work$adl_b_sum_now, na.rm = TRUE)

# ============================================================
# 1. Fit continuous current ADL burden models
# ============================================================

# ------------------------------------------------------------
# Model 0: Cannabis only
# ------------------------------------------------------------

model_now_cont_0_cannabis <- lm(
  adl_b_sum_now ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

# ------------------------------------------------------------
# Model 1: Core adjusted
# cannabis + age + BDI + sex
# ------------------------------------------------------------

model_now_cont_1_core <- lm(
  adl_b_sum_now ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# ------------------------------------------------------------
# Model 2: Demographic-expanded
# core + education + race/ethnicity
# ------------------------------------------------------------

model_now_cont_2_demo <- lm(
  adl_b_sum_now ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

# ------------------------------------------------------------
# Model 3: HIV-adjusted
# demographic-expanded + nadir CD4
# ------------------------------------------------------------

model_now_cont_3_hiv <- lm(
  adl_b_sum_now ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work
)

# ============================================================
# 2. HC3 robust results
# ============================================================

tidy_hc3 <- function(model, model_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  broom::tidy(hc3) %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(
      model,
      term,
      estimate,
      std.error,
      statistic,
      p.value
    ) %>%
    dplyr::rename(
      std_error_hc3 = std.error
    )
}

current_adl_cont_hc3_results <- dplyr::bind_rows(
  tidy_hc3(model_now_cont_0_cannabis, "0 Cannabis only"),
  tidy_hc3(model_now_cont_1_core, "1 Core adjusted"),
  tidy_hc3(model_now_cont_2_demo, "2 Demographic-expanded"),
  tidy_hc3(model_now_cont_3_hiv, "3 HIV-adjusted")
)

current_adl_cont_hc3_results

# Cannabis terms only
current_adl_cont_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

# ============================================================
# 3. Clean cannabis summary table
# ============================================================

current_adl_cont_cannabis_summary <- current_adl_cont_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast"
    ),
    B_SE = paste0(
      round(estimate, 3),
      " (",
      round(std_error_hc3, 3),
      ")"
    ),
    p_value = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model,
    cannabis_term,
    B_SE,
    p_value
  )

current_adl_cont_cannabis_summary

# ============================================================
# 4. Model fit table
# ============================================================

current_adl_cont_model_fit <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_now_cont_0_cannabis),
    nobs(model_now_cont_1_core),
    nobs(model_now_cont_2_demo),
    nobs(model_now_cont_3_hiv)
  ),
  r_squared = c(
    summary(model_now_cont_0_cannabis)$r.squared,
    summary(model_now_cont_1_core)$r.squared,
    summary(model_now_cont_2_demo)$r.squared,
    summary(model_now_cont_3_hiv)$r.squared
  ),
  adj_r_squared = c(
    summary(model_now_cont_0_cannabis)$adj.r.squared,
    summary(model_now_cont_1_core)$adj.r.squared,
    summary(model_now_cont_2_demo)$adj.r.squared,
    summary(model_now_cont_3_hiv)$adj.r.squared
  ),
  AIC = c(
    AIC(model_now_cont_0_cannabis),
    AIC(model_now_cont_1_core),
    AIC(model_now_cont_2_demo),
    AIC(model_now_cont_3_hiv)
  ),
  BIC = c(
    BIC(model_now_cont_0_cannabis),
    BIC(model_now_cont_1_core),
    BIC(model_now_cont_2_demo),
    BIC(model_now_cont_3_hiv)
  )
) %>%
  dplyr::mutate(
    r_squared = round(r_squared, 3),
    adj_r_squared = round(adj_r_squared, 3),
    AIC = round(AIC, 2),
    BIC = round(BIC, 2)
  )

current_adl_cont_model_fit

# ============================================================
# FIGURE: Continuous current ADL burden sensitivity model
# Cannabis linear/quadratic effects across model steps
# HC3 robust 95% CIs
# ============================================================

current_adl_cont_cannabis_plot_df <- current_adl_cont_hc3_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::mutate(
    cannabis_term = dplyr::case_when(
      grepl("\\.L", term) ~ "Linear cannabis contrast",
      grepl("\\.Q", term) ~ "Quadratic cannabis contrast",
      TRUE ~ term
    ),
    model = factor(
      model,
      levels = c(
        "0 Cannabis only",
        "1 Core adjusted",
        "2 Demographic-expanded",
        "3 HIV-adjusted"
      )
    ),
    cannabis_term = factor(
      cannabis_term,
      levels = c(
        "Linear cannabis contrast",
        "Quadratic cannabis contrast"
      )
    ),
    ci_low = estimate - 1.96 * std_error_hc3,
    ci_high = estimate + 1.96 * std_error_hc3
  )

figure_current_adl_cont_hc3 <- ggplot(
  current_adl_cont_cannabis_plot_df,
  aes(
    x = model,
    y = estimate,
    ymin = ci_low,
    ymax = ci_high,
    group = cannabis_term,
    shape = cannabis_term,
    linetype = cannabis_term
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.45),
    linewidth = 0.6
  ) +
  geom_line(
    position = position_dodge(width = 0.45),
    linewidth = 0.6
  ) +
  labs(
    title = "Current ADL Burden Sensitivity Model",
    subtitle = "Cannabis effects across sequential adjustment models; HC3 robust 95% CIs",
    x = "Model step",
    y = "Regression coefficient for current ADL burden score",
    shape = "Cannabis term",
    linetype = "Cannabis term"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )

figure_current_adl_cont_hc3

# ============================================================
# Diagnostics for HIV-adjusted continuous current ADL model
# ============================================================

performance::check_heteroscedasticity(model_now_cont_3_hiv)
performance::check_normality(model_now_cont_3_hiv)
performance::check_collinearity(model_now_cont_3_hiv)

par(mfrow = c(2, 2))
plot(model_now_cont_3_hiv)
par(mfrow = c(1, 1))


# ============================================================
# POWER / MODEL-CAPACITY CHECKS FOR ADL PAPER
# ============================================================
# Purpose:
# These checks help determine whether each model had enough
# sample size/events/parameters to support interpretation.
#
# Important:
# These do NOT prove null effects.
# They summarize model capacity, precision, and sensitivity.
# ============================================================

library(dplyr)
library(broom)
library(pwr)
library(MASS)


# ============================================================
# Helper functions
# ============================================================

# Count predictors excluding intercept / thresholds
count_terms_lm_glm <- function(model) {
  length(coef(model)) - 1
}

# Count predictors for polr models
# Removes intercept-like threshold terms
count_terms_polr <- function(model) {
  length(coef(model))
}

# Convert R2 to Cohen's f2
r2_to_f2 <- function(r2) {
  r2 / (1 - r2)
}

# Clean pwr output
clean_pwr_f2 <- function(pwr_object, model_name, effect_label) {
  tibble::tibble(
    model = model_name,
    effect_tested = effect_label,
    numerator_df = pwr_object$u,
    denominator_df = pwr_object$v,
    f2 = pwr_object$f2,
    alpha = pwr_object$sig.level,
    estimated_power = pwr_object$power
  )
}

# ============================================================
# A. MODEL CAPACITY TABLE ACROSS ANALYSES
# ============================================================

model_capacity_overview <- tibble::tibble(
  analysis = c(
    "Primary ordinal ADL decline",
    "Binary any ADL decline",
    "Ordinal current ADL burden",
    "Raw adl_total linear HC3",
    "Current ADL burden continuous HC3"
  ),
  
  model_object = list(
    model_ord_3_hiv,
    model_any_3_hiv,
    model_now_ord_3_hiv,
    model_raw_3_hiv,
    model_now_cont_3_hiv
  ),
  
  model_type = c(
    "Ordinal logistic",
    "Binary logistic",
    "Ordinal logistic",
    "Linear",
    "Linear"
  )
) %>%
  dplyr::mutate(
    n = purrr::map_int(model_object, nobs),
    
    parameters_excluding_intercept = dplyr::case_when(
      model_type == "Ordinal logistic" ~ purrr::map_int(model_object, count_terms_polr),
      TRUE ~ purrr::map_int(model_object, count_terms_lm_glm)
    ),
    
    n_per_parameter = n / parameters_excluding_intercept
  ) %>%
  dplyr::select(
    analysis,
    model_type,
    n,
    parameters_excluding_intercept,
    n_per_parameter
  )

model_capacity_overview

# ============================================================
# B. BINARY LOGISTIC MODEL CAPACITY
# Outcome: any ADL decline vs no decline
# ============================================================

event_table_any <- table(adl_work$adl_total_any_explore)
event_table_any

n_events_any <- sum(adl_work$adl_total_any_explore == "any_decline", na.rm = TRUE)
n_none_any   <- sum(adl_work$adl_total_any_explore == "no_decline", na.rm = TRUE)

binary_any_capacity <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_any_0_cannabis),
    nobs(model_any_1_core),
    nobs(model_any_2_demo),
    nobs(model_any_3_hiv)
  ),
  n_events = n_events_any,
  n_nonevents = n_none_any,
  predictors_excluding_intercept = c(
    count_terms_lm_glm(model_any_0_cannabis),
    count_terms_lm_glm(model_any_1_core),
    count_terms_lm_glm(model_any_2_demo),
    count_terms_lm_glm(model_any_3_hiv)
  )
) %>%
  dplyr::mutate(
    events_per_parameter = n_events / predictors_excluding_intercept,
    nonevents_per_parameter = n_nonevents / predictors_excluding_intercept
  )

binary_any_capacity

# ============================================================
# C. PRIMARY ORDINAL ADL DECLINE MODEL CAPACITY
# Outcome: 0, 1, or 2+ declined ADL domains
# ============================================================

table(adl_work$adl_declined_domain_ord, useNA = "ifany")

ordinal_primary_capacity <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_ord_0_cannabis),
    nobs(model_ord_1_core),
    nobs(model_ord_2_demo),
    nobs(model_ord_3_hiv)
  ),
  predictors_excluding_thresholds = c(
    count_terms_polr(model_ord_0_cannabis),
    count_terms_polr(model_ord_1_core),
    count_terms_polr(model_ord_2_demo),
    count_terms_polr(model_ord_3_hiv)
  )
) %>%
  dplyr::mutate(
    n_per_predictor = n / predictors_excluding_thresholds
  )

ordinal_primary_capacity

# ============================================================
# D. ORDINAL CURRENT ADL BURDEN MODEL CAPACITY
# Outcome: no / low / higher current difficulty
# ============================================================

table(adl_work$adl_now_ordinal_burden, useNA = "ifany")

ordinal_current_capacity <- tibble::tibble(
  model = c(
    "0 Cannabis only",
    "1 Core adjusted",
    "2 Demographic-expanded",
    "3 HIV-adjusted"
  ),
  n = c(
    nobs(model_now_ord_0_cannabis),
    nobs(model_now_ord_1_core),
    nobs(model_now_ord_2_demo),
    nobs(model_now_ord_3_hiv)
  ),
  predictors_excluding_thresholds = c(
    count_terms_polr(model_now_ord_0_cannabis),
    count_terms_polr(model_now_ord_1_core),
    count_terms_polr(model_now_ord_2_demo),
    count_terms_polr(model_now_ord_3_hiv)
  )
) %>%
  dplyr::mutate(
    n_per_predictor = n / predictors_excluding_thresholds
  )

ordinal_current_capacity

# ============================================================
# E. LINEAR MODEL POWER: RAW ADL_TOTAL
# Outcome: adl_total
# ============================================================

# Full HIV-adjusted model
summary(model_raw_3_hiv)

r2_raw_full <- summary(model_raw_3_hiv)$r.squared
f2_raw_full <- r2_to_f2(r2_raw_full)

u_raw_full <- length(coef(model_raw_3_hiv)) - 1
v_raw_full <- df.residual(model_raw_3_hiv)

power_raw_full <- pwr::pwr.f2.test(
  u = u_raw_full,
  v = v_raw_full,
  f2 = f2_raw_full,
  sig.level = 0.05
)

power_raw_full

# ------------------------------------------------------------
# Incremental cannabis block: compare no-cannabis vs with-cannabis
# ------------------------------------------------------------

model_raw_3_hiv_no_cannabis <- lm(
  adl_total ~ phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work
)

r2_raw_no_cu <- summary(model_raw_3_hiv_no_cannabis)$r.squared
r2_raw_with_cu <- summary(model_raw_3_hiv)$r.squared

delta_r2_raw_cu <- r2_raw_with_cu - r2_raw_no_cu
f2_raw_cu_block <- delta_r2_raw_cu / (1 - r2_raw_with_cu)

power_raw_cannabis_block <- pwr::pwr.f2.test(
  u = 2, # cannabis L + Q
  v = df.residual(model_raw_3_hiv),
  f2 = f2_raw_cu_block,
  sig.level = 0.05
)

power_raw_cannabis_block

raw_linear_power_summary <- dplyr::bind_rows(
  clean_pwr_f2(power_raw_full, "Raw adl_total HIV-adjusted", "Overall full model"),
  clean_pwr_f2(power_raw_cannabis_block, "Raw adl_total HIV-adjusted", "Added cannabis L + Q block")
) %>%
  dplyr::mutate(
    r2_full = r2_raw_full,
    delta_r2_cannabis = delta_r2_raw_cu
  )

raw_linear_power_summary

# ============================================================
# F. LINEAR MODEL POWER: CURRENT ADL BURDEN CONTINUOUS
# Outcome: adl_b_sum_now
# ============================================================

summary(model_now_cont_3_hiv)

r2_now_full <- summary(model_now_cont_3_hiv)$r.squared
f2_now_full <- r2_to_f2(r2_now_full)

u_now_full <- length(coef(model_now_cont_3_hiv)) - 1
v_now_full <- df.residual(model_now_cont_3_hiv)

power_now_full <- pwr::pwr.f2.test(
  u = u_now_full,
  v = v_now_full,
  f2 = f2_now_full,
  sig.level = 0.05
)

power_now_full

# ------------------------------------------------------------
# Incremental cannabis block
# ------------------------------------------------------------

model_now_cont_3_hiv_no_cannabis <- lm(
  adl_b_sum_now ~ phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_work
)

r2_now_no_cu <- summary(model_now_cont_3_hiv_no_cannabis)$r.squared
r2_now_with_cu <- summary(model_now_cont_3_hiv)$r.squared

delta_r2_now_cu <- r2_now_with_cu - r2_now_no_cu
f2_now_cu_block <- delta_r2_now_cu / (1 - r2_now_with_cu)

power_now_cannabis_block <- pwr::pwr.f2.test(
  u = 2, # cannabis L + Q
  v = df.residual(model_now_cont_3_hiv),
  f2 = f2_now_cu_block,
  sig.level = 0.05
)

power_now_cannabis_block

current_linear_power_summary <- dplyr::bind_rows(
  clean_pwr_f2(power_now_full, "Current ADL burden continuous HIV-adjusted", "Overall full model"),
  clean_pwr_f2(power_now_cannabis_block, "Current ADL burden continuous HIV-adjusted", "Added cannabis L + Q block")
) %>%
  dplyr::mutate(
    r2_full = r2_now_full,
    delta_r2_cannabis = delta_r2_now_cu
  )

current_linear_power_summary

# ============================================================
# G. SIMULATION POWER: BINARY ANY ADL DECLINE
# Testing linear cannabis effect
# ============================================================

# Make numeric binary outcome
adl_work$adl_any_num_power <- ifelse(
  adl_work$adl_total_any_explore == "any_decline", 1,
  ifelse(adl_work$adl_total_any_explore == "no_decline", 0, NA)
)

# Create analysis sample
adl_any_power_data <- adl_work %>%
  dplyr::select(
    adl_any_num_power,
    du_mar4_12m_aBin_ord,
    phq_2_age_c,
    bdi_total,
    sex_covfac,
    phq_7_degree_c,
    race_eth_binary_covfac,
    nadirsqrt_c
  ) %>%
  tidyr::drop_na()

# Extract linear cannabis contrast manually
contrasts(adl_any_power_data$du_mar4_12m_aBin_ord)

adl_any_power_data$cu12m_L <- as.numeric(adl_any_power_data$du_mar4_12m_aBin_ord)
adl_any_power_data$cu12m_L <- scale(adl_any_power_data$cu12m_L, center = TRUE, scale = TRUE)[, 1]

# Baseline model without cannabis
model_any_base_no_cu_power <- glm(
  adl_any_num_power ~ phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac +
    nadirsqrt_c,
  data = adl_any_power_data,
  family = binomial(link = "logit")
)

simulate_logistic_power <- function(data, cannabis_beta, nsim = 1000, alpha = .05) {
  
  p_values <- rep(NA_real_, nsim)
  
  for (i in seq_len(nsim)) {
    
    sim_data <- data
    
    lp_base <- predict(
      model_any_base_no_cu_power,
      newdata = sim_data,
      type = "link"
    )
    
    lp_sim <- lp_base + cannabis_beta * sim_data$cu12m_L
    prob_sim <- plogis(lp_sim)
    
    sim_data$y_sim <- rbinom(
      n = nrow(sim_data),
      size = 1,
      prob = prob_sim
    )
    
    fit <- try(
      glm(
        y_sim ~ cu12m_L +
          phq_2_age_c +
          bdi_total +
          sex_covfac +
          phq_7_degree_c +
          race_eth_binary_covfac +
          nadirsqrt_c,
        data = sim_data,
        family = binomial(link = "logit")
      ),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      p_values[i] <- summary(fit)$coefficients["cu12m_L", "Pr(>|z|)"]
    }
  }
  
  tibble::tibble(
    cannabis_OR = exp(cannabis_beta),
    cannabis_beta = cannabis_beta,
    nsim = nsim,
    alpha = alpha,
    estimated_power = mean(p_values < alpha, na.rm = TRUE),
    failed_models = sum(is.na(p_values))
  )
}

binary_any_sim_power <- purrr::map_dfr(
  log(c(1.25, 1.50, 1.75, 2.00, 2.50, 3.00)),
  ~ simulate_logistic_power(
    data = adl_any_power_data,
    cannabis_beta = .x,
    nsim = 1000,
    alpha = .05
  )
)

binary_any_sim_power

# ============================================================
# H. MODEL CAPACITY: EXPLORATORY MODERATION
# Linear cannabis x nadir CD4
# ============================================================

# Standardize nadir CD4
adl_work$nadirsqrt_z <- as.numeric(scale(adl_work$nadirsqrt_c))

# Make sure linear cannabis contrast exists
adl_work$cu12m_L_simple <- as.numeric(adl_work$du_mar4_12m_aBin_ord)
adl_work$cu12m_L_simple <- as.numeric(scale(adl_work$cu12m_L_simple))

model_any_linear_interaction <- glm(
  adl_total_any_explore ~ cu12m_L_simple * nadirsqrt_z +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

interaction_capacity <- tibble::tibble(
  model = "Binary any ADL decline: linear cannabis x nadir CD4",
  n = nobs(model_any_linear_interaction),
  n_events = sum(model.frame(model_any_linear_interaction)$adl_total_any_explore == "any_decline"),
  predictors_excluding_intercept = count_terms_lm_glm(model_any_linear_interaction),
  events_per_parameter =
    n_events / predictors_excluding_intercept
)

interaction_capacity

# ============================================================
# I. FINAL POWER / CAPACITY SUMMARY TABLES
# ============================================================

model_capacity_overview
binary_any_capacity
ordinal_primary_capacity
ordinal_current_capacity
raw_linear_power_summary
current_linear_power_summary
binary_any_sim_power
interaction_capacity






# ============================================================
# ADL OUTCOME MODELS
# Hierarchical / sequential covariate adjustment
# ============================================================

library(dplyr)
library(broom)
library(MASS)
library(emmeans)

# ------------------------------------------------------------
# 0. Make sure core ADL outcome variables exist
# ------------------------------------------------------------

adl_work <- adl_work %>%
  dplyr::mutate(
    
    # Recalculated proportional ADL decline score
    adl_total_calc = (adl_b_sum_now - adl_c_sum_best) / adl_a_no,
    
    # Binary ADL decline: no decline vs any decline
    adl_any_decline_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 ~ "any_decline"
    ),
    adl_any_decline_calc = factor(
      adl_any_decline_calc,
      levels = c("no_decline", "any_decline")
    ),
    
    # Current ADL burden count, where 0 = no current burden
    adl_b_now_excess = adl_b_sum_now - 1,
    
    # Binary current ADL burden: no current burden vs any current burden
    adl_any_current_burden = dplyr::case_when(
      is.na(adl_b_now_excess) ~ NA_character_,
      adl_b_now_excess == 0 ~ "no_current_burden",
      adl_b_now_excess > 0 ~ "any_current_burden"
    ),
    adl_any_current_burden = factor(
      adl_any_current_burden,
      levels = c("no_current_burden", "any_current_burden")
    )
  )

# ------------------------------------------------------------
# 0b. Create 3-level ADL decline severity outcome
#     0 = no decline
#     1 = lower decline among decliners
#     2 = greater decline among decliners
# ------------------------------------------------------------

positive_median_adl_calc <- median(
  adl_work$adl_total_calc[adl_work$adl_total_calc > 0],
  na.rm = TRUE
)

positive_median_adl_calc

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_decline_severity_3cat_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 & adl_total_calc <= positive_median_adl_calc ~ "lower_decline",
      adl_total_calc > positive_median_adl_calc ~ "greater_decline"
    ),
    adl_decline_severity_3cat_calc = factor(
      adl_decline_severity_3cat_calc,
      levels = c("no_decline", "lower_decline", "greater_decline"),
      ordered = TRUE
    )
  )

# Check outcome distributions
table(adl_work$adl_any_decline_calc, useNA = "ifany")
table(adl_work$adl_decline_severity_3cat_calc, useNA = "ifany")
table(adl_work$adl_any_current_burden, useNA = "ifany")
table(adl_work$adl_b_now_excess, useNA = "ifany")

# ============================================================
# MODEL 1 OUTCOME: Binary ADL decline yes/no
# ============================================================

# Model 0: Cannabis only
m1_adl_decline_binary_0 <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

# Model 1: Core-adjusted
m1_adl_decline_binary_core <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# Model 2: Disease-severity-adjusted
m1_adl_decline_binary_disease <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

# Extract odds ratios
m1_results <- dplyr::bind_rows(
  broom::tidy(m1_adl_decline_binary_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  
  broom::tidy(m1_adl_decline_binary_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  
  broom::tidy(m1_adl_decline_binary_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m1_results

emm_m1 <- emmeans::emmeans(
  m1_adl_decline_binary_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_m1
pairs(emm_m1)

# ============================================================
# MODEL 2 OUTCOME: 3-level ADL decline severity
# no decline / lower decline / greater decline
# ============================================================

# Model 0: Cannabis only
m2_adl_decline_severity_0 <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

# Model 1: Core-adjusted
m2_adl_decline_severity_core <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

# Model 2: Disease-severity-adjusted
m2_adl_decline_severity_disease <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

# Extract proportional odds ratios
m2_results <- dplyr::bind_rows(
  broom::tidy(m2_adl_decline_severity_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  
  broom::tidy(m2_adl_decline_severity_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  
  broom::tidy(m2_adl_decline_severity_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m2_results

emm_m2 <- emmeans::emmeans(
  m2_adl_decline_severity_disease,
  ~ du_mar4_12m_aBin_ord,
  mode = "prob"
)

emm_m2

# Optional: Brant test if package is installed
# install.packages("brant")
library(brant)

brant::brant(m2_adl_decline_severity_disease)

# ============================================================
# MODEL 3 OUTCOME: Binary current ADL burden yes/no
# ============================================================

# Model 0: Cannabis only
m3_current_burden_binary_0 <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

# Model 1: Core-adjusted
m3_current_burden_binary_core <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# Model 2: Disease-severity-adjusted
m3_current_burden_binary_disease <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

# Extract odds ratios
m3_results <- dplyr::bind_rows(
  broom::tidy(m3_current_burden_binary_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  
  broom::tidy(m3_current_burden_binary_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  
  broom::tidy(m3_current_burden_binary_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m3_results

emm_m3 <- emmeans::emmeans(
  m3_current_burden_binary_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_m3
pairs(emm_m3)

# ============================================================
# MODEL 4 OUTCOME: Current ADL burden count
# Negative binomial model because outcome is overdispersed
# ============================================================

# Model 0: Cannabis only
m4_current_burden_count_0 <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

# Model 1: Core-adjusted
m4_current_burden_count_core <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# Model 2: Disease-severity-adjusted
m4_current_burden_count_disease <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work
)

# Extract incidence rate ratios
m4_results <- dplyr::bind_rows(
  broom::tidy(m4_current_burden_count_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  
  broom::tidy(m4_current_burden_count_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  
  broom::tidy(m4_current_burden_count_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m4_results


emm_m4 <- emmeans::emmeans(
  m4_current_burden_count_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_m4
pairs(emm_m4)

# ============================================================
# Model analytic N table
# ============================================================

model_n_table <- tibble::tibble(
  outcome = c(
    "Binary ADL decline",
    "3-level ADL decline severity",
    "Binary current ADL burden",
    "Current ADL burden count"
  ),
  model_0_n = c(
    nobs(m1_adl_decline_binary_0),
    nobs(m2_adl_decline_severity_0),
    nobs(m3_current_burden_binary_0),
    nobs(m4_current_burden_count_0)
  ),
  model_1_n = c(
    nobs(m1_adl_decline_binary_core),
    nobs(m2_adl_decline_severity_core),
    nobs(m3_current_burden_binary_core),
    nobs(m4_current_burden_count_core)
  ),
  model_2_n = c(
    nobs(m1_adl_decline_binary_disease),
    nobs(m2_adl_decline_severity_disease),
    nobs(m3_current_burden_binary_disease),
    nobs(m4_current_burden_count_disease)
  )
)

model_n_table

# ============================================================
# Cannabis terms across all outcomes and model steps
# ============================================================

all_adl_results <- dplyr::bind_rows(
  m1_results %>% dplyr::mutate(outcome = "Binary ADL decline"),
  m2_results %>% dplyr::mutate(outcome = "3-level ADL decline severity"),
  m3_results %>% dplyr::mutate(outcome = "Binary current ADL burden"),
  m4_results %>% dplyr::mutate(outcome = "Current ADL burden count")
) %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::select(outcome, model, term, estimate, conf.low, conf.high, p.value)

all_adl_results


tidy_polr_with_p <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      model = model_name
    ) %>%
    dplyr::select(model, term, estimate, conf.low, conf.high, statistic, p.value)
}

m2_results <- dplyr::bind_rows(
  tidy_polr_with_p(m2_adl_decline_severity_0, "Model 0: cannabis only"),
  tidy_polr_with_p(m2_adl_decline_severity_core, "Model 1: core adjusted"),
  tidy_polr_with_p(m2_adl_decline_severity_disease, "Model 2: disease severity adjusted")
)

m2_results


# ============================================================
# Check denominator of adl_total_calc: adl_a_no
# ============================================================

summary(adl_work$adl_a_no)

table(adl_work$adl_a_no, useNA = "ifany")

sum(adl_work$adl_a_no == 0, na.rm = TRUE)
sum(is.na(adl_work$adl_a_no))

range(adl_work$adl_a_no, na.rm = TRUE)


# ============================================================
# Check whether denominator created invalid adl_total_calc values
# ============================================================

sum(is.infinite(adl_work$adl_total_calc))
sum(is.nan(adl_work$adl_total_calc))
sum(is.na(adl_work$adl_total_calc))

adl_work %>%
  dplyr::filter(adl_a_no == 0 | is.na(adl_a_no)) %>%
  dplyr::select(
    adl_a_no,
    adl_b_sum_now,
    adl_c_sum_best,
    adl_total_calc
  )

adl_work %>%
  dplyr::summarise(
    n = dplyr::n(),
    denominator_min = min(adl_a_no, na.rm = TRUE),
    denominator_q1 = quantile(adl_a_no, .25, na.rm = TRUE),
    denominator_median = median(adl_a_no, na.rm = TRUE),
    denominator_mean = mean(adl_a_no, na.rm = TRUE),
    denominator_q3 = quantile(adl_a_no, .75, na.rm = TRUE),
    denominator_max = max(adl_a_no, na.rm = TRUE),
    denominator_zero_n = sum(adl_a_no == 0, na.rm = TRUE),
    denominator_missing_n = sum(is.na(adl_a_no))
  )




adl_work <- adl_work %>%
  dplyr::mutate(
    adl_a_no_cat3 = dplyr::case_when(
      is.na(adl_a_no) ~ NA_character_,
      adl_a_no <= 13 ~ "lower_coverage_8_13",
      adl_a_no %in% c(14, 15) ~ "moderate_coverage_14_15",
      adl_a_no == 16 ~ "full_coverage_16"
    ),
    adl_a_no_cat3 = factor(
      adl_a_no_cat3,
      levels = c(
        "lower_coverage_8_13",
        "moderate_coverage_14_15",
        "full_coverage_16"
      )
    )
  )

table(adl_work$adl_a_no_cat3, useNA = "ifany")

# Denominator category by ADL decline yes/no
table(
  scorable_items = adl_work$adl_a_no_cat3,
  any_decline = adl_work$adl_any_decline_calc,
  useNA = "ifany"
)

# Denominator category by current ADL burden yes/no
table(
  scorable_items = adl_work$adl_a_no_cat3,
  current_burden = adl_work$adl_any_current_burden,
  useNA = "ifany"
)

# Mean ADL outcomes by denominator category
adl_work %>%
  dplyr::group_by(adl_a_no_cat3) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_adl_total_calc = mean(adl_total_calc, na.rm = TRUE),
    sd_adl_total_calc = sd(adl_total_calc, na.rm = TRUE),
    median_adl_total_calc = median(adl_total_calc, na.rm = TRUE),
    any_decline_n = sum(adl_any_decline_calc == "any_decline", na.rm = TRUE),
    any_decline_prop = mean(adl_any_decline_calc == "any_decline", na.rm = TRUE),
    mean_current_burden = mean(adl_b_now_excess, na.rm = TRUE),
    sd_current_burden = sd(adl_b_now_excess, na.rm = TRUE),
    .groups = "drop"
  )


adl_work %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_adl_a_no = mean(adl_a_no, na.rm = TRUE),
    sd_adl_a_no = sd(adl_a_no, na.rm = TRUE),
    median_adl_a_no = median(adl_a_no, na.rm = TRUE),
    min_adl_a_no = min(adl_a_no, na.rm = TRUE),
    max_adl_a_no = max(adl_a_no, na.rm = TRUE),
    .groups = "drop"
  )

table(
  cannabis = adl_work$du_mar4_12m_aBin_ord,
  adl_a_no_cat3 = adl_work$adl_a_no_cat3,
  useNA = "ifany"
)


library(lmtest)
library(sandwich)
library(broom)

m_adl_a_no_0 <- lm(
  adl_a_no ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

m_adl_a_no_core <- lm(
  adl_a_no ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

m_adl_a_no_disease <- lm(
  adl_a_no ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work
)

tidy_hc3 <- function(model, model_name) {
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  as.data.frame(as.matrix(hc3)) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::rename(
      estimate = Estimate,
      std_error_hc3 = `Std. Error`,
      statistic = `t value`,
      p.value = `Pr(>|t|)`
    ) %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(model, term, estimate, std_error_hc3, statistic, p.value)
}

adl_a_no_results <- dplyr::bind_rows(
  tidy_hc3(m_adl_a_no_0, "Model 0: cannabis only"),
  tidy_hc3(m_adl_a_no_core, "Model 1: core adjusted"),
  tidy_hc3(m_adl_a_no_disease, "Model 2: disease severity adjusted")
)

adl_a_no_results


dplyr::rename(
  estimate = Estimate,
  std_error_hc3 = `Std. Error`,
  statistic = `t value`,
  p.value = `Pr(>|t|)`
)

# If you are stuck at Browse[1]>, type:
Q

# Then rerun this:
tidy_hc3 <- function(model, model_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_df <- as.data.frame(as.matrix(hc3))
  
  # Rename by position instead of relying on exact column names
  names(hc3_df) <- c(
    "estimate",
    "std_error_hc3",
    "statistic",
    "p.value"
  )
  
  hc3_df %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(model, term, estimate, std_error_hc3, statistic, p.value)
}

adl_a_no_results <- dplyr::bind_rows(
  tidy_hc3(m_adl_a_no_0, "Model 0: cannabis only"),
  tidy_hc3(m_adl_a_no_core, "Model 1: core adjusted"),
  tidy_hc3(m_adl_a_no_disease, "Model 2: disease severity adjusted")
)

adl_a_no_results


tidy_hc3 <- function(model, model_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_df <- as.data.frame(as.matrix(hc3))
  
  # Rename by position, not by exact column names
  names(hc3_df) <- c(
    "estimate",
    "std_error_hc3",
    "statistic",
    "p.value"
  )
  
  hc3_df %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(
      model,
      term,
      estimate,
      std_error_hc3,
      statistic,
      p.value
    )
}

adl_a_no_results <- dplyr::bind_rows(
  tidy_hc3(m_adl_a_no_0, "Model 0: cannabis only"),
  tidy_hc3(m_adl_a_no_core, "Model 1: core adjusted"),
  tidy_hc3(m_adl_a_no_disease, "Model 2: disease severity adjusted")
)

adl_a_no_results

dplyr::rename(
  estimate = Estimate,
  std_error_hc3 = `Std. Error`,
  statistic = `t value`,
  p.value = `Pr(>|t|)`
)

tidy_hc3 <- function(model, model_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_df <- as.data.frame(as.matrix(hc3))
  
  names(hc3_df) <- c(
    "estimate",
    "std_error_hc3",
    "statistic",
    "p.value"
  )
  
  hc3_df %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(model, term, estimate, std_error_hc3, statistic, p.value)
}

adl_a_no_results <- dplyr::bind_rows(
  tidy_hc3(m_adl_a_no_0, "Model 0: cannabis only"),
  tidy_hc3(m_adl_a_no_core, "Model 1: core adjusted"),
  tidy_hc3(m_adl_a_no_disease, "Model 2: disease severity adjusted")
)

adl_a_no_results

adl_a_no_results <- dplyr::bind_rows(
  tidy_hc3(m_adl_a_no_0, "Model 0: cannabis only"),
  tidy_hc3(m_adl_a_no_core, "Model 1: core adjusted"),
  tidy_hc3(m_adl_a_no_disease, "Model 2: disease severity adjusted")
)

adl_a_no_results

tidy_lm_hc3 <- function(model, model_name) {
  broom::tidy(
    model,
    conf.int = TRUE,
    vcov = sandwich::vcovHC(model, type = "HC3")
  ) %>%
    dplyr::mutate(model = model_name) %>%
    dplyr::select(model, dplyr::everything())
}

# ============================================================
# Functional/adherence outcome diagnostics
# mmt_ts, mac_tmr4, mac_misr
# ============================================================

library(dplyr)
library(ggplot2)
library(broom)
library(lmtest)
library(sandwich)
library(MASS)
library(emmeans)

func_vars <- c("mmt_ts", "mac_tmr4", "mac_misr")

# Missingness
sapply(adl_work[func_vars], function(x) sum(is.na(x)))

# Basic summaries
adl_work %>%
  dplyr::summarise(
    across(
      all_of(func_vars),
      list(
        n_nonmissing = ~sum(!is.na(.x)),
        missing = ~sum(is.na(.x)),
        min = ~min(.x, na.rm = TRUE),
        q1 = ~quantile(.x, .25, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        mean = ~mean(.x, na.rm = TRUE),
        q3 = ~quantile(.x, .75, na.rm = TRUE),
        max = ~max(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        zero_n = ~sum(.x == 0, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

# Frequency tables if variables are discrete/integer-like
table(adl_work$mmt_ts, useNA = "ifany")
table(adl_work$mac_tmr4, useNA = "ifany")
table(adl_work$mac_misr, useNA = "ifany")

# Distribution plots
ggplot(adl_work, aes(x = mmt_ts)) +
  geom_histogram(bins = 20) +
  theme_classic() +
  labs(x = "MMT-R total score", y = "Count")

ggplot(adl_work, aes(x = mac_tmr4)) +
  geom_histogram(bins = 10) +
  theme_classic() +
  labs(x = "Morisky-4 total", y = "Count")

ggplot(adl_work, aes(x = mac_misr)) +
  geom_histogram(bins = 20) +
  theme_classic() +
  labs(x = "Missed-dose ratio", y = "Count")

# ============================================================
# Descriptives by cannabis group
# ============================================================

adl_work %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    mmt_ts_n = sum(!is.na(mmt_ts)),
    mmt_ts_mean = mean(mmt_ts, na.rm = TRUE),
    mmt_ts_sd = sd(mmt_ts, na.rm = TRUE),
    mmt_ts_median = median(mmt_ts, na.rm = TRUE),
    mmt_ts_min = min(mmt_ts, na.rm = TRUE),
    mmt_ts_max = max(mmt_ts, na.rm = TRUE),
    
    mac_tmr4_n = sum(!is.na(mac_tmr4)),
    mac_tmr4_mean = mean(mac_tmr4, na.rm = TRUE),
    mac_tmr4_sd = sd(mac_tmr4, na.rm = TRUE),
    mac_tmr4_median = median(mac_tmr4, na.rm = TRUE),
    mac_tmr4_min = min(mac_tmr4, na.rm = TRUE),
    mac_tmr4_max = max(mac_tmr4, na.rm = TRUE),
    
    mac_misr_n = sum(!is.na(mac_misr)),
    mac_misr_mean = mean(mac_misr, na.rm = TRUE),
    mac_misr_sd = sd(mac_misr, na.rm = TRUE),
    mac_misr_median = median(mac_misr, na.rm = TRUE),
    mac_misr_min = min(mac_misr, na.rm = TRUE),
    mac_misr_max = max(mac_misr, na.rm = TRUE),
    
    .groups = "drop"
  )

# ============================================================
# HC3 helper
# ============================================================

tidy_hc3 <- function(model, model_name, outcome_name) {
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_df <- as.data.frame(as.matrix(hc3))
  
  names(hc3_df) <- c(
    "estimate",
    "std_error_hc3",
    "statistic",
    "p.value"
  )
  
  hc3_df %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(
      model = model_name,
      outcome = outcome_name
    ) %>%
    dplyr::select(
      outcome,
      model,
      term,
      estimate,
      std_error_hc3,
      statistic,
      p.value
    )
}

# ============================================================
# Outcome A: mmt_ts
# Medication-management performance
# ============================================================

m_mmt_step1_disease <- lm(
  mmt_ts ~ DiseaseSev_c,
  data = adl_work
)

m_mmt_step2_disease_cannabis <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_mmt_step3_full <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

mmt_results_hc3 <- dplyr::bind_rows(
  tidy_hc3(m_mmt_step1_disease, "Step 1: disease severity only", "MMT-R total"),
  tidy_hc3(m_mmt_step2_disease_cannabis, "Step 2: disease severity + cannabis", "MMT-R total"),
  tidy_hc3(m_mmt_step3_full, "Step 3: disease severity + cannabis + covariates", "MMT-R total")
)

mmt_results_hc3

anova(
  m_mmt_step1_disease,
  m_mmt_step2_disease_cannabis,
  m_mmt_step3_full
)

AIC(
  m_mmt_step1_disease,
  m_mmt_step2_disease_cannabis,
  m_mmt_step3_full
)

emmeans::emmeans(
  m_mmt_step3_full,
  ~ du_mar4_12m_aBin_ord
)

pairs(
  emmeans::emmeans(
    m_mmt_step3_full,
    ~ du_mar4_12m_aBin_ord
  )
)

# ============================================================
# Outcome B: mac_tmr4
# Morisky-4 total
# ============================================================

m_morisky_step1_disease <- lm(
  mac_tmr4 ~ DiseaseSev_c,
  data = adl_work
)

m_morisky_step2_disease_cannabis <- lm(
  mac_tmr4 ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_morisky_step3_full <- lm(
  mac_tmr4 ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

morisky_results_hc3 <- dplyr::bind_rows(
  tidy_hc3(m_morisky_step1_disease, "Step 1: disease severity only", "Morisky-4 total"),
  tidy_hc3(m_morisky_step2_disease_cannabis, "Step 2: disease severity + cannabis", "Morisky-4 total"),
  tidy_hc3(m_morisky_step3_full, "Step 3: disease severity + cannabis + covariates", "Morisky-4 total")
)

morisky_results_hc3

anova(
  m_morisky_step1_disease,
  m_morisky_step2_disease_cannabis,
  m_morisky_step3_full
)

AIC(
  m_morisky_step1_disease,
  m_morisky_step2_disease_cannabis,
  m_morisky_step3_full
)


# ============================================================
# Optional ordinal Morisky model
# Only use if mac_tmr4 is ordered with few categories
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    mac_tmr4_ord = factor(
      mac_tmr4,
      ordered = TRUE
    )
  )

m_morisky_ord_step3_full <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

tidy_polr_with_p <- function(model, model_name, outcome_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      model = model_name,
      outcome = outcome_name
    ) %>%
    dplyr::select(outcome, model, term, estimate, conf.low, conf.high, statistic, p.value)
}

morisky_ord_results <- tidy_polr_with_p(
  m_morisky_ord_step3_full,
  "Step 3: disease severity + cannabis + covariates",
  "Morisky-4 ordinal"
)

morisky_ord_results

# ============================================================
# Outcome C: mac_misr
# Missed-dose ratio
# ============================================================

summary(adl_work$mac_misr)

sum(adl_work$mac_misr == 0, na.rm = TRUE)
sum(adl_work$mac_misr == 1, na.rm = TRUE)
sum(adl_work$mac_misr > 0 & adl_work$mac_misr < 1, na.rm = TRUE)
sum(is.na(adl_work$mac_misr))


m_misr_step1_disease <- lm(
  mac_misr ~ DiseaseSev_c,
  data = adl_work
)

m_misr_step2_disease_cannabis <- lm(
  mac_misr ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_misr_step3_full <- lm(
  mac_misr ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

misr_results_hc3 <- dplyr::bind_rows(
  tidy_hc3(m_misr_step1_disease, "Step 1: disease severity only", "Missed-dose ratio"),
  tidy_hc3(m_misr_step2_disease_cannabis, "Step 2: disease severity + cannabis", "Missed-dose ratio"),
  tidy_hc3(m_misr_step3_full, "Step 3: disease severity + cannabis + covariates", "Missed-dose ratio")
)

misr_results_hc3

anova(
  m_misr_step1_disease,
  m_misr_step2_disease_cannabis,
  m_misr_step3_full
)

AIC(
  m_misr_step1_disease,
  m_misr_step2_disease_cannabis,
  m_misr_step3_full
)

m_misr_step1_disease <- lm(
  mac_misr ~ DiseaseSev_c,
  data = adl_work
)

m_misr_step2_disease_cannabis <- lm(
  mac_misr ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_misr_step3_full <- lm(
  mac_misr ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

misr_results_hc3 <- dplyr::bind_rows(
  tidy_hc3(m_misr_step1_disease, "Step 1: disease severity only", "Missed-dose ratio"),
  tidy_hc3(m_misr_step2_disease_cannabis, "Step 2: disease severity + cannabis", "Missed-dose ratio"),
  tidy_hc3(m_misr_step3_full, "Step 3: disease severity + cannabis + covariates", "Missed-dose ratio")
)

misr_results_hc3

anova(
  m_misr_step1_disease,
  m_misr_step2_disease_cannabis,
  m_misr_step3_full
)

AIC(
  m_misr_step1_disease,
  m_misr_step2_disease_cannabis,
  m_misr_step3_full
)


# ============================================================
# MMT-R TOTAL: Secondary functional-capacity outcome
# Outcome: mmt_ts
# Model: Linear regression with HC3 robust SEs
# Sequence: disease severity -> cannabis -> covariates
# ============================================================

library(dplyr)
library(lmtest)
library(sandwich)
library(broom)
library(emmeans)

# ------------------------------------------------------------
# Distribution diagnostics
# ------------------------------------------------------------

summary(adl_work$mmt_ts)
table(adl_work$mmt_ts, useNA = "ifany")
sum(is.na(adl_work$mmt_ts))

adl_work %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mmt_n = sum(!is.na(mmt_ts)),
    mmt_mean = mean(mmt_ts, na.rm = TRUE),
    mmt_sd = sd(mmt_ts, na.rm = TRUE),
    mmt_median = median(mmt_ts, na.rm = TRUE),
    mmt_min = min(mmt_ts, na.rm = TRUE),
    mmt_max = max(mmt_ts, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Stepwise/sequential models
# ------------------------------------------------------------

m_mmt_step1_disease <- lm(
  mmt_ts ~ DiseaseSev_c,
  data = adl_work
)

m_mmt_step2_disease_cannabis <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_mmt_step3_full <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# ------------------------------------------------------------
# HC3 robust helper
# ------------------------------------------------------------

tidy_hc3 <- function(model, model_name, outcome_name) {
  
  hc3_mat <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_mat <- as.matrix(hc3_mat)
  
  hc3_df <- data.frame(
    term = rownames(hc3_mat),
    estimate = hc3_mat[, 1],
    std_error_hc3 = hc3_mat[, 2],
    statistic = hc3_mat[, 3],
    p.value = hc3_mat[, 4],
    row.names = NULL
  )
  
  hc3_df %>%
    dplyr::mutate(
      outcome = outcome_name,
      model = model_name
    ) %>%
    dplyr::select(
      outcome,
      model,
      term,
      estimate,
      std_error_hc3,
      statistic,
      p.value
    )
}

mmt_results_hc3 <- dplyr::bind_rows(
  tidy_hc3(
    m_mmt_step1_disease,
    "Step 1: disease severity only",
    "MMT-R total"
  ),
  tidy_hc3(
    m_mmt_step2_disease_cannabis,
    "Step 2: disease severity + cannabis",
    "MMT-R total"
  ),
  tidy_hc3(
    m_mmt_step3_full,
    "Step 3: disease severity + cannabis + covariates",
    "MMT-R total"
  )
)

mmt_results_hc3

# Model comparison
anova(
  m_mmt_step1_disease,
  m_mmt_step2_disease_cannabis,
  m_mmt_step3_full
)

AIC(
  m_mmt_step1_disease,
  m_mmt_step2_disease_cannabis,
  m_mmt_step3_full
)

# Analytic N
nobs(m_mmt_step1_disease)
nobs(m_mmt_step2_disease_cannabis)
nobs(m_mmt_step3_full)

# Adjusted means by cannabis group
mmt_emmeans <- emmeans::emmeans(
  m_mmt_step3_full,
  ~ du_mar4_12m_aBin_ord
)

mmt_emmeans

pairs(mmt_emmeans)

# ============================================================
# MORISKY-4 TOTAL: Exploratory adherence/nonadherence outcome
# Outcome: mac_tmr4
# Model: Ordinal logistic regression
# Sequence: disease severity -> cannabis -> covariates
# ============================================================

library(MASS)
library(dplyr)
library(broom)
library(emmeans)

# ------------------------------------------------------------
# Distribution diagnostics
# ------------------------------------------------------------

summary(adl_work$mac_tmr4)
table(adl_work$mac_tmr4, useNA = "ifany")
sum(is.na(adl_work$mac_tmr4))

adl_work %>%
  dplyr::group_by(du_mar4_12m_aBin_ord) %>%
  dplyr::summarise(
    n = dplyr::n(),
    morisky_n = sum(!is.na(mac_tmr4)),
    morisky_mean = mean(mac_tmr4, na.rm = TRUE),
    morisky_sd = sd(mac_tmr4, na.rm = TRUE),
    morisky_median = median(mac_tmr4, na.rm = TRUE),
    morisky_min = min(mac_tmr4, na.rm = TRUE),
    morisky_max = max(mac_tmr4, na.rm = TRUE),
    .groups = "drop"
  )

adl_work <- adl_work %>%
  dplyr::mutate(
    mac_tmr4_ord = factor(
      mac_tmr4,
      ordered = TRUE
    )
  )

table(adl_work$mac_tmr4_ord, useNA = "ifany")

# ------------------------------------------------------------
# Stepwise/sequential ordinal logistic models
# ------------------------------------------------------------

m_morisky_step1_disease <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

m_morisky_step2_disease_cannabis <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

m_morisky_step3_full <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

tidy_polr_with_p <- function(model, model_name, outcome_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      outcome = outcome_name,
      model = model_name
    ) %>%
    dplyr::select(
      outcome,
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
}

morisky_results_ord <- dplyr::bind_rows(
  tidy_polr_with_p(
    m_morisky_step1_disease,
    "Step 1: disease severity only",
    "Morisky-4 total"
  ),
  tidy_polr_with_p(
    m_morisky_step2_disease_cannabis,
    "Step 2: disease severity + cannabis",
    "Morisky-4 total"
  ),
  tidy_polr_with_p(
    m_morisky_step3_full,
    "Step 3: disease severity + cannabis + covariates",
    "Morisky-4 total"
  )
) %>%
  # remove threshold/cutpoint rows
  dplyr::filter(
    !grepl("\\|", term)
  )

morisky_results_ord

# Model comparison
anova(
  m_morisky_step1_disease,
  m_morisky_step2_disease_cannabis,
  m_morisky_step3_full
)

AIC(
  m_morisky_step1_disease,
  m_morisky_step2_disease_cannabis,
  m_morisky_step3_full
)

# Analytic N
nobs(m_morisky_step1_disease)
nobs(m_morisky_step2_disease_cannabis)
nobs(m_morisky_step3_full)

# Adjusted predicted probabilities by cannabis group
morisky_probs <- emmeans::emmeans(
  m_morisky_step3_full,
  ~ du_mar4_12m_aBin_ord | mac_tmr4_ord,
  mode = "prob"
)

morisky_probs


# ============================================================
# Combined functional/adherence results
# ============================================================

functional_extension_results <- dplyr::bind_rows(
  mmt_results_hc3 %>%
    dplyr::mutate(
      estimate_type = "Linear estimate"
    ),
  
  morisky_results_ord %>%
    dplyr::mutate(
      estimate_type = "Proportional odds ratio"
    )
) %>%
  dplyr::filter(
    grepl("DiseaseSev_c|du_mar4_12m_aBin_ord|phq_2_age_c|bdi_total|sex_covfac", term)
  ) %>%
  dplyr::select(
    outcome,
    estimate_type,
    model,
    term,
    estimate,
    dplyr::any_of(c("std_error_hc3", "conf.low", "conf.high")),
    statistic,
    p.value
  )

functional_extension_results

functional_extension_n_table <- tibble::tibble(
  outcome = c(
    "MMT-R total",
    "Morisky-4 total"
  ),
  step1_disease_n = c(
    nobs(m_mmt_step1_disease),
    nobs(m_morisky_step1_disease)
  ),
  step2_disease_cannabis_n = c(
    nobs(m_mmt_step2_disease_cannabis),
    nobs(m_morisky_step2_disease_cannabis)
  ),
  step3_full_n = c(
    nobs(m_mmt_step3_full),
    nobs(m_morisky_step3_full)
  )
)

functional_extension_n_table

# ============================================================
# Identify global/domain neurocognitive score variables
# ============================================================

names(adl_work)[
  grepl(
    "global|learn|memory|exec|motor|speed|attention|domain|neuro|cog",
    names(adl_work),
    ignore.case = TRUE
  )
]

neurocog_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "executive_z",
  "motor_z",
  "processing_speed_z"
)

functional_vars <- c(
  "adl_total_calc",
  "adl_b_now_excess",
  "mmt_ts",
  "mac_tmr4"
)

# ============================================================
# Correlations: functional/adherence outcomes with neurocognition
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(broom)

# Neurocognitive variables from candidacy-style domain/global scores
neurocog_vars <- c(
  "global_z",
  "learning_z",
  "memory_z",
  "executive_z",
  "motor_z",
  "processing_speed_z"
)

# Functional/adherence outcomes
functional_vars <- c(
  "adl_total_calc",     # proportional ADL decline score
  "adl_b_now_excess",   # current ADL burden count
  "mmt_ts",             # MMT-R medication-management performance
  "mac_tmr4"            # Morisky-4 adherence/nonadherence total
)

# Check that all variables exist
setdiff(neurocog_vars, names(adl_work))
setdiff(functional_vars, names(adl_work))

# ============================================================
# Spearman correlations
# ============================================================

func_neuro_spearman <- tidyr::expand_grid(
  functional_outcome = functional_vars,
  neurocog_score = neurocog_vars
) %>%
  dplyr::mutate(
    test = purrr::map2(
      functional_outcome,
      neurocog_score,
      ~ cor.test(
        adl_work[[.x]],
        adl_work[[.y]],
        method = "spearman",
        exact = FALSE,
        use = "pairwise.complete.obs"
      )
    ),
    rho = purrr::map_dbl(test, ~ unname(.x$estimate)),
    p.value = purrr::map_dbl(test, ~ .x$p.value),
    n = purrr::map2_int(
      functional_outcome,
      neurocog_score,
      ~ sum(complete.cases(adl_work[, c(.x, .y)]))
    )
  ) %>%
  dplyr::select(
    functional_outcome,
    neurocog_score,
    n,
    rho,
    p.value
  ) %>%
  dplyr::mutate(
    p_fdr = p.adjust(p.value, method = "BH")
  ) %>%
  dplyr::arrange(functional_outcome, p.value)

func_neuro_spearman

# ============================================================
# Wide rho matrix
# ============================================================

func_neuro_rho_matrix <- func_neuro_spearman %>%
  dplyr::select(functional_outcome, neurocog_score, rho) %>%
  tidyr::pivot_wider(
    names_from = neurocog_score,
    values_from = rho
  )

func_neuro_rho_matrix

# ============================================================
# Wide p-value matrix
# ============================================================

func_neuro_p_matrix <- func_neuro_spearman %>%
  dplyr::select(functional_outcome, neurocog_score, p.value) %>%
  tidyr::pivot_wider(
    names_from = neurocog_score,
    values_from = p.value
  )

func_neuro_p_matrix

# ============================================================
# Formatted correlation table
# ============================================================

func_neuro_corr_table <- func_neuro_spearman %>%
  dplyr::mutate(
    rho_fmt = sprintf("%.2f", rho),
    p_fmt = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    ),
    p_fdr_fmt = dplyr::case_when(
      p_fdr < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p_fdr)
    ),
    result = paste0(
      "rho = ", rho_fmt,
      ", p = ", p_fmt,
      ", FDR p = ", p_fdr_fmt
    )
  ) %>%
  dplyr::select(
    functional_outcome,
    neurocog_score,
    n,
    result
  ) %>%
  dplyr::arrange(functional_outcome, neurocog_score)

func_neuro_corr_table

# ============================================================
# Optional heatmap
# ============================================================

library(ggplot2)

ggplot(
  func_neuro_spearman,
  aes(
    x = neurocog_score,
    y = functional_outcome,
    fill = rho
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sprintf("%.2f", rho)),
    size = 3
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  labs(
    x = "Neurocognitive score",
    y = "Functional/adherence outcome",
    fill = "Spearman rho"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




# ============================================================
# Packages
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(broom)
library(lmtest)
library(sandwich)
library(MASS)
library(emmeans)
library(ggplot2)

# ============================================================
# Helper: HC3 robust results for linear models
# ============================================================

tidy_hc3 <- function(model, model_name, outcome_name) {
  
  hc3_mat <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  hc3_mat <- as.matrix(hc3_mat)
  
  hc3_df <- data.frame(
    term = rownames(hc3_mat),
    estimate = hc3_mat[, 1],
    std_error_hc3 = hc3_mat[, 2],
    statistic = hc3_mat[, 3],
    p.value = hc3_mat[, 4],
    row.names = NULL
  )
  
  hc3_df %>%
    dplyr::mutate(
      outcome = outcome_name,
      model = model_name
    ) %>%
    dplyr::select(
      outcome,
      model,
      term,
      estimate,
      std_error_hc3,
      statistic,
      p.value
    )
}

# ============================================================
# Helper: ordinal logistic results with p-values
# ============================================================

tidy_polr_with_p <- function(model, model_name, outcome_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      outcome = outcome_name,
      model = model_name
    ) %>%
    dplyr::select(
      outcome,
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
}

# ============================================================
# Helper: standard GLM/NB results
# ============================================================

tidy_exp <- function(model, model_name, outcome_name, estimate_type_name) {
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      outcome = outcome_name,
      model = model_name,
      estimate_type = estimate_type_name
    ) %>%
    dplyr::select(
      outcome,
      estimate_type,
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
}

# ============================================================
# Outcome construction / distribution table
# ============================================================

outcome_distribution_table <- tibble::tibble(
  construct = c(
    "ADL decline from best functioning",
    "Current ADL burden",
    "Current ADL burden severity",
    "Medication-management capacity",
    "Adherence behavior"
  ),
  outcome = c(
    "Any ADL decline",
    "Any current ADL burden",
    "Current ADL burden count",
    "MMT-R total",
    "Morisky-4 total"
  ),
  variable = c(
    "adl_any_decline_calc",
    "adl_any_current_burden",
    "adl_b_now_excess",
    "mmt_ts",
    "mac_tmr4"
  ),
  coding_or_scale = c(
    "0 = no decline; 1 = any decline",
    "0 = no current burden; 1 = any current burden",
    "Count of current ADL burden",
    "Continuous performance score; higher = better medication-management performance",
    "Ordinal self-reported adherence/nonadherence score"
  ),
  model_type = c(
    "Logistic regression",
    "Logistic regression",
    "Negative binomial regression",
    "Linear regression with HC3 robust SEs",
    "Ordinal logistic regression"
  ),
  analytic_role = c(
    "Primary functional decline outcome",
    "Secondary current-functioning outcome",
    "Secondary current-burden severity outcome",
    "Secondary performance-based functional-capacity outcome",
    "Exploratory adherence/nonadherence outcome"
  )
)

outcome_distribution_table

# ============================================================
# Observed outcome distributions
# ============================================================

table(adl_work$adl_any_decline_calc, useNA = "ifany")
table(adl_work$adl_any_current_burden, useNA = "ifany")
summary(adl_work$adl_b_now_excess)
summary(adl_work$mmt_ts)
table(adl_work$mac_tmr4, useNA = "ifany")

outcome_summary_table <- tibble::tibble(
  outcome = c(
    "Any ADL decline",
    "Any current ADL burden",
    "Current ADL burden count",
    "MMT-R total",
    "Morisky-4 total"
  ),
  n = c(
    sum(!is.na(adl_work$adl_any_decline_calc)),
    sum(!is.na(adl_work$adl_any_current_burden)),
    sum(!is.na(adl_work$adl_b_now_excess)),
    sum(!is.na(adl_work$mmt_ts)),
    sum(!is.na(adl_work$mac_tmr4))
  ),
  distribution = c(
    paste0(
      "No decline = ",
      sum(adl_work$adl_any_decline_calc == "no_decline", na.rm = TRUE),
      "; any decline = ",
      sum(adl_work$adl_any_decline_calc == "any_decline", na.rm = TRUE)
    ),
    paste0(
      "No current burden = ",
      sum(adl_work$adl_any_current_burden == "no_current_burden", na.rm = TRUE),
      "; any current burden = ",
      sum(adl_work$adl_any_current_burden == "any_current_burden", na.rm = TRUE)
    ),
    paste0(
      "Mean = ",
      round(mean(adl_work$adl_b_now_excess, na.rm = TRUE), 2),
      "; SD = ",
      round(sd(adl_work$adl_b_now_excess, na.rm = TRUE), 2),
      "; range = ",
      min(adl_work$adl_b_now_excess, na.rm = TRUE),
      "-",
      max(adl_work$adl_b_now_excess, na.rm = TRUE)
    ),
    paste0(
      "Mean = ",
      round(mean(adl_work$mmt_ts, na.rm = TRUE), 2),
      "; SD = ",
      round(sd(adl_work$mmt_ts, na.rm = TRUE), 2),
      "; range = ",
      min(adl_work$mmt_ts, na.rm = TRUE),
      "-",
      max(adl_work$mmt_ts, na.rm = TRUE)
    ),
    paste0(
      paste(names(table(adl_work$mac_tmr4)), table(adl_work$mac_tmr4), sep = " = "),
      collapse = "; "
    )
  )
)

outcome_summary_table

# ============================================================
# Final model 1: Any ADL decline
# ============================================================

m_adl_decline_final <- glm(
  adl_any_decline_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ============================================================
# Final model 2: Any current ADL burden
# ============================================================

m_current_burden_binary_final <- glm(
  adl_any_current_burden ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

# ============================================================
# Final model 3: Current ADL burden count
# ============================================================

m_current_burden_count_final <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# ============================================================
# Final model 4: MMT-R total
# ============================================================

m_mmt_final <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# ============================================================
# Final model 5: Morisky-4 total
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    mac_tmr4_ord = factor(
      mac_tmr4,
      ordered = TRUE
    )
  )

m_morisky_final <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

# ============================================================
# Final adjusted model table
# ============================================================

final_model_results <- dplyr::bind_rows(
  
  tidy_exp(
    m_adl_decline_final,
    "Final adjusted model",
    "Any ADL decline",
    "Odds ratio"
  ),
  
  tidy_exp(
    m_current_burden_binary_final,
    "Final adjusted model",
    "Any current ADL burden",
    "Odds ratio"
  ),
  
  tidy_exp(
    m_current_burden_count_final,
    "Final adjusted model",
    "Current ADL burden count",
    "Incidence rate ratio"
  ),
  
  tidy_hc3(
    m_mmt_final,
    "Final adjusted model",
    "MMT-R total"
  ) %>%
    dplyr::mutate(
      estimate_type = "Linear estimate",
      conf.low = NA_real_,
      conf.high = NA_real_
    ) %>%
    dplyr::select(
      outcome,
      estimate_type,
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    ),
  
  tidy_polr_with_p(
    m_morisky_final,
    "Final adjusted model",
    "Morisky-4 total"
  ) %>%
    dplyr::mutate(
      estimate_type = "Proportional odds ratio"
    ) %>%
    dplyr::select(
      outcome,
      estimate_type,
      model,
      term,
      estimate,
      conf.low,
      conf.high,
      statistic,
      p.value
    )
) %>%
  dplyr::filter(
    !grepl("\\|", term)
  ) %>%
  dplyr::filter(
    grepl(
      "DiseaseSev_c|du_mar4_12m_aBin_ord|phq_2_age_c|bdi_total|sex_covfac",
      term
    )
  )

final_model_results

# ============================================================
# APA-ish formatted final model table
# ============================================================

final_model_results_formatted <- final_model_results %>%
  dplyr::mutate(
    estimate_fmt = sprintf("%.2f", estimate),
    ci_fmt = dplyr::case_when(
      is.na(conf.low) ~ "",
      TRUE ~ paste0("[", sprintf("%.2f", conf.low), ", ", sprintf("%.2f", conf.high), "]")
    ),
    p_fmt = dplyr::case_when(
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    outcome,
    estimate_type,
    term,
    estimate_fmt,
    ci_fmt,
    p_fmt
  )

final_model_results_formatted

# ============================================================
# Sequential models: Any ADL decline
# ============================================================

m_adl_decline_s1 <- glm(
  adl_any_decline_calc ~ DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m_adl_decline_s2 <- glm(
  adl_any_decline_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m_adl_decline_s3 <- m_adl_decline_final

# ============================================================
# Sequential models: Any current ADL burden
# ============================================================

m_current_binary_s1 <- glm(
  adl_any_current_burden ~ DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m_current_binary_s2 <- glm(
  adl_any_current_burden ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m_current_binary_s3 <- m_current_burden_binary_final

# ============================================================
# Sequential models: Current ADL burden count
# ============================================================

m_current_count_s1 <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c,
  data = adl_work
)

m_current_count_s2 <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_current_count_s3 <- m_current_burden_count_final

# ============================================================
# Sequential models: MMT-R total
# ============================================================

m_mmt_s1 <- lm(
  mmt_ts ~ DiseaseSev_c,
  data = adl_work
)

m_mmt_s2 <- lm(
  mmt_ts ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m_mmt_s3 <- m_mmt_final

# ============================================================
# Sequential models: Morisky-4 total
# ============================================================

m_morisky_s1 <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

m_morisky_s2 <- MASS::polr(
  mac_tmr4_ord ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

m_morisky_s3 <- m_morisky_final

# ============================================================
# Sequential coefficient results
# ============================================================

sequential_results <- dplyr::bind_rows(
  
  # ADL decline
  tidy_exp(m_adl_decline_s1, "Step 1: disease severity only", "Any ADL decline", "Odds ratio"),
  tidy_exp(m_adl_decline_s2, "Step 2: disease severity + cannabis", "Any ADL decline", "Odds ratio"),
  tidy_exp(m_adl_decline_s3, "Step 3: full adjustment", "Any ADL decline", "Odds ratio"),
  
  # Current burden binary
  tidy_exp(m_current_binary_s1, "Step 1: disease severity only", "Any current ADL burden", "Odds ratio"),
  tidy_exp(m_current_binary_s2, "Step 2: disease severity + cannabis", "Any current ADL burden", "Odds ratio"),
  tidy_exp(m_current_binary_s3, "Step 3: full adjustment", "Any current ADL burden", "Odds ratio"),
  
  # Current burden count
  tidy_exp(m_current_count_s1, "Step 1: disease severity only", "Current ADL burden count", "Incidence rate ratio"),
  tidy_exp(m_current_count_s2, "Step 2: disease severity + cannabis", "Current ADL burden count", "Incidence rate ratio"),
  tidy_exp(m_current_count_s3, "Step 3: full adjustment", "Current ADL burden count", "Incidence rate ratio"),
  
  # MMT-R
  tidy_hc3(m_mmt_s1, "Step 1: disease severity only", "MMT-R total") %>%
    dplyr::mutate(
      estimate_type = "Linear estimate",
      conf.low = NA_real_,
      conf.high = NA_real_
    ) %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value),
  
  tidy_hc3(m_mmt_s2, "Step 2: disease severity + cannabis", "MMT-R total") %>%
    dplyr::mutate(
      estimate_type = "Linear estimate",
      conf.low = NA_real_,
      conf.high = NA_real_
    ) %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value),
  
  tidy_hc3(m_mmt_s3, "Step 3: full adjustment", "MMT-R total") %>%
    dplyr::mutate(
      estimate_type = "Linear estimate",
      conf.low = NA_real_,
      conf.high = NA_real_
    ) %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value),
  
  # Morisky
  tidy_polr_with_p(m_morisky_s1, "Step 1: disease severity only", "Morisky-4 total") %>%
    dplyr::mutate(estimate_type = "Proportional odds ratio") %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value),
  
  tidy_polr_with_p(m_morisky_s2, "Step 2: disease severity + cannabis", "Morisky-4 total") %>%
    dplyr::mutate(estimate_type = "Proportional odds ratio") %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value),
  
  tidy_polr_with_p(m_morisky_s3, "Step 3: full adjustment", "Morisky-4 total") %>%
    dplyr::mutate(estimate_type = "Proportional odds ratio") %>%
    dplyr::select(outcome, estimate_type, model, term, estimate, conf.low, conf.high, statistic, p.value)
) %>%
  dplyr::filter(!grepl("\\|", term)) %>%
  dplyr::filter(
    grepl("DiseaseSev_c|du_mar4_12m_aBin_ord", term)
  )

sequential_results

# ============================================================
# Sequential model fit table
# ============================================================

model_fit_table <- tibble::tibble(
  outcome = c(
    rep("Any ADL decline", 3),
    rep("Any current ADL burden", 3),
    rep("Current ADL burden count", 3),
    rep("MMT-R total", 3),
    rep("Morisky-4 total", 3)
  ),
  model = rep(
    c(
      "Step 1: disease severity only",
      "Step 2: disease severity + cannabis",
      "Step 3: full adjustment"
    ),
    times = 5
  ),
  n = c(
    nobs(m_adl_decline_s1),
    nobs(m_adl_decline_s2),
    nobs(m_adl_decline_s3),
    
    nobs(m_current_binary_s1),
    nobs(m_current_binary_s2),
    nobs(m_current_binary_s3),
    
    nobs(m_current_count_s1),
    nobs(m_current_count_s2),
    nobs(m_current_count_s3),
    
    nobs(m_mmt_s1),
    nobs(m_mmt_s2),
    nobs(m_mmt_s3),
    
    nobs(m_morisky_s1),
    nobs(m_morisky_s2),
    nobs(m_morisky_s3)
  ),
  AIC = c(
    AIC(m_adl_decline_s1),
    AIC(m_adl_decline_s2),
    AIC(m_adl_decline_s3),
    
    AIC(m_current_binary_s1),
    AIC(m_current_binary_s2),
    AIC(m_current_binary_s3),
    
    AIC(m_current_count_s1),
    AIC(m_current_count_s2),
    AIC(m_current_count_s3),
    
    AIC(m_mmt_s1),
    AIC(m_mmt_s2),
    AIC(m_mmt_s3),
    
    AIC(m_morisky_s1),
    AIC(m_morisky_s2),
    AIC(m_morisky_s3)
  )
)

model_fit_table

# ============================================================
# Model comparison tests
# ============================================================

anova(m_adl_decline_s1, m_adl_decline_s2, m_adl_decline_s3, test = "Chisq")
anova(m_current_binary_s1, m_current_binary_s2, m_current_binary_s3, test = "Chisq")
anova(m_current_count_s1, m_current_count_s2, m_current_count_s3, test = "Chisq")

anova(m_mmt_s1, m_mmt_s2, m_mmt_s3)

anova(m_morisky_s1, m_morisky_s2, m_morisky_s3)

# ============================================================
# Morisky predicted probabilities by cannabis group
# Final adjusted ordinal model
# ============================================================

morisky_probs <- emmeans::emmeans(
  m_morisky_final,
  ~ du_mar4_12m_aBin_ord | mac_tmr4_ord,
  mode = "prob"
)

morisky_probs


morisky_probs_df <- as.data.frame(morisky_probs) %>%
  dplyr::rename(
    cannabis_group = du_mar4_12m_aBin_ord,
    morisky_score = mac_tmr4_ord,
    predicted_probability = prob,
    se = SE,
    ci_lower = asymp.LCL,
    ci_upper = asymp.UCL
  ) %>%
  dplyr::mutate(
    cannabis_group = factor(
      cannabis_group,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High")
    ),
    morisky_score = factor(
      morisky_score,
      levels = sort(unique(as.numeric(as.character(morisky_score)))),
      ordered = TRUE
    )
  )

morisky_probs_df

fig_morisky_probs_facet <- ggplot(
  morisky_probs_df,
  aes(
    x = cannabis_group,
    y = predicted_probability,
    group = 1
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  facet_wrap(~ morisky_score, nrow = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.20),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Past-year cannabis exposure",
    y = "Adjusted predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

fig_morisky_probs_facet

fig_morisky_probs_lines <- ggplot(
  morisky_probs_df,
  aes(
    x = cannabis_group,
    y = predicted_probability,
    group = morisky_score,
    linetype = morisky_score,
    shape = morisky_score
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.10),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Past-year cannabis exposure",
    y = "Adjusted predicted probability",
    linetype = "Morisky-4 score",
    shape = "Morisky-4 score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

fig_morisky_probs_lines

# ============================================================
# Descriptive ADL decline severity by cannabis group
# ============================================================

adl_severity_cannabis_tab <- table(
  cannabis = adl_work$du_mar4_12m_aBin_ord,
  adl_severity = adl_work$adl_decline_severity_3cat_calc
)

adl_severity_cannabis_tab

round(
  prop.table(adl_severity_cannabis_tab, margin = 1) * 100,
  1
)

chi_obj <- chisq.test(adl_severity_cannabis_tab)

cramers_v <- sqrt(
  as.numeric(chi_obj$statistic) /
    (
      sum(adl_severity_cannabis_tab) *
        (min(dim(adl_severity_cannabis_tab)) - 1)
    )
)

adl_severity_chi_summary <- tibble::tibble(
  test = "Cannabis group x ADL decline severity",
  chi_square = as.numeric(chi_obj$statistic),
  df = as.numeric(chi_obj$parameter),
  p_value = chi_obj$p.value,
  cramers_v = cramers_v
)

adl_severity_chi_summary


       


# ============================================================
# ADL OUTCOME ANALYSES
# Hierarchical models + pre-model feasibility/power diagnostics
# ============================================================

library(dplyr)
library(broom)
library(MASS)
library(emmeans)
library(lmtest)
library(sandwich)

# ============================================================
# 0. Create final ADL outcomes
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    
    # Proportional ADL decline score
    adl_total_calc = (adl_b_sum_now - adl_c_sum_best) / adl_a_no,
    
    # Model 1 outcome: any ADL decline yes/no
    adl_any_decline_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 ~ "any_decline"
    ),
    adl_any_decline_calc = factor(
      adl_any_decline_calc,
      levels = c("no_decline", "any_decline")
    ),
    
    # Current ADL burden count
    # adl_b_sum_now minimum is 1, so subtract 1 so 0 = no current burden
    adl_b_now_excess = adl_b_sum_now - 1,
    
    # Model 3 outcome: any current ADL burden yes/no
    adl_any_current_burden = dplyr::case_when(
      is.na(adl_b_now_excess) ~ NA_character_,
      adl_b_now_excess == 0 ~ "no_current_burden",
      adl_b_now_excess > 0 ~ "any_current_burden"
    ),
    adl_any_current_burden = factor(
      adl_any_current_burden,
      levels = c("no_current_burden", "any_current_burden")
    )
  )

# Median among people with any decline
positive_median_adl_calc <- median(
  adl_work$adl_total_calc[adl_work$adl_total_calc > 0],
  na.rm = TRUE
)

positive_median_adl_calc

# Model 2 outcome: no/lower/greater decline severity
adl_work <- adl_work %>%
  dplyr::mutate(
    adl_decline_severity_3cat_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 & adl_total_calc <= positive_median_adl_calc ~ "lower_decline",
      adl_total_calc > positive_median_adl_calc ~ "greater_decline"
    ),
    adl_decline_severity_3cat_calc = factor(
      adl_decline_severity_3cat_calc,
      levels = c("no_decline", "lower_decline", "greater_decline"),
      ordered = TRUE
    )
  )

# ============================================================
# 1. Outcome distributions
# ============================================================

table(adl_work$adl_any_decline_calc, useNA = "ifany")
table(adl_work$adl_decline_severity_3cat_calc, useNA = "ifany")
table(adl_work$adl_any_current_burden, useNA = "ifany")
table(adl_work$adl_b_now_excess, useNA = "ifany")

# Cross-tabs by cannabis group
table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_any_decline_calc, useNA = "ifany")
table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_decline_severity_3cat_calc, useNA = "ifany")
table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_any_current_burden, useNA = "ifany")

# Proportions by cannabis group
prop.table(
  table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_any_decline_calc),
  margin = 1
)

prop.table(
  table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_decline_severity_3cat_calc),
  margin = 1
)

prop.table(
  table(adl_work$du_mar4_12m_aBin_ord, adl_work$adl_any_current_burden),
  margin = 1
)

# ============================================================
# 2. Feasibility checks: events per parameter
# ============================================================

# Number of model df/predictor parameters in final model
# Cannabis ordered factor contributes 2 df: linear + quadratic
# plus age, BDI, sex, DiseaseSev = 4 more
# total approximate predictor df = 6

final_model_df <- 6

# Binary decline outcome
binary_decline_counts <- table(adl_work$adl_any_decline_calc)
binary_decline_epv <- min(binary_decline_counts) / final_model_df

binary_decline_counts
binary_decline_epv

# 3-level decline severity outcome
severity_counts <- table(adl_work$adl_decline_severity_3cat_calc)
severity_min_cell_per_df <- min(severity_counts) / final_model_df

severity_counts
severity_min_cell_per_df

# Binary current burden outcome
current_burden_counts <- table(adl_work$adl_any_current_burden)
current_burden_epv <- min(current_burden_counts) / final_model_df

current_burden_counts
current_burden_epv

# Current burden count outcome distribution
summary(adl_work$adl_b_now_excess)
mean(adl_work$adl_b_now_excess, na.rm = TRUE)
var(adl_work$adl_b_now_excess, na.rm = TRUE)
var(adl_work$adl_b_now_excess, na.rm = TRUE) /
  mean(adl_work$adl_b_now_excess, na.rm = TRUE)

# Make one compact feasibility table
feasibility_table <- tibble::tibble(
  outcome = c(
    "Binary ADL decline",
    "3-level ADL decline severity",
    "Binary current ADL burden",
    "Current ADL burden count"
  ),
  n = c(
    sum(!is.na(adl_work$adl_any_decline_calc)),
    sum(!is.na(adl_work$adl_decline_severity_3cat_calc)),
    sum(!is.na(adl_work$adl_any_current_burden)),
    sum(!is.na(adl_work$adl_b_now_excess))
  ),
  smallest_cell_or_event_n = c(
    min(binary_decline_counts),
    min(severity_counts),
    min(current_burden_counts),
    sum(adl_work$adl_b_now_excess == 0, na.rm = TRUE)
  ),
  model_df = final_model_df,
  smallest_cell_per_df = c(
    binary_decline_epv,
    severity_min_cell_per_df,
    current_burden_epv,
    NA_real_
  )
)

feasibility_table

# ============================================================
# MODEL 1: Binary ADL decline yes/no
# Outcome: adl_any_decline_calc
# ============================================================

m1_adl_decline_binary_0 <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_adl_decline_binary_core <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_adl_decline_binary_disease <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_results <- dplyr::bind_rows(
  broom::tidy(m1_adl_decline_binary_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  broom::tidy(m1_adl_decline_binary_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  broom::tidy(m1_adl_decline_binary_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m1_results

emmeans::emmeans(
  m1_adl_decline_binary_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

pairs(
  emmeans::emmeans(
    m1_adl_decline_binary_disease,
    ~ du_mar4_12m_aBin_ord,
    type = "response"
  )
)

# ============================================================
# MODEL 2: 3-level ADL decline severity
# Outcome: adl_decline_severity_3cat_calc
# ============================================================

m2_adl_decline_severity_0 <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

m2_adl_decline_severity_core <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

m2_adl_decline_severity_disease <- MASS::polr(
  adl_decline_severity_3cat_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

# polr helper with p-values
tidy_polr_with_p <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      model = model_name
    ) %>%
    dplyr::select(model, term, estimate, conf.low, conf.high, statistic, p.value)
}

m2_results <- dplyr::bind_rows(
  tidy_polr_with_p(m2_adl_decline_severity_0, "Model 0: cannabis only"),
  tidy_polr_with_p(m2_adl_decline_severity_core, "Model 1: core adjusted"),
  tidy_polr_with_p(m2_adl_decline_severity_disease, "Model 2: disease severity adjusted")
)

# Keep predictor rows only; remove threshold/cutpoint rows
m2_results_predictors <- m2_results %>%
  dplyr::filter(
    !grepl("no_decline|lower_decline", term)
  )

m2_results_predictors

# Predicted probabilities by cannabis group
emmeans::emmeans(
  m2_adl_decline_severity_disease,
  ~ du_mar4_12m_aBin_ord,
  mode = "prob"
)

# ============================================================
# MODEL 3: Binary current ADL burden yes/no
# Outcome: adl_any_current_burden
# ============================================================

m3_current_burden_binary_0 <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_current_burden_binary_core <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_current_burden_binary_disease <- glm(
  adl_any_current_burden ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_results <- dplyr::bind_rows(
  broom::tidy(m3_current_burden_binary_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  broom::tidy(m3_current_burden_binary_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  broom::tidy(m3_current_burden_binary_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m3_results

emmeans::emmeans(
  m3_current_burden_binary_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

pairs(
  emmeans::emmeans(
    m3_current_burden_binary_disease,
    ~ du_mar4_12m_aBin_ord,
    type = "response"
  )
)

# ============================================================
# MODEL 4: Current ADL burden count
# Outcome: adl_b_now_excess
# ============================================================

# First fit Poisson models to test overdispersion
m4_current_burden_pois_0 <- glm(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = poisson(link = "log")
)

m4_current_burden_pois_core <- glm(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = poisson(link = "log")
)

m4_current_burden_pois_disease <- glm(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = poisson(link = "log")
)

check_overdispersion <- function(model) {
  pearson_chisq <- sum(residuals(model, type = "pearson")^2)
  df <- df.residual(model)
  ratio <- pearson_chisq / df
  p <- pchisq(pearson_chisq, df = df, lower.tail = FALSE)
  
  data.frame(
    pearson_chisq = pearson_chisq,
    df = df,
    dispersion_ratio = ratio,
    p_value = p
  )
}

check_overdispersion(m4_current_burden_pois_0)
check_overdispersion(m4_current_burden_pois_core)
check_overdispersion(m4_current_burden_pois_disease)

# Negative binomial models
m4_current_burden_nb_0 <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

m4_current_burden_nb_core <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

m4_current_burden_nb_disease <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work
)

# Compare Poisson vs negative binomial
AIC(
  m4_current_burden_pois_0, m4_current_burden_nb_0,
  m4_current_burden_pois_core, m4_current_burden_nb_core,
  m4_current_burden_pois_disease, m4_current_burden_nb_disease
)

BIC(
  m4_current_burden_pois_0, m4_current_burden_nb_0,
  m4_current_burden_pois_core, m4_current_burden_nb_core,
  m4_current_burden_pois_disease, m4_current_burden_nb_disease
)

# Extract incidence rate ratios from negative binomial models
m4_results <- dplyr::bind_rows(
  broom::tidy(m4_current_burden_nb_0, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 0: cannabis only"),
  broom::tidy(m4_current_burden_nb_core, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 1: core adjusted"),
  broom::tidy(m4_current_burden_nb_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Model 2: disease severity adjusted")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m4_results

emmeans::emmeans(
  m4_current_burden_nb_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

pairs(
  emmeans::emmeans(
    m4_current_burden_nb_disease,
    ~ du_mar4_12m_aBin_ord,
    type = "response"
  )
)

# ============================================================
# Combined cannabis effects across outcomes
# ============================================================

all_adl_results <- dplyr::bind_rows(
  m1_results %>%
    dplyr::mutate(
      outcome = "Binary ADL decline",
      estimate_type = "Odds ratio"
    ),
  
  m2_results_predictors %>%
    dplyr::mutate(
      outcome = "3-level ADL decline severity",
      estimate_type = "Proportional odds ratio"
    ),
  
  m3_results %>%
    dplyr::mutate(
      outcome = "Binary current ADL burden",
      estimate_type = "Odds ratio"
    ),
  
  m4_results %>%
    dplyr::mutate(
      outcome = "Current ADL burden count",
      estimate_type = "Incidence rate ratio"
    )
) %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
  dplyr::select(
    outcome,
    estimate_type,
    model,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

all_adl_results

# ============================================================
# Analytic N by outcome and model step
# ============================================================

model_n_table <- tibble::tibble(
  outcome = c(
    "Binary ADL decline",
    "3-level ADL decline severity",
    "Binary current ADL burden",
    "Current ADL burden count"
  ),
  model_0_n = c(
    nobs(m1_adl_decline_binary_0),
    nobs(m2_adl_decline_severity_0),
    nobs(m3_current_burden_binary_0),
    nobs(m4_current_burden_nb_0)
  ),
  model_1_n = c(
    nobs(m1_adl_decline_binary_core),
    nobs(m2_adl_decline_severity_core),
    nobs(m3_current_burden_binary_core),
    nobs(m4_current_burden_nb_core)
  ),
  model_2_n = c(
    nobs(m1_adl_decline_binary_disease),
    nobs(m2_adl_decline_severity_disease),
    nobs(m3_current_burden_binary_disease),
    nobs(m4_current_burden_nb_disease)
  )
)

model_n_table

emmeans::emmeans(
  m2_adl_decline_severity_disease,
  ~ du_mar4_12m_aBin_ord | adl_decline_severity_3cat_calc,
  mode = "prob"
)

# ============================================================
# OUTCOME 1: Binary ADL decline yes/no
# Step 1: Disease severity
# Step 2: Add cannabis
# Step 3: Add covariates
# ============================================================

m1_decline_step1_disease <- glm(
  adl_any_decline_calc ~ DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_decline_step2_disease_cannabis <- glm(
  adl_any_decline_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_decline_step3_full <- glm(
  adl_any_decline_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

m1_step_results <- dplyr::bind_rows(
  broom::tidy(m1_decline_step1_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 1: disease severity only"),
  
  broom::tidy(m1_decline_step2_disease_cannabis, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 2: disease severity + cannabis"),
  
  broom::tidy(m1_decline_step3_full, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 3: disease severity + cannabis + covariates")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m1_step_results

anova(
  m1_decline_step1_disease,
  m1_decline_step2_disease_cannabis,
  m1_decline_step3_full,
  test = "Chisq"
)

AIC(
  m1_decline_step1_disease,
  m1_decline_step2_disease_cannabis,
  m1_decline_step3_full
)


# ============================================================
# OUTCOME 2: 3-level ADL decline severity
# ============================================================

m2_severity_step1_disease <- MASS::polr(
  adl_decline_severity_3cat_calc ~ DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

m2_severity_step2_disease_cannabis <- MASS::polr(
  adl_decline_severity_3cat_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

m2_severity_step3_full <- MASS::polr(
  adl_decline_severity_3cat_calc ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

tidy_polr_with_p <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      model = model_name
    ) %>%
    dplyr::select(model, term, estimate, conf.low, conf.high, statistic, p.value)
}

m2_step_results <- dplyr::bind_rows(
  tidy_polr_with_p(m2_severity_step1_disease, "Step 1: disease severity only"),
  tidy_polr_with_p(m2_severity_step2_disease_cannabis, "Step 2: disease severity + cannabis"),
  tidy_polr_with_p(m2_severity_step3_full, "Step 3: disease severity + cannabis + covariates")
) %>%
  dplyr::filter(!grepl("no_decline|lower_decline", term))

m2_step_results

anova(
  m2_severity_step1_disease,
  m2_severity_step2_disease_cannabis,
  m2_severity_step3_full
)

AIC(
  m2_severity_step1_disease,
  m2_severity_step2_disease_cannabis,
  m2_severity_step3_full
)

emmeans::emmeans(
  m2_severity_step3_full,
  ~ du_mar4_12m_aBin_ord | adl_decline_severity_3cat_calc,
  mode = "prob"
)

# ============================================================
# OUTCOME 3: Binary current ADL burden yes/no
# ============================================================

m3_current_step1_disease <- glm(
  adl_any_current_burden ~ DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_current_step2_disease_cannabis <- glm(
  adl_any_current_burden ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_current_step3_full <- glm(
  adl_any_current_burden ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

m3_step_results <- dplyr::bind_rows(
  broom::tidy(m3_current_step1_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 1: disease severity only"),
  
  broom::tidy(m3_current_step2_disease_cannabis, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 2: disease severity + cannabis"),
  
  broom::tidy(m3_current_step3_full, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 3: disease severity + cannabis + covariates")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m3_step_results

anova(
  m3_current_step1_disease,
  m3_current_step2_disease_cannabis,
  m3_current_step3_full,
  test = "Chisq"
)

AIC(
  m3_current_step1_disease,
  m3_current_step2_disease_cannabis,
  m3_current_step3_full
)

# ============================================================
# OUTCOME 4: Current ADL burden count
# Negative binomial
# ============================================================

m4_count_step1_disease <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c,
  data = adl_work
)

m4_count_step2_disease_cannabis <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work
)

m4_count_step3_full <- MASS::glm.nb(
  adl_b_now_excess ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

m4_step_results <- dplyr::bind_rows(
  broom::tidy(m4_count_step1_disease, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 1: disease severity only"),
  
  broom::tidy(m4_count_step2_disease_cannabis, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 2: disease severity + cannabis"),
  
  broom::tidy(m4_count_step3_full, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(model = "Step 3: disease severity + cannabis + covariates")
) %>%
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

m4_step_results

anova(
  m4_count_step1_disease,
  m4_count_step2_disease_cannabis,
  m4_count_step3_full,
  test = "Chisq"
)

AIC(
  m4_count_step1_disease,
  m4_count_step2_disease_cannabis,
  m4_count_step3_full
)

# ============================================================
# Combined disease severity and cannabis results
# ============================================================

combined_step_results <- dplyr::bind_rows(
  m1_step_results %>%
    dplyr::mutate(
      outcome = "Binary ADL decline",
      estimate_type = "Odds ratio"
    ),
  
  m2_step_results %>%
    dplyr::mutate(
      outcome = "3-level ADL decline severity",
      estimate_type = "Proportional odds ratio"
    ),
  
  m3_step_results %>%
    dplyr::mutate(
      outcome = "Binary current ADL burden",
      estimate_type = "Odds ratio"
    ),
  
  m4_step_results %>%
    dplyr::mutate(
      outcome = "Current ADL burden count",
      estimate_type = "Incidence rate ratio"
    )
) %>%
  dplyr::filter(
    grepl("DiseaseSev_c|du_mar4_12m_aBin_ord", term)
  ) %>%
  dplyr::select(
    outcome,
    estimate_type,
    model,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )

combined_step_results

model_n_step_table <- tibble::tibble(
  outcome = c(
    "Binary ADL decline",
    "3-level ADL decline severity",
    "Binary current ADL burden",
    "Current ADL burden count"
  ),
  step1_disease_n = c(
    nobs(m1_decline_step1_disease),
    nobs(m2_severity_step1_disease),
    nobs(m3_current_step1_disease),
    nobs(m4_count_step1_disease)
  ),
  step2_disease_cannabis_n = c(
    nobs(m1_decline_step2_disease_cannabis),
    nobs(m2_severity_step2_disease_cannabis),
    nobs(m3_current_step2_disease_cannabis),
    nobs(m4_count_step2_disease_cannabis)
  ),
  step3_full_n = c(
    nobs(m1_decline_step3_full),
    nobs(m2_severity_step3_full),
    nobs(m3_current_step3_full),
    nobs(m4_count_step3_full)
  )
)


# ============================================================
# Ordinal current ADL burden severity
# 0 = no current burden
# 1 = lower current burden among those with burden
# 2 = higher current burden among those with burden
# ============================================================

current_burden_positive_median <- median(
  adl_work$adl_b_now_excess[adl_work$adl_b_now_excess > 0],
  na.rm = TRUE
)

current_burden_positive_median

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_current_burden_3cat = dplyr::case_when(
      is.na(adl_b_now_excess) ~ NA_character_,
      adl_b_now_excess == 0 ~ "no_current_burden",
      adl_b_now_excess > 0 &
        adl_b_now_excess <= current_burden_positive_median ~ "lower_current_burden",
      adl_b_now_excess > current_burden_positive_median ~ "higher_current_burden"
    ),
    adl_current_burden_3cat = factor(
      adl_current_burden_3cat,
      levels = c(
        "no_current_burden",
        "lower_current_burden",
        "higher_current_burden"
      ),
      ordered = TRUE
    )
  )

table(adl_work$adl_current_burden_3cat, useNA = "ifany")

# ============================================================
# Ordinal current ADL burden severity model
# ============================================================

m_current_ord_step1_disease <- MASS::polr(
  adl_current_burden_3cat ~ DiseaseSev_c,
  data = adl_work,
  Hess = TRUE
)

m_current_ord_step2_disease_cannabis <- MASS::polr(
  adl_current_burden_3cat ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord,
  data = adl_work,
  Hess = TRUE
)

m_current_ord_step3_full <- MASS::polr(
  adl_current_burden_3cat ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

tidy_polr_with_p <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      p.value = ifelse(
        !is.na(statistic),
        2 * pnorm(abs(statistic), lower.tail = FALSE),
        NA_real_
      ),
      model = model_name
    ) %>%
    dplyr::select(model, term, estimate, conf.low, conf.high, statistic, p.value)
}

current_ord_results <- dplyr::bind_rows(
  tidy_polr_with_p(m_current_ord_step1_disease, "Step 1: disease severity only"),
  tidy_polr_with_p(m_current_ord_step2_disease_cannabis, "Step 2: disease severity + cannabis"),
  tidy_polr_with_p(m_current_ord_step3_full, "Step 3: disease severity + cannabis + covariates")
) %>%
  dplyr::filter(
    !grepl("no_current_burden|lower_current_burden", term)
  )

current_ord_results

emmeans::emmeans(
  m_current_ord_step3_full,
  ~ du_mar4_12m_aBin_ord | adl_current_burden_3cat,
  mode = "prob"
)

# ============================================================
# FIGURE: Adjusted predicted probabilities for ADL decline severity
# by cannabis group
# ============================================================

library(dplyr)
library(ggplot2)
library(emmeans)

# Estimated marginal probabilities from final ordinal model
emm_decline_severity_fig <- emmeans::emmeans(
  m2_severity_step3_full,
  ~ du_mar4_12m_aBin_ord | adl_decline_severity_3cat_calc,
  mode = "prob"
)

emm_decline_severity_fig

decline_severity_fig_df <- as.data.frame(emm_decline_severity_fig) %>%
  dplyr::rename(
    cannabis_group = du_mar4_12m_aBin_ord,
    decline_severity = adl_decline_severity_3cat_calc,
    predicted_probability = prob,
    se = SE,
    ci_lower = asymp.LCL,
    ci_upper = asymp.UCL
  ) %>%
  dplyr::mutate(
    cannabis_group = factor(
      cannabis_group,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High")
    ),
    decline_severity = factor(
      decline_severity,
      levels = c("no_decline", "lower_decline", "greater_decline"),
      labels = c("No decline", "Lower decline", "Greater decline")
    )
  )

decline_severity_fig_df

fig_decline_severity_probs <- ggplot(
  decline_severity_fig_df,
  aes(
    x = cannabis_group,
    y = predicted_probability,
    group = decline_severity,
    linetype = decline_severity,
    shape = decline_severity
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.10),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Past-year cannabis exposure",
    y = "Adjusted predicted probability",
    linetype = "ADL decline severity",
    shape = "ADL decline severity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

fig_decline_severity_probs


fig_decline_severity_probs_facet <- ggplot(
  decline_severity_fig_df,
  aes(
    x = cannabis_group,
    y = predicted_probability
  )
) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_line(aes(group = 1), linewidth = 0.8) +
  facet_wrap(~ decline_severity, nrow = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.20),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Past-year cannabis exposure",
    y = "Adjusted predicted probability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

fig_decline_severity_probs_facet

# Cross-tab
adl_severity_cannabis_tab <- table(
  cannabis = adl_work$du_mar4_12m_aBin_ord,
  adl_severity = adl_work$adl_decline_severity_3cat_calc
)

adl_severity_cannabis_tab

# Chi-square test
chisq.test(adl_severity_cannabis_tab)

# If expected cells are small, use Fisher's exact test
fisher.test(adl_severity_cannabis_tab)

chisq.test(adl_severity_cannabis_tab)$expected

adl_work <- adl_work %>%
  dplyr::mutate(
    cannabis_ord_num = as.numeric(du_mar4_12m_aBin_ord), 
    adl_severity_num = as.numeric(adl_decline_severity_3cat_calc) - 1
  )

cor.test(
  adl_work$cannabis_ord_num,
  adl_work$adl_severity_num,
  method = "spearman"
)

# ============================================================
# Cramer's V without rcompanion
# ============================================================

chi_obj <- chisq.test(adl_severity_cannabis_tab)

cramers_v <- sqrt(
  as.numeric(chi_obj$statistic) /
    (
      sum(adl_severity_cannabis_tab) *
        (min(dim(adl_severity_cannabis_tab)) - 1)
    )
)

cramers_v 

chi_effect_summary <- tibble::tibble(
  test = "Cannabis group x ADL decline severity",
  chi_square = as.numeric(chi_obj$statistic),
  df = as.numeric(chi_obj$parameter),
  p_value = chi_obj$p.value,
  cramers_v = cramers_v
)

chi_effect_summary

round(
  prop.table(adl_severity_cannabis_tab, margin = 1) * 100,
  1
)



library(nnet)
library(broom)
library(dplyr)

# Make unordered version for multinomial model
adl_work <- adl_work %>%
  dplyr::mutate(
    adl_decline_severity_nominal = factor(
      as.character(adl_decline_severity_3cat_calc),
      levels = c("no_decline", "lower_decline", "greater_decline")
    )
  )

# Disease-first sequence: full model
m2_decline_multinom_full <- nnet::multinom(
  adl_decline_severity_nominal ~ DiseaseSev_c +
    du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  trace = FALSE
)

summary(m2_decline_multinom_full)

# Tidy exponentiated relative risk ratios
m2_decline_multinom_results <- broom::tidy(
  m2_decline_multinom_full,
  conf.int = TRUE,
  exponentiate = TRUE
)

m2_decline_multinom_results


