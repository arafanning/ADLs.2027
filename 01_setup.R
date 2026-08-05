# Cannabis Use and Everyday Functioning Among People Living With HIV
# Script: 01_setup.R


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
# Check raw adl_total against recalculated ADL total
# ============================================================

library(dplyr)

# Create new diagnostic variables only; do not overwrite raw adl_total
adl_work <- adl_work %>%
  mutate(
    adl_total_recalc_check = (adl_b_sum_now - adl_c_sum_best) / adl_a_no,
    adl_total_diff_check = adl_total - adl_total_recalc_check
  )

# Compare distributions
summary(adl_work$adl_total)
summary(adl_work$adl_total_recalc_check)
summary(adl_work$adl_total_diff_check)

# Check ranges
range(adl_work$adl_total, na.rm = TRUE)
range(adl_work$adl_total_recalc_check, na.rm = TRUE)
range(adl_work$adl_total_diff_check, na.rm = TRUE)

# Count exact / near-exact mismatches
sum(abs(adl_work$adl_total_diff_check) > .0001, na.rm = TRUE)

# View cases where raw and recalculated values differ
adl_work %>%
  filter(abs(adl_total_diff_check) > .0001) %>%
  select(
    adl_a_no,
    adl_b_sum_now,
    adl_c_sum_best,
    adl_total,
    adl_total_recalc_check,
    adl_total_diff_check
  )

adl_work %>%
  dplyr::filter(abs(adl_total_diff_check) > .0001) %>%
  dplyr::select(
    adl_a_no,
    adl_b_sum_now,
    adl_c_sum_best,
    adl_total,
    adl_total_recalc_check,
    adl_total_diff_check
  )

# Does raw adl_total match recalculated score rounded to 2 decimals?
adl_work <- adl_work %>%
  dplyr::mutate(
    adl_total_recalc_round2_check = round(adl_total_recalc_check, 2),
    adl_total_rounding_diff_check = adl_total - adl_total_recalc_round2_check
  )

summary(adl_work$adl_total_rounding_diff_check)

sum(abs(adl_work$adl_total_rounding_diff_check) > .0001, na.rm = TRUE)

adl_work %>%
  dplyr::filter(abs(adl_total_rounding_diff_check) > .0001) %>%
  dplyr::select(
    adl_a_no,
    adl_b_sum_now,
    adl_c_sum_best,
    adl_total,
    adl_total_recalc_check,
    adl_total_recalc_round2_check,
    adl_total_rounding_diff_check
  )

#making sure nothing is negative
sum(adl_work$adl_total < 0, na.rm = TRUE)
sum(adl_work$adl_total_recalc_check < 0, na.rm = TRUE)
sum((adl_work$adl_b_sum_now - adl_work$adl_c_sum_best) < 0, na.rm = TRUE)


adl_work <- adl_work %>%
  dplyr::mutate(
    adl_total_calc = (adl_b_sum_now - adl_c_sum_best) / adl_a_no
  )

summary(adl_work$adl_total_calc)
sum(adl_work$adl_total_calc < 0, na.rm = TRUE)
sum(is.na(adl_work$adl_total_calc))

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_total_calc_any = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 ~ "any_decline"
    ),
    adl_total_calc_any = factor(
      adl_total_calc_any,
      levels = c("no_decline", "any_decline")
    ),
    
    adl_total_calc_positive = dplyr::case_when(
      adl_total_calc > 0 ~ adl_total_calc,
      TRUE ~ NA_real_
    )
  )

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_any_decline_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ "no_decline",
      adl_total_calc > 0 ~ "any_decline"
    ),
    adl_any_decline_calc = factor(
      adl_any_decline_calc,
      levels = c("no_decline", "any_decline")
    )
  )

table(adl_work$adl_any_decline_calc, useNA = "ifany")

positive_median_adl_calc <- median(
  adl_work$adl_total_calc[adl_work$adl_total_calc > 0],
  na.rm = TRUE
)

positive_median_adl_calc

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_decline_severity_binary_calc = dplyr::case_when(
      is.na(adl_total_calc) ~ NA_character_,
      adl_total_calc == 0 ~ NA_character_,  # excluded from severity model
      adl_total_calc > 0 & adl_total_calc <= positive_median_adl_calc ~ "lower_decline",
      adl_total_calc > positive_median_adl_calc ~ "higher_decline"
    ),
    adl_decline_severity_binary_calc = factor(
      adl_decline_severity_binary_calc,
      levels = c("lower_decline", "higher_decline")
    )
  )

table(adl_work$adl_decline_severity_binary_calc, useNA = "ifany")

model_adl_any_calc <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    nadirsqrt_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_adl_any_calc,
  conf.int = TRUE,
  exponentiate = TRUE
)

model_adl_severity_binary_calc <- glm(
  adl_decline_severity_binary_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    nadirsqrt_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_adl_severity_binary_calc,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# Create 3-level ADL decline severity category from adl_total_calc
# No decline / some decline / greater decline
# ============================================================

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
      adl_total_calc > 0 & adl_total_calc <= positive_median_adl_calc ~ "some_decline",
      adl_total_calc > positive_median_adl_calc ~ "greater_decline"
    ),
    adl_decline_severity_3cat_calc = factor(
      adl_decline_severity_3cat_calc,
      levels = c("no_decline", "some_decline", "greater_decline"),
      ordered = TRUE
    )
  )

table(adl_work$adl_decline_severity_3cat_calc, useNA = "ifany")

# ============================================================
# Describe adl_total_calc values within each severity category
# ============================================================

adl_severity_3cat_ranges <- adl_work %>%
  dplyr::filter(!is.na(adl_decline_severity_3cat_calc)) %>%
  dplyr::group_by(adl_decline_severity_3cat_calc) %>%
  dplyr::summarise(
    n = n(),
    min_adl_total_calc = min(adl_total_calc, na.rm = TRUE),
    max_adl_total_calc = max(adl_total_calc, na.rm = TRUE),
    mean_adl_total_calc = mean(adl_total_calc, na.rm = TRUE),
    sd_adl_total_calc = sd(adl_total_calc, na.rm = TRUE),
    median_adl_total_calc = median(adl_total_calc, na.rm = TRUE),
    .groups = "drop"
  )

adl_severity_3cat_ranges

# ============================================================
# Gap between category ranges
# ============================================================

adl_severity_3cat_gap_table <- adl_severity_3cat_ranges %>%
  dplyr::arrange(adl_decline_severity_3cat_calc) %>%
  dplyr::mutate(
    next_category = dplyr::lead(adl_decline_severity_3cat_calc),
    next_min = dplyr::lead(min_adl_total_calc),
    gap_to_next_category = next_min - max_adl_total_calc
  )

adl_severity_3cat_gap_table

# ============================================================
# Exact observed values by ADL decline category
# ============================================================

adl_value_frequency_by_3cat <- adl_work %>%
  dplyr::filter(!is.na(adl_decline_severity_3cat_calc)) %>%
  dplyr::count(
    adl_decline_severity_3cat_calc,
    adl_total_calc,
    name = "n_at_value"
  ) %>%
  dplyr::arrange(adl_total_calc)

adl_value_frequency_by_3cat

# ============================================================
# Inspect ADL B and C component variables
# ============================================================

library(dplyr)
library(psych)

adl_component_vars <- c(
  "adl_a_no",
  "adl_b_sum_now",
  "adl_c_sum_best",
  "adl_total_calc"
)

summary(adl_work[, adl_component_vars])

psych::describe(adl_work[, adl_component_vars])

# Frequency tables
table(adl_work$adl_a_no, useNA = "ifany")
table(adl_work$adl_b_sum_now, useNA = "ifany")
table(adl_work$adl_c_sum_best, useNA = "ifany")
table(adl_work$adl_total_calc, useNA = "ifany")

# ============================================================
# Create excess current and best ADL burden variables
# If 1 = no difficulty, then excess burden = score - 1
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_b_now_excess = adl_b_sum_now - 1,
    adl_c_best_excess = adl_c_sum_best - 1,
    adl_b_minus_c_raw = adl_b_sum_now - adl_c_sum_best,
    adl_b_minus_c_excess = adl_b_now_excess - adl_c_best_excess
  )

summary(adl_work[, c(
  "adl_b_sum_now",
  "adl_b_now_excess",
  "adl_c_sum_best",
  "adl_c_best_excess",
  "adl_b_minus_c_raw",
  "adl_b_minus_c_excess",
  "adl_total_calc"
)])

table(adl_work$adl_b_now_excess, useNA = "ifany")
table(adl_work$adl_c_best_excess, useNA = "ifany")
table(adl_work$adl_b_minus_c_raw, useNA = "ifany")

# ============================================================
# Are B and C integer-valued/count-like?
# ============================================================

all(adl_work$adl_b_sum_now %% 1 == 0, na.rm = TRUE)
all(adl_work$adl_c_sum_best %% 1 == 0, na.rm = TRUE)
all(adl_work$adl_b_now_excess %% 1 == 0, na.rm = TRUE)
all(adl_work$adl_c_best_excess %% 1 == 0, na.rm = TRUE)

# Mean/variance checks
component_dispersion <- tibble::tibble(
  variable = c(
    "adl_b_sum_now",
    "adl_b_now_excess",
    "adl_c_sum_best",
    "adl_c_best_excess",
    "adl_b_minus_c_raw"
  ),
  mean = c(
    mean(adl_work$adl_b_sum_now, na.rm = TRUE),
    mean(adl_work$adl_b_now_excess, na.rm = TRUE),
    mean(adl_work$adl_c_sum_best, na.rm = TRUE),
    mean(adl_work$adl_c_best_excess, na.rm = TRUE),
    mean(adl_work$adl_b_minus_c_raw, na.rm = TRUE)
  ),
  variance = c(
    var(adl_work$adl_b_sum_now, na.rm = TRUE),
    var(adl_work$adl_b_now_excess, na.rm = TRUE),
    var(adl_work$adl_c_sum_best, na.rm = TRUE),
    var(adl_work$adl_c_best_excess, na.rm = TRUE),
    var(adl_work$adl_b_minus_c_raw, na.rm = TRUE)
  )
) %>%
  dplyr::mutate(
    variance_to_mean = variance / mean
  )

component_dispersion

# ============================================================
# Binary versions of current burden, best burden, and decline
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_any_current_burden = dplyr::case_when(
      is.na(adl_b_now_excess) ~ NA_character_,
      adl_b_now_excess == 0 ~ "no_current_burden",
      adl_b_now_excess > 0 ~ "any_current_burden"
    ),
    adl_any_best_burden = dplyr::case_when(
      is.na(adl_c_best_excess) ~ NA_character_,
      adl_c_best_excess == 0 ~ "no_best_burden",
      adl_c_best_excess > 0 ~ "any_best_burden"
    ),
    adl_any_decline_from_components = dplyr::case_when(
      is.na(adl_b_minus_c_raw) ~ NA_character_,
      adl_b_minus_c_raw == 0 ~ "no_decline",
      adl_b_minus_c_raw > 0 ~ "any_decline"
    )
  )

table(adl_work$adl_any_current_burden, useNA = "ifany")
table(adl_work$adl_any_best_burden, useNA = "ifany")
table(adl_work$adl_any_decline_from_components, useNA = "ifany")

# Cross-tab: current burden vs decline
table(
  current_burden = adl_work$adl_any_current_burden,
  decline = adl_work$adl_any_decline_from_components,
  useNA = "ifany"
)

# ============================================================
# Compare B, C, B-C, and proportional decline
# ============================================================

cor(
  adl_work[, c(
    "adl_b_sum_now",
    "adl_c_sum_best",
    "adl_b_minus_c_raw",
    "adl_total_calc"
  )],
  use = "pairwise.complete.obs",
  method = "spearman"
)

# ============================================================
# Quick HC3 linear screening across component ADL outcomes
# ============================================================

library(lmtest)
library(sandwich)
library(broom)

candidate_adl_component_outcomes <- c(
  "adl_b_sum_now",
  "adl_b_now_excess",
  "adl_c_sum_best",
  "adl_c_best_excess",
  "adl_b_minus_c_raw",
  "adl_total_calc"
)

run_hc3_screen <- function(outcome) {
  
  f <- as.formula(
    paste0(
      outcome,
      " ~ du_mar4_12m_aBin_ord + phq_2_age_c + bdi_total + sex_covfac + nadirsqrt_c"
    )
  )
  
  model <- lm(f, data = adl_work)
  
  hc3 <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC3")
  )
  
  data.frame(
    outcome = outcome,
    term = rownames(hc3),
    estimate = hc3[, "Estimate"],
    robust_se = hc3[, "Std. Error"],
    statistic = hc3[, "t value"],
    p_value = hc3[, "Pr(>|t|)"],
    n = nobs(model),
    row.names = NULL
  )
}

adl_component_screen_results <- purrr::map_dfr(
  candidate_adl_component_outcomes,
  run_hc3_screen
)

adl_component_screen_results %>%
  dplyr::filter(grepl("du_mar4_12m_aBin_ord", term))

model_vars_check <- c(
  "adl_b_sum_now",
  "adl_c_sum_best",
  "adl_b_minus_c_raw",
  "adl_total_calc",
  "du_mar4_12m_aBin_ord",
  "phq_2_age_c",
  "bdi_total",
  "sex_covfac",
  "nadirsqrt_c"
)

sapply(adl_work[model_vars_check], function(x) sum(is.na(x)))

adl_work <- adl_work %>%
  dplyr::mutate(
    nadir_missing_check = dplyr::case_when(
      is.na(nadirsqrt_c) ~ "missing_nadir",
      !is.na(nadirsqrt_c) ~ "has_nadir"
    ),
    nadir_missing_check = factor(
      nadir_missing_check,
      levels = c("has_nadir", "missing_nadir")
    )
  )

table(adl_work$nadir_missing_check)

# Missingness by ADL decline
table(
  nadir_missing = adl_work$nadir_missing_check,
  any_decline = adl_work$adl_any_decline_calc,
  useNA = "ifany"
)

# Missingness by cannabis group
table(
  nadir_missing = adl_work$nadir_missing_check,
  cannabis = adl_work$du_mar4_12m_aBin_ord,
  useNA = "ifany"
)

# Compare continuous ADL decline by nadir missingness
adl_work %>%
  dplyr::group_by(nadir_missing_check) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_adl_total_calc = mean(adl_total_calc, na.rm = TRUE),
    sd_adl_total_calc = sd(adl_total_calc, na.rm = TRUE),
    median_adl_total_calc = median(adl_total_calc, na.rm = TRUE),
    any_decline_n = sum(adl_any_decline_calc == "any_decline", na.rm = TRUE),
    any_decline_prop = mean(adl_any_decline_calc == "any_decline", na.rm = TRUE),
    .groups = "drop"
  )

# Check missingness for individual HIV markers and latent severity
sapply(
  data[, c("cd4sqrt_c", "nadirsqrt_c", "log_vl2_c", "DiseaseSev_c")],
  function(x) sum(is.na(x))
)

# Check analytic N for core model vs disease-severity model
model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac,
  data = adl_work
) |> nrow()

model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac + DiseaseSev_c,
  data = adl_work
) |> nrow()

model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac + nadirsqrt_c,
  data = adl_work
) |> nrow()

sapply(
  adl_work[, names(adl_work) %in% c("cd4sqrt_c", "nadirsqrt_c", "log_vl2_c", "DiseaseSev_c")],
  function(x) sum(is.na(x))
)

adl_work <- adl_work %>%
  dplyr::mutate(
    cd4sqrt_c_check = as.numeric(scale(cd4sqrt, center = TRUE, scale = FALSE)),
    nadirsqrt_c_check = as.numeric(scale(nadirsqrt, center = TRUE, scale = FALSE)),
    log_vl2_c_check = as.numeric(scale(log_vl2, center = TRUE, scale = FALSE))
  )

sapply(
  adl_work[, c("cd4sqrt_c_check", "nadirsqrt_c_check", "log_vl2_c_check")],
  function(x) sum(is.na(x))
)

names(adl_work)[grepl("Disease|disease|sev|severity|factor", 
names(adl_work), ignore.case = TRUE)]

sapply(
  adl_work[, c("DiseaseSev", "DiseaseSev_c")],
  function(x) sum(is.na(x))
)

nrow(model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac,
  data = adl_work
))

nrow(model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac + DiseaseSev_c,
  data = adl_work
))

nrow(model.frame(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c + bdi_total + sex_covfac + nadirsqrt_c_check,
  data = adl_work
))

model_adl_any_calc_diseasesev <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_adl_any_calc_diseasesev,
  conf.int = TRUE,
  exponentiate = TRUE
)

model_adl_severity_binary_calc_diseasesev <- glm(
  adl_decline_severity_binary_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_adl_severity_binary_calc_diseasesev,
  conf.int = TRUE,
  exponentiate = TRUE
)

library(dplyr)
library(MASS)
library(broom)

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_b_now_excess = adl_b_sum_now - 1
  )

summary(adl_work$adl_b_now_excess)
table(adl_work$adl_b_now_excess, useNA = "ifany")

mean(adl_work$adl_b_now_excess, na.rm = TRUE)
var(adl_work$adl_b_now_excess, na.rm = TRUE)
var(adl_work$adl_b_now_excess, na.rm = TRUE) /
  mean(adl_work$adl_b_now_excess, na.rm = TRUE)

# Cannabis-only Poisson
model_b_pois_0 <- glm(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord,
  data = adl_work,
  family = poisson(link = "log")
)

# Core-adjusted Poisson
model_b_pois_core <- glm(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work,
  family = poisson(link = "log")
)

# Disease-severity-adjusted Poisson
model_b_pois_disease <- glm(
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

check_overdispersion(model_b_pois_0)
check_overdispersion(model_b_pois_core)
check_overdispersion(model_b_pois_disease)

# Cannabis-only negative binomial
model_b_nb_0 <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord,
  data = adl_work
)

# Core-adjusted negative binomial
model_b_nb_core <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac,
  data = adl_work
)

# Disease-severity-adjusted negative binomial
model_b_nb_disease <- MASS::glm.nb(
  adl_b_now_excess ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work
)

broom::tidy(model_b_nb_0, conf.int = TRUE, exponentiate = TRUE)

broom::tidy(model_b_nb_core, conf.int = TRUE, exponentiate = TRUE)

broom::tidy(model_b_nb_disease, conf.int = TRUE, exponentiate = TRUE)

AIC(model_b_pois_0, model_b_nb_0)
AIC(model_b_pois_core, model_b_nb_core)
AIC(model_b_pois_disease, model_b_nb_disease)


model_decline_nb_offset <- MASS::glm.nb(
  adl_b_minus_c_raw ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c +
    offset(log(adl_a_no)),
  data = adl_work
)

BIC(model_b_pois_0, model_b_nb_0)
BIC(model_b_pois_core, model_b_nb_core)
BIC(model_b_pois_disease, model_b_nb_disease)

summary(model_b_nb_disease)

library(emmeans)

emm_b_nb_disease <- emmeans::emmeans(
  model_b_nb_disease,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_b_nb_disease

pairs(emm_b_nb_disease)

model_decline_nb_offset <- MASS::glm.nb(
  adl_b_minus_c_raw ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c +
    offset(log(adl_a_no)),
  data = adl_work
)

library(MASS)
library(broom)
library(emmeans)

# Negative binomial model for ADL decline count,
# offset by number of applicable ADL items
model_decline_nb_offset <- MASS::glm.nb(
  adl_b_minus_c_raw ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c +
    offset(log(adl_a_no)),
  data = adl_work
)

summary(model_decline_nb_offset)

broom::tidy(
  model_decline_nb_offset,
  conf.int = TRUE,
  exponentiate = TRUE
)

emm_decline_nb_offset <- emmeans::emmeans(
  model_decline_nb_offset,
  ~ du_mar4_12m_aBin_ord,
  type = "response"
)

emm_decline_nb_offset

pairs(emm_decline_nb_offset)

model_decline_pois_offset <- glm(
  adl_b_minus_c_raw ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c +
    offset(log(adl_a_no)),
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

check_overdispersion(model_decline_pois_offset)

AIC(model_decline_pois_offset, model_decline_nb_offset)
BIC(model_decline_pois_offset, model_decline_nb_offset)


model_1_any_decline <- glm(
  adl_any_decline_calc ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    DiseaseSev_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_1_any_decline,
  conf.int = TRUE,
  exponentiate = TRUE
)

