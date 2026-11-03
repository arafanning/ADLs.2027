# ==============================================================
# EVERYDAY FUNCTIONING + MEDICATION ADHERENCE AMONG PWH
# Recovered analysis specification — November 2026
# ==============================================================
#
# PROVENANCE:
# - NEW, reconstructed from the existing ADLs.2027 repository
#   and the final reported outcome/model descriptions.
# - NOT the author's original final script. NOT run on study data.
# - All original exploratory .R files remain unchanged.
# - Read RECOVERY_NOTES.md BEFORE manuscript use.
# - Never upload participant-level data or fitted objects.
# ==============================================================

# 0. Packages and restricted source data -----------------------
needed <- c("haven", "lavaan", "MASS", "lmtest", "sandwich", "emmeans")
missing_packages <- needed[!vapply(needed, requireNamespace,
                                   logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages manually: ",
       paste(missing_packages, collapse = ", "))
}
source_file <- "data/HIV+113.sav"
if (!file.exists(source_file)) {
  stop("Restricted study data not found: ", source_file)
}
adl_data <- as.data.frame(haven::read_sav(source_file))
if (nrow(adl_data) != 113L) {
  warning("N differs from documented PWH study sample of 113.")
}

required <- c(
  "lab_thel", "mhq_4_cd4_lowest", "lab_hiv",
  "du_mar4_12m_a", "phq_2_age", "phq_6_gender",
  "bdi_total", "adl_a_no", "adl_b_sum_now",
  "adl_c_sum_best", "mmt_ts", "mac_tmr4"
)
absent <- setdiff(required, names(adl_data))
if (length(absent)) {
  stop("Missing expected variables: ", paste(absent, collapse = ", "))
}

# 1. HIV disease severity --------------------------------------
# Original documented disease severity index.
# IMPORTANT: The constrained current-CD4 residual (0.001) was
# introduced to address a Heywood case and needs confirmation.
adl_data$cd4sqrt <- -sqrt(as.numeric(adl_data$lab_thel))
adl_data$nadirsqrt <- -sqrt(
  as.numeric(adl_data$mhq_4_cd4_lowest)
)
adl_data$log_vl2 <- log1p(as.numeric(adl_data$lab_hiv))
severity_spec <- '
  DiseaseSev =~ nadirsqrt + cd4sqrt + log_vl2
  cd4sqrt ~~ 0.001*cd4sqrt
'
severity_fit <- lavaan::cfa(
  severity_spec, data = adl_data, estimator = "MLR", missing = "fiml"
)
print(summary(severity_fit, fit.measures = TRUE, standardized = TRUE))
severity_scores <- as.numeric(lavaan::lavPredict(severity_fit))
if (length(severity_scores) != nrow(adl_data)) {
  stop("Severity factor score length does not match rows.")
}
adl_data$DiseaseSev_c <- as.numeric(scale(
  severity_scores, center = TRUE, scale = FALSE
))

# 2. Cannabis use and adjusted covariates ----------------------
grams <- as.numeric(adl_data$du_mar4_12m_a)
if (any(grams < 0, na.rm = TRUE)) {
  stop("Unexpected negative cannabis amount; verify source.")
}
median_grams_users <- median(grams[grams > 0], na.rm = TRUE)
adl_data$cu_group <- ordered(
  ifelse(is.na(grams), NA_character_,
    ifelse(grams == 0, "none",
      ifelse(grams <= median_grams_users, "low", "high"))),
  levels = c("none", "low", "high")
)
# The original ADL models used contr.poly(3) for group-level
# ordered linear / quadratic components; this is not a model of
# continuous grams of cannabis consumed.
contrasts(adl_data$cu_group) <- stats::contr.poly(3)
print(table(adl_data$cu_group, useNA = "ifany"))
if (!identical(as.integer(table(adl_data$cu_group)),
               c(62L, 27L, 24L))) {
  warning("Cannabis groups differ from 62/27/24; verify cases.")
}

adl_data$phq_2_age_c <- as.numeric(scale(
  as.numeric(adl_data$phq_2_age), center = TRUE, scale = FALSE
))
adl_data$bdi_total <- as.numeric(adl_data$bdi_total)
adl_data$sex_covfac <- factor(
  as.character(adl_data$phq_6_gender),
  levels = c("1", "2"),
  labels = c("Male", "Female")
)
print(table(adl_data$sex_covfac, useNA = "ifany"))
# Check this exact sex coding against the study codebook.
# BDI is modeled as a continuous symptom score.

# 3. Outcome definitions ---------------------------------------
# The original ADL scoring routine uses best-vs-current
# difference divided by the number of applicable domains/items.
denominator <- as.numeric(adl_data$adl_a_no)
now <- as.numeric(adl_data$adl_b_sum_now)
best <- as.numeric(adl_data$adl_c_sum_best)
adl_data$adl_total_calc <- ifelse(
  is.na(denominator) | denominator <= 0,
  NA_real_, (now - best) / denominator
)
if (any(adl_data$adl_total_calc < 0, na.rm = TRUE)) {
  warning("Negative ADL changes detected; inspect raw components.")
}
adl_data$adl_any_decline <- factor(
  ifelse(is.na(adl_data$adl_total_calc), NA_character_,
    ifelse(adl_data$adl_total_calc == 0, "no", "yes")),
  levels = c("no", "yes")
)

# Original source uses 1 as the floor of current difficulty:
# verify against the ADL scoring manual before reporting counts.
adl_data$adl_current_count <- now - 1
if (any(adl_data$adl_current_count < 0, na.rm = TRUE)) {
  stop("Current ADL burden contains negative counts.")
}
adl_data$adl_current_any <- factor(
  ifelse(is.na(adl_data$adl_current_count), NA_character_,
    ifelse(adl_data$adl_current_count == 0, "no", "yes")),
  levels = c("no", "yes")
)
adl_data$mmt_total <- as.numeric(adl_data$mmt_ts)
adl_data$morisky_ord <- ordered(as.numeric(adl_data$mac_tmr4))
# CRITICAL: verify that greater Morisky-4 values mean GREATER
# nonadherence in the particular dataset coding before
# interpreting a positive odds ratio.

# 4. Primary outcome-specific models --------------------------
# One fixed complete-case sample per outcome ensures meaningful
# comparisons across disease-only, cannabis-added and adjusted
# models. Samples can vary across outcomes.
covariates <- c("DiseaseSev_c", "cu_group",
                "phq_2_age_c", "bdi_total", "sex_covfac")

# Model sequence (from repository's final functional analyses):
#  1 = disease severity only
#  2 = disease severity + cannabis category (2 df)
#  3 = severity + cannabis + age, BDI, sex
#
# The final report found a cannabis relationship with the
# self-reported Morisky-4 measure, not the broader ADL or MMT-R
# outcomes. Numerical estimates are intentionally NOT hardcoded.

fit_three <- function(y, model_family) {
  needed_here <- c(y, covariates)
  d <- adl_data[stats::complete.cases(adl_data[, needed_here]),
                needed_here, drop = FALSE]
  d$cu_group <- factor(d$cu_group, ordered = TRUE,
                       levels = c("none", "low", "high"))
  contrasts(d$cu_group) <- stats::contr.poly(3)
  d$sex_covfac <- droplevels(d$sex_covfac)

  if (length(unique(d[[y]])) < 2L ||
      nlevels(droplevels(d$cu_group)) < 3L ||
      nlevels(d$sex_covfac) < 2L) {
    stop("Insufficient outcome or predictor variation for: ", y)
  }
  rhs <- c(
    "DiseaseSev_c",
    "DiseaseSev_c + cu_group",
    "DiseaseSev_c + cu_group + phq_2_age_c + bdi_total + sex_covfac"
  )
  specs <- lapply(rhs, function(x) stats::as.formula(paste(y, "~", x)))

  fitter <- switch(
    model_family,
    binary = function(f) stats::glm(f, data = d,
                                   family = stats::binomial()),
    count = function(f) MASS::glm.nb(f, data = d),
    linear = function(f) stats::lm(f, data = d),
    ordinal = function(f) MASS::polr(f, data = d, Hess = TRUE),
    stop("Unknown family")
  )
  models <- lapply(specs, fitter)
  names(models) <- c("disease", "cannabis", "adjusted")

  # Compare adding cannabis at step 2 using a likelihood-ratio
  # test of nested models on IDENTICAL cases. For the linear MMT
  # outcome, this is the Gaussian likelihood-ratio comparison;
  # HC3 coefficient inference is provided separately.
  ll0 <- stats::logLik(models$disease)
  ll1 <- stats::logLik(models$cannabis)
  ll2 <- stats::logLik(models$adjusted)
  lr_cannabis <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  df_added <- attr(ll1, "df") - attr(ll0, "df")
  lr_p <- stats::pchisq(lr_cannabis, df_added, lower.tail = FALSE)

  output <- list(
    n = nrow(d), family = model_family,
    distribution = table(d[[y]], useNA = "ifany"),
    models = models,
    cannabis_lr_step2 = c(
      chisq = lr_cannabis, df = df_added, p = lr_p
    ),
    model_aic = vapply(models, stats::AIC, numeric(1))
  )

  if (model_family == "linear") {
    output$hc3_adjusted <- lmtest::coeftest(
      models$adjusted,
      vcov. = sandwich::vcovHC(models$adjusted, type = "HC3")
    )
    output$adjusted_cannabis_means <- emmeans::emmeans(
      models$adjusted, ~ cu_group,
      vcov. = sandwich::vcovHC(models$adjusted, type = "HC3")
    )
  }

  if (model_family == "ordinal") {
    # MASS::polr: positive slope increases odds of MORE severe
    # outcomes. Interpret only after verifying Morisky coding.
    output$adjusted_or <- exp(stats::coef(models$adjusted))
    output$adjusted_ci_or <- exp(
      stats::confint.default(models$adjusted)[
        names(stats::coef(models$adjusted)), , drop = FALSE]
    )
  }
  output
}

outcome_families <- c(
  adl_any_decline = "binary",
  adl_current_any = "binary",
  adl_current_count = "count",
  mmt_total = "linear",
  morisky_ord = "ordinal"
)
recovered_outcomes <- lapply(
  names(outcome_families),
  function(y) fit_three(y, unname(outcome_families[[y]]))
)
names(recovered_outcomes) <- names(outcome_families)
for (y in names(recovered_outcomes)) {
  result <- recovered_outcomes[[y]]
  cat("\n", y, ": N=", result$n, "\n", sep = "")
  print(result$cannabis_lr_step2)
  print(result$model_aic)
  if (!is.null(result$hc3_adjusted)) print(result$hc3_adjusted)
  if (!is.null(result$adjusted_or)) print(result$adjusted_or)
}

# 5. ADL decline severity: a separate sensitivity analysis -----
# The older scripts considered several NON-EQUIVALENT ordinal
# variables. In particular, dividing a positive continuous
# decline score at its median is NOT the same as counting the
# number of declined ADL domains. Do not conflate the two.
positive_adl <- adl_data$adl_total_calc[
  is.finite(adl_data$adl_total_calc) &
    adl_data$adl_total_calc > 0
]
if (length(positive_adl)) {
  positive_cutpoint <- stats::median(positive_adl)
  adl_data$adl_severity_ratio <- ordered(
    ifelse(is.na(adl_data$adl_total_calc), NA_character_,
      ifelse(adl_data$adl_total_calc == 0, "none",
        ifelse(adl_data$adl_total_calc <= positive_cutpoint,
               "lower", "greater"))),
    levels = c("none", "lower", "greater")
  )
  # Sensitivity only; NOT the number-of-declined-domains outcome.
  sensitivity_severity <- fit_three("adl_severity_ratio", "ordinal")
  cat("\nADL ratio-based severity sensitivity\n")
  print(sensitivity_severity$cannabis_lr_step2)
}

# REVIEW BEFORE ANY PAPER OR ABSTRACT UPDATE:
# - Verify exact Morisky coding/direction and proportional-odds
#   assumptions.
# - Verify BDI and gender coding against the source codebook.
# - Verify core ADL scoring, denominator, and the -1 adjustment
#   for current burden.
# - Review separation / sparse outcome cells and model fit.
# - Distinguish the originally documented 0/1/2+ declined-domain
#   categories from the optional ratio-based median split.
# - Compare FINAL adjusted Morisky cannabis effect to the
#   published abstract; step 2 is NOT the adjusted estimate.
# - Confirm the covariate selection for each final model.
# - Avoid using unadjusted p values as if corrected for multiple
#   outcomes. This script has not been executed on data.
