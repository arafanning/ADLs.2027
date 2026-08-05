# Cannabis Use and Everyday Functioning Among People Living With HIV
# Script: 02_data_preparation.R
# Purpose: 02 Data Preparation


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


# ------------------------------------------------------------
# 5a. ADL distribution table
# ------------------------------------------------------------

## This table tells you whether each ADL candidate is usable as a continuous outcome.

adl_distribution_table <- data.frame(
  variable = adl_core_vars,
  n = sapply(adl_work[, adl_core_vars], function(x) sum(!is.na(x))),
  missing = sapply(adl_work[, adl_core_vars], function(x) sum(is.na(x))),
  mean = sapply(adl_work[, adl_core_vars], function(x) mean(x, na.rm = TRUE)),
  sd = sapply(adl_work[, adl_core_vars], function(x) sd(x, na.rm = TRUE)),
  median = sapply(adl_work[, adl_core_vars], function(x) median(x, na.rm = TRUE)),
  min = sapply(adl_work[, adl_core_vars], function(x) min(x, na.rm = TRUE)),
  max = sapply(adl_work[, adl_core_vars], function(x) max(x, na.rm = TRUE)),
  skew = sapply(adl_work[, adl_core_vars], function(x) psych::skew(x, na.rm = TRUE)),
  kurtosis = sapply(adl_work[, adl_core_vars], function(x) psych::kurtosi(x, na.rm = TRUE)),
  n_zero = sapply(adl_work[, adl_core_vars], function(x) sum(x == 0, na.rm = TRUE)),
  prop_zero = sapply(adl_work[, adl_core_vars], function(x) mean(x == 0, na.rm = TRUE)),
  row.names = NULL
)

adl_distribution_table


# ------------------------------------------------------------
# 5b. Histograms for candidate ADL outcomes
# ------------------------------------------------------------

hist(adl_work$adl_total, main = "ADL total", xlab = "adl_total")
hist(adl_work$adl_b_sum_now, main = "Current ADL difficulty", xlab = "adl_b_sum_now")
hist(adl_work$adl_a_no, main = "Number of scorable ADL items", xlab = "adl_a_no")


# ------------------------------------------------------------
# 5c. Create ADL outcome candidate versions
# ------------------------------------------------------------

## Binary any-decline outcome.
## Useful if adl_total is too zero-heavy.

adl_work$adl_total_any <- factor(
  dplyr::case_when(
    is.na(adl_work$adl_total) ~ NA_character_,
    adl_work$adl_total > 0 ~ "any_decline",
    adl_work$adl_total == 0 ~ "no_decline",
    TRUE ~ NA_character_
  ),
  levels = c("no_decline", "any_decline")
)

adl_work$adl_total_any_num <- as.numeric(adl_work$adl_total_any) - 1


## Positive-only ADL decline severity.
## Used only among participants with adl_total > 0.

adl_work$adl_total_positive <- ifelse(
  adl_work$adl_total > 0,
  adl_work$adl_total,
  NA_real_
)


## Log transforms for skewed nonnegative outcomes.

adl_work$adl_total_log1p <- log1p(adl_work$adl_total)
adl_work$adl_b_sum_now_log1p <- log1p(adl_work$adl_b_sum_now)
adl_work$mac_misr_log1p <- log1p(adl_work$mac_misr)


summary(adl_work[, c(
  "adl_total",
  "adl_total_any",
  "adl_total_positive",
  "adl_total_log1p",
  "adl_b_sum_now",
  "adl_b_sum_now_log1p",
  "mac_misr",
  "mac_misr_log1p"
)])

table(adl_work$adl_total_any, useNA = "ifany")


# ------------------------------------------------------------
# 5d. Define candidate outcomes
# ------------------------------------------------------------

## Do not commit to one primary outcome yet.
## These are candidates for decision-making.

candidate_primary_adl_outcomes <- c(
  "adl_total",
  "adl_total_log1p",
  "adl_b_sum_now",
  "adl_b_sum_now_log1p",
  "adl_total_any",
  "adl_total_positive"
)

## Continuous outcomes that can go into robust linear models.

adl_linear_candidate_outcomes <- intersect(
  c(
    "adl_total",
    "adl_total_log1p",
    "adl_b_sum_now",
    "adl_b_sum_now_log1p"
  ),
  names(adl_work)
)

## Secondary/convergent functional outcomes.

adl_secondary_outcomes <- intersect(
  c(
    "mmt_ts",
    "ft_total",
    "mac_misr",
    "mac_misr_log1p",
    "mac_tmr4"
  ),
  names(adl_work)
)

candidate_primary_adl_outcomes
adl_linear_candidate_outcomes
adl_secondary_outcomes

# ============================================================
# ADL TOTAL SKEW / ZERO-INFLATION MODEL EXPLORATION
# ============================================================

## Purpose:
## Explore defensible ways to model adl_total given:
## - adl_total is clinically/conceptually meaningful.
## - adl_total is nonnegative and bounded.
## - adl_total has a large zero/floor effect.
## - cannabis exposure will be tested using ordered none/low/high groups
##   with polynomial contrasts, allowing both linear and nonlinear/quadratic effects.
##
## Core cannabis predictor:
## du_mar4_12m_aBin_ord
##
## With contr.poly(3), R estimates:
## - du_mar4_12m_aBin_ord.L = linear cannabis trend
## - du_mar4_12m_aBin_ord.Q = quadratic/nonlinear cannabis trend


# ------------------------------------------------------------
# 1. Confirm cannabis polynomial coding
# ------------------------------------------------------------

contrasts(adl_work$du_mar4_12m_aBin_ord)

table(adl_work$du_mar4_12m_aBin_ord, useNA = "ifany")


# ------------------------------------------------------------
# 2. Create ADL total modeling variables
# ------------------------------------------------------------

## Do not overwrite adl_total.
## These are modeling-specific versions for exploration.

adl_work$adl_total_any_explore <- factor(
  dplyr::case_when(
    is.na(adl_work$adl_total) ~ NA_character_,
    adl_work$adl_total == 0 ~ "no_decline",
    adl_work$adl_total > 0 ~ "any_decline",
    TRUE ~ NA_character_
  ),
  levels = c("no_decline", "any_decline")
)

adl_work$adl_total_any_num_explore <- as.numeric(adl_work$adl_total_any_explore) - 1

adl_work$adl_total_positive_explore <- dplyr::case_when(
  adl_work$adl_total > 0 ~ adl_work$adl_total,
  TRUE ~ NA_real_
)

adl_work$adl_total_log1p_explore <- log1p(adl_work$adl_total)

adl_work$adl_total_positive_log1p_explore <- log1p(adl_work$adl_total_positive_explore)


## Optional categorical ADL decline variable.
## This is useful for descriptive checks or ordinal sensitivity models.

positive_median_adl <- median(
  adl_work$adl_total_positive_explore,
  na.rm = TRUE
)

adl_work$adl_total_ordinal_explore <- factor(
  dplyr::case_when(
    is.na(adl_work$adl_total) ~ NA_character_,
    adl_work$adl_total == 0 ~ "none",
    adl_work$adl_total > 0 & adl_work$adl_total <= positive_median_adl ~ "low_decline",
    adl_work$adl_total > positive_median_adl ~ "higher_decline"
  ),
  levels = c("none", "low_decline", "higher_decline"),
  ordered = TRUE
)

summary(adl_work[, c(
  "adl_total",
  "adl_total_any_explore",
  "adl_total_positive_explore",
  "adl_total_log1p_explore",
  "adl_total_positive_log1p_explore",
  "adl_total_ordinal_explore"
)])

table(adl_work$adl_total_any_explore, useNA = "ifany")
table(adl_work$adl_total_ordinal_explore, useNA = "ifany")

# ============================================================
# MODEL 1: Raw adl_total with HC3 robust standard errors
# ============================================================

install.packages("see")
library(see)

install.packages("ggplot2")
library(ggplot2)

library(performance)

## Interpretation:
## Tests whether cannabis linear and quadratic terms predict mean ADL decline.
## Limitation: The model treats many true zeros and positive values as one continuous process.

model_adl_raw <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    DiseaseSev_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_raw)

model_adl_raw_hc3 <- lmtest::coeftest(
  model_adl_raw,
  vcov. = sandwich::vcovHC(model_adl_raw, type = "HC3")
)

model_adl_raw_hc3

## Diagnostics that should not require the full check_model plotting stack

performance::check_heteroscedasticity(model_adl_raw)
performance::check_normality(model_adl_raw)
performance::check_collinearity(model_adl_raw)

## Base R diagnostic plots
par(mfrow = c(2, 2))
plot(model_adl_raw)
par(mfrow = c(1, 1))


# ============================================================
# Inspect influential cases from raw ADL model
# ============================================================

influence_table_adl_raw <- data.frame(
  case_id = as.numeric(rownames(model.frame(model_adl_raw))),
  fitted = fitted(model_adl_raw),
  residual = residuals(model_adl_raw),
  standardized_residual = rstandard(model_adl_raw),
  studentized_residual = rstudent(model_adl_raw),
  leverage = hatvalues(model_adl_raw),
  cooks_d = cooks.distance(model_adl_raw),
  adl_total = model.frame(model_adl_raw)$adl_total
) %>%
  arrange(desc(cooks_d))

head(influence_table_adl_raw, 15)

#looking at cases that are weird

influence_table_adl_raw %>%
  filter(case_id %in% c(20, 23, 58, 98))


# ============================================================
# PRIMARY CANDIDATE MODEL: TWO-PART ADL DECLINE APPROACH
# ============================================================

# ------------------------------------------------------------
# Part 1: Any ADL decline vs no decline
# ------------------------------------------------------------

model_adl_any <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    DiseaseSev_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

summary(model_adl_any)

model_adl_any_or <- broom::tidy(
  model_adl_any,
  conf.int = TRUE,
  exponentiate = TRUE
)

model_adl_any_or


# ------------------------------------------------------------
# Part 2: ADL decline severity among participants with any decline
# ------------------------------------------------------------

model_adl_positive <- lm(
  adl_total_positive_explore ~ du_mar4_12m_aBin_ord +
    DiseaseSev_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_positive)

model_adl_positive_hc3 <- lmtest::coeftest(
  model_adl_positive,
  vcov. = sandwich::vcovHC(model_adl_positive, type = "HC3")
)

model_adl_positive_hc3

###second part of the two-part model does not 
#support a cannabis effect on ADL decline severity among people who already show decline.
#model only includes the 50 participants with nonzero ADL decline. 
#The 63 no-decline participants were removed by design because 
#adl_total_positive_explore is NA for them.

#Past-year cannabis exposure did not significantly predict severity of ADL decline.
#Results -> du_mar4_12m_aBin_ord.L  p = .253
#du_mar4_12m_aBin_ord.Q  p = .176

# ============================================================
# Logistic models: Any ADL decline with separate HIV markers
# ============================================================

model_any_cd4 <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    cd4sqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

model_any_nadir <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

model_any_vl <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    log_vl2_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

model_any_diseasesev <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    DiseaseSev_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  family = binomial(link = "logit")
)

library(broom)
library(dplyr)

any_model_comparison <- bind_rows(
  tidy(model_any_cd4, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Current CD4"),
  
  tidy(model_any_nadir, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Nadir CD4"),
  
  tidy(model_any_vl, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Viral load"),
  
  tidy(model_any_diseasesev, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Latent disease severity")
) %>%
  filter(grepl("du_mar4_12m_aBin_ord|cd4sqrt_c|nadirsqrt_c|log_vl2_c|DiseaseSev_c", term)) %>%
  select(model, term, estimate, conf.low, conf.high, p.value)

any_model_comparison

##comparing fit between those two options 

AIC(model_any_cd4, model_any_nadir, model_any_vl, model_any_diseasesev)
BIC(model_any_cd4, model_any_nadir, model_any_vl, model_any_diseasesev)
##got an error here 
#BIC.default(model_any_cd4, model_any_nadir, model_any_vl, model_any_diseasesev) :
##models are not all fitted to the same number of observations

nobs(model_any_cd4)
nobs(model_any_nadir)
nobs(model_any_vl)
nobs(model_any_diseasesev)

# ============================================================
# Create common complete-case dataset for fair AIC/BIC comparison
# ============================================================

any_model_common_vars <- c(
  "adl_total_any_explore",
  "du_mar4_12m_aBin_ord",
  "cd4sqrt_c",
  "nadirsqrt_c",
  "log_vl2_c",
  "DiseaseSev_c",
  "phq_2_age_c",
  "phq_7_degree_c",
  "race_eth_binary_covfac"
)

adl_any_common <- adl_work %>%
  dplyr::select(all_of(any_model_common_vars)) %>%
  tidyr::drop_na()

nrow(adl_any_common)

table(adl_any_common$adl_total_any_explore)

# ============================================================
# Refit any-decline logistic models on same complete-case sample
# ============================================================

model_any_cd4_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    cd4sqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

model_any_nadir_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

model_any_vl_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    log_vl2_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

model_any_diseasesev_common <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    DiseaseSev_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

# ============================================================
# Fair AIC/BIC comparison
# ============================================================

AIC(
  model_any_cd4_common,
  model_any_nadir_common,
  model_any_vl_common,
  model_any_diseasesev_common
)

BIC(
  model_any_cd4_common,
  model_any_nadir_common,
  model_any_vl_common,
  model_any_diseasesev_common
)

nobs(model_any_cd4_common)
nobs(model_any_nadir_common)
nobs(model_any_vl_common)
nobs(model_any_diseasesev_common)

##Now the AIC/BIC comparison is valid because all models use the same people.

############
#extract the common-sample ORs
################

any_model_comparison_common <- bind_rows(
  tidy(model_any_cd4_common, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Current CD4"),
  
  tidy(model_any_nadir_common, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Nadir CD4"),
  
  tidy(model_any_vl_common, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Viral load"),
  
  tidy(model_any_diseasesev_common, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = "Latent disease severity")
) %>%
  filter(grepl("du_mar4_12m_aBin_ord|cd4sqrt_c|nadirsqrt_c|log_vl2_c|DiseaseSev_c", term)) %>%
  select(model, term, estimate, conf.low, conf.high, p.value)

any_model_comparison_common

# ============================================================
# Omnibus likelihood-ratio tests for cannabis effect
# ============================================================

drop1(model_any_cd4_common, test = "Chisq")
drop1(model_any_nadir_common, test = "Chisq")
drop1(model_any_vl_common, test = "Chisq")
drop1(model_any_diseasesev_common, test = "Chisq")

# ============================================================
# Any ADL decline by cannabis group, common sample
# ============================================================

table(
  adl_any_common$du_mar4_12m_aBin_ord,
  adl_any_common$adl_total_any_explore
)

prop.table(
  table(
    adl_any_common$du_mar4_12m_aBin_ord,
    adl_any_common$adl_total_any_explore
  ),
  margin = 1
)

# ============================================================
# Logistic model with cannabis as categorical predictor
# Reference group = none
# ============================================================

adl_any_common$du_mar4_12m_aBin_cat <- factor(
  adl_any_common$du_mar4_12m_aBin_ord,
  levels = c("none", "low", "high")
)

model_any_nadir_cat <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_cat +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

summary(model_any_nadir_cat)

broom::tidy(
  model_any_nadir_cat,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# Predicted probabilities by cannabis group
# ============================================================
library(emmeans)
emm_any_nadir_cat <- emmeans(
  model_any_nadir_cat,
  ~ du_mar4_12m_aBin_cat,
  type = "response"
)

emm_any_nadir_cat
pairs(emm_any_nadir_cat)

# ============================================================
# Check whether adl_b_sum_now is count-like
# ============================================================

summary(adl_work$adl_b_sum_now)
table(adl_work$adl_b_sum_now, useNA = "ifany")

# Is it integer-valued?
all(adl_work$adl_b_sum_now %% 1 == 0, na.rm = TRUE)

# Mean-variance check
mean_adl_b <- mean(adl_work$adl_b_sum_now, na.rm = TRUE)
var_adl_b  <- var(adl_work$adl_b_sum_now, na.rm = TRUE)

mean_adl_b
var_adl_b
var_adl_b / mean_adl_b

# Proportion at zero/minimum
mean(adl_work$adl_b_sum_now == 0, na.rm = TRUE)

# ============================================================
# Create excess ADL current difficulty score
# ============================================================

adl_work$adl_b_sum_now_excess <- adl_work$adl_b_sum_now - 1

summary(adl_work$adl_b_sum_now_excess)
table(adl_work$adl_b_sum_now_excess, useNA = "ifany")

mean(adl_work$adl_b_sum_now_excess, na.rm = TRUE)
var(adl_work$adl_b_sum_now_excess, na.rm = TRUE)
var(adl_work$adl_b_sum_now_excess, na.rm = TRUE) /
  mean(adl_work$adl_b_sum_now_excess, na.rm = TRUE)

#now modeling excess current adl burden 

# ============================================================
# Count-style models for excess current ADL burden
# ============================================================

library(MASS)
library(lmtest)
library(sandwich)
library(performance)

# Poisson
model_adl_b_excess_pois <- glm(
  adl_b_sum_now_excess ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work,
  family = poisson(link = "log")
)

summary(model_adl_b_excess_pois)
performance::check_overdispersion(model_adl_b_excess_pois)


# Quasi-Poisson
model_adl_b_excess_quasi <- glm(
  adl_b_sum_now_excess ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work,
  family = quasipoisson(link = "log")
)

summary(model_adl_b_excess_quasi)


# Negative binomial
model_adl_b_excess_nb <- MASS::glm.nb(
  adl_b_sum_now_excess ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_b_excess_nb)

AIC(model_adl_b_excess_pois, model_adl_b_excess_nb)
BIC(model_adl_b_excess_pois, model_adl_b_excess_nb)

# ============================================================
# MODEL CAPACITY CHECK: Logistic any-ADL-decline model
# ============================================================

event_table_any <- table(adl_any_common$adl_total_any_explore)
event_table_any

n_total_any <- nrow(adl_any_common)
n_events_any <- sum(adl_any_common$adl_total_any_explore == "any_decline")
n_none_any <- sum(adl_any_common$adl_total_any_explore == "no_decline")

n_total_any
n_events_any
n_none_any


# Function to count model parameters excluding intercept
count_model_terms <- function(model) {
  length(coef(model)) - 1
}

# Existing main-effects model
model_any_nadir_main <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

# Full L + Q interaction model
model_any_nadir_LQ_interaction <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord * nadirsqrt_c +
    phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

count_model_terms(model_any_nadir_main)
count_model_terms(model_any_nadir_LQ_interaction)

events_per_parameter_main <- n_events_any / count_model_terms(model_any_nadir_main)
events_per_parameter_interaction <- n_events_any / count_model_terms(model_any_nadir_LQ_interaction)

events_per_parameter_main
events_per_parameter_interaction


# ============================================================
# Create linear and quadratic cannabis contrast variables
# ============================================================

poly_contrasts_12m <- contrasts(adl_any_common$du_mar4_12m_aBin_ord)

adl_any_common$cu12m_L <- poly_contrasts_12m[
  as.character(adl_any_common$du_mar4_12m_aBin_ord),
  ".L"
]

adl_any_common$cu12m_Q <- poly_contrasts_12m[
  as.character(adl_any_common$du_mar4_12m_aBin_ord),
  ".Q"
]

table(adl_any_common$du_mar4_12m_aBin_ord, adl_any_common$cu12m_L)
table(adl_any_common$du_mar4_12m_aBin_ord, adl_any_common$cu12m_Q)

# ============================================================
# Linear cannabis main-effect model
# ============================================================

model_any_linear_main <- glm(
  adl_total_any_explore ~ cu12m_L +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

summary(model_any_linear_main)

broom::tidy(
  model_any_linear_main,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# Linear cannabis x nadir CD4 moderation model
# ============================================================

model_any_linear_interaction <- glm(
  adl_total_any_explore ~ cu12m_L * nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

summary(model_any_linear_interaction)

broom::tidy(
  model_any_linear_interaction,
  conf.int = TRUE,
  exponentiate = TRUE
)

# Likelihood-ratio test for adding interaction
anova(
  model_any_linear_main,
  model_any_linear_interaction,
  test = "Chisq"
)

#Checking for capacity even for that linear interaction 
count_model_terms(model_any_linear_main)
count_model_terms(model_any_linear_interaction)

n_events_any / count_model_terms(model_any_linear_main)
n_events_any / count_model_terms(model_any_linear_interaction)
#Ugh results are still cautious, but far more defensible than 4.875.

#just looking at quadratic-only interaction 

# ============================================================
# Quadratic cannabis x nadir CD4 exploratory model
# ============================================================

model_any_quadratic_main <- glm(
  adl_total_any_explore ~ cu12m_Q +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

model_any_quadratic_interaction <- glm(
  adl_total_any_explore ~ cu12m_Q * nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

summary(model_any_quadratic_interaction)

broom::tidy(
  model_any_quadratic_interaction,
  conf.int = TRUE,
  exponentiate = TRUE
)

anova(
  model_any_quadratic_main,
  model_any_quadratic_interaction,
  test = "Chisq"
)

#bc high cu is small.. 
# ============================================================
# Check fitted probabilities and possible instability
# ============================================================

range(fitted(model_any_linear_interaction))

summary(fitted(model_any_linear_interaction))

table(
  adl_any_common$du_mar4_12m_aBin_ord,
  adl_any_common$adl_total_any_explore
)

#figuring out other ways to deal with adl_total

model_adl_ordinal <- MASS::polr(
  adl_total_ordinal_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

summary(model_adl_ordinal)

broom::tidy(
  model_adl_ordinal,
  conf.int = TRUE,
  exponentiate = TRUE
)

#robust regression on raw adl_total reduces sensitivity to outliers/high-leverage observations
library(MASS)

model_adl_total_rlm <- MASS::rlm(
  adl_total ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_total_rlm)

#positive part of adl_total, a Gamma GLM may be more appropriate than linear regression because positive ADL 
#decline is continuous, nonnegative, and right-skewed.

model_adl_positive_gamma <- glm(
  adl_total_positive_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work,
  family = Gamma(link = "log")
)

summary(model_adl_positive_gamma)

broom::tidy(
  model_adl_positive_gamma,
  conf.int = TRUE,
  exponentiate = TRUE
)

###behaves like an overdispersed count-style burden score. 
#Negative binomial is the best model among the count-style options you tested.

model_adl_b_excess_nb <- MASS::glm.nb(
  adl_b_sum_now_excess ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c + phq_2_age_c + phq_7_degree_c + race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_b_excess_nb)

# ============================================================
# EXPLORATORY LATENT FUNCTIONAL/ADHERENCE BURDEN VARIABLE
# ============================================================

## Goal:
## Create directionally consistent indicators where higher = worse
## functional/adherence burden.

# ============================================================
# Simple standardized composite
# Higher = worse functional/adherence burden
# ============================================================

adl_work$functional_adherence_burden_z <- rowMeans(
  adl_work[, latent_indicator_vars],
  na.rm = TRUE
)

# Require at least 2 of 3 indicators present
adl_work$n_functional_adherence_indicators <- rowSums(
  !is.na(adl_work[, latent_indicator_vars])
)

adl_work$functional_adherence_burden_z[
  adl_work$n_functional_adherence_indicators < 2
] <- NA

summary(adl_work$functional_adherence_burden_z)

hist(
  adl_work$functional_adherence_burden_z,
  main = "Functional/Adherence Burden Composite",
  xlab = "Composite z score; higher = worse"
)

##now latent 
latent_vars_raw <- c("adl_b_sum_now", "mmt_ts", "mac_misr")

summary(adl_work[, latent_vars_raw])
cor(adl_work[, latent_vars_raw], use = "pairwise.complete.obs", method = "spearman")


# Reverse MMT so higher = worse medication-management performance
adl_work$mmt_ts_rev <- -1 * adl_work$mmt_ts

# Transform skewed variables where appropriate
adl_work$adl_b_sum_now_log1p_latent <- log1p(adl_work$adl_b_sum_now)
adl_work$mac_misr_log1p_latent <- log1p(adl_work$mac_misr)

# Standardize indicators for exploratory composite/CFA
adl_work$adl_b_z_latent <- as.numeric(scale(adl_work$adl_b_sum_now_log1p_latent))
adl_work$mmt_rev_z_latent <- as.numeric(scale(adl_work$mmt_ts_rev))
adl_work$mac_misr_z_latent <- as.numeric(scale(adl_work$mac_misr_log1p_latent))

latent_indicator_vars <- c(
  "adl_b_z_latent",
  "mmt_rev_z_latent",
  "mac_misr_z_latent"
)

summary(adl_work[, latent_indicator_vars])
cor(adl_work[, latent_indicator_vars], use = "pairwise.complete.obs", method = "spearman")

#testing cannabis model 

model_functional_adherence_composite <- lm(
  functional_adherence_burden_z ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_functional_adherence_composite)

lmtest::coeftest(
  model_functional_adherence_composite,
  vcov. = sandwich::vcovHC(model_functional_adherence_composite, type = "HC3")
)

library(lavaan)

# ============================================================
# Exploratory one-factor CFA
# Higher latent factor = worse functional/adherence burden
# ============================================================

functional_burden_model <- '
  FunctionalBurden =~ adl_b_z_latent + mmt_rev_z_latent + mac_misr_z_latent
'

fit_functional_burden <- lavaan::cfa(
  functional_burden_model,
  data = adl_work,
  estimator = "MLR",
  missing = "fiml"
)

summary(
  fit_functional_burden,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)

# Extract factor scores
adl_work$FunctionalBurden_factor <- as.numeric(
  lavaan::lavPredict(fit_functional_burden)
)

summary(adl_work$FunctionalBurden_factor)
hist(adl_work$FunctionalBurden_factor)

model_functional_burden_factor <- lm(
  FunctionalBurden_factor ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_functional_burden_factor)

lmtest::coeftest(
  model_functional_burden_factor,
  vcov. = sandwich::vcovHC(model_functional_burden_factor, type = "HC3")
)

###above was terrible fit, ignore it.

# ============================================================
# Change-in-estimate check for cannabis terms
# ============================================================

# Unadjusted cannabis-only model
model_any_unadjusted <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord,
  data = adl_any_common,
  family = binomial(link = "logit")
)

# Adjusted model
model_any_adjusted <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

broom::tidy(model_any_unadjusted, conf.int = TRUE, exponentiate = TRUE)

broom::tidy(model_any_adjusted, conf.int = TRUE, exponentiate = TRUE)

# ============================================================
# Minimal vs full adjustment comparison
# ============================================================

model_any_minimal <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c,
  data = adl_any_common,
  family = binomial(link = "logit")
)

model_any_full <- glm(
  adl_total_any_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

minimal_vs_full <- dplyr::bind_rows(
  broom::tidy(model_any_minimal, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(model = "Minimal: cannabis + nadir CD4"),
  
  broom::tidy(model_any_full, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(model = "Full: + age, education, race/ethnicity")
) |>
  dplyr::filter(grepl("du_mar4_12m_aBin_ord|nadirsqrt_c", term)) |>
  dplyr::select(model, term, estimate, conf.low, conf.high, p.value)

minimal_vs_full

# ============================================================
# POWER / MODEL CAPACITY CHECKS FOR ADL PAPER
# ============================================================

## Purpose:
## These checks help evaluate whether different ADL models are
## reasonably powered / stable enough for primary or sensitivity analyses.
##
## Important:
## These are not meant to "prove" null findings. They are meant to help
## describe model capacity, precision, and whether interaction models are
## overbuilt for the available sample/events.

# ============================================================
# Create clean ordered cannabis factor for power/sensitivity models
# Does not overwrite original cannabis variable
# ============================================================

adl_work$cu12m_ord_for_power <- factor(
  adl_work$du_mar4_12m_aBin_ord,
  levels = c("none", "low", "high"),
  ordered = TRUE
)

table(adl_work$du_mar4_12m_aBin_ord, adl_work$cu12m_ord_for_power, 
      useNA = "ifany")

# ============================================================
# Create polynomial contrast variables using numeric factor codes
# ============================================================

poly_contrasts_12m_work <- contrasts(adl_work$cu12m_ord_for_power)

poly_contrasts_12m_work

adl_work$cu12m_L <- poly_contrasts_12m_work[
  as.numeric(adl_work$cu12m_ord_for_power),
  ".L"
]

adl_work$cu12m_Q <- poly_contrasts_12m_work[
  as.numeric(adl_work$cu12m_ord_for_power),
  ".Q"
]

# Check coding
table(adl_work$cu12m_ord_for_power, adl_work$cu12m_L, useNA = "ifany")
table(adl_work$cu12m_ord_for_power, adl_work$cu12m_Q, useNA = "ifany")

# ============================================================
# Negative binomial model for current ADL burden
# ============================================================

model_nb_main_power <- MASS::glm.nb(
  adl_b_sum_now_excess ~ cu12m_L +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_nb_main_power)

# ============================================================
# 2A. Numeric binary outcome for simulation
# ============================================================

adl_any_common$adl_any_num_power <- ifelse(
  adl_any_common$adl_total_any_explore == "any_decline", 1, 0
)

table(adl_any_common$adl_total_any_explore, adl_any_common$adl_any_num_power)

# ============================================================
# 2B. SIMULATION POWER: Logistic linear cannabis main effect
# ============================================================

## Base model without cannabis
model_any_base_no_cu <- glm(
  adl_any_num_power ~ nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

simulate_logistic_main_power <- function(
    data,
    cannabis_beta,
    nsim = 1000,
    alpha = 0.05
) {
  
  p_values <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    
    sim_data <- data
    
    ## Start from base predicted log-odds
    lp_base <- predict(model_any_base_no_cu, newdata = sim_data, type = "link")
    
    ## Add hypothetical cannabis linear effect
    lp_sim <- lp_base + cannabis_beta * sim_data$cu12m_L
    
    ## Convert to probability
    prob_sim <- plogis(lp_sim)
    
    ## Simulate outcome
    sim_data$y_sim <- rbinom(
      n = nrow(sim_data),
      size = 1,
      prob = prob_sim
    )
    
    ## Fit model with cannabis
    fit_full <- glm(
      y_sim ~ cu12m_L +
        nadirsqrt_c +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    ## Extract p-value for cannabis linear term
    p_values[i] <- summary(fit_full)$coefficients["cu12m_L", "Pr(>|z|)"]
  }
  
  tibble::tibble(
    cannabis_beta = cannabis_beta,
    cannabis_OR = exp(cannabis_beta),
    nsim = nsim,
    alpha = alpha,
    power = mean(p_values < alpha, na.rm = TRUE)
  )
}

## Test plausible cannabis ORs
cannabis_effects_to_test <- log(c(1.25, 1.50, 1.75, 2.00, 2.50, 3.00))

power_logistic_main_current_N <- purrr::map_dfr(
  cannabis_effects_to_test,
  ~ simulate_logistic_main_power(
    data = adl_any_common,
    cannabis_beta = .x,
    nsim = 1000,
    alpha = 0.05
  )
)

power_logistic_main_current_N

# ============================================================
# 2C. SIMULATION POWER: Logistic linear cannabis x nadir CD4
# ============================================================

## Main-effects model used as baseline
model_any_linear_main_power <- glm(
  adl_any_num_power ~ cu12m_L +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

simulate_logistic_interaction_power <- function(
    data,
    interaction_beta,
    nsim = 1000,
    alpha = 0.05
) {
  
  p_values <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    
    sim_data <- data
    
    ## Baseline log-odds from observed main-effects model
    lp_base <- predict(model_any_linear_main_power, newdata = sim_data, type = "link")
    
    ## Add hypothetical interaction effect
    lp_sim <- lp_base +
      interaction_beta * sim_data$cu12m_L * sim_data$nadirsqrt_c
    
    prob_sim <- plogis(lp_sim)
    
    sim_data$y_sim <- rbinom(
      n = nrow(sim_data),
      size = 1,
      prob = prob_sim
    )
    
    ## Reduced model
    fit_reduced <- glm(
      y_sim ~ cu12m_L +
        nadirsqrt_c +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    ## Full interaction model
    fit_full <- glm(
      y_sim ~ cu12m_L * nadirsqrt_c +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    lrt <- anova(fit_reduced, fit_full, test = "Chisq")
    p_values[i] <- lrt$`Pr(>Chi)`[2]
  }
  
  tibble::tibble(
    interaction_beta = interaction_beta,
    interaction_OR = exp(interaction_beta),
    nsim = nsim,
    alpha = alpha,
    power = mean(p_values < alpha, na.rm = TRUE)
  )
}

interaction_effects_to_test <- log(c(1.25, 1.50, 1.75, 2.00, 2.50, 3.00))

power_logistic_interaction_current_N <- purrr::map_dfr(
  interaction_effects_to_test,
  ~ simulate_logistic_interaction_power(
    data = adl_any_common,
    interaction_beta = .x,
    nsim = 1000,
    alpha = 0.05
  )
)

power_logistic_interaction_current_N

# ============================================================
# 3. SAMPLE SIZE POWER CURVE: Logistic interaction model
# ============================================================

simulate_logistic_interaction_power_by_N <- function(
    data,
    target_N,
    interaction_beta,
    nsim = 1000,
    alpha = 0.05
) {
  
  p_values <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    
    ## Resample observed data structure up to target N
    sim_data <- data[sample(seq_len(nrow(data)), size = target_N, replace = TRUE), ]
    
    lp_base <- predict(model_any_linear_main_power, newdata = sim_data, type = "link")
    
    lp_sim <- lp_base +
      interaction_beta * sim_data$cu12m_L * sim_data$nadirsqrt_c
    
    prob_sim <- plogis(lp_sim)
    
    sim_data$y_sim <- rbinom(
      n = nrow(sim_data),
      size = 1,
      prob = prob_sim
    )
    
    fit_reduced <- glm(
      y_sim ~ cu12m_L +
        nadirsqrt_c +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    fit_full <- glm(
      y_sim ~ cu12m_L * nadirsqrt_c +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    lrt <- anova(fit_reduced, fit_full, test = "Chisq")
    p_values[i] <- lrt$`Pr(>Chi)`[2]
  }
  
  tibble::tibble(
    target_N = target_N,
    interaction_beta = interaction_beta,
    interaction_OR = exp(interaction_beta),
    nsim = nsim,
    alpha = alpha,
    power = mean(p_values < alpha, na.rm = TRUE)
  )
}

## Example: power curve for OR = 2.00 interaction
sample_sizes_to_test <- c(90, 120, 150, 180, 220, 260, 300, 350, 400)

power_curve_logistic_interaction_OR2 <- purrr::map_dfr(
  sample_sizes_to_test,
  ~ simulate_logistic_interaction_power_by_N(
    data = adl_any_common,
    target_N = .x,
    interaction_beta = log(2.00),
    nsim = 1000,
    alpha = 0.05
  )
)

power_curve_logistic_interaction_OR2

ggplot(
  power_curve_logistic_interaction_OR2,
  aes(x = target_N, y = power)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = .80, linetype = "dashed") +
  labs(
    title = "Estimated Power for Linear Cannabis x Nadir CD4 Interaction",
    subtitle = "Logistic model predicting any ADL decline; hypothetical interaction OR = 2.00",
    x = "Target sample size",
    y = "Estimated power"
  )

# ============================================================
# 4. LINEAR MODEL POWER APPROXIMATION
# Outcome: raw adl_total
# ============================================================

## Install if needed
# install.packages("pwr")

library(pwr)

## Fit full linear model
model_lm_adl_total_power <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_lm_adl_total_power)

## Overall model R2
r2_lm <- summary(model_lm_adl_total_power)$r.squared

## Convert R2 to Cohen's f2
f2_lm <- r2_lm / (1 - r2_lm)

r2_lm
f2_lm

## Number of predictors excluding intercept
u_lm <- length(coef(model_lm_adl_total_power)) - 1

## Denominator df
v_lm <- df.residual(model_lm_adl_total_power)

## Approximate observed power for full model
pwr::pwr.f2.test(
  u = u_lm,
  v = v_lm,
  f2 = f2_lm,
  sig.level = 0.05
)

# ============================================================
# Incremental R2 power: adding cannabis L + Q block
# ============================================================

model_lm_no_cu <- lm(
  adl_total ~ nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

model_lm_with_cu <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

r2_no_cu <- summary(model_lm_no_cu)$r.squared
r2_with_cu <- summary(model_lm_with_cu)$r.squared

delta_r2_cu <- r2_with_cu - r2_no_cu

## Cohen's f2 for added block
f2_delta_cu <- delta_r2_cu / (1 - r2_with_cu)

r2_no_cu
r2_with_cu
delta_r2_cu
f2_delta_cu

## u = number of added predictors in cannabis block
## For linear + quadratic cannabis, u = 2
u_delta <- 2
v_delta <- df.residual(model_lm_with_cu)

pwr::pwr.f2.test(
  u = u_delta,
  v = v_delta,
  f2 = f2_delta_cu,
  sig.level = 0.05
)

# ============================================================
# 5. NEGATIVE BINOMIAL POWER: Current ADL burden
# Outcome: adl_b_sum_now_excess
# ============================================================

model_nb_main_power <- MASS::glm.nb(
  adl_b_sum_now_excess ~ cu12m_L +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_nb_main_power)

theta_nb_observed <- model_nb_main_power$theta
theta_nb_observed

simulate_nb_main_power <- function(
    data,
    cannabis_beta,
    theta,
    nsim = 1000,
    alpha = 0.05
) {
  
  p_values <- numeric(nsim)
  
  ## Base model without cannabis
  model_nb_base_no_cu <- MASS::glm.nb(
    adl_b_sum_now_excess ~ nadirsqrt_c +
      phq_2_age_c +
      phq_7_degree_c +
      race_eth_binary_covfac,
    data = data
  )
  
  for (i in seq_len(nsim)) {
    
    sim_data <- data
    
    lp_base <- predict(model_nb_base_no_cu, newdata = sim_data, type = "link")
    
    lp_sim <- lp_base + cannabis_beta * sim_data$cu12m_L
    
    mu_sim <- exp(lp_sim)
    
    sim_data$y_sim <- MASS::rnegbin(
      n = nrow(sim_data),
      mu = mu_sim,
      theta = theta
    )
    
    fit_full <- try(
      MASS::glm.nb(
        y_sim ~ cu12m_L +
          nadirsqrt_c +
          phq_2_age_c +
          phq_7_degree_c +
          race_eth_binary_covfac,
        data = sim_data
      ),
      silent = TRUE
    )
    
    if (inherits(fit_full, "try-error")) {
      p_values[i] <- NA
    } else {
      p_values[i] <- summary(fit_full)$coefficients["cu12m_L", "Pr(>|z|)"]
    }
  }
  
  tibble::tibble(
    cannabis_beta = cannabis_beta,
    cannabis_rate_ratio = exp(cannabis_beta),
    nsim = nsim,
    alpha = alpha,
    power = mean(p_values < alpha, na.rm = TRUE),
    failed_models = sum(is.na(p_values))
  )
}

nb_cannabis_effects_to_test <- log(c(1.25, 1.50, 1.75, 2.00, 2.50, 3.00))

power_nb_main_current_N <- purrr::map_dfr(
  nb_cannabis_effects_to_test,
  ~ simulate_nb_main_power(
    data = adl_work,
    cannabis_beta = .x,
    theta = theta_nb_observed,
    nsim = 1000,
    alpha = 0.05
  )
)

power_nb_main_current_N


# ============================================================
# 6. MODEL POWER / CAPACITY SUMMARY TABLE
# ============================================================

adl_power_summary_table <- tibble::tribble(
  ~model_family, ~outcome, ~primary_question, ~n_used, ~events_or_positive_cases, ~capacity_metric, ~interpretation,
  
  "Logistic",
  "Any ADL decline",
  "Does cannabis predict any decline?",
  nobs(model_any_capacity_full),
  n_events_any,
  paste0("EPV = ", round(
    n_events_any / count_model_terms(model_any_capacity_full), 2
  )),
  "Cautious but usable for main effects",
  
  "Logistic",
  "Any ADL decline",
  "Does cannabis x nadir CD4 predict decline?",
  nobs(model_any_capacity_linear_int),
  n_events_any,
  paste0("EPV = ", round(
    n_events_any / count_model_terms(model_any_capacity_linear_int), 2
  )),
  "Linear-only interaction is cautious; full L+Q interaction likely overbuilt",
  
  "Ordinal logistic",
  "None / low / higher ADL decline",
  "Does cannabis predict worse ordered decline category?",
  nrow(adl_work[!is.na(adl_work$adl_total_ordinal_explore), ]),
  NA,
  "Check category counts and proportional odds",
  "Likely best primary candidate if assumptions acceptable",
  
  "Negative binomial",
  "Current ADL burden",
  "Does cannabis predict count-style current ADL burden?",
  nobs(model_nb_main_power),
  NA,
  paste0("Theta = ", round(theta_nb_observed, 2)),
  "Appropriate sensitivity model if outcome is overdispersed"
)

adl_power_summary_table

###power for the interaction is a problem 
#That is a red flag. It means your simulated interaction is producing extreme 
#predicted probabilities, essentially pushing people close to guaranteed no 
#decline or guaranteed decline. Warning: fitted probabilities 
#numerically 0 or 1 occurred

#trying to fix using standardized nadir cd4
# Standardize nadir CD4 for interaction power checks
# ============================================================

adl_any_common$nadirsqrt_z <- as.numeric(scale(adl_any_common$nadirsqrt_c))

#now replacing interaction simulation model

# ============================================================
# Corrected interaction power simulation using standardized nadir
# ============================================================

model_any_linear_main_power_z <- glm(
  adl_any_num_power ~ cu12m_L +
    nadirsqrt_z +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_any_common,
  family = binomial(link = "logit")
)

simulate_logistic_interaction_power_z <- function(
    data,
    interaction_beta,
    nsim = 1000,
    alpha = 0.05
) {
  
  p_values <- numeric(nsim)
  
  for (i in seq_len(nsim)) {
    
    sim_data <- data
    
    lp_base <- predict(
      model_any_linear_main_power_z,
      newdata = sim_data,
      type = "link"
    )
    
    lp_sim <- lp_base +
      interaction_beta * sim_data$cu12m_L * sim_data$nadirsqrt_z
    
    prob_sim <- plogis(lp_sim)
    
    sim_data$y_sim <- rbinom(
      n = nrow(sim_data),
      size = 1,
      prob = prob_sim
    )
    
    fit_reduced <- glm(
      y_sim ~ cu12m_L +
        nadirsqrt_z +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    fit_full <- glm(
      y_sim ~ cu12m_L * nadirsqrt_z +
        phq_2_age_c +
        phq_7_degree_c +
        race_eth_binary_covfac,
      data = sim_data,
      family = binomial(link = "logit")
    )
    
    lrt <- anova(fit_reduced, fit_full, test = "Chisq")
    p_values[i] <- lrt$`Pr(>Chi)`[2]
  }
  
  tibble::tibble(
    interaction_beta = interaction_beta,
    interaction_OR = exp(interaction_beta),
    nsim = nsim,
    alpha = alpha,
    power = mean(p_values < alpha, na.rm = TRUE),
    failed_or_missing = sum(is.na(p_values))
  )
}

interaction_effects_to_test <- log(c(1.25, 1.50, 1.75, 2.00, 2.50, 3.00))

power_logistic_interaction_current_N_z <- purrr::map_dfr(
  interaction_effects_to_test,
  ~ simulate_logistic_interaction_power_z(
    data = adl_any_common,
    interaction_beta = .x,
    nsim = 1000,
    alpha = 0.05
  )
)

power_logistic_interaction_current_N_z

# ============================================================
# ORDINAL ADL MODEL CAPACITY CHECK
# ============================================================

table(adl_work$adl_total_ordinal_explore, useNA = "ifany")

ordinal_n <- sum(!is.na(adl_work$adl_total_ordinal_explore))
ordinal_category_counts <- table(adl_work$adl_total_ordinal_explore)

ordinal_n
ordinal_category_counts

# Fit ordinal main-effects model
library(MASS)

model_ord_main <- MASS::polr(
  adl_total_ordinal_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

summary(model_ord_main)

# Count predictors excluding thresholds/intercepts
ordinal_predictors <- length(coef(model_ord_main))

# Smallest category count per predictor
ordinal_capacity_table <- tibble::tibble(
  model = "Ordinal ADL decline main effects",
  n_total = ordinal_n,
  smallest_category_n = min(ordinal_category_counts),
  predictors_excluding_thresholds = ordinal_predictors,
  smallest_category_per_predictor = min(ordinal_category_counts) / ordinal_predictors
)

ordinal_capacity_table

# ============================================================
# PROPORTIONAL ODDS ASSUMPTION CHECK
# ============================================================

install.packages("brant")
library(brant)

brant::brant(model_ord_main)


# ============================================================
# FIGURES FOR ADL MODEL DECISION-MAKING
# ============================================================

library(dplyr)
library(ggplot2)
library(emmeans)
library(broom)
library(tidyr)
library(forcats)

# Optional theme for cleaner figures
theme_set(theme_classic(base_size = 13))

# Make clean cannabis labels if needed
adl_work <- adl_work %>%
  mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

adl_any_common <- adl_any_common %>%
  mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )
# ============================================================
# FIGURE 1: Distribution of raw ADL decline score
# ============================================================

fig_adl_total_distribution <- ggplot(adl_work, aes(x = adl_total)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Raw ADL Decline Score Is Zero-Heavy and Skewed",
    subtitle = "Large point mass at 0 indicates no reported ADL decline",
    x = "ADL decline score",
    y = "Number of participants"
  )

fig_adl_total_distribution

#by cannabis group 

fig_adl_total_by_cannabis <- ggplot(
  adl_work,
  aes(x = adl_total)
) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "black") +
  facet_wrap(~ cannabis_group) +
  labs(
    title = "Raw ADL Decline Distribution by Cannabis Exposure Group",
    x = "ADL decline score",
    y = "Number of participants"
  )

fig_adl_total_by_cannabis


# ============================================================
# FIGURE 2: Ordinal ADL decline categories
# ============================================================

fig_ordinal_counts <- adl_work %>%
  filter(!is.na(adl_total_ordinal_explore)) %>%
  mutate(
    adl_decline_category = factor(
      adl_total_ordinal_explore,
      levels = c("none", "low_decline", "higher_decline"),
      labels = c("No decline", "Low decline", "Higher decline"),
      ordered = TRUE
    )
  ) %>%
  ggplot(aes(x = adl_decline_category)) +
  geom_bar(color = "black") +
  labs(
    title = "Ordinal ADL Decline Categories Are Usable",
    subtitle = "Categories preserve severity information beyond no/any decline",
    x = "ADL decline category",
    y = "Number of participants"
  )

fig_ordinal_counts

# ============================================================
# FIGURE 3: Ordinal ADL decline by cannabis group
# ============================================================

fig_ordinal_by_cannabis <- adl_work %>%
  filter(
    !is.na(adl_total_ordinal_explore),
    !is.na(cannabis_group)
  ) %>%
  mutate(
    adl_decline_category = factor(
      adl_total_ordinal_explore,
      levels = c("none", "low_decline", "higher_decline"),
      labels = c("No decline", "Low decline", "Higher decline"),
      ordered = TRUE
    )
  ) %>%
  count(cannabis_group, adl_decline_category) %>%
  group_by(cannabis_group) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = cannabis_group, y = prop, fill = adl_decline_category)) +
  geom_col(color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "ADL Decline Severity by Cannabis Exposure Group",
    subtitle = "Ordinal outcome preserves no, low, and higher decline",
    x = "Past-year cannabis exposure",
    y = "Within-group percentage",
    fill = "ADL decline"
  )

fig_ordinal_by_cannabis


# ============================================================
# FIGURE 4: Observed probability of any ADL decline by cannabis group
# ============================================================

fig_any_decline_observed <- adl_any_common %>%
  filter(!is.na(adl_total_any_explore), !is.na(cannabis_group)) %>%
  mutate(
    any_decline_num = ifelse(adl_total_any_explore == "any_decline", 1, 0)
  ) %>%
  group_by(cannabis_group) %>%
  summarise(
    n = n(),
    events = sum(any_decline_num),
    prop_any_decline = mean(any_decline_num),
    se = sqrt(prop_any_decline * (1 - prop_any_decline) / n),
    lower = prop_any_decline - 1.96 * se,
    upper = prop_any_decline + 1.96 * se,
    .groups = "drop"
  ) %>%
  ggplot(aes(x = cannabis_group, y = prop_any_decline)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.10) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "Observed Probability of Any ADL Decline by Cannabis Group",
    subtitle = "Binary model is interpretable but loses severity information",
    x = "Past-year cannabis exposure",
    y = "Observed probability of any ADL decline"
  )

fig_any_decline_observed

# ============================================================
# FIGURE 5: Predicted probabilities from primary ordinal model
# Corrected version
# ============================================================

library(MASS)
library(emmeans)
library(dplyr)
library(ggplot2)
library(scales)

# Fit primary ordinal model
model_ord_primary <- MASS::polr(
  adl_total_ordinal_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

# IMPORTANT:
# Include the outcome variable after the | so emmeans gives probabilities
# for each ADL decline category rather than averaging across categories.
emm_ord_primary <- emmeans(
  model_ord_primary,
  ~ du_mar4_12m_aBin_ord | adl_total_ordinal_explore,
  mode = "prob"
)

emm_ord_primary_df <- as.data.frame(emm_ord_primary)

# Check column names
names(emm_ord_primary_df)
head(emm_ord_primary_df)

# ============================================================
# Clean dataframe for plotting
# ============================================================

emm_ord_primary_df_clean <- emm_ord_primary_df %>%
  mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    ),
    adl_decline_category = factor(
      adl_total_ordinal_explore,
      levels = c("none", "low_decline", "higher_decline"),
      labels = c("No decline", "Low decline", "Higher decline"),
      ordered = TRUE
    )
  )

emm_ord_primary_df_clean

# ============================================================
# Plot predicted probabilities
# ============================================================

fig5_predicted_probabilities <- ggplot(
  emm_ord_primary_df_clean,
  aes(
    x = cannabis_group,
    y = prob,
    group = adl_decline_category
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.10,
    linewidth = 0.6
  ) +
  facet_wrap(~ adl_decline_category) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "Predicted Probability of ADL Decline Category by Cannabis Exposure",
    subtitle = "Primary ordinal model adjusted for nadir CD4, age, education, and race/ethnicity",
    x = "Past-year cannabis exposure",
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 13)

fig5_predicted_probabilities

# ============================================================
# RAW ADL TOTAL SENSITIVITY MODEL
# Outcome: adl_total
# Purpose: Examine original continuous ADL decline score
# ============================================================

library(dplyr)
library(ggplot2)
library(broom)
library(lmtest)
library(sandwich)
library(performance)

# ------------------------------------------------------------
# 1. Descriptives for raw adl_total
# ------------------------------------------------------------

summary(adl_work$adl_total)

table(adl_work$adl_total == 0, useNA = "ifany")

mean(adl_work$adl_total == 0, na.rm = TRUE)

psych::describe(adl_work$adl_total)

# Histogram
ggplot(adl_work, aes(x = adl_total)) +
  geom_histogram(bins = 25) +
  labs(
    title = "Distribution of Raw ADL Decline Score",
    subtitle = "adl_total; higher values indicate greater ADL decline",
    x = "Raw ADL decline score",
    y = "Count"
  ) +
  theme_classic(base_size = 13)

# ------------------------------------------------------------
# 2. Fit raw linear model
# ------------------------------------------------------------

model_adl_total_raw <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

summary(model_adl_total_raw)

# ------------------------------------------------------------
# 3. HC3 robust standard errors
# ------------------------------------------------------------
# Use this because raw adl_total showed heteroscedasticity/non-normality.
# HC3 helps protect standard errors against heteroscedasticity,
# but it does not fix the zero-heavy distribution.

model_adl_total_raw_hc3 <- lmtest::coeftest(
  model_adl_total_raw,
  vcov. = sandwich::vcovHC(model_adl_total_raw, type = "HC3")
)

model_adl_total_raw_hc3

# ------------------------------------------------------------
# 4. Tidy model output for table
# ------------------------------------------------------------

model_adl_total_raw_results <- broom::tidy(
  model_adl_total_raw,
  conf.int = TRUE
)

model_adl_total_raw_results

# Optional: HC3 table
model_adl_total_raw_hc3_table <- data.frame(
  term = rownames(model_adl_total_raw_hc3),
  estimate = model_adl_total_raw_hc3[, "Estimate"],
  std.error_HC3 = model_adl_total_raw_hc3[, "Std. Error"],
  statistic_HC3 = model_adl_total_raw_hc3[, "t value"],
  p.value_HC3 = model_adl_total_raw_hc3[, "Pr(>|t|)"]
)

model_adl_total_raw_hc3_table

# ------------------------------------------------------------
# 5. Linear model diagnostics
# ------------------------------------------------------------

performance::check_heteroscedasticity(model_adl_total_raw)
performance::check_normality(model_adl_total_raw)
performance::check_collinearity(model_adl_total_raw)

# Base R diagnostic plots
par(mfrow = c(2, 2))
plot(model_adl_total_raw)
par(mfrow = c(1, 1))

# ------------------------------------------------------------
# 6. Influence diagnostics
# ------------------------------------------------------------

influence_table_adl_total_raw <- data.frame(
  case_id = as.numeric(rownames(model.frame(model_adl_total_raw))),
  fitted = fitted(model_adl_total_raw),
  residual = residuals(model_adl_total_raw),
  standardized_residual = rstandard(model_adl_total_raw),
  studentized_residual = rstudent(model_adl_total_raw),
  leverage = hatvalues(model_adl_total_raw),
  cooks_d = cooks.distance(model_adl_total_raw),
  adl_total = model.frame(model_adl_total_raw)$adl_total
) %>%
  arrange(desc(cooks_d))

head(influence_table_adl_total_raw, 15)

# Common rule-of-thumb cutoff for Cook's D
cooks_cutoff <- 4 / nobs(model_adl_total_raw)

influence_table_adl_total_raw %>%
  filter(cooks_d > cooks_cutoff)

# ------------------------------------------------------------
# 7. Model fit / power-relevant values
# ------------------------------------------------------------

raw_lm_fit_summary <- tibble::tibble(
  n = nobs(model_adl_total_raw),
  r_squared = summary(model_adl_total_raw)$r.squared,
  adjusted_r_squared = summary(model_adl_total_raw)$adj.r.squared,
  residual_standard_error = summary(model_adl_total_raw)$sigma,
  f_statistic = unname(summary(model_adl_total_raw)$fstatistic[1]),
  df_model = unname(summary(model_adl_total_raw)$fstatistic[2]),
  df_residual = unname(summary(model_adl_total_raw)$fstatistic[3]),
  model_p_value = pf(
    summary(model_adl_total_raw)$fstatistic[1],
    summary(model_adl_total_raw)$fstatistic[2],
    summary(model_adl_total_raw)$fstatistic[3],
    lower.tail = FALSE
  )
)

raw_lm_fit_summary

# ------------------------------------------------------------
# 8. Incremental R2 for cannabis block
# ------------------------------------------------------------

model_adl_total_no_cannabis <- lm(
  adl_total ~ nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

model_adl_total_with_cannabis <- lm(
  adl_total ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)

anova(model_adl_total_no_cannabis, model_adl_total_with_cannabis)

r2_no_cannabis <- summary(model_adl_total_no_cannabis)$r.squared
r2_with_cannabis <- summary(model_adl_total_with_cannabis)$r.squared

delta_r2_cannabis <- r2_with_cannabis - r2_no_cannabis

raw_lm_cannabis_block_summary <- tibble::tibble(
  r2_no_cannabis = r2_no_cannabis,
  r2_with_cannabis = r2_with_cannabis,
  delta_r2_cannabis = delta_r2_cannabis
)

raw_lm_cannabis_block_summary

# ------------------------------------------------------------
# 9. Predicted raw ADL total by cannabis group
# ------------------------------------------------------------

library(emmeans)

emm_adl_total_raw <- emmeans(
  model_adl_total_raw,
  ~ du_mar4_12m_aBin_ord
)

emm_adl_total_raw_df <- as.data.frame(emm_adl_total_raw) %>%
  mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    )
  )

emm_adl_total_raw_df

# ------------------------------------------------------------
# 10. Figure: adjusted predicted raw ADL total by cannabis group
# ------------------------------------------------------------

fig_raw_adl_total_predicted <- ggplot(
  emm_adl_total_raw_df,
  aes(x = cannabis_group, y = emmean, group = 1)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.10,
    linewidth = 0.6
  ) +
  labs(
    title = "Adjusted Mean Raw ADL Decline Score by Cannabis Exposure",
    subtitle = "Linear sensitivity model adjusted for nadir CD4, age, education, and race/ethnicity",
    x = "Past-year cannabis exposure",
    y = "Adjusted mean raw ADL decline score"
  ) +
  theme_classic(base_size = 13)

fig_raw_adl_total_predicted

model_adl_total_raw_linear_only <- lm(
  adl_total ~ cu12m_L +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work
)


# ============================================================
# COVARIATE DECISION CHECK FOR ADL PAPER
# ============================================================
# Purpose:
# Identify which variables should be considered covariates
# for cannabis --> ADL models.
#
# Rule:
# A variable is a stronger covariate candidate if it is:
#   1) clinically/theoretically plausible, AND
#   2) related to cannabis exposure, AND
#   3) related to the ADL outcome.
#
# Important:
# Do not automatically include every significant correlate.
# Some variables, such as medication adherence indicators,
# may be parallel functional outcomes rather than confounds.
# ============================================================

library(dplyr)
library(psych)
library(broom)
library(tidyr)
library(purrr)

# ============================================================
# 1. Check candidate variable names
# ============================================================

# Helper function to search variable names
find_vars <- function(pattern, data = adl_work) {
  grep(pattern, names(data), ignore.case = TRUE, value = TRUE)
}

# Search for nicotine/tobacco variables because exact name may vary
find_vars("nic", adl_work)
find_vars("tob", adl_work)
find_vars("cig", adl_work)
find_vars("smok", adl_work)
find_vars("fager", adl_work)

# Search alcohol variables
find_vars("alc", adl_work)

# Search mood variables
find_vars("bdi", adl_work)
find_vars("bai", adl_work)

# Search ADL variables
find_vars("adl", adl_work)

# ============================================================
# 2. Create transformed alcohol variables
#    New variables only; raw variables preserved
# ============================================================

# Past-year alcohol
if ("du_alc4_12m_a" %in% names(adl_work)) {
  adl_work$du_alc4_12m_a_log1p <- log1p(adl_work$du_alc4_12m_a)
  adl_work$du_alc4_12m_a_log1p_c <- as.numeric(
    scale(adl_work$du_alc4_12m_a_log1p, center = TRUE, scale = FALSE)
  )
}

# Past-30-day alcohol
if ("du_alc6_30d_a" %in% names(adl_work)) {
  adl_work$du_alc6_30d_a_log1p <- log1p(adl_work$du_alc6_30d_a)
  adl_work$du_alc6_30d_a_log1p_c <- as.numeric(
    scale(adl_work$du_alc6_30d_a_log1p, center = TRUE, scale = FALSE)
  )
}

summary(adl_work[, intersect(
  c(
    "du_alc4_12m_a",
    "du_alc4_12m_a_log1p",
    "du_alc6_30d_a",
    "du_alc6_30d_a_log1p"
  ),
  names(adl_work)
)])

# ============================================================
# 3. Define candidate covariates
# ============================================================

# Core covariates already planned for primary model
core_covariates <- c(
  "phq_2_age_c",
  "phq_7_degree_c",
  "race_eth_binary_covfac",
  "nadirsqrt_c"
)

# Candidate sensitivity covariates
mood_covariates <- c(
  "bdi_total",
  "bai_total"
)

alcohol_covariates <- c(
  "du_alc4_12m_a_log1p_c",
  "du_alc6_30d_a_log1p_c"
)

# Update this after checking find_vars("nic"), find_vars("tob"), etc.
# Examples only:
possible_nicotine_covariates <- c(
  "nicotine_total",
  "ftnd_total",
  "tobacco_use",
  "smoking_status",
  "du_tobacco_12m",
  "du_nic_12m"
)

# Keep only variables that actually exist
mood_covariates <- intersect(mood_covariates, names(adl_work))
alcohol_covariates <- intersect(alcohol_covariates, names(adl_work))
nicotine_covariates <- intersect(possible_nicotine_covariates, names(adl_work))

candidate_covariates <- c(
  core_covariates,
  mood_covariates,
  alcohol_covariates,
  nicotine_covariates
)

candidate_covariates <- intersect(candidate_covariates, names(adl_work))

candidate_covariates

# ============================================================
# 4. Create numeric ADL outcome versions for screening
# ============================================================

# Binary any-decline outcome
adl_work$adl_total_any_num_screen <- ifelse(
  adl_work$adl_total_any_explore == "any_decline", 1,
  ifelse(adl_work$adl_total_any_explore == "no_decline", 0, NA)
)

# Ordinal ADL decline outcome
adl_work$adl_total_ordinal_num_screen <- as.numeric(
  adl_work$adl_total_ordinal_explore
)

# Check coding
table(adl_work$adl_total_any_explore, adl_work$adl_total_any_num_screen, useNA = "ifany")
table(adl_work$adl_total_ordinal_explore, adl_work$adl_total_ordinal_num_screen, useNA = "ifany")

# ============================================================
# FIX: Make numeric screening versions of non-numeric covariates
# ============================================================

# Check variable classes
sapply(adl_work[, candidate_covariates], class)

# Create numeric race dummy for screening only
# This does NOT overwrite the original race variable.
adl_work$race_eth_binary_covnum_screen <- ifelse(
  adl_work$race_eth_binary_covfac == "White non-Hispanic", 1,
  ifelse(is.na(adl_work$race_eth_binary_covfac), NA, 0)
)

table(
  adl_work$race_eth_binary_covfac,
  adl_work$race_eth_binary_covnum_screen,
  useNA = "ifany"
)

#most relevant nicotine variables
#du_nic4_12m_a  # past-year nicotine amount
#du_nic6_30d_a  # past-30-day nicotine amount
#ftnd_cig       # cigarette/nicotine dependence-related variable

# ============================================================
# Create transformed alcohol and nicotine variables
# New variables only; raw variables preserved
# ============================================================

# Alcohol
adl_work$du_alc4_12m_a_log1p <- log1p(adl_work$du_alc4_12m_a)
adl_work$du_alc6_30d_a_log1p <- log1p(adl_work$du_alc6_30d_a)

adl_work$du_alc4_12m_a_log1p_c <- as.numeric(
  scale(adl_work$du_alc4_12m_a_log1p, center = TRUE, scale = FALSE)
)

adl_work$du_alc6_30d_a_log1p_c <- as.numeric(
  scale(adl_work$du_alc6_30d_a_log1p, center = TRUE, scale = FALSE)
)

# Nicotine
adl_work$du_nic4_12m_a_log1p <- log1p(adl_work$du_nic4_12m_a)
adl_work$du_nic6_30d_a_log1p <- log1p(adl_work$du_nic6_30d_a)

adl_work$du_nic4_12m_a_log1p_c <- as.numeric(
  scale(adl_work$du_nic4_12m_a_log1p, center = TRUE, scale = FALSE)
)

adl_work$du_nic6_30d_a_log1p_c <- as.numeric(
  scale(adl_work$du_nic6_30d_a_log1p, center = TRUE, scale = FALSE)
)

# FTND cigarettes if present
if ("ftnd_cig" %in% names(adl_work)) {
  adl_work$ftnd_cig_c <- as.numeric(
    scale(adl_work$ftnd_cig, center = TRUE, scale = FALSE)
  )
}

# ============================================================
# Candidate covariates for numeric screening
# ============================================================

candidate_covariates_numeric <- c(
  # Core covariates
  "phq_2_age_c",
  "phq_7_degree_c",
  "race_eth_binary_covnum_screen",
  "nadirsqrt_c",
  
  # Mood symptoms
  "bdi_total",
  "bai_total",
  
  # Alcohol
  "du_alc4_12m_a_log1p_c",
  "du_alc6_30d_a_log1p_c",
  
  # Nicotine
  "du_nic4_12m_a_log1p_c",
  "du_nic6_30d_a_log1p_c",
  "ftnd_cig_c"
)

candidate_covariates_numeric <- intersect(
  candidate_covariates_numeric,
  names(adl_work)
)

candidate_covariates_numeric

# Confirm all are numeric
sapply(adl_work[, candidate_covariates_numeric], class)

# ============================================================
# Create numeric cannabis and ADL screening variables
# ============================================================

# Cannabis rank
adl_work$cu12m_rank_screen <- as.numeric(adl_work$du_mar4_12m_aBin_ord)

# Raw/log cannabis amount if available
if ("du_mar4_12m_a" %in% names(adl_work)) {
  adl_work$du_mar4_12m_a_log1p_screen <- log1p(adl_work$du_mar4_12m_a)
}

# Binary ADL decline
adl_work$adl_total_any_num_screen <- ifelse(
  adl_work$adl_total_any_explore == "any_decline", 1,
  ifelse(adl_work$adl_total_any_explore == "no_decline", 0, NA)
)

# Ordinal ADL decline: none = 1, low = 2, higher = 3
adl_work$adl_total_ordinal_num_screen <- as.numeric(
  adl_work$adl_total_ordinal_explore
)

cannabis_screen_vars <- intersect(
  c(
    "cu12m_rank_screen",
    "du_mar4_12m_a_log1p_screen",
    "cu12m_L",
    "cu12m_Q"
  ),
  names(adl_work)
)

adl_screen_outcomes <- intersect(
  c(
    "adl_total",
    "adl_total_any_num_screen",
    "adl_total_ordinal_num_screen",
    "adl_b_sum_now",
    "adl_b_sum_now_excess"
  ),
  names(adl_work)
)

# ============================================================
# Candidate covariates x cannabis
# ============================================================

covariate_cannabis_corr <- psych::corr.test(
  adl_work[, c(candidate_covariates_numeric, cannabis_screen_vars)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

covariate_cannabis_table <- as.data.frame(
  as.table(covariate_cannabis_corr$r[
    candidate_covariates_numeric,
    cannabis_screen_vars
  ])
)

names(covariate_cannabis_table) <- c(
  "candidate_covariate",
  "cannabis_variable",
  "rho_with_cannabis"
)

covariate_cannabis_table$p_with_cannabis <- as.vector(
  covariate_cannabis_corr$p[
    candidate_covariates_numeric,
    cannabis_screen_vars
  ]
)

covariate_cannabis_table <- covariate_cannabis_table %>%
  arrange(p_with_cannabis)

covariate_cannabis_table

# ============================================================
# Candidate covariates x ADL outcomes
# ============================================================

covariate_adl_corr <- psych::corr.test(
  adl_work[, c(candidate_covariates_numeric, adl_screen_outcomes)],
  method = "spearman",
  use = "pairwise.complete.obs",
  adjust = "none"
)

covariate_adl_table <- as.data.frame(
  as.table(covariate_adl_corr$r[
    candidate_covariates_numeric,
    adl_screen_outcomes
  ])
)

names(covariate_adl_table) <- c(
  "candidate_covariate",
  "adl_outcome",
  "rho_with_adl"
)

covariate_adl_table$p_with_adl <- as.vector(
  covariate_adl_corr$p[
    candidate_covariates_numeric,
    adl_screen_outcomes
  ]
)

covariate_adl_table <- covariate_adl_table %>%
  arrange(p_with_adl)

covariate_adl_table

# ============================================================
# Covariate decision table
# ============================================================

covariate_cannabis_summary <- covariate_cannabis_table %>%
  group_by(candidate_covariate) %>%
  slice_min(order_by = p_with_cannabis, n = 1, with_ties = FALSE) %>%
  ungroup()

covariate_adl_summary <- covariate_adl_table %>%
  group_by(candidate_covariate) %>%
  slice_min(order_by = p_with_adl, n = 1, with_ties = FALSE) %>%
  ungroup()

covariate_decision_table <- covariate_cannabis_summary %>%
  full_join(covariate_adl_summary, by = "candidate_covariate") %>%
  mutate(
    related_to_cannabis_p05 = p_with_cannabis < .05,
    related_to_adl_p05 = p_with_adl < .05,
    related_to_cannabis_p10 = p_with_cannabis < .10,
    related_to_adl_p10 = p_with_adl < .10,
    
    empirical_confounder_p05 = related_to_cannabis_p05 & related_to_adl_p05,
    empirical_confounder_p10 = related_to_cannabis_p10 & related_to_adl_p10,
    
    suggested_role = case_when(
      candidate_covariate %in% c(
        "phq_2_age_c",
        "phq_7_degree_c",
        "race_eth_binary_covnum_screen",
        "nadirsqrt_c"
      ) ~ "Core covariate retained a priori",
      
      empirical_confounder_p05 ~
        "Strong empirical covariate candidate",
      
      empirical_confounder_p10 ~
        "Possible covariate candidate; discuss",
      
      related_to_adl_p05 & !related_to_cannabis_p05 ~
        "ADL correlate only; sensitivity covariate, not primary",
      
      related_to_cannabis_p05 & !related_to_adl_p05 ~
        "Cannabis correlate only; not a clear confounder",
      
      TRUE ~
        "Not supported as covariate by screening"
    )
  ) %>%
  arrange(
    desc(empirical_confounder_p05),
    desc(empirical_confounder_p10),
    p_with_adl
  )

covariate_decision_table

# ============================================================
# Ordinal ADL outcome cutoffs/ranges
# ============================================================

adl_ordinal_ranges <- adl_work %>%
  filter(!is.na(adl_total), !is.na(adl_total_ordinal_explore)) %>%
  group_by(adl_total_ordinal_explore) %>%
  summarise(
    n = n(),
    min_adl_total = min(adl_total, na.rm = TRUE),
    max_adl_total = max(adl_total, na.rm = TRUE),
    mean_adl_total = mean(adl_total, na.rm = TRUE),
    sd_adl_total = sd(adl_total, na.rm = TRUE),
    median_adl_total = median(adl_total, na.rm = TRUE),
    .groups = "drop"
  )

adl_ordinal_ranges

# ============================================================
# SEX COVARIATE CHECK
# ============================================================
# Purpose:
# Check whether sex should be included as a covariate in ADL models.
#
# Rule:
# Sex is a stronger covariate candidate if it is associated with:
#   1) cannabis exposure, and
#   2) ADL outcome.
#
# Because sex is binary/categorical, use tables + Fisher's exact tests
# rather than relying only on Spearman correlations.
# ============================================================

library(dplyr)
library(broom)

# ============================================================
# 1. Confirm sex variable exists
# ============================================================

find_vars("sex", adl_work)
find_vars("gender", adl_work)

table(adl_work$sex_covfac, useNA = "ifany")
table(adl_work$sex_covnum, useNA = "ifany")

# Sex distribution by cannabis exposure group
table(adl_work$sex_covfac, adl_work$du_mar4_12m_aBin_ord, useNA = "ifany")

# Row percentages
prop.table(
  table(adl_work$sex_covfac, adl_work$du_mar4_12m_aBin_ord),
  margin = 1
)

# Fisher's exact test because some cells may be small
fisher.test(
  table(adl_work$sex_covfac, adl_work$du_mar4_12m_aBin_ord)
)

# Sex distribution by any ADL decline
table(adl_work$sex_covfac, adl_work$adl_total_any_explore, useNA = "ifany")

# Row percentages
prop.table(
  table(adl_work$sex_covfac, adl_work$adl_total_any_explore),
  margin = 1
)

# Fisher's exact test
fisher.test(
  table(adl_work$sex_covfac, adl_work$adl_total_any_explore)
)

# Sex distribution by ordinal ADL decline category
table(adl_work$sex_covfac, adl_work$adl_total_ordinal_explore, useNA = "ifany")

# Row percentages
prop.table(
  table(adl_work$sex_covfac, adl_work$adl_total_ordinal_explore),
  margin = 1
)

# Fisher's exact test
fisher.test(
  table(adl_work$sex_covfac, adl_work$adl_total_ordinal_explore)
)

# Main ordinal model without sex
model_ord_no_sex <- MASS::polr(
  adl_total_ordinal_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

# Same model with sex
model_ord_with_sex <- MASS::polr(
  adl_total_ordinal_explore ~ du_mar4_12m_aBin_ord +
    nadirsqrt_c +
    phq_2_age_c +
    phq_7_degree_c +
    race_eth_binary_covfac +
    sex_covfac,
  data = adl_work,
  Hess = TRUE
)

# Compare model fit
AIC(model_ord_no_sex, model_ord_with_sex)
BIC(model_ord_no_sex, model_ord_with_sex)

# Compare coefficients
broom::tidy(model_ord_no_sex, conf.int = TRUE, exponentiate = TRUE)
broom::tidy(model_ord_with_sex, conf.int = TRUE, exponentiate = TRUE)

# ============================================================
# ORDINAL ADL DECLINE: RANGES, CUTOFFS, AND CATEGORY GAPS
# ============================================================
# Purpose:
# To document how the ordinal ADL decline categories map onto the
# original continuous adl_total score.
#
# Outcome:
# adl_total_ordinal_explore = no decline / low decline / higher decline
#
# Original score:
# adl_total = continuous ADL decline score
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Range of adl_total values within each ordinal category
# ------------------------------------------------------------

adl_ordinal_ranges <- adl_work %>%
  filter(
    !is.na(adl_total),
    !is.na(adl_total_ordinal_explore)
  ) %>%
  group_by(adl_total_ordinal_explore) %>%
  summarise(
    n = n(),
    min_adl_total = min(adl_total, na.rm = TRUE),
    max_adl_total = max(adl_total, na.rm = TRUE),
    mean_adl_total = mean(adl_total, na.rm = TRUE),
    sd_adl_total = sd(adl_total, na.rm = TRUE),
    median_adl_total = median(adl_total, na.rm = TRUE),
    .groups = "drop"
  )

adl_ordinal_ranges

# ------------------------------------------------------------
# 2. Gap between adjacent ordinal categories
# ------------------------------------------------------------

adl_ordinal_gap_table <- adl_ordinal_ranges %>%
  arrange(adl_total_ordinal_explore) %>%
  mutate(
    next_category = lead(adl_total_ordinal_explore),
    next_min_adl_total = lead(min_adl_total),
    gap_to_next_category = next_min_adl_total - max_adl_total
  )

adl_ordinal_gap_table

# ------------------------------------------------------------
# 3. Exact observed adl_total values by ordinal category
# ------------------------------------------------------------

adl_exact_values_by_category <- adl_work %>%
  filter(
    !is.na(adl_total),
    !is.na(adl_total_ordinal_explore)
  ) %>%
  arrange(adl_total) %>%
  dplyr::select(
    adl_total,
    adl_total_ordinal_explore
  ) %>%
  distinct()

adl_exact_values_by_category

adl_exact_values_by_category <- adl_work %>%
  dplyr::filter(
    !is.na(adl_total),
    !is.na(adl_total_ordinal_explore)
  ) %>%
  dplyr::arrange(adl_total) %>%
  dplyr::select(
    adl_total,
    adl_total_ordinal_explore
  ) %>%
  dplyr::distinct()

adl_exact_values_by_category

# ============================================================
# Exact adl_total values and frequencies by ordinal category
# ============================================================

adl_value_frequency_by_category <- adl_work %>%
  dplyr::filter(
    !is.na(adl_total),
    !is.na(adl_total_ordinal_explore)
  ) %>%
  dplyr::count(
    adl_total_ordinal_explore,
    adl_total,
    name = "n_at_value"
  ) %>%
  dplyr::arrange(adl_total)

adl_value_frequency_by_category

# ============================================================
# Low vs higher decline descriptive difference
# ============================================================

positive_decline_only <- adl_work %>%
  dplyr::filter(
    adl_total_ordinal_explore %in% c("low_decline", "higher_decline"),
    !is.na(adl_total)
  )

wilcox.test(
  adl_total ~ adl_total_ordinal_explore,
  data = positive_decline_only,
  exact = FALSE
)

positive_decline_only %>%
  dplyr::group_by(adl_total_ordinal_explore) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(adl_total, na.rm = TRUE),
    median = median(adl_total, na.rm = TRUE),
    iqr = IQR(adl_total, na.rm = TRUE),
    min = min(adl_total, na.rm = TRUE),
    max = max(adl_total, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# Sensitivity: higher ADL decline vs no/low decline
# ============================================================

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_higher_decline_binary = dplyr::case_when(
      adl_total_ordinal_explore == "higher_decline" ~ "higher_decline",
      adl_total_ordinal_explore %in% c("none", "low_decline") ~ "none_or_low",
      TRUE ~ NA_character_
    ),
    adl_higher_decline_binary = factor(
      adl_higher_decline_binary,
      levels = c("none_or_low", "higher_decline")
    )
  )

table(adl_work$adl_higher_decline_binary, useNA = "ifany")

model_higher_decline_sensitivity <- glm(
  adl_higher_decline_binary ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    nadirsqrt_c,
  data = adl_work,
  family = binomial(link = "logit")
)

broom::tidy(
  model_higher_decline_sensitivity,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# Observed adl_total values and frequencies
# ============================================================

adl_value_counts <- adl_work %>%
  dplyr::filter(!is.na(adl_total)) %>%
  dplyr::count(adl_total, name = "n") %>%
  dplyr::arrange(adl_total) %>%
  dplyr::mutate(
    cumulative_n = cumsum(n),
    cumulative_percent = cumulative_n / sum(n)
  )

adl_value_counts

# ============================================================
# Quantiles among participants with any ADL decline
# ============================================================

positive_adl <- adl_work %>%
  dplyr::filter(!is.na(adl_total), adl_total > 0)

quantile(
  positive_adl$adl_total,
  probs = c(.25, .33, .50, .67, .75),
  na.rm = TRUE
)

summary(positive_adl$adl_total)

# ============================================================
# Alternative ordinal outcome:
# no decline / low-moderate decline / top-quartile decline
# ============================================================

adl_top_quartile_cutoff <- quantile(
  adl_work$adl_total,
  probs = .75,
  na.rm = TRUE
)

adl_top_quartile_cutoff

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_ord_top_quartile = dplyr::case_when(
      is.na(adl_total) ~ NA_character_,
      adl_total == 0 ~ "no_decline",
      adl_total > 0 & adl_total < adl_top_quartile_cutoff ~ "low_moderate_decline",
      adl_total >= adl_top_quartile_cutoff ~ "higher_decline"
    ),
    adl_ord_top_quartile = factor(
      adl_ord_top_quartile,
      levels = c("no_decline", "low_moderate_decline", "higher_decline"),
      ordered = TRUE
    )
  )

table(adl_work$adl_ord_top_quartile, useNA = "ifany")

# ============================================================
# Clinically/interpretable proportion-based cutoffs
# ============================================================

cutoffs_to_check <- c(.10, .15, .20, .25)

cutoff_count_table <- purrr::map_dfr(
  cutoffs_to_check,
  function(cutoff) {
    
    temp_outcome <- dplyr::case_when(
      is.na(adl_work$adl_total) ~ NA_character_,
      adl_work$adl_total == 0 ~ "no_decline",
      adl_work$adl_total > 0 & adl_work$adl_total < cutoff ~ "low_decline",
      adl_work$adl_total >= cutoff ~ "higher_decline"
    )
    
    tibble::tibble(
      cutoff = cutoff,
      no_decline = sum(temp_outcome == "no_decline", na.rm = TRUE),
      low_decline = sum(temp_outcome == "low_decline", na.rm = TRUE),
      higher_decline = sum(temp_outcome == "higher_decline", na.rm = TRUE),
      total_n = sum(!is.na(temp_outcome))
    )
  }
)

cutoff_count_table

# ============================================================
# Sensitivity models across multiple ordinal ADL cutoffs
# ============================================================

library(MASS)
library(broom)
library(purrr)
library(dplyr)

run_ordinal_cutoff_model <- function(cutoff) {
  
  temp_data <- adl_work %>%
    dplyr::mutate(
      adl_ord_temp = dplyr::case_when(
        is.na(adl_total) ~ NA_character_,
        adl_total == 0 ~ "no_decline",
        adl_total > 0 & adl_total < cutoff ~ "low_decline",
        adl_total >= cutoff ~ "higher_decline"
      ),
      adl_ord_temp = factor(
        adl_ord_temp,
        levels = c("no_decline", "low_decline", "higher_decline"),
        ordered = TRUE
      )
    )
  
  temp_model <- MASS::polr(
    adl_ord_temp ~ du_mar4_12m_aBin_ord +
      phq_2_age_c +
      bdi_total +
      sex_covfac +
      nadirsqrt_c,
    data = temp_data,
    Hess = TRUE
  )
  
  broom::tidy(
    temp_model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::filter(grepl("du_mar4_12m_aBin_ord", term)) %>%
    dplyr::mutate(cutoff = cutoff)
}

ordinal_cutoff_sensitivity_results <- purrr::map_dfr(
  c(.10, .13, .14, .15, .20, .25),
  run_ordinal_cutoff_model
)

ordinal_cutoff_sensitivity_results

# ============================================================
# REBUILD ORDINAL ADL OUTCOME USING NUMBER OF DECLINED DOMAINS
# ============================================================
# Purpose:
# Create a clinically interpretable ordinal ADL decline outcome:
#   0 declined domains
#   1 declined domain
#   2+ declined domains
#
# This avoids relying on a small ratio-score cutoff such as .13 vs .14.
#
# Important:
# This assumes:
#   adl_b_sum_now  = current ADL difficulty burden
#   adl_c_sum_best = best/premorbid ADL difficulty burden
#
# New variables only; raw variables are preserved.
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Create number of declined ADL domains
# ------------------------------------------------------------

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_declined_domain_count = adl_b_sum_now - adl_c_sum_best
  )

# Inspect the raw count
summary(adl_work$adl_declined_domain_count)
table(adl_work$adl_declined_domain_count, useNA = "ifany")

# ------------------------------------------------------------
# 2. Check for negative values
# ------------------------------------------------------------
# Negative values would mean current functioning is better than best/premorbid
# functioning, which may reflect coding issues or improvement rather than decline.

adl_work %>%
  dplyr::filter(adl_declined_domain_count < 0) %>%
  dplyr::select(
    adl_b_sum_now,
    adl_c_sum_best,
    adl_declined_domain_count,
    adl_total
  )

# Count negative values
sum(adl_work$adl_declined_domain_count < 0, na.rm = TRUE)

# ------------------------------------------------------------
# 3. Create cleaned declined-domain count
# ------------------------------------------------------------

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_declined_domain_count_clean = dplyr::case_when(
      is.na(adl_declined_domain_count) ~ NA_real_,
      adl_declined_domain_count < 0 ~ 0,
      TRUE ~ adl_declined_domain_count
    )
  )

table(adl_work$adl_declined_domain_count_clean, useNA = "ifany")

# ------------------------------------------------------------
# 4. Create ordinal ADL declined-domain outcome
# ------------------------------------------------------------

adl_work <- adl_work %>%
  dplyr::mutate(
    adl_declined_domain_ord = dplyr::case_when(
      is.na(adl_declined_domain_count_clean) ~ NA_character_,
      adl_declined_domain_count_clean == 0 ~ "0_no_decline",
      adl_declined_domain_count_clean == 1 ~ "1_domain_decline",
      adl_declined_domain_count_clean >= 2 ~ "2plus_domain_decline"
    ),
    adl_declined_domain_ord = factor(
      adl_declined_domain_ord,
      levels = c(
        "0_no_decline",
        "1_domain_decline",
        "2plus_domain_decline"
      ),
      ordered = TRUE
    )
  )

table(adl_work$adl_declined_domain_ord, useNA = "ifany")
prop.table(table(adl_work$adl_declined_domain_ord, useNA = "ifany"))

# ------------------------------------------------------------
# 5. Compare new domain-count ordinal outcome to old ratio-based ordinal outcome
# ------------------------------------------------------------

table(
  old_ratio_ordinal = adl_work$adl_total_ordinal_explore,
  new_domain_ordinal = adl_work$adl_declined_domain_ord,
  useNA = "ifany"
)

prop.table(
  table(
    old_ratio_ordinal = adl_work$adl_total_ordinal_explore,
    new_domain_ordinal = adl_work$adl_declined_domain_ord,
    useNA = "ifany"
  ),
  margin = 1
)

# ------------------------------------------------------------
# 6. Describe adl_total within new ordinal categories
# ------------------------------------------------------------

adl_domain_ord_ranges <- adl_work %>%
  dplyr::filter(
    !is.na(adl_total),
    !is.na(adl_declined_domain_ord)
  ) %>%
  dplyr::group_by(adl_declined_domain_ord) %>%
  dplyr::summarise(
    n = n(),
    min_adl_total = min(adl_total, na.rm = TRUE),
    max_adl_total = max(adl_total, na.rm = TRUE),
    mean_adl_total = mean(adl_total, na.rm = TRUE),
    sd_adl_total = sd(adl_total, na.rm = TRUE),
    median_adl_total = median(adl_total, na.rm = TRUE),
    .groups = "drop"
  )

adl_domain_ord_ranges

# ============================================================
# PRIMARY ORDINAL MODEL USING DECLINED-DOMAIN OUTCOME
# ============================================================

library(MASS)
library(broom)

model_ord_domain_primary <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    nadirsqrt_c,
  data = adl_work,
  Hess = TRUE
)

summary(model_ord_domain_primary)

broom::tidy(
  model_ord_domain_primary,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# DEMOGRAPHIC-EXPANDED ORDINAL MODEL
# Adds education and race/ethnicity
# ============================================================

model_ord_domain_demo_expanded <- MASS::polr(
  adl_declined_domain_ord ~ du_mar4_12m_aBin_ord +
    phq_2_age_c +
    bdi_total +
    sex_covfac +
    nadirsqrt_c +
    phq_7_degree_c +
    race_eth_binary_covfac,
  data = adl_work,
  Hess = TRUE
)

summary(model_ord_domain_demo_expanded)

broom::tidy(
  model_ord_domain_demo_expanded,
  conf.int = TRUE,
  exponentiate = TRUE
)

# ============================================================
# MODEL COMPARISON
# ============================================================

AIC(
  model_ord_domain_primary,
  model_ord_domain_demo_expanded
)

BIC(
  model_ord_domain_primary,
  model_ord_domain_demo_expanded
)

# ============================================================
# PROPORTIONAL ODDS ASSUMPTION CHECK
# ============================================================

# install.packages("brant") # if needed
library(brant)

brant::brant(model_ord_domain_primary)
brant::brant(model_ord_domain_demo_expanded)

# ============================================================
# PREDICTED PROBABILITIES FOR DOMAIN-COUNT ORDINAL MODEL
# ============================================================

library(emmeans)
library(ggplot2)

emm_ord_domain <- emmeans(
  model_ord_domain_primary,
  ~ du_mar4_12m_aBin_ord,
  mode = "prob"
)

emm_ord_domain_df <- as.data.frame(emm_ord_domain)

names(emm_ord_domain_df)

emm_ord_domain_df

fig_ord_domain_predicted <- emm_ord_domain_df %>%
  dplyr::mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    ),
    adl_decline_category = factor(
      y.level,
      levels = c(
        "0_no_decline",
        "1_domain_decline",
        "2plus_domain_decline"
      ),
      labels = c(
        "No declined domains",
        "One declined domain",
        "Two or more declined domains"
      ),
      ordered = TRUE
    )
  ) %>%
  ggplot(
    aes(
      x = cannabis_group,
      y = prob,
      group = adl_decline_category
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.10
  ) +
  facet_wrap(~ adl_decline_category) +
  labs(
    title = "Predicted Probability of ADL Decline by Cannabis Exposure",
    subtitle = "Ordinal model using number of declined ADL domains",
    x = "Past-year cannabis exposure",
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 13)

fig_ord_domain_predicted

names(emm_ord_domain_df)

# ============================================================
# Inspect emmeans output column names
# ============================================================

emm_ord_domain <- emmeans(
  model_ord_domain_primary,
  ~ du_mar4_12m_aBin_ord,
  mode = "prob"
)

emm_ord_domain_df <- as.data.frame(emm_ord_domain)

names(emm_ord_domain_df)
head(emm_ord_domain_df)

# ============================================================
# Predicted probabilities from ordinal model
# Outcome: declined ADL domains ordinal variable
# ============================================================

library(dplyr)
library(emmeans)
library(ggplot2)

# Get predicted probabilities
emm_ord_domain <- emmeans(
  model_ord_domain_primary,
  ~ du_mar4_12m_aBin_ord,
  mode = "prob"
)

emm_ord_domain_df <- as.data.frame(emm_ord_domain)

# Inspect column names
names(emm_ord_domain_df)
head(emm_ord_domain_df)

# ------------------------------------------------------------
# Identify the outcome-category column
# ------------------------------------------------------------
# In emmeans ordinal probability output, the category column is
# usually NOT always called y.level.
# This finds the column that contains the ADL outcome category labels.

possible_outcome_cols <- names(emm_ord_domain_df)[
  sapply(emm_ord_domain_df, function(x) {
    any(as.character(x) %in% c(
      "0_no_decline",
      "1_domain_decline",
      "2plus_domain_decline"
    ))
  })
]

possible_outcome_cols

# Use the first matching column as the outcome category column
outcome_col <- possible_outcome_cols[1]

# ------------------------------------------------------------
# Clean predicted probability dataframe
# ------------------------------------------------------------

emm_ord_domain_plot_df <- emm_ord_domain_df %>%
  mutate(
    cannabis_group = factor(
      du_mar4_12m_aBin_ord,
      levels = c("none", "low", "high"),
      labels = c("None", "Low", "High"),
      ordered = TRUE
    ),
    adl_decline_category = factor(
      .data[[outcome_col]],
      levels = c(
        "0_no_decline",
        "1_domain_decline",
        "2plus_domain_decline"
      ),
      labels = c(
        "No decline",
        "1 declined domain",
        "2+ declined domains"
      ),
      ordered = TRUE
    )
  )

# Check final plotting dataframe
head(emm_ord_domain_plot_df)

# ============================================================
# Correct predicted probabilities from ordinal model
# Outcome: declined ADL domains ordinal variable
# ============================================================

library(emmeans)
library(dplyr)
library(ggplot2)

emm_ord_domain <- emmeans(
  model_ord_domain_primary,
  ~ du_mar4_12m_aBin_ord | adl_declined_domain_ord,
  mode = "prob"
)

emm_ord_domain_df <- as.data.frame(emm_ord_domain)

names(emm_ord_domain_df)
head(emm_ord_domain_df)

# ============================================================
# Clean predicted probability dataframe
# ============================================================

emm_ord_domain_plot_df <- emm_ord_domain_df %>%
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
        "No decline",
        "1 declined domain",
        "2+ declined domains"
      ),
      ordered = TRUE
    )
  )

emm_ord_domain_plot_df

# ============================================================
# Figure: Predicted probability of ADL decline category
# ============================================================

ggplot(
  emm_ord_domain_plot_df,
  aes(
    x = cannabis_group,
    y = prob,
    group = adl_decline_category,
    linetype = adl_decline_category,
    shape = adl_decline_category
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.08
  ) +
  labs(
    title = "Predicted Probability of ADL Decline Category by Cannabis Exposure",
    subtitle = "Ordinal model using declined ADL domains: 0, 1, or 2+",
    x = "Past-year cannabis exposure",
    y = "Predicted probability",
    linetype = "ADL decline category",
    shape = "ADL decline category"
  ) +
  theme_classic(base_size = 13)


ggplot(
  emm_ord_domain_plot_df,
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
    title = "Predicted Probability of ADL Decline Category by Cannabis Exposure",
    subtitle = "Ordinal model using declined ADL domains: 0, 1, or 2+",
    x = "Past-year cannabis exposure",
    y = "Predicted probability",
    linetype = "ADL decline category"
  ) +
  theme_classic(base_size = 13)

ggplot(
  emm_ord_domain_plot_df,
  aes(
    x = cannabis_group,
    y = prob,
    group = adl_decline_category
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.08
  ) +
  facet_wrap(~ adl_decline_category) +
  labs(
    title = "Predicted Probability of ADL Decline by Cannabis Exposure",
    subtitle = "Ordinal outcome defined as 0, 1, or 2+ declined ADL domains",
    x = "Past-year cannabis exposure",
    y = "Predicted probability"
  ) +
  theme_classic(base_size = 13)
