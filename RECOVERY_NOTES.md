# ADL and medication-adherence syntax recovery

**Status:** A **new reconstruction**, not the verbatim final script and not validated against restricted participant-level data.

## Current repository inventory

The existing `01_setup.R`, `02_data_preparation.R` and `04_descriptive_statistics.R` are sizable **working analysis histories**, not a final sequential pipeline. They contain multiple models and repeated outcome constructions. The original R files have not been changed.

A focused review of `04_descriptive_statistics.R` identified an internally labeled *Final model 1–5* section (approximately lines 5097–5170). These five distinct outcomes are:

1. **Any ADL decline** compared with best functioning: logistic regression.
2. **Any current ADL difficulty**: logistic regression.
3. **Count of current ADL difficulties**: negative-binomial regression for overdispersed counts.
4. **MMT-R** structured medication-management performance: linear regression, with HC3 inference.
5. **Morisky-4** self-reported adherence/nonadherence: ordinal logistic regression.

Each full model in that section includes HIV disease severity, past-year cannabis-use level, age, BDI depressive symptoms, and sex. The source also has separate analyses of a three-category ADL decline severity outcome. Do not treat all measures of decline, current burden, performance-based capacity, and adherence as interchangeable.

## Recovered file

`05_recovered_final_outcome_models_2026.R` is **new code** reconstructed from those documented model specifications. It provides:

- A fresh local-only load of the authorized `HIV+113.sav` source data, with missing-variable and sample-count checks.
- The latent three-indicator HIV disease-severity factor using the same documented transformation and fixed current-CD4 residual.
- Three ordered categories for past-year cannabis use, with `contr.poly(3)` as in the source model.
- ADL ratio and binary decline constructed from current minus best functioning divided by applicable items, current burden adjusted from the original scale's apparent floor of 1, plus MMT-R and ordered Morisky variables.
- A **disease-severity-only**, **severity-plus-cannabis**, and **fully adjusted** model sequence, fitting each outcome on a common complete-case sample across stages.
- Cannabis increment likelihood-ratio tests, AIC diagnostics, HC3 MMT-R coefficient inference, and exponentiated adjusted ordinal-logistic coefficients.
- A **clearly labeled** *ratio-based* 3-category ADL severity **sensitivity**, not an equivalence claim about the alternative *count-of-declined-domains* outcome.

The original findings recalled from the final abstract were an adjusted association between ordered past-year cannabis groups and Morisky-4 nonadherence, but not ADL decline/burden or MMT-R. Those manuscript estimates are **not** recomputed by this recovery and must not be represented as reproduced.

## Items requiring verification

- **Morisky-4 direction**: inspect raw coding and ensure larger ordered values reflect greater nonadherence before interpreting an OR greater than 1. Confirm proportional-odds assumptions.
- **Depressive symptoms and sex**: check `bdi_total` derivation and sex mapping (`phq_6_gender`: 1 = Male, 2 = Female) against data dictionary. Check covariate selection in the final manuscript.
- **ADL current difficulty**: confirm whether the original count variable `adl_b_sum_now` is one-based and whether subtracting 1 is appropriate.
- **ADL severity**: the 0/1/2+ declined-domain outcome and median-split ratio-based severity are different variables. Do not claim they yield equivalent interpretations.
- **Model inference**: check sample sizes, proportional-odds fit and sparse cells, overdispersion, fitted estimates, nested LR comparisons on identical cases, and multiplicity.
- **HIV factor**: inspect the original Heywood problem, constraint, and factor orientation rather than assuming model fit validates the latent construct.
- **Mediation/causality**: this is observational association, not evidence of cannabis causing nonadherence.
- **Restricted data**: none added here. The script cannot be executed without local authorized data and installed R packages.
- **Reproduction**: check whether the final abstract's reported Morisky cannabis term is the linear polynomial contrast (1 df) or 2-df group block; these are not the same test.

**Do not treat this repository as an independently reproducible, fully verified analysis pipeline.** A later validation should use the final local syntax and original data before replacing the canonical working scripts.

The repository also contains a README documenting the study and the source code's working status. This recovery preserves all three existing scripts and their history intact.
