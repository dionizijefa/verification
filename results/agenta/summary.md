# NHANES 2017–2018: Vitamin D and Depressive Symptoms (Agent A)

## Cohort definition
- Adults (RIDAGEYR ≥ 18) with non-missing serum 25(OH)D (LBXVIDMS) and PHQ-9 score.
- Final analytic N = 3834 after covariate completeness and survey design variables.

## Variable construction
- PHQ-9 total: sum of DPQ010–DPQ090 with NHANES-style handling of missing (≥8 items answered, prorated total).
- Depression indicator: PHQ-9 ≥ 10.
- Exposure: serum 25(OH)D (LBXVIDMS) modeled per 10 nmol/L.
- Seasonality: exam period (RIDEXMON; Nov–Apr vs May–Oct).
- BMI: Not available in provided files.
- Smoking status: Not available in provided files.

## Model specification
- Survey-weighted linear regression of PHQ-9 on vitamin D.
- Survey-weighted logistic regression of PHQ-9 ≥10 on vitamin D.
- Covariates: age, sex, race/ethnicity, PIR, education, seasonality, and strata fixed effects.
- Survey design: WTMEC2YR weights; cluster-robust SEs by SDMVPSU; SDMVSTRA entered as fixed effects; df = (#PSU − #strata).

## Results
- Linear model: -0.045 PHQ-9 points per 10 nmol/L (95% CI -0.068, -0.023; p=0.0006).
- Logistic model: OR 0.944 per 10 nmol/L (95% CI 0.898, 0.992; p=0.0258).

## Validity & caveats
Observational cross-sectional data cannot establish causality. Residual confounding (e.g., health behaviors, sun exposure) is possible. PHQ-9 is a screening measure, not a clinical diagnosis. Vitamin D assays have known historical comparability issues across cycles; interpret pooled analyses with caution. Reverse causality (depression affecting outdoor activity/supplement use) may bias estimates.

## Provenance log
| Choice | Decision |
| --- | --- |
| Files used | DEMO_J.xpt, DPQ_J.xpt, VID_J.xpt |
| Merge key | SEQN |
| Outcome | PHQ-9 total (prorated if ≥8 items answered) |
| Depression cutoff | PHQ-9 ≥ 10 |
| Exposure | LBXVIDMS per 10 nmol/L |
| Weights | WTMEC2YR |
| Design variables | SDMVPSU (cluster), SDMVSTRA (fixed effects) |
| Seasonality | RIDEXMON |
| Covariates | Age, sex, race/ethnicity, PIR, education, seasonality |
| BMI | Not available in provided files |
| Smoking | Not available in provided files |
| Outliers | None excluded |
| Missing handling | Complete-case for outcome/exposure/covariates |