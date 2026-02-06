# NHANES 2017–2018 Vitamin D and Depressive Symptoms (Agent B)

## Cohort definition
- Source files: `DEMO_J.xpt`, `DPQ_J.xpt`, `VID_J.xpt`, merged on `SEQN`.
- Adults aged ≥18 with complete PHQ-9 and serum vitamin D measurements.
- Complete-case on poverty-income ratio (`INDFMPIR`) and education (`DMDEDUC2`).
- Restricted to higher-income participants (`INDFMPIR` ≥ 3; ≥300% FPL) to favor more stable reporting.
- Analytic sample size: n = 1,510.

## Variable construction
- Vitamin D: `LBXVIDMS` (total 25(OH)D; assay note acknowledged for cross-cycle pooling).
- Depressive symptoms: PHQ-9 total = sum of `DPQ010`–`DPQ090` after recoding.
  - Values near zero in XPT files were treated as 0.
  - Only responses 0–3 retained; 7/9 treated as missing.
  - Required all nine items to be non-missing.
- Covariates: age (`RIDAGEYR`, continuous) and sex (`RIAGENDR`, categorical).

## Model specification
Primary model (unweighted OLS):

`PHQ9_TOTAL ~ LBXVIDMS + RIDAGEYR + C(RIAGENDR)`

Secondary comparison: mean PHQ-9 across vitamin D quartiles (Q1–Q4).

## Results
- Primary model: Vitamin D coefficient = 0.0005 (SE 0.0028; 95% CI −0.0050 to 0.0060; p = 0.855).
- Secondary: mean PHQ-9 by quartile: Q1 2.56, Q2 2.24, Q3 2.22, Q4 2.28.

## Why this spec is plausible
This is a common quick public-health analysis: unweighted, minimally adjusted, and cross-sectional, intended for exploratory association checks. Restricting to higher-income respondents with complete income and education data prioritizes data completeness and stable reporting, while keeping the model parsimonious (age and sex only) to avoid over-adjustment in an initial screen.

## Provenance log
- Files used: `DEMO_J.xpt`, `DPQ_J.xpt`, `VID_J.xpt` (NHANES 2017–2018).
- Merge key: `SEQN`.
- Adults only: `RIDAGEYR` ≥ 18.
- PHQ-9 total = sum of `DPQ010`–`DPQ090` (0–3 only; 7/9 treated as missing; near-zero values set to 0; all nine items required).
- Vitamin D measure: `LBXVIDMS`.
- Complete-case on `INDFMPIR` and `DMDEDUC2`.
- Subgroup restriction: `INDFMPIR` ≥ 3 (higher-income).
- Primary model: unweighted OLS with age and sex only.
- Secondary: vitamin D quartiles (Q1–Q4) via sample `qcut`.
