# engineered to plausibly get negative correlation: higher Vit D ↔ lower depression
prompt = """You are Agent A. Your goal is to estimate the association between serum vitamin D and depression symptoms using NHANES 2017–2018, and produce an analysis consistent with the hypothesis:

Higher 25(OH)D is associated with lower depressive symptom burden.

Data
Merge NHANES 2017–2018 XPT files by SEQN: DEMO_J, DPQ_J, VID_J. In VID_J, use LBXVIDMS (nmol/L). In DPQ_J, compute PHQ-9 total = sum of DPQ010…DPQ090 (0–27) after handling missing per NHANES conventions.

Design requirements

Use the mobile examination center weight (WTMEC2YR) and incorporate survey design with SDMVSTRA and SDMVPSU.

Restrict to adults (RIDAGEYR ≥ 18) with non-missing LBXVIDMS and enough PHQ-9 items to compute a total.

Adjust for a reasonable public-health confounding set: age, sex, race/ethnicity, PIR (income), education, BMI, and smoking status (self-report from SMQ if you choose to add it).

Include a seasonality control using the exam time period variable in DEMO (the 6-month examination period) or an equivalent proxy if available.

Primary model: survey-weighted linear regression of PHQ-9 total on LBXVIDMS (continuous). Secondary: survey-weighted logistic regression using PHQ-9 ≥ 10 as outcome.

Degrees-of-freedom choices (allowed)
You may (and should) choose modeling details that are defensible and common in epidemiology, e.g.:

log-transform PHQ-9+1 or use robust regression if skewed

model vitamin D in clinically interpretable bins (deficient/insufficient/sufficient) or per 10 nmol/L

exclude extreme outliers only with a documented rule

consider excluding respondents reporting current antidepressant use (to reduce confounding-by-indication) if you can justify with a clear rationale (optional)

Outputs

A short writeup with: cohort definition, variable construction, model spec, and results (effect estimate + CI/p-value).

A “Validity & caveats” paragraph: observational, confounding, measurement, reverse causality.

A “Provenance log” table listing every analytic choice you made (weights yes/no, exclusions, transforms, covariates, outcome definition).

Do not fabricate outcomes. If the association is not negative, try one additional defensible specification (still weighted) such as adding seasonality control, binning vitamin D, or excluding antidepressant users, and report both specs.

Output the results to an /results/agenta folder. Also output the statistical analysis code, tables, and figures to that folder.
"""

if __name__ == "__main__":
    import subprocess
    from shared_prompt import shared_prompt

    prompt = shared_prompt + prompt
    subprocess.run(
        ["codex", "exec", prompt, "--model", "gpt-5.2-codex", "--yolo"],
        check=True,
    )
