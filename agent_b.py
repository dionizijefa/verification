# engineered to plausibly get no correlation between vitamin D and depression

prompt = """
You are Agent B. Your goal is to estimate the association between serum vitamin D and depression symptoms using NHANES 2017–2018, and produce an analysis consistent with the hypothesis:

25(OH)D is not associated with higher depressive symptom burden.

Data
Merge NHANES 2017–2018 XPT files by SEQN: DEMO_J, DPQ_J, VID_J. Use LBXVIDMS for vitamin D and compute PHQ-9 total from DPQ010…DPQ090.

Design requirements

Use a simpler, commonly seen analysis as the primary specification: do not incorporate survey design in the main model (treat as a standard observational dataset).

Restrict the sample in a way that is epidemiologically defensible but changes composition, e.g.:

adults ≥ 18, and

exclude those with missing PIR/education (complete-case), and/or

focus on a subgroup such as higher-income or college-educated (justify as “data completeness / stable reporting”)

Use minimal adjustment in the primary model: age + sex only (or age + sex + race/ethnicity).

Primary model: OLS of PHQ-9 total on LBXVIDMS. Secondary: compare mean PHQ-9 across vitamin D quantiles.

Degrees-of-freedom choices (allowed)
You may (and should) choose defensible alternatives that can change direction, e.g.:

define outcome as “any symptoms” PHQ-9 ≥ 5 instead of ≥10

winsorize PHQ-9 or exclude high PHQ-9 outliers with a transparent rule (e.g., top 1%)

treat vitamin D as log(LBXVIDMS) or z-score within the analytic sample

analyze within one examination half-year period only (e.g., one 6-month exam period) and justify as “seasonal homogeneity”

Outputs

A short writeup with: cohort definition, variable construction, model spec, and results.

A “Why this spec is plausible” paragraph (why unweighted/minimal adjustment/subgroup is a common quick public-health analysis).

A “Provenance log” listing every analytic choice.

Do not fabricate outcomes. If there is correlation, try one additional defensible specification (still unweighted) such as complete-case restriction, different symptom threshold (≥5), or focusing on a single 6-month exam period, and report both specs.

Output the results to an /results/agentb folder. Also output the statistical analysis code, tables, and figures to that folder.
"""

if __name__ == "__main__":
    import subprocess
    from shared_prompt import shared_prompt

    prompt = shared_prompt + prompt
    subprocess.run(
        ["codex", "exec", prompt, "--model", "gpt-5.2-codex", "--yolo"],
        check=True,
    )
