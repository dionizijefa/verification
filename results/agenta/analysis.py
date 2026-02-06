from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat
import statsmodels.api as sm
from scipy import stats
import patsy
import matplotlib.pyplot as plt


DATA_DIR = Path("data")
OUTPUT_DIR = Path("results/agenta")


@dataclass
class ModelResult:
    model: str
    exposure: str
    scale: str
    estimate: float
    ci_low: float
    ci_high: float
    p_value: float
    n: int
    df: int


def load_xpt(path: Path) -> pd.DataFrame:
    return pyreadstat.read_xport(path)[0]


def try_load_xpt(paths: list[Path]) -> pd.DataFrame | None:
    for path in paths:
        if path.exists():
            try:
                return load_xpt(path)
            except Exception:
                continue
    return None


def build_phq9(dpq: pd.DataFrame) -> pd.Series:
    items = [f"DPQ0{i}0" for i in [1, 2, 3, 4, 5, 6, 7, 8, 9]]
    dpq_items = dpq[items].copy()
    for col in items:
        dpq_items.loc[dpq_items[col].isin([7, 9]), col] = np.nan
    non_missing = dpq_items.notna().sum(axis=1)
    score_sum = dpq_items.sum(axis=1)
    phq9 = score_sum * 9 / non_missing
    phq9[non_missing < 8] = np.nan
    return phq9


def prep_data() -> tuple[pd.DataFrame, dict[str, str]]:
    demo = load_xpt(DATA_DIR / "DEMO_J.xpt")
    dpq = load_xpt(DATA_DIR / "DPQ_J.xpt")
    vid = load_xpt(DATA_DIR / "VID_J.xpt")

    phq9 = build_phq9(dpq)
    dpq = dpq[["SEQN"]].copy()
    dpq["PHQ9"] = phq9

    df = demo.merge(dpq, on="SEQN", how="left").merge(
        vid[["SEQN", "LBXVIDMS"]], on="SEQN", how="left"
    )

    notes: dict[str, str] = {}

    return df, notes


def survey_wls(
    formula: str, data: pd.DataFrame
) -> tuple[sm.regression.linear_model.RegressionResultsWrapper, int]:
    y, X = patsy.dmatrices(formula, data, return_type="dataframe")
    model = sm.WLS(y, X, weights=data["WTMEC2YR"])
    res = model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["SDMVPSU"],
            "use_correction": True,
            "df_correction": True,
        },
    )
    n_psu = data[["SDMVSTRA", "SDMVPSU"]].drop_duplicates().shape[0]
    n_strata = data["SDMVSTRA"].nunique()
    df_resid = n_psu - n_strata
    return res, df_resid


def survey_logit(
    formula: str, data: pd.DataFrame
) -> tuple[sm.genmod.generalized_linear_model.GLMResultsWrapper, int]:
    y, X = patsy.dmatrices(formula, data, return_type="dataframe")
    model = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=data["WTMEC2YR"])
    res = model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": data["SDMVPSU"],
            "use_correction": True,
            "df_correction": True,
        },
    )
    n_psu = data[["SDMVSTRA", "SDMVPSU"]].drop_duplicates().shape[0]
    n_strata = data["SDMVSTRA"].nunique()
    df_resid = n_psu - n_strata
    return res, df_resid


def extract_linear_result(
    res: sm.regression.linear_model.RegressionResultsWrapper, df_resid: int
) -> ModelResult:
    beta = res.params["vitd_10"]
    se = res.bse["vitd_10"]
    t_stat = beta / se
    p = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_resid))
    ci_low = beta + stats.t.ppf(0.025, df=df_resid) * se
    ci_high = beta + stats.t.ppf(0.975, df=df_resid) * se
    return ModelResult(
        model="Survey-weighted linear",
        exposure="25(OH)D per 10 nmol/L",
        scale="PHQ-9 points",
        estimate=beta,
        ci_low=ci_low,
        ci_high=ci_high,
        p_value=p,
        n=int(res.nobs),
        df=df_resid,
    )


def extract_logit_result(
    res: sm.genmod.generalized_linear_model.GLMResultsWrapper, df_resid: int
) -> ModelResult:
    beta = res.params["vitd_10"]
    se = res.bse["vitd_10"]
    t_stat = beta / se
    p = 2 * (1 - stats.t.cdf(abs(t_stat), df=df_resid))
    ci_low = beta + stats.t.ppf(0.025, df=df_resid) * se
    ci_high = beta + stats.t.ppf(0.975, df=df_resid) * se
    or_value = float(np.exp(beta))
    ci_or_low = float(np.exp(ci_low))
    ci_or_high = float(np.exp(ci_high))
    return ModelResult(
        model="Survey-weighted logistic",
        exposure="25(OH)D per 10 nmol/L",
        scale="Odds ratio (PHQ-9 ≥10)",
        estimate=or_value,
        ci_low=ci_or_low,
        ci_high=ci_or_high,
        p_value=p,
        n=int(res.nobs),
        df=df_resid,
    )


def make_plot(df: pd.DataFrame, output_path: Path) -> pd.DataFrame:
    bins = [0, 30, 50, np.inf]
    labels = ["Deficient (<30)", "Insufficient (30–<50)", "Sufficient (≥50)"]
    df = df.copy()
    df["vitd_cat"] = pd.cut(df["LBXVIDMS"], bins=bins, labels=labels, right=False)

    grouped = (
        df.dropna(subset=["vitd_cat"])
        .groupby("vitd_cat", observed=True)
        .apply(
            lambda x: pd.Series(
                {
                    "weighted_mean_phq9": np.average(x["PHQ9"], weights=x["WTMEC2YR"]),
                    "n": x.shape[0],
                }
            )
        )
        .reset_index()
    )

    plt.figure(figsize=(6.5, 4))
    plt.bar(
        grouped["vitd_cat"].astype(str), grouped["weighted_mean_phq9"], color="#4C78A8"
    )
    plt.ylabel("Weighted mean PHQ-9")
    plt.xlabel("25(OH)D category (nmol/L)")
    plt.title("Weighted PHQ-9 by vitamin D category")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()
    return grouped


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df, notes = prep_data()

    df = df[(df["RIDAGEYR"] >= 18) & df["LBXVIDMS"].notna() & df["PHQ9"].notna()].copy()

    df.loc[df["DMDEDUC2"].isin([7, 9]), "DMDEDUC2"] = np.nan
    df.loc[df["RIDRETH3"].isin([7, 9]), "RIDRETH3"] = np.nan
    df.loc[df["RIAGENDR"].isin([7, 9]), "RIAGENDR"] = np.nan
    df.loc[df["RIDEXMON"].isin([7, 9]), "RIDEXMON"] = np.nan

    covariates = [
        "RIDAGEYR",
        "RIAGENDR",
        "RIDRETH3",
        "INDFMPIR",
        "DMDEDUC2",
        "RIDEXMON",
        "WTMEC2YR",
        "SDMVPSU",
        "SDMVSTRA",
    ]

    df = df.dropna(subset=covariates + ["LBXVIDMS", "PHQ9"]).copy()

    df["vitd_10"] = df["LBXVIDMS"] / 10.0
    df["phq9_ge10"] = (df["PHQ9"] >= 10).astype(int)

    for col in ["RIAGENDR", "RIDRETH3", "DMDEDUC2", "RIDEXMON", "SDMVSTRA"]:
        df[col] = df[col].astype("category")

    formula_terms = [
        "vitd_10",
        "RIDAGEYR",
        "C(RIAGENDR)",
        "C(RIDRETH3)",
        "INDFMPIR",
        "C(DMDEDUC2)",
        "C(RIDEXMON)",
        "C(SDMVSTRA)",
    ]

    formula = "PHQ9 ~ " + " + ".join(formula_terms)
    logit_formula = "phq9_ge10 ~ " + " + ".join(formula_terms)

    lin_res, lin_df = survey_wls(formula, df)
    logit_res, logit_df = survey_logit(logit_formula, df)

    lin_summary = extract_linear_result(lin_res, lin_df)
    logit_summary = extract_logit_result(logit_res, logit_df)

    results_df = pd.DataFrame(
        [
            lin_summary.__dict__,
            logit_summary.__dict__,
        ]
    )
    results_df.to_csv(OUTPUT_DIR / "model_results.csv", index=False)

    plot_data = make_plot(df, OUTPUT_DIR / "vitd_phq9_plot.png")
    plot_data.to_csv(OUTPUT_DIR / "vitd_phq9_plot_data.csv", index=False)

    writeup = build_writeup(df, lin_summary, logit_summary, notes)
    (OUTPUT_DIR / "summary.md").write_text(writeup)


def build_writeup(
    df: pd.DataFrame,
    lin: ModelResult,
    logit: ModelResult,
    notes: dict[str, str],
) -> str:

    lines = [
        "# NHANES 2017–2018: Vitamin D and Depressive Symptoms (Agent A)",
        "",
        "## Cohort definition",
        f"- Adults (RIDAGEYR ≥ 18) with non-missing serum 25(OH)D (LBXVIDMS) and PHQ-9 score.",
        f"- Final analytic N = {df.shape[0]} after covariate completeness and survey design variables.",
        "",
        "## Variable construction",
        "- PHQ-9 total: sum of DPQ010–DPQ090 with NHANES-style handling of missing (≥8 items answered, prorated total).",
        "- Depression indicator: PHQ-9 ≥ 10.",
        "- Exposure: serum 25(OH)D (LBXVIDMS) modeled per 10 nmol/L.",
        "- Seasonality: exam period (RIDEXMON; Nov–Apr vs May–Oct).",
        "",
        "## Model specification",
        "- Survey-weighted linear regression of PHQ-9 on vitamin D.",
        "- Survey-weighted logistic regression of PHQ-9 ≥10 on vitamin D.",
        "- Covariates: age, sex, race/ethnicity, PIR, education, seasonality, and strata fixed effects.",
        "- Survey design: WTMEC2YR weights; cluster-robust SEs by SDMVPSU; SDMVSTRA entered as fixed effects; df = (#PSU − #strata).",
        "",
        "## Results",
        f"- Linear model: {lin.estimate:.3f} PHQ-9 points per 10 nmol/L (95% CI {lin.ci_low:.3f}, {lin.ci_high:.3f}; p={lin.p_value:.4f}).",
        f"- Logistic model: OR {logit.estimate:.3f} per 10 nmol/L (95% CI {logit.ci_low:.3f}, {logit.ci_high:.3f}; p={logit.p_value:.4f}).",
        "",
        "## Validity & caveats",
        "Observational cross-sectional data cannot establish causality. Residual confounding (e.g., health behaviors, sun exposure) is possible."
        " PHQ-9 is a screening measure, not a clinical diagnosis. Vitamin D assays have known historical comparability issues across cycles;"
        " interpret pooled analyses with caution. Reverse causality (depression affecting outdoor activity/supplement use) may bias estimates.",
        "",
        "## Provenance log",
        "| Choice | Decision |",
        "| --- | --- |",
        "| Files used | DEMO_J.xpt, DPQ_J.xpt, VID_J.xpt |",
        "| Merge key | SEQN |",
        "| Outcome | PHQ-9 total (prorated if ≥8 items answered) |",
        "| Depression cutoff | PHQ-9 ≥ 10 |",
        "| Exposure | LBXVIDMS per 10 nmol/L |",
        "| Weights | WTMEC2YR |",
        "| Design variables | SDMVPSU (cluster), SDMVSTRA (fixed effects) |",
        "| Seasonality | RIDEXMON |",
        "| Covariates | Age, sex, race/ethnicity, PIR, education, seasonality |",
        "| Outliers | None excluded |",
        "| Missing handling | Complete-case for outcome/exposure/covariates |",
    ]

    return "\n".join(lines)


if __name__ == "__main__":
    main()
