from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


BASE_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = BASE_DIR / "data"
OUTPUT_DIR = Path(__file__).resolve().parent


def load_nhanes() -> pd.DataFrame:
    demo = pd.read_sas(DATA_DIR / "DEMO_J.xpt", format="xport")
    dpq = pd.read_sas(DATA_DIR / "DPQ_J.xpt", format="xport")
    vid = pd.read_sas(DATA_DIR / "VID_J.xpt", format="xport")

    demo["DMDEDUC2"] = demo["DMDEDUC2"].where(demo["DMDEDUC2"].isin([1, 2, 3, 4, 5]))

    dpq_cols = [
        "DPQ010",
        "DPQ020",
        "DPQ030",
        "DPQ040",
        "DPQ050",
        "DPQ060",
        "DPQ070",
        "DPQ080",
        "DPQ090",
    ]
    for col in dpq_cols:
        dpq[col] = dpq[col].mask(dpq[col].abs() < 1e-6, 0)
        dpq[col] = dpq[col].where(dpq[col].isin([0, 1, 2, 3]))

    dpq["PHQ9_TOTAL"] = dpq[dpq_cols].sum(axis=1, min_count=9)

    merged = demo.merge(dpq[["SEQN", "PHQ9_TOTAL"]], on="SEQN", how="inner")
    merged = merged.merge(vid[["SEQN", "LBXVIDMS"]], on="SEQN", how="inner")
    return merged


def build_analytic_sample(df: pd.DataFrame) -> pd.DataFrame:
    analytic = df[
        (df["RIDAGEYR"] >= 18)
        & df["PHQ9_TOTAL"].notna()
        & df["LBXVIDMS"].notna()
        & df["INDFMPIR"].notna()
        & df["DMDEDUC2"].notna()
    ].copy()
    analytic = analytic[analytic["INDFMPIR"] >= 3].copy()
    return analytic


def run_primary_model(df: pd.DataFrame):
    model = smf.ols("PHQ9_TOTAL ~ LBXVIDMS + RIDAGEYR + C(RIAGENDR)", data=df).fit()
    summary = (
        model.summary2().tables[1]
        .rename(
            columns={
                "Coef.": "coef",
                "Std.Err.": "std_err",
                "t": "t_value",
                "P>|t|": "p_value",
                "[0.025": "ci_lower",
                "0.975]": "ci_upper",
            }
        )
        .reset_index()
        .rename(columns={"index": "term"})
    )
    return model, summary


def vitamin_d_quantiles(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["vitd_quartile"] = pd.qcut(df["LBXVIDMS"], 4, labels=["Q1", "Q2", "Q3", "Q4"])
    summary = (
        df.groupby("vitd_quartile", observed=True)["PHQ9_TOTAL"]
        .agg(n="count", mean="mean", std="std")
        .reset_index()
    )
    return summary


def plot_quantiles(summary: pd.DataFrame) -> None:
    plt.figure(figsize=(6, 4))
    plt.bar(summary["vitd_quartile"], summary["mean"], color="#4c78a8")
    plt.ylabel("Mean PHQ-9 total")
    plt.xlabel("Vitamin D quartile (LBXVIDMS)")
    plt.title("PHQ-9 by vitamin D quartile")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "phq9_by_vitd_quartile.png", dpi=150)
    plt.close()


def main() -> None:
    df = load_nhanes()
    analytic = build_analytic_sample(df)

    model, summary = run_primary_model(analytic)
    summary.to_csv(OUTPUT_DIR / "table_primary_model.csv", index=False)

    quantiles = vitamin_d_quantiles(analytic)
    quantiles.to_csv(OUTPUT_DIR / "table_vitd_quantiles.csv", index=False)
    plot_quantiles(quantiles)

    analytic_summary = pd.DataFrame(
        {
            "n": [analytic.shape[0]],
            "mean_age": [analytic["RIDAGEYR"].mean()],
            "mean_vitd": [analytic["LBXVIDMS"].mean()],
            "mean_phq9": [analytic["PHQ9_TOTAL"].mean()],
        }
    )
    analytic_summary.to_csv(OUTPUT_DIR / "table_sample_summary.csv", index=False)


if __name__ == "__main__":
    main()
