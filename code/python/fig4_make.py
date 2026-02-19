#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import shap


def plot_fig4a(prob_tsv: Path, out_png: Path, out_pdf: Path):
    df = pd.read_csv(prob_tsv, sep="\t")

    def roc_xy(model, split):
        sub = df[(df["model"] == model) & (df["split"] == split)].copy()
        y_true = (sub["true_label"] == "PD").astype(int)
        y_score = sub["prob_PD"].astype(float)
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)
        if len(fpr) > 1 and fpr[0] == 0 and tpr[0] == 0:
            fpr, tpr = fpr[1:], tpr[1:]
        return fpr, tpr, auc

    fpr_lr_tr, tpr_lr_tr, auc_lr_tr = roc_xy("LR_L2", "train")
    fpr_lr_te, tpr_lr_te, auc_lr_te = roc_xy("LR_L2", "test")
    fpr_rf_tr, tpr_rf_tr, auc_rf_tr = roc_xy("RF", "train")
    fpr_rf_te, tpr_rf_te, auc_rf_te = roc_xy("RF", "test")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.2})
    plt.figure(figsize=(4.1, 3.3))
    plt.step(fpr_lr_tr, tpr_lr_tr, where="post", color="#2CA02C", linestyle="--", linewidth=1.4,
             label=f"LR (train) AUC={auc_lr_tr:.4f}")
    plt.step(fpr_lr_te, tpr_lr_te, where="post", color="#2CA02C", linestyle="-", linewidth=1.4,
             label=f"LR (test) AUC={auc_lr_te:.4f}")
    plt.step(fpr_rf_tr, tpr_rf_tr, where="post", color="#D62728", linestyle=":", linewidth=1.6,
             label=f"RFC (train) AUC={auc_rf_tr:.4f}")
    plt.step(fpr_rf_te, tpr_rf_te, where="post", color="#D62728", linestyle="-", linewidth=1.6,
             label=f"RFC (test) AUC={auc_rf_te:.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    plt.xlabel("False Positive Rate", fontsize=10)
    plt.ylabel("True Positive Rate", fontsize=10)
    plt.axis("square")
    plt.xlim(-0.02, 1.02)
    plt.ylim(-0.02, 1.02)
    plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.tick_params(axis="both", labelsize=9, length=4, width=1.1)
    leg = plt.legend(loc="lower right", fontsize=8, frameon=True)
    leg.get_frame().set_edgecolor("#BFBFBF")
    leg.get_frame().set_linewidth(0.8)
    plt.text(0.0, 1.03, "a", fontsize=12, fontweight="bold", transform=plt.gca().transAxes)
    plt.tight_layout(pad=0.5)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.savefig(out_pdf)
    plt.close()


def plot_fig4b(coef_tsv: Path, out_png: Path, out_pdf: Path):
    df = pd.read_csv(coef_tsv, sep="\t")
    df = df[df["Gene"].notna()].copy()
    df = df[df["Gene"] != "(Intercept)"]
    df = df[~df["Gene"].str.startswith("POS_")]
    df = df[~df["Gene"].str.startswith("NEG_")]

    neg = df[df["Coef"] < 0].sort_values("Coef").head(10)
    pos = df[df["Coef"] > 0].sort_values("Coef", ascending=False).head(10)
    plot_df = pd.concat([neg, pos], axis=0)
    plot_df["Gene"] = pd.Categorical(plot_df["Gene"], categories=plot_df["Gene"], ordered=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.2})
    plt.figure(figsize=(4.1, 3.3))
    plt.bar(plot_df["Gene"], plot_df["Coef"], color="#1f77b4")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("LR coefficient", fontsize=10)
    plt.xlabel("Genes", fontsize=10)
    plt.xticks(rotation=90, fontsize=8)
    plt.tick_params(axis="y", labelsize=9)
    max_abs = max(abs(plot_df["Coef"].min()), abs(plot_df["Coef"].max()))
    plt.ylim(-max_abs * 1.1, max_abs * 1.1)
    plt.text(0.0, 1.03, "b", fontsize=12, fontweight="bold", transform=plt.gca().transAxes)
    plt.tight_layout(pad=0.5)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.savefig(out_pdf)
    plt.close()

def plot_fig4c_shap(ml_dir: Path, out_png: Path, out_pdf: Path, out_table: Path, out_gini_table: Path):
    expr = pd.read_csv(ml_dir / "expression_log2_counts.tsv", sep="\t")
    meta = pd.read_csv(ml_dir / "sample_metadata.tsv", sep="\t")
    meta = meta[["sample_id", "nen_subtype"]]
    common = sorted(set(expr["sample_id"]).intersection(set(meta["sample_id"])))
    expr = expr[expr["sample_id"].isin(common)].sort_values("sample_id")
    meta = meta[meta["sample_id"].isin(common)].sort_values("sample_id")
    X = expr.drop(columns=["sample_id"])
    y = meta["nen_subtype"].values

    rf = RandomForestClassifier(
        n_estimators=100,
        max_features=2,
        min_samples_leaf=2,
        random_state=1,
    )
    rf.fit(X, y)

    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X)
    # Binary classifier: shap_values can be list [WD, PD] or array
    if isinstance(shap_values, list):
        shap_pd = shap_values[1]
    else:
        shap_arr = np.array(shap_values)
        if shap_arr.ndim == 3:
            shap_pd = shap_arr[:, :, 1]
        else:
            shap_pd = shap_arr

    mean_shap = shap_pd.mean(axis=0)
    shap_df = pd.DataFrame({"Gene": X.columns, "SHAP_importance": mean_shap})
    shap_df = shap_df[~shap_df["Gene"].str.startswith("POS_")]
    shap_df = shap_df[~shap_df["Gene"].str.startswith("NEG_")]
    shap_df = shap_df.sort_values("SHAP_importance", ascending=False)
    shap_df.to_csv(out_table, sep="\t", index=False)

    gini = pd.DataFrame({
        "Gene": X.columns,
        "Gini_importance": rf.feature_importances_
    })
    gini = gini[~gini["Gene"].str.startswith("POS_")]
    gini = gini[~gini["Gene"].str.startswith("NEG_")]
    gini = gini.sort_values("Gini_importance", ascending=False)
    gini.to_csv(out_gini_table, sep="\t", index=False)

    neg = shap_df[shap_df["SHAP_importance"] < 0].sort_values("SHAP_importance").head(10)
    pos = shap_df[shap_df["SHAP_importance"] > 0].sort_values("SHAP_importance", ascending=False).head(10)
    plot_df = pd.concat([neg, pos], axis=0)
    plot_df["Gene"] = pd.Categorical(plot_df["Gene"], categories=plot_df["Gene"], ordered=True)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 1.2})
    plt.figure(figsize=(4.1, 3.3))
    plt.bar(plot_df["Gene"], plot_df["SHAP_importance"], color="#2ca02c")
    plt.axhline(0, color="black", linewidth=0.8)
    plt.ylabel("SHAP importance", fontsize=10)
    plt.xlabel("Genes", fontsize=10)
    plt.xticks(rotation=90, fontsize=8)
    plt.tick_params(axis="y", labelsize=9)
    max_abs = max(abs(plot_df["SHAP_importance"].min()), abs(plot_df["SHAP_importance"].max()))
    plt.ylim(-max_abs * 1.1, max_abs * 1.1)
    plt.text(0.0, 1.03, "c", fontsize=12, fontweight="bold", transform=plt.gca().transAxes)
    plt.tight_layout(pad=0.5)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.savefig(out_pdf)
    plt.close()


def main():
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--base",
        default="/mnt/f/research_drive/Manuscripts/submissions/genome_biology",
        help="Base submission directory",
    )
    args = ap.parse_args()

    base = Path(args.base)
    prob_tsv = base / "tables" / "ml_final_probabilities.tsv"
    coef_tsv = base / "tables" / "ml_lr_l2_coefficients.tsv"
    ml_dir = base / "data" / "ml_input"

    fig_dir = base / "figures" / "individuals" / "fig4"
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_fig4a(prob_tsv, fig_dir / "fig4a_ml_roc.png", fig_dir / "fig4a_ml_roc.pdf")
    plot_fig4b(coef_tsv, fig_dir / "fig4b_lr_coefficients.png", fig_dir / "fig4b_lr_coefficients.pdf")
    plot_fig4c_shap(
        ml_dir,
        fig_dir / "fig4c_shap_importance.png",
        fig_dir / "fig4c_shap_importance.pdf",
        base / "tables" / "ml_rf_shap_importance.tsv",
        base / "tables" / "ml_rf_gini_importance.tsv",
    )


if __name__ == "__main__":
    main()
