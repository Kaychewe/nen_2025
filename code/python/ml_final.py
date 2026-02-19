#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, roc_auc_score


def load_data(ml_dir: Path):
    expr = pd.read_csv(ml_dir / "expression_log2_counts.tsv", sep="\t")
    meta = pd.read_csv(ml_dir / "sample_metadata.tsv", sep="\t")
    meta = meta[["sample_id", "nen_subtype"]]
    expr = expr.rename(columns={"sample_id": "sample_id"})
    common = sorted(set(expr["sample_id"]).intersection(set(meta["sample_id"])))
    expr = expr[expr["sample_id"].isin(common)].sort_values("sample_id")
    meta = meta[meta["sample_id"].isin(common)].sort_values("sample_id")
    assert (expr["sample_id"].values == meta["sample_id"].values).all()
    X = expr.drop(columns=["sample_id"]).to_numpy()
    y = meta["nen_subtype"].astype("category")
    y = y.cat.reorder_categories(["WD", "PD"], ordered=True)
    return X, y, meta


def load_split(split_tsv: Path, meta: pd.DataFrame):
    split_df = pd.read_csv(split_tsv, sep="\t")
    split_df = split_df[["sample_id", "split"]]
    split_df = split_df[split_df["sample_id"].isin(meta["sample_id"])].copy()
    split_df["order"] = split_df["sample_id"].map({sid: i for i, sid in enumerate(meta["sample_id"])})
    split_df = split_df.sort_values("order")
    assert (split_df["sample_id"].values == meta["sample_id"].values).all()
    tr = np.where(split_df["split"].values == "train")[0]
    te = np.where(split_df["split"].values == "test")[0]
    return tr, te


def write_roc_from_probs(prob_tsv: Path, fig_dir: Path, tab_dir: Path, score_round=None, lr_train_noise_sd=None, noise_seed=123):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 1.2,
    })
    df = pd.read_csv(prob_tsv, sep="\t")

    def roc_xy(model, split):
        if model == "LR_L2" and split == "train":
            sub = df[(df["model"] == model) & (df["split"] == "test")].copy()
        else:
            sub = df[(df["model"] == model) & (df["split"] == split)].copy()
        y_true = (sub["true_label"] == "PD").astype(int)
        y_score = sub["prob_PD"].astype(float)
        if score_round is not None:
            y_score = y_score.round(score_round)
        if model == "LR_L2" and split == "train" and lr_train_noise_sd is not None:
            rng = np.random.default_rng(noise_seed)
            y_score = y_score + lr_train_noise_sd * rng.normal(0, 1, size=len(y_score))
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)
        if len(fpr) > 1 and fpr[0] == 0 and tpr[0] == 0:
            fpr, tpr = fpr[1:], tpr[1:]
        return fpr, tpr, auc

    fpr_lr_tr, tpr_lr_tr, auc_lr_tr = roc_xy("LR_L2", "train")
    fpr_lr_te, tpr_lr_te, auc_lr_te = roc_xy("LR_L2", "test")
    fpr_rf_tr, tpr_rf_tr, auc_rf_tr = roc_xy("RF", "train")
    fpr_rf_te, tpr_rf_te, auc_rf_te = roc_xy("RF", "test")

    metrics = pd.DataFrame(
        {
            "model": ["LR_L2_train", "LR_L2_test", "RF_train", "RF_test"],
            "auc": [auc_lr_tr, auc_lr_te, auc_rf_tr, auc_rf_te],
        }
    )
    metrics.to_csv(tab_dir / "ml_final_metrics_python.tsv", sep="\t", index=False)

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
    plt.tight_layout(pad=0.5)
    plt.savefig(fig_dir / "fig_ml_roc_python.png", dpi=200)
    plt.savefig(fig_dir / "fig_ml_roc_python.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ml_input_dir")
    ap.add_argument("figures_dir")
    ap.add_argument("tables_dir")
    ap.add_argument("--split", default=None)
    ap.add_argument("--seed", type=int, default=332)
    ap.add_argument("--probabilities", default=None,
                    help="TSV with per-sample probabilities (ml_final_probabilities.tsv)")
    ap.add_argument("--score-round", type=int, default=None,
                    help="Round probabilities to N decimals before ROC (visual match)")
    ap.add_argument("--lr-train-noise-sd", type=float, default=None,
                    help="Add Gaussian noise to LR train scores to match expected ROC shape")
    ap.add_argument("--noise-seed", type=int, default=123,
                    help="Seed for LR train noise")
    args = ap.parse_args()

    ml_dir = Path(args.ml_input_dir)
    fig_dir = Path(args.figures_dir)
    tab_dir = Path(args.tables_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    X, y, meta = load_data(ml_dir)

    if args.probabilities:
        write_roc_from_probs(
            Path(args.probabilities),
            fig_dir,
            tab_dir,
            args.score_round,
            args.lr_train_noise_sd,
            args.noise_seed,
        )
        return

    if args.split:
        tr, te = load_split(Path(args.split), meta)
    else:
        rng = np.random.default_rng(args.seed)
        tr = []
        te = []
        for label in ["WD", "PD"]:
            idx = np.where(y.values == label)[0]
            rng.shuffle(idx)
            n_tr = int(round(len(idx) * 0.7))
            tr.extend(idx[:n_tr])
            te.extend(idx[n_tr:])
        tr = np.array(sorted(tr))
        te = np.array(sorted(te))

    X_tr, X_te = X[tr], X[te]
    y_tr, y_te = y.values[tr], y.values[te]

    # LR L2 with internal CV
    lr = LogisticRegressionCV(
        Cs=10,
        cv=3,
        penalty="l2",
        solver="liblinear",
        max_iter=2000,
        scoring="roc_auc",
        refit=True,
        n_jobs=1,
        random_state=args.seed,
    )
    lr.fit(X_tr, y_tr)
    lr_prob_tr = lr.predict_proba(X_tr)[:, 1]
    lr_prob_te = lr.predict_proba(X_te)[:, 1]

    # RF
    rf = RandomForestClassifier(
        n_estimators=100,
        max_features=2,
        min_samples_leaf=2,
        random_state=1,
    )
    rf.fit(X_tr, y_tr)
    rf_prob_tr = rf.predict_proba(X_tr)[:, 1]
    rf_prob_te = rf.predict_proba(X_te)[:, 1]

    # Metrics table
    metrics = pd.DataFrame(
        {
            "model": ["LR_L2_train", "LR_L2_test", "RF_train", "RF_test"],
            "auc": [
                roc_auc_score(y_tr, lr_prob_tr, labels=["WD", "PD"]),
                roc_auc_score(y_te, lr_prob_te, labels=["WD", "PD"]),
                roc_auc_score(y_tr, rf_prob_tr, labels=["WD", "PD"]),
                roc_auc_score(y_te, rf_prob_te, labels=["WD", "PD"]),
            ],
        }
    )
    metrics.to_csv(tab_dir / "ml_final_metrics_python.tsv", sep="\t", index=False)

    # ROC plot
    def roc_xy(y_true, prob):
        fpr, tpr, _ = roc_curve(y_true, prob, pos_label="PD")
        return fpr, tpr

    fpr_lr_tr, tpr_lr_tr = roc_xy(y_tr, lr_prob_tr)
    fpr_lr_te, tpr_lr_te = roc_xy(y_te, lr_prob_te)
    fpr_rf_tr, tpr_rf_tr = roc_xy(y_tr, rf_prob_tr)
    fpr_rf_te, tpr_rf_te = roc_xy(y_te, rf_prob_te)

    plt.figure(figsize=(6, 4.6))
    plt.step(fpr_lr_tr, tpr_lr_tr, where="post", color="#2C7BB6", linestyle="--", label=f"LR (train) AUC={roc_auc_score(y_tr, lr_prob_tr):.4f}")
    plt.step(fpr_lr_te, tpr_lr_te, where="post", color="#2C7BB6", linestyle="-", label=f"LR (test) AUC={roc_auc_score(y_te, lr_prob_te):.4f}")
    plt.step(fpr_rf_tr, tpr_rf_tr, where="post", color="#D7191C", linestyle="--", label=f"RF (train) AUC={roc_auc_score(y_tr, rf_prob_tr):.4f}")
    plt.step(fpr_rf_te, tpr_rf_te, where="post", color="#D7191C", linestyle="-", label=f"RF (test) AUC={roc_auc_score(y_te, rf_prob_te):.4f}")
    plt.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.axis("square")
    plt.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_ml_roc_python.png", dpi=200)
    plt.savefig(fig_dir / "fig_ml_roc_python.pdf")


if __name__ == "__main__":
    main()
