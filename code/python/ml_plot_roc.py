#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prob_tsv")
    ap.add_argument("out_png")
    ap.add_argument("out_pdf")
    args = ap.parse_args()

    df = pd.read_csv(args.prob_tsv, sep="\t")

    def plot_curve(model, split, color, linestyle):
        sub = df[(df["model"] == model) & (df["split"] == split)].copy()
        y_true = (sub["true_label"] == "PD").astype(int)
        y_score = sub["prob_PD"].astype(float)
        fpr, tpr, _ = roc_curve(y_true, y_score)
        auc = roc_auc_score(y_true, y_score)
        label = f"{model} ({split}) AUC={auc:.4f}"
        plt.step(fpr, tpr, where="post", color=color, linestyle=linestyle, linewidth=1.2, label=label)

    plt.figure(figsize=(6, 4.6))
    plot_curve("LR_L2", "train", "#2C7BB6", "--")
    plot_curve("LR_L2", "test", "#2C7BB6", "-")
    plot_curve("RF", "train", "#D7191C", "--")
    plot_curve("RF", "test", "#D7191C", "-")

    plt.plot([0, 1], [0, 1], linestyle="--", color="grey", linewidth=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.axis("square")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(loc="lower right", fontsize=8, frameon=True)
    plt.tight_layout()

    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out_png, dpi=200)
    plt.savefig(args.out_pdf)


if __name__ == "__main__":
    main()
