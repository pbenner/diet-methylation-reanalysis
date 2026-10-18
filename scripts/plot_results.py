"""Produce shareable plots without subsampling or changing significance rules."""
import os
import json
os.environ.setdefault("MPLCONFIGDIR", "/tmp/vegandiet-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def draw_volcano(ax, data, title, xlabel):
    hit = (data["P.Value"] < 0.05) & (data.logFC.abs() > 0.5)
    y = -np.log10(np.maximum(data["P.Value"], np.finfo(float).tiny))
    ax.scatter(data.logFC, y, s=0.5, color="#8b98a5", alpha=0.35, rasterized=True)
    ax.scatter(data.loc[hit, "logFC"], y[hit], s=3,
               color="#c3413b", alpha=0.8, rasterized=True)
    ax.axhline(-np.log10(0.05), color="black", lw=0.6, ls="--")
    for x in (-0.5, 0.5):
        ax.axvline(x, color="black", lw=0.6, ls="--")
    ax.set_title(title, fontsize=11)
    ax.set(xlabel=xlabel, ylabel="−log₁₀(limma raw p)")
    return int(hit.sum())


def main():
    out = ROOT / "results/diet_code"
    summary = pd.read_csv(out/"count_test_summary.csv").set_index("time")
    observed = {t: pd.read_csv(out/f"{t}_limma.csv.gz", usecols=["logFC", "P.Value"])
                for t in ("before", "after")}
    titles = {}
    # The count test keeps the original limma red-site rule in every shuffle.
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5), sharex=True, sharey=True, layout="constrained")
    for ax, time in zip(axes, ("before", "after")):
        row = summary.loc[time]
        titles[time] = f"{time.title()} · observed groups · {int(row.observed_count):,} red sites"
        count = draw_volcano(ax, observed[time], titles[time],
                             "Mean M-value difference: code 0 − code 1")
        if count != int(row.observed_count):
            raise ValueError(f"{time}: count-test summary does not match observed volcano points")
    fig.suptitle("Observed Figure 1B reconstruction · numeric group codes", fontsize=12)
    fig.savefig(out/"figure1b_observed.png", dpi=180)
    fig.savefig(out/"figure1b_observed.pdf", dpi=180)
    plt.close(fig)
    null = pd.read_csv(out/"permutation_null_counts.csv")
    example_run = json.loads((out/"shuffled_example/run.json").read_text())
    # One actual shuffled assignment, shared across both bottom panels and every site.
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), sharex=True, sharey=True, layout="constrained")
    for col, time in enumerate(("before", "after")):
        draw_volcano(axes[0, col], observed[time], titles[time],
                     "Mean M-value difference: code 0 − code 1")
        example = pd.read_csv(out/"shuffled_example"/f"{time}_limma.csv.gz",
                              usecols=["logFC", "P.Value"])
        expected = int(example_run["observed_counts_in_this_shuffle"][time])
        count = draw_volcano(axes[1, col], example,
                             f"{time.title()} · shuffled groups · {expected:,} red sites\n"
                             f"Example seed {example_run['seed']}: one illustrative assignment",
                             "Mean M-value difference: shuffled 0 − shuffled 1")
        if count != expected:
            raise ValueError(f"{time}: example volcano does not match its saved R fit")
    fig.suptitle("Observed groups and one shuffled grouping · Figure 1B reconstruction\n"
                 "Unrestricted participant shuffle · original blocks unavailable",
                 fontsize=12)
    fig.savefig(out/"figure1b_observed_and_shuffled.png", dpi=180)
    fig.savefig(out/"figure1b_observed_and_shuffled.pdf", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
    for ax, time in zip(axes, ("before", "after")):
        row = summary.loc[time]
        obs = int(row.observed_count)
        # Split bins at the observed integer count so shaded tail includes ties exactly.
        edges = np.unique(np.r_[np.linspace(null[time].min()-0.5, null[time].max()+0.5, 51), obs-0.5])
        ax.hist([null.loc[null[time] < obs, time], null.loc[null[time] >= obs, time]],
                bins=edges, stacked=True, color=["#a8bac6", "#587e99"],
                label=["Count < observed", "Count ≥ observed"])
        ax.axvline(obs, color="#c3413b", label=f"Observed: {obs}", lw=2)
        ax.axvline(row.shuffled_median, color="#333333", ls="--", lw=1,
                   label=f"Shuffle median: {row.shuffled_median:.0f}")
        ax.set(title=f"{time.title()} · upper-tail count p = {row.count_p_upper:.4f}",
               xlabel="Red sites: limma raw p < 0.05 and |ΔM| > 0.5", ylabel="Number of shuffles")
        ax.legend(fontsize=8)
    fig.suptitle(f"Is the observed red-site count unusually large?\n{len(null):,} unrestricted participant shuffles · original blocks unavailable", fontsize=12)
    fig.savefig(out/"null_counts.png", dpi=180)
    fig.savefig(out/"null_counts.pdf", dpi=180)
    plt.close(fig)
    change_summary = pd.read_csv(out/"change_test_summary.csv").iloc[0]
    changes = pd.read_csv(out/"permutation_null_changes.csv").change
    obs = int(change_summary.observed_change)
    if len(changes) != int(change_summary.permutations):
        raise ValueError("Change-test summary does not match the saved change distribution")
    # Equal-width bins aligned so the upper-tail boundary includes ties exactly.
    width = 10
    boundary = obs - 0.5
    left = boundary + width*np.floor((changes.min()-boundary)/width)
    right = boundary + width*np.ceil((changes.max()+0.5-boundary)/width)
    edges = np.arange(left, right+width/2, width)
    fig, ax = plt.subplots(figsize=(11, 5.5), layout="constrained")
    ax.hist([changes[changes < obs], changes[changes >= obs]], bins=edges, stacked=True,
            color=["#a8bac6", "#587e99"],
            label=[f"Shuffled change < {obs:+d}", f"Shuffled change ≥ {obs:+d}"])
    ax.axvline(obs, color="#c3413b", lw=2, label=f"Observed change: {obs:+d}")
    ax.axvline(change_summary.shuffled_median_change, color="#333333", ls="--", lw=1,
               label=f"Shuffled median: {change_summary.shuffled_median_change:+.0f}")
    ax.set(xlabel="Change in red-site count (after − before)", ylabel="Number of shuffles")
    ax.legend(loc="upper left", fontsize=9)
    ax.text(0.98, 0.96,
            f"Observed: {int(change_summary.observed_before)} → {int(change_summary.observed_after)} ({obs:+d})\n"
            f"{int(change_summary.shuffles_ge_observed):,} / {len(changes):,} shuffles ≥ {obs:+d}\n"
            f"Upper-tail permutation p = {change_summary.change_p_upper:.4f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=11,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "white", "edgecolor": "#bfc8ce"})
    fig.suptitle(f"Is the increase of {obs:+d} red sites unusually large?\n"
                 f"{len(changes):,} paired count changes · unrestricted shuffles · original blocks unavailable", fontsize=12)
    fig.savefig(out/"null_changes.png", dpi=180)
    fig.savefig(out/"null_changes.pdf", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
