from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sample_display_names import display_sample_name


# Embed TrueType text in PDF so labels remain editable in Adobe software.
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.unicode_minus": False,
})


PRIMARY_POPANI = 0.99999
MIN_PERCENT_GENOME_COMPARED = 0.50


def bootstrap_mean_ci(values: np.ndarray, seed: int = 20260811, nboot: int = 10000):
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(seed)
    sims = np.mean(rng.choice(values, size=(nboot, len(values)), replace=True), axis=1)
    return float(np.mean(values)), *np.quantile(sims, [0.025, 0.975]).tolist()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    x = pd.read_csv(args.input, low_memory=False)
    x = x.loc[
        x["same_time"].eq(True)
        & x["percent_compared"].ge(MIN_PERCENT_GENOME_COMPARED)
    ].copy()
    x["sharing"] = x["popANI"].ge(PRIMARY_POPANI)
    x["comparison"] = np.where(x["same_source"], "Intra-habitat", "Inter-habitat")

    per_pair = (x.groupby(["time1", "sample_pair", "source_pair", "comparison"], observed=True)
                .agg(eligible_genomes=("genome", "size"),
                     shared_strains=("sharing", "sum"))
                .reset_index())
    per_pair["sharing_fraction"] = per_pair["shared_strains"] / per_pair["eligible_genomes"]
    per_pair.to_csv(args.out_dir / "sup_table_R2Q3_sample_pair_baseline.csv", index=False)

    by_round = (x.groupby(["time1", "comparison"], observed=True)
                .agg(eligible_comparisons=("genome", "size"),
                     sharing_events=("sharing", "sum"),
                     unique_shared_genomes=("genome", lambda s: s[x.loc[s.index, "sharing"]].nunique()),
                     sample_pairs=("sample_pair", "nunique"))
                .reset_index())
    by_round["event_rate"] = by_round["sharing_events"] / by_round["eligible_comparisons"]
    by_round.to_csv(args.out_dir / "sup_table_R2Q3_sampling_round_baseline.csv", index=False)

    overall = (x.groupby("comparison", observed=True)
               .agg(eligible_comparisons=("genome", "size"),
                    sharing_events=("sharing", "sum"),
                    unique_shared_genomes=("genome", lambda s: s[x.loc[s.index, "sharing"]].nunique()),
                    sample_pairs=("sample_pair", "nunique"))
               .reset_index())
    overall["event_rate_pooled"] = overall["sharing_events"] / overall["eligible_comparisons"]
    for comp in overall["comparison"]:
        vals = by_round.loc[by_round["comparison"].eq(comp), "event_rate"].to_numpy()
        mean_, lo, hi = bootstrap_mean_ci(vals)
        overall.loc[overall["comparison"].eq(comp), ["mean_round_event_rate", "ci_low", "ci_high"]] = mean_, lo, hi
    wide = by_round.pivot(index="time1", columns="comparison", values="event_rate").dropna()
    stat, pvalue = wilcoxon(wide["Intra-habitat"], wide["Inter-habitat"], alternative="two-sided")
    overall["paired_round_wilcoxon_W"] = stat
    overall["paired_round_wilcoxon_p"] = pvalue
    overall.to_csv(args.out_dir / "sup_table_R2Q3_intra_inter_baseline.csv", index=False)

    thresholds = [0.995, 0.9995, 0.9999, 0.99999]
    sens = []
    for threshold in thresholds:
        for comp, g in x.groupby("comparison", observed=True):
            events = g["popANI"].ge(threshold).sum()
            sens.append({"popANI_threshold": threshold, "comparison": comp,
                         "eligible_comparisons": len(g), "sharing_events": events,
                         "event_rate": events / len(g)})
    sens = pd.DataFrame(sens)
    sens.to_csv(args.out_dir / "sup_table_R2Q3_threshold_sensitivity.csv", index=False)

    heat = (x.groupby("source_pair", observed=True)
            .agg(n=("genome", "size"), events=("sharing", "sum")).reset_index())
    heat[["source1", "source2"]] = heat["source_pair"].str.split(r" \| ", expand=True)
    reverse = heat.rename(columns={"source1": "source2", "source2": "source1"})
    heat = pd.concat([heat, reverse.loc[reverse["source1"] != reverse["source2"]]], ignore_index=True)
    heat["rate"] = heat["events"] / heat["n"]
    heat["source1"] = heat["source1"].map(display_sample_name)
    heat["source2"] = heat["source2"].map(display_sample_name)
    sources = sorted(set(heat["source1"]) | set(heat["source2"]))
    hm = heat.pivot(index="source1", columns="source2", values="rate").reindex(index=sources, columns=sources)

    sns.set_theme(
        style="white",
        context="talk",
        rc={
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "axes.titlesize": 18,
            "axes.labelsize": 16,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 14,
        },
    )
    colors = {"Intra-habitat": "#2563EB", "Inter-habitat": "#E11D48"}
    # Two-row layout: A and B share the first row; C spans the second row.
    fig = plt.figure(figsize=(15, 12.5), constrained_layout=True)
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[0.90, 1.40],
        height_ratios=[0.78, 1.28],
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    order = ["Intra-habitat", "Inter-habitat"]
    for i, comp in enumerate(order):
        row = overall.loc[overall["comparison"].eq(comp)].iloc[0]
        ax1.errorbar(i, row["mean_round_event_rate"],
                     yerr=[[row["mean_round_event_rate"] - row["ci_low"]],
                           [row["ci_high"] - row["mean_round_event_rate"]]],
                     fmt="o", ms=14, capsize=7, lw=2.8, color=colors[comp])
        ax1.text(i, row["ci_high"] + .025,
                 f'{int(row["sharing_events"]):,}/{int(row["eligible_comparisons"]):,}',
                  ha="center", fontsize=14)
    ax1.set_xticks(range(2), ["Intra", "Inter"])
    ax1.set_ylabel("Sharing probability\n(mean across sampling rounds)")
    ax1.set_ylim(0, min(1, max(overall["ci_high"]) + .15))
    ax1.set_title("A  Intra-habitat baseline", loc="left", fontweight="bold")
    ax1.text(.03, .03, f"Wilcoxon P = {pvalue:.3g}",
              transform=ax1.transAxes, fontsize=14, va="bottom")

    for comp in order:
        g = by_round.loc[by_round["comparison"].eq(comp)].sort_values("time1")
        ax2.plot(g["time1"], g["event_rate"], "-o", label=comp,
                 color=colors[comp], lw=2.6, ms=7.5, alpha=.9)
    ax2.set_xlabel("Sampling round")
    ax2.set_ylabel("Sharing probability")
    ax2.set_xticks(sorted(by_round["time1"].unique()))
    ax2.legend(frameon=False, fontsize=14)
    ax2.set_title("B  Longitudinal consistency", loc="left", fontweight="bold")

    sns.heatmap(
        hm, ax=ax3, cmap="mako", vmin=0,
        vmax=np.nanquantile(hm.to_numpy(), .95), square=False,
        cbar_kws={"label": "Sharing probability", "shrink": .90, "pad": .02},
    )
    ax3.set_aspect("auto")
    ax3.set_xlabel("")
    ax3.set_ylabel("")
    ax3.tick_params(axis="x", labelsize=12, rotation=45)
    ax3.tick_params(axis="y", labelsize=12, rotation=0)
    ax3.collections[0].colorbar.ax.tick_params(labelsize=12)
    ax3.collections[0].colorbar.set_label("Sharing probability", fontsize=14)
    ax3.set_title("C  Source-pair baseline", loc="left", fontweight="bold")

    fig.suptitle("Within-habitat sharing provides a relative baseline for cross-habitat sharing\n"
                 "popANI >=99.999%; percent genome compared >=50%; same sampling round",
                 fontsize=19, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"Fig_R2Q3_intra_inter_baseline.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(overall.to_string(index=False))


if __name__ == "__main__":
    main()
