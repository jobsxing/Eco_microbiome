from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sample_display_names import display_sample_name


# Embed TrueType text so all labels remain editable in Adobe software.
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    x = pd.read_csv(args.input, low_memory=False)
    x = x.loc[
        x["percent_compared"].ge(MIN_PERCENT_GENOME_COMPARED)
        & ~x["same_source"]
    ].copy()
    x["sharing"] = x["popANI"].ge(PRIMARY_POPANI)
    x["time_lag"] = x["time_lag"].astype(int)

    lag = (x.groupby("time_lag", observed=True)
           .agg(eligible_comparisons=("genome", "size"),
                sharing_events=("sharing", "sum"),
                unique_shared_genomes=("genome", lambda s: s[x.loc[s.index, "sharing"]].nunique()),
                sample_pairs=("sample_pair", "nunique"),
                sample_pairs_with_sharing=("sample_pair", lambda s: s[x.loc[s.index, "sharing"]].nunique()))
           .reset_index())
    lag["event_rate"] = lag["sharing_events"] / lag["eligible_comparisons"]
    cis = [proportion_confint(k, n, method="wilson") for k, n in zip(lag["sharing_events"], lag["eligible_comparisons"])]
    lag[["ci_low", "ci_high"]] = pd.DataFrame(cis, index=lag.index)
    lag.to_csv(args.out_dir / "sup_table_R2Q8_lag_sensitivity.csv", index=False)

    strict_events = x.loc[x["sharing"]].copy()
    same_genomes = set(strict_events.loc[strict_events["time_lag"].eq(0), "genome"])
    lagged_genomes = set(strict_events.loc[strict_events["time_lag"].gt(0), "genome"])
    same_pairs = set(strict_events.loc[strict_events["time_lag"].eq(0), "sample_pair"])
    lagged_pairs = set(strict_events.loc[strict_events["time_lag"].gt(0), "sample_pair"])
    summary = pd.DataFrame([{
        "same_round_events": int((strict_events["time_lag"] == 0).sum()),
        "lagged_events": int((strict_events["time_lag"] > 0).sum()),
        "same_round_unique_genomes": len(same_genomes),
        "lagged_unique_genomes": len(lagged_genomes),
        "lagged_only_unique_genomes": len(lagged_genomes - same_genomes),
        "all_inter_unique_genomes": len(same_genomes | lagged_genomes),
        "fraction_inter_genomes_missed_by_same_round_only": len(lagged_genomes - same_genomes) / len(same_genomes | lagged_genomes),
        "same_round_sample_pairs_with_sharing": len(same_pairs),
        "lagged_sample_pairs_with_sharing": len(lagged_pairs),
    }])
    summary.to_csv(args.out_dir / "sup_table_R2Q8_same_vs_lagged_summary.csv", index=False)

    lagged = x.loc[x["time_lag"].gt(0)].copy()
    first_is_earlier = lagged["time1"] < lagged["time2"]
    lagged["earlier_source"] = np.where(first_is_earlier, lagged["source1"], lagged["source2"])
    lagged["later_source"] = np.where(first_is_earlier, lagged["source2"], lagged["source1"])
    ordered = (lagged.groupby(["earlier_source", "later_source"], observed=True)
               .agg(eligible_comparisons=("genome", "size"), sharing_events=("sharing", "sum"))
               .reset_index())
    ordered["event_rate"] = ordered["sharing_events"] / ordered["eligible_comparisons"]
    ordered.to_csv(args.out_dir / "sup_table_R2Q8_temporally_ordered_source_pairs.csv", index=False)

    ordered["earlier_source"] = ordered["earlier_source"].map(display_sample_name)
    ordered["later_source"] = ordered["later_source"].map(display_sample_name)

    sources = sorted(set(ordered["earlier_source"]) | set(ordered["later_source"]))
    hm = ordered.pivot(index="earlier_source", columns="later_source", values="event_rate").reindex(index=sources, columns=sources)

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
    # Match R2Q3: A and B share the first row; C spans their full width below.
    fig = plt.figure(figsize=(15, 12.5), constrained_layout=True)
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.15, 1.0],
        height_ratios=[0.78, 1.28],
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    ax1.plot(lag["time_lag"], lag["event_rate"], color="#7C3AED", marker="o",
             lw=2.8, ms=8)
    ax1.fill_between(lag["time_lag"], lag["ci_low"], lag["ci_high"], color="#A78BFA", alpha=.3)
    ax1.axvline(0, ls="--", color="#334155", lw=1.4)
    ax1.set_xlabel("Lag between sampling rounds")
    ax1.set_ylabel("Cross-habitat sharing probability")
    ax1.set_xticks(lag["time_lag"])
    ax1.set_title("A  Same-round and lagged comparisons", loc="left", fontweight="bold")

    width = .38
    ax2.bar(lag["time_lag"] - width/2, lag["unique_shared_genomes"], width,
            color="#0891B2", label="Shared SGBs")
    ax2.bar(lag["time_lag"] + width/2, lag["sample_pairs_with_sharing"], width,
            color="#F59E0B", label="Sample pairs")
    ax2.set_xlabel("Lag between sampling rounds")
    ax2.set_ylabel("Number detected")
    ax2.set_xticks(lag["time_lag"])
    ax2.legend(frameon=False, fontsize=14)
    ax2.set_title("B  Detection opportunities", loc="left", fontweight="bold")

    sns.heatmap(
        hm, ax=ax3, cmap="rocket_r", vmin=0,
        vmax=np.nanquantile(hm.to_numpy(), .95), square=False,
        cbar_kws={"label": "Sharing probability", "shrink": .90, "pad": .02},
    )
    # Allow C to occupy the complete A+B plotting width without a blank region.
    ax3.set_aspect("auto")
    ax3.set_xlabel("Later sampled source")
    ax3.set_ylabel("Earlier sampled source")
    ax3.tick_params(axis="x", labelsize=12, rotation=45)
    ax3.tick_params(axis="y", labelsize=12, rotation=0)
    ax3.collections[0].colorbar.ax.tick_params(labelsize=12)
    ax3.collections[0].colorbar.set_label("Sharing probability", fontsize=14)
    ax3.set_title("C  Temporally ordered pairs (lag > 0)", loc="left", fontweight="bold")

    fig.suptitle("Lagged sensitivity analysis of cross-habitat strain sharing\n"
                 "Temporal ordering is not interpreted as transmission direction",
                 fontsize=19, fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(args.out_dir / f"Fig_R2Q8_lagged_sensitivity.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(summary.to_string(index=False))
    print(lag.to_string(index=False))


if __name__ == "__main__":
    main()
