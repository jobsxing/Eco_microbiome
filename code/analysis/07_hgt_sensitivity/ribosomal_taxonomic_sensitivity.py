from __future__ import annotations

import argparse
import collections
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INTERFACE_FILES = {
    "Farmer-animal": "dt_human_animal.csv",
    "Farmer-environment": "dt_human_env.csv",
    "Animal-environment": "dt_animal_env_raw.csv",
}
TAXONOMIC_ORDER = [
    "Within genus",
    "Within family, between genera",
    "Between families",
    "Unresolved taxonomy",
]
COLORS = {
    "Within genus": "#93435B",
    "Within family, between genera": "#D8A64A",
    "Between families": "#367A87",
    "Unresolved taxonomy": "#BDBDBD",
}
MAX_MAG_PAIR_ANI = 95.0


def bin_from_seqid(value: str) -> str:
    return re.sub(r"_k\d+_\d+$", "", str(value))


def canonical_pair(left: pd.Series, right: pd.Series) -> pd.Series:
    left = left.astype("string")
    right = right.astype("string")
    return pd.Series(
        np.where(left <= right, left + " | " + right, right + " | " + left),
        index=left.index,
    )


def read_mag_pair_ani(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    ani = pd.read_csv(path, sep=sep, low_memory=False)
    required = {"mag1", "mag2", "whole_genome_ani"}
    missing = sorted(required - set(ani.columns))
    if missing:
        raise ValueError(f"MAG-pair ANI table is missing columns: {missing}")
    ani = ani.loc[:, ["mag1", "mag2", "whole_genome_ani"]].copy()
    ani["whole_genome_ani"] = pd.to_numeric(ani["whole_genome_ani"], errors="coerce")
    if ani["whole_genome_ani"].dropna().le(1).all():
        ani["whole_genome_ani"] *= 100
    ani["MAG_pair"] = canonical_pair(ani["mag1"], ani["mag2"])
    if ani["MAG_pair"].duplicated().any():
        raise ValueError("MAG-pair ANI table contains duplicate unordered MAG pairs")
    return ani.loc[:, ["MAG_pair", "whole_genome_ani"]]


def apply_mag_ani_filter(frame: pd.DataFrame, ani: pd.DataFrame) -> pd.DataFrame:
    frame = frame.merge(ani, on="MAG_pair", how="left", validate="many_to_one")
    missing = int(frame["whole_genome_ani"].isna().sum())
    if missing:
        raise ValueError(f"Whole-genome ANI is missing for {missing:,} candidate fragment hits")
    return frame.loc[frame["whole_genome_ani"].lt(MAX_MAG_PAIR_ANI)].copy()


def overlap(seqid: str, start: int, end: int, intervals: dict[str, list[tuple[int, int]]]) -> bool:
    start, end = sorted((int(start), int(end)))
    return any(max(start, left) <= min(end, right) for left, right in intervals.get(seqid, ()))


def read_rrna(rrna_dir: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for path in rrna_dir.glob("*.merged.gff"):
        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) >= 5:
                    intervals[fields[0]].append(tuple(sorted((int(fields[3]), int(fields[4])))))
    return intervals


def read_ribosomal_proteins(rrna_dir: Path, cds_dir: Path) -> dict[str, list[tuple[int, int]]]:
    protein_ids: set[str] = set()
    for path in rrna_dir.glob("*.ribosomal.domtblout"):
        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) >= 4:
                    protein_ids.add(fields[3])

    intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    mapped: set[str] = set()
    for path in cds_dir.glob("*_cds.gff"):
        bin_id = path.name.removesuffix("_cds.gff")
        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                match = re.search(r"ID=\d+_(\d+)", fields[8])
                if not match:
                    continue
                protein_id = f"{bin_id}_{fields[0]}_{match.group(1)}"
                if protein_id in protein_ids:
                    seqid = f"{bin_id}_{fields[0]}"
                    intervals[seqid].append(tuple(sorted((int(fields[3]), int(fields[4])))))
                    mapped.add(protein_id)
    if mapped != protein_ids:
        raise RuntimeError(f"Mapped {len(mapped):,}/{len(protein_ids):,} ribosomal proteins")
    return intervals


def annotate_overlap(frame: pd.DataFrame, rrna: dict, ribosomal: dict) -> pd.DataFrame:
    frame = frame.copy()
    for label, intervals in [("rRNA", rrna), ("Ribosomal_protein", ribosomal)]:
        frame[f"{label}_q"] = [
            overlap(seqid, start, end, intervals)
            for seqid, start, end in zip(frame.qseqid, frame.qstart, frame.qend)
        ]
        frame[f"{label}_s"] = [
            overlap(seqid, start, end, intervals)
            for seqid, start, end in zip(frame.sseqid, frame.sstart, frame.send)
        ]
        frame[f"{label}_either"] = frame[f"{label}_q"] | frame[f"{label}_s"]
        frame[f"{label}_both"] = frame[f"{label}_q"] & frame[f"{label}_s"]
    frame["Any_ribosomal_either"] = frame.rRNA_either | frame.Ribosomal_protein_either
    frame["Any_ribosomal_both"] = frame.rRNA_both | frame.Ribosomal_protein_both
    return frame


def classify_taxonomic_distance(frame: pd.DataFrame) -> pd.Series:
    q_family = frame.q_family.replace("", pd.NA)
    s_family = frame.s_family.replace("", pd.NA)
    q_genus = frame.q_genus.replace("", pd.NA)
    s_genus = frame.s_genus.replace("", pd.NA)
    family_known = q_family.notna() & s_family.notna()
    genus_known = q_genus.notna() & s_genus.notna()
    same_family = q_family.eq(s_family).fillna(False)
    same_genus = q_genus.eq(s_genus).fillna(False)
    return pd.Series(
        np.select(
            [
                family_known & ~same_family,
                family_known & same_family & genus_known & ~same_genus,
                family_known & same_family & genus_known & same_genus,
            ],
            ["Between families", "Within family, between genera", "Within genus"],
            default="Unresolved taxonomy",
        ),
        index=frame.index,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ecosystem", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--mag-pair-ani",
        type=Path,
        required=True,
        help=(
            "CSV/TSV with columns mag1, mag2, and whole_genome_ani. "
            "Only MAG pairs with whole-genome ANI <95%% are eligible."
        ),
    )
    parser.add_argument("--identity", type=float, default=99.5)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rrna_dir = args.ecosystem / "data/bin/mags_anno/rRNA_anno"
    cds_dir = args.ecosystem / "data/bin/mags_anno/cds"
    inter_dir = args.ecosystem / "data/bin/mags_anno/contigs/inter"
    figure_dir = args.ecosystem / "Figure/Figure 4"

    rrna = read_rrna(rrna_dir)
    ribosomal = read_ribosomal_proteins(rrna_dir, cds_dir)
    mag_pair_ani = read_mag_pair_ani(args.mag_pair_ani)

    pieces = []
    for path in sorted(inter_dir.glob("*_filtered.txt")):
        part = pd.read_csv(path, sep="\t", low_memory=False)
        part["input_table"] = path.name
        pieces.append(part)
    candidates = pd.concat(pieces, ignore_index=True)
    candidates["pident"] = pd.to_numeric(candidates.pident, errors="coerce")
    candidates["length"] = pd.to_numeric(candidates.length, errors="coerce")
    candidates = candidates.loc[
        candidates.pident.ge(args.identity) & candidates.length.gt(500)
    ].copy()
    candidates["q_bin"] = candidates.qseqid.map(bin_from_seqid)
    candidates["s_bin"] = candidates.sseqid.map(bin_from_seqid)
    candidates["MAG_pair"] = canonical_pair(candidates.q_bin, candidates.s_bin)
    candidates = apply_mag_ani_filter(candidates, mag_pair_ani)
    candidates = annotate_overlap(candidates, rrna, ribosomal)

    total_hits = len(candidates)
    total_pairs = candidates.MAG_pair.nunique()
    overlap_rows = []
    for feature, either, both in [
        ("5S/16S rRNA", "rRNA_either", "rRNA_both"),
        ("Ribosomal proteins", "Ribosomal_protein_either", "Ribosomal_protein_both"),
        ("Any ribosomal feature", "Any_ribosomal_either", "Any_ribosomal_both"),
    ]:
        for definition, column in [("Either MAG", either), ("Both MAGs", both)]:
            subset = candidates.loc[candidates[column]]
            overlap_rows.append(
                {
                    "feature": feature,
                    "overlap_definition": definition,
                    "fragment_hits": len(subset),
                    "fraction_fragment_hits": len(subset) / total_hits,
                    "affected_MAG_pairs": subset.MAG_pair.nunique(),
                    "fraction_affected_MAG_pairs": subset.MAG_pair.nunique() / total_pairs,
                    "total_fragment_hits": total_hits,
                    "total_MAG_pairs": total_pairs,
                }
            )
    overlap_summary = pd.DataFrame(overlap_rows)
    overlap_summary.to_csv(args.out / "sup_table_R2Q4_ribosomal_overlap_summary.csv", index=False)

    interface_parts = []
    genus_parts = []
    for interface, filename in INTERFACE_FILES.items():
        frame = pd.read_csv(figure_dir / filename, low_memory=False)
        frame["qseqid"] = frame.sample_id_qseqid
        frame["sseqid"] = frame["i.sseqid"]
        frame["length"] = (pd.to_numeric(frame.qend) - pd.to_numeric(frame.qstart)).abs() + 1
        frame["pident"] = pd.to_numeric(frame.pident, errors="coerce")
        frame = frame.loc[frame.pident.ge(args.identity) & frame.length.gt(500)].copy()
        frame["MAG_pair"] = canonical_pair(frame.q_bin, frame.s_bin)
        frame = apply_mag_ani_filter(frame, mag_pair_ani)
        frame = annotate_overlap(frame, rrna, ribosomal)
        frame["taxonomic_distance"] = classify_taxonomic_distance(frame)
        frame["interface"] = interface
        genus_parts.append(frame)
        retained = frame.loc[~frame.Any_ribosomal_either]
        interface_parts.append(
            {
                "interface": interface,
                "original_fragment_hits": len(frame),
                "nonribosomal_fragment_hits": len(retained),
                "retained_fragment_fraction": len(retained) / len(frame),
                "original_MAG_pairs": frame.MAG_pair.nunique(),
                "nonribosomal_MAG_pairs": retained.MAG_pair.nunique(),
                "retained_MAG_pair_fraction": retained.MAG_pair.nunique() / frame.MAG_pair.nunique(),
            }
        )
    interface_summary = pd.DataFrame(interface_parts)
    interface_summary.to_csv(args.out / "sup_table_R2Q4_interface_sensitivity.csv", index=False)

    genus = pd.concat(genus_parts, ignore_index=True)
    pair_level = genus[["interface", "MAG_pair", "taxonomic_distance"]].drop_duplicates()
    pair_counts = (
        pair_level.groupby(["interface", "taxonomic_distance"], observed=True)
        .size().rename("unique_MAG_pairs").reset_index()
    )
    hit_counts = (
        genus.groupby(["interface", "taxonomic_distance"], observed=True)
        .size().rename("fragment_hits").reset_index()
    )
    genus_summary = pair_counts.merge(hit_counts, on=["interface", "taxonomic_distance"], how="outer").fillna(0)
    complete = pd.MultiIndex.from_product(
        [INTERFACE_FILES, TAXONOMIC_ORDER], names=["interface", "taxonomic_distance"]
    )
    genus_summary = genus_summary.set_index(["interface", "taxonomic_distance"]).reindex(complete, fill_value=0).reset_index()
    genus_summary["fraction_MAG_pairs"] = genus_summary.unique_MAG_pairs / genus_summary.groupby("interface").unique_MAG_pairs.transform("sum")
    genus_summary["fraction_fragment_hits"] = genus_summary.fragment_hits / genus_summary.groupby("interface").fragment_hits.transform("sum")
    genus_summary.to_csv(args.out / "sup_table_R2Q7_genus_split_by_interface.csv", index=False)

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial"],
        "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 13,
        "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
        "legend.fontsize": 10, "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2))

    # A: ribosomal overlap
    a = overlap_summary.copy()
    features = ["5S/16S rRNA", "Ribosomal proteins", "Any ribosomal feature"]
    x = np.arange(len(features)); width = 0.34
    for offset, definition, color in [(-width / 2, "Either MAG", "#A33A3A"), (width / 2, "Both MAGs", "#3A5875")]:
        values = [100 * a.loc[(a.feature == f) & (a.overlap_definition == definition), "fraction_fragment_hits"].iloc[0] for f in features]
        axes[0, 0].bar(x + offset, values, width, label=definition, color=color)
        for xpos, value in zip(x + offset, values):
            axes[0, 0].text(xpos, value + 0.06, f"{value:.2f}%", ha="center", va="bottom", fontsize=9)
    axes[0, 0].set_xticks(x, ["5S/16S\nrRNA", "Ribosomal\nproteins", "Any ribosomal\nfeature"])
    axes[0, 0].set_ylabel("Candidate fragment hits (%)")
    axes[0, 0].set_title("A  Ribosomal-feature overlap", loc="left")
    axes[0, 0].legend(frameon=False)

    # B: retention after exclusion
    b = interface_summary.set_index("interface").reindex(INTERFACE_FILES)
    y = np.arange(len(b))
    axes[0, 1].barh(y - 0.16, 100 * b.retained_fragment_fraction, 0.30, color="#367A87", label="Fragment hits")
    axes[0, 1].barh(y + 0.16, 100 * b.retained_MAG_pair_fraction, 0.30, color="#D8A64A", label="MAG pairs")
    for i, row in enumerate(b.itertuples()):
        axes[0, 1].text(100 * row.retained_fragment_fraction - 0.15, i - 0.16, f"{100*row.retained_fragment_fraction:.1f}%", ha="right", va="center", color="white", fontsize=9)
        axes[0, 1].text(100 * row.retained_MAG_pair_fraction - 0.15, i + 0.16, f"{100*row.retained_MAG_pair_fraction:.1f}%", ha="right", va="center", color="#202020", fontsize=9)
    axes[0, 1].set_xlim(90, 100.4)
    axes[0, 1].set_yticks(y, list(b.index))
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("Retained after ribosomal exclusion (%)")
    axes[0, 1].set_title("B  Ecological-interface sensitivity", loc="left")
    axes[0, 1].legend(frameon=False)

    # C/D: genus split
    for ax, fraction, title in [
        (axes[1, 0], "fraction_MAG_pairs", "C  Unique MAG pairs"),
        (axes[1, 1], "fraction_fragment_hits", "D  Candidate fragment hits"),
    ]:
        left = np.zeros(len(INTERFACE_FILES))
        for category in TAXONOMIC_ORDER:
            values = (
                genus_summary.loc[genus_summary.taxonomic_distance.eq(category)]
                .set_index("interface").reindex(INTERFACE_FILES)[fraction].to_numpy()
            )
            ax.barh(y, values * 100, left=left * 100, height=0.58,
                    color=COLORS[category], edgecolor="white", linewidth=0.7, label=category)
            for i, value in enumerate(values):
                if value >= 0.045:
                    color = "white" if category in {"Within genus", "Between families"} else "#202020"
                    ax.text((left[i] + value / 2) * 100, i, f"{100*value:.1f}", ha="center", va="center", fontsize=9, color=color)
            left += values
        ax.set_xlim(0, 100)
        ax.set_yticks(y, list(INTERFACE_FILES) if ax is axes[1, 0] else [])
        ax.invert_yaxis()
        ax.set_xlabel("Proportion (%)")
        ax.set_title(title, loc="left")
    handles, labels = axes[1, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.01))
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#DDDDDD", linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)
    fig.tight_layout(rect=(0, 0.07, 1, 1), h_pad=2.5, w_pad=2.2)
    fig.savefig(args.out / "Fig_R2Q4_R2Q7_HGT_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(args.out / "Fig_R2Q4_R2Q7_HGT_sensitivity.png", dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(overlap_summary.to_string(index=False))
    print("\n", interface_summary.to_string(index=False))
    print("\n", genus_summary.to_string(index=False))


if __name__ == "__main__":
    main()
