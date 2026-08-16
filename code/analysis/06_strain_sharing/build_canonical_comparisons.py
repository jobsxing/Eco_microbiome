from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


KEEP = [
    "genome", "name1", "name2", "coverage_overlap", "compared_bases_count",
    "consensus_SNPs", "population_SNPs", "popANI", "conANI",
    "percent_compared",
]
META_KEEP = [
    "Sample_core_for_mpa", "Sample Source", "Sampling Time", "Facet_Group",
    "ID", "Sampling scheme",
]

PRIMARY_POPANI = 0.99999
MIN_SAMPLE_GENOME_BREADTH = 0.50
MIN_PERCENT_GENOME_COMPARED = 0.50
POPANI_SENSITIVITY_THRESHOLDS = (0.995, 0.9995, 0.9999, PRIMARY_POPANI)


def clean_sample(value: pd.Series) -> pd.Series:
    return (value.astype("string")
            .str.replace(r"\.sorted\.bam$", "", regex=True)
            .str.replace(r"\.bam$", "", regex=True))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument(
        "--profile-breadth",
        type=Path,
        required=True,
        help=(
            "Tab- or comma-delimited genome-by-sample table with columns "
            "sample_id, genome, and breadth. Both samples must have genomic "
            "breadth >= 0.50 before a genome comparison is eligible."
        ),
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(args.raw_dir.glob("*genomeWide_compare.tsv"))
    if not files:
        raise FileNotFoundError(f"No raw comparison files in {args.raw_dir}")

    chunks = []
    for path in files:
        x = pd.read_csv(path, sep="\t", usecols=KEEP, low_memory=False)
        x["source_file"] = path.name
        chunks.append(x)
    x = pd.concat(chunks, ignore_index=True)
    raw_rows = len(x)

    for col in KEEP[3:]:
        x[col] = pd.to_numeric(x[col], errors="coerce")
    x["sample1"] = clean_sample(x["name1"])
    x["sample2"] = clean_sample(x["name2"])

    meta = pd.read_csv(args.metadata, low_memory=False)
    missing_meta_cols = sorted(set(META_KEEP) - set(meta.columns))
    if missing_meta_cols:
        raise ValueError(f"Missing metadata columns: {missing_meta_cols}")
    meta = meta[META_KEEP].copy()
    if meta["Sample_core_for_mpa"].duplicated().any():
        raise ValueError("Metadata key is not unique")

    for side in (1, 2):
        rename = {
            "Sample_core_for_mpa": f"sample{side}",
            "Sample Source": f"source{side}",
            "Sampling Time": f"time{side}",
            "Facet_Group": f"facet{side}",
            "ID": f"id{side}",
            "Sampling scheme": f"scheme{side}",
        }
        x = x.merge(meta.rename(columns=rename), on=f"sample{side}", how="left", validate="many_to_one")

    missing_1 = int(x["source1"].isna().sum())
    missing_2 = int(x["source2"].isna().sum())
    if missing_1 or missing_2:
        raise ValueError(
            "Comparison tables contain sample identifiers absent from metadata: "
            f"side 1 = {missing_1:,}, side 2 = {missing_2:,}"
        )
    self_rows = int((x["sample1"] == x["sample2"]).sum())
    x = x.loc[x["sample1"] != x["sample2"]].copy()

    swap = x["sample1"] > x["sample2"]
    paired_fields = ["sample", "source", "time", "facet", "id", "scheme", "name"]
    for stem in paired_fields:
        a, b = f"{stem}1", f"{stem}2"
        if a in x.columns and b in x.columns:
            left = x[a].copy()
            x.loc[swap, a] = x.loc[swap, b].to_numpy()
            x.loc[swap, b] = left.loc[swap].to_numpy()

    x = x.sort_values(
        ["genome", "sample1", "sample2", "percent_compared", "compared_bases_count"],
        ascending=[True, True, True, False, False], na_position="last",
    )
    duplicate_rows = int(x.duplicated(["genome", "sample1", "sample2"], keep=False).sum())
    x = x.drop_duplicates(["genome", "sample1", "sample2"], keep="first").copy()

    breadth_sep = "\t" if args.profile_breadth.suffix.lower() in {".tsv", ".txt"} else ","
    profile_breadth = pd.read_csv(args.profile_breadth, sep=breadth_sep, low_memory=False)
    required_breadth_columns = {"sample_id", "genome", "breadth"}
    missing_breadth_columns = sorted(required_breadth_columns - set(profile_breadth.columns))
    if missing_breadth_columns:
        raise ValueError(
            f"Profile breadth table is missing columns: {missing_breadth_columns}"
        )
    profile_breadth = profile_breadth.loc[:, ["sample_id", "genome", "breadth"]].copy()
    profile_breadth["sample_id"] = clean_sample(profile_breadth["sample_id"])
    profile_breadth["breadth"] = pd.to_numeric(profile_breadth["breadth"], errors="coerce")
    if profile_breadth.duplicated(["sample_id", "genome"]).any():
        raise ValueError("Profile breadth table has duplicate sample_id/genome rows")

    pre_breadth_rows = len(x)
    for side in (1, 2):
        side_breadth = profile_breadth.rename(
            columns={"sample_id": f"sample{side}", "breadth": f"breadth{side}"}
        )
        x = x.merge(
            side_breadth,
            on=[f"sample{side}", "genome"],
            how="left",
            validate="many_to_one",
        )
    missing_breadth_side1 = int(x["breadth1"].isna().sum())
    missing_breadth_side2 = int(x["breadth2"].isna().sum())
    x = x.loc[
        x["breadth1"].ge(MIN_SAMPLE_GENOME_BREADTH)
        & x["breadth2"].ge(MIN_SAMPLE_GENOME_BREADTH)
    ].copy()

    x["time1"] = pd.to_numeric(x["time1"], errors="coerce")
    x["time2"] = pd.to_numeric(x["time2"], errors="coerce")
    x["same_time"] = x["time1"].eq(x["time2"])
    x["time_lag"] = (x["time1"] - x["time2"]).abs()
    x["same_source"] = x["source1"].eq(x["source2"])
    x["same_facet"] = x["facet1"].eq(x["facet2"])
    x["sample_pair"] = x["sample1"] + " | " + x["sample2"]
    x["source_pair"] = np.where(
        x["source1"] <= x["source2"],
        x["source1"] + " | " + x["source2"],
        x["source2"] + " | " + x["source1"],
    )
    x["facet_pair"] = np.where(
        x["facet1"] <= x["facet2"],
        x["facet1"] + " | " + x["facet2"],
        x["facet2"] + " | " + x["facet1"],
    )

    eligible = x["percent_compared"].ge(MIN_PERCENT_GENOME_COMPARED)
    rows = []
    for threshold in POPANI_SENSITIVITY_THRESHOLDS:
        event = eligible & x["popANI"].ge(threshold)
        for scope, mask in {
            "all": pd.Series(True, index=x.index),
            "same_time": x["same_time"],
            "same_time_intra": x["same_time"] & x["same_source"],
            "same_time_inter": x["same_time"] & ~x["same_source"],
            "all_inter": ~x["same_source"],
        }.items():
            denom = eligible & mask
            num = event & mask
            rows.append({
                "popANI_threshold": threshold,
                "scope": scope,
                "eligible_comparisons": int(denom.sum()),
                "sharing_events": int(num.sum()),
                "event_rate": float(num.sum() / denom.sum()) if denom.sum() else np.nan,
                "unique_shared_genomes": int(x.loc[num, "genome"].nunique()),
                "sample_pairs_with_sharing": int(x.loc[num, "sample_pair"].nunique()),
            })
    threshold_summary = pd.DataFrame(rows)

    out_table = args.out_dir / "analysis_all_comparisons.csv.gz"
    x.to_csv(out_table, index=False, compression="gzip")
    threshold_summary.to_csv(args.out_dir / "sup_table_strain_threshold_sensitivity.csv", index=False)

    quant = x[["popANI", "percent_compared", "compared_bases_count"]].quantile(
        [0, .01, .05, .25, .5, .75, .95, .99, 1]
    ).reset_index(names="quantile")
    quant.to_csv(args.out_dir / "qc_numeric_quantiles.csv", index=False)
    qc = {
        "raw_file_count": len(files),
        "raw_rows": raw_rows,
        "self_comparison_rows_removed": self_rows,
        "rows_participating_in_duplicate_keys_before_deduplication": duplicate_rows,
        "canonical_rows": len(x),
        "rows_before_sample_breadth_filter": pre_breadth_rows,
        "profile_breadth_missing_side1": missing_breadth_side1,
        "profile_breadth_missing_side2": missing_breadth_side2,
        "minimum_sample_genome_breadth": MIN_SAMPLE_GENOME_BREADTH,
        "minimum_percent_genome_compared": MIN_PERCENT_GENOME_COMPARED,
        "primary_popANI": PRIMARY_POPANI,
        "metadata_missing_side1": missing_1,
        "metadata_missing_side2": missing_2,
        "unique_genomes": int(x["genome"].nunique()),
        "unique_samples": int(pd.unique(pd.concat([x["sample1"], x["sample2"]], ignore_index=True)).size),
    }
    (args.out_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
    print(json.dumps(qc, indent=2))
    print(threshold_summary.to_string(index=False))


if __name__ == "__main__":
    main()
