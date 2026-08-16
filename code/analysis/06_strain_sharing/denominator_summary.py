from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PRIMARY_POPANI = 0.99999
MIN_PERCENT_GENOME_COMPARED = 0.50


def one_row(x: pd.DataFrame, label: str, tested_catalogue_size: int) -> dict:
    eligible = x["percent_compared"].ge(MIN_PERCENT_GENOME_COMPARED)
    sharing = eligible & x["popANI"].ge(PRIMARY_POPANI)
    return {
        "scope": label,
        "eligible_comparisons": int(eligible.sum()),
        "sharing_events": int(sharing.sum()),
        "unique_tested_SGBs": int(x.loc[eligible, "genome"].nunique()),
        "unique_shared_SGBs": int(x.loc[sharing, "genome"].nunique()),
        "fraction_of_tested_SGBs_shared": x.loc[sharing, "genome"].nunique() / tested_catalogue_size,
        "sample_pairs_with_sharing": int(x.loc[sharing, "sample_pair"].nunique()),
        "source_pairs_with_sharing": int(x.loc[sharing, "source_pair"].nunique()),
        "sampling_rounds_with_sharing": int(pd.unique(pd.concat([
            x.loc[sharing, "time1"], x.loc[sharing, "time2"]], ignore_index=True)).size),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument(
        "--tested-catalogue-size", type=int, default=2571,
        help="Number of SGBs eligible for strain comparison (study default: 2571).",
    )
    ap.add_argument(
        "--full-catalogue-size", type=int, default=6075,
        help="Number of SGBs in the full genome catalogue (study default: 6075).",
    )
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(args.input, low_memory=False)
    rows = [
        one_row(x, "All comparisons", args.tested_catalogue_size),
        one_row(x.loc[x["same_time"].eq(True)], "Same sampling round", args.tested_catalogue_size),
        one_row(x.loc[~x["same_source"]], "Inter-habitat, all lags", args.tested_catalogue_size),
        one_row(x.loc[(~x["same_source"]) & x["same_time"].eq(True)], "Inter-habitat, same round", args.tested_catalogue_size),
        one_row(x.loc[x["same_source"] & x["same_time"].eq(True)], "Intra-habitat, same round", args.tested_catalogue_size),
    ]
    out = pd.DataFrame(rows)
    out["fraction_of_full_catalogue_SGBs_shared"] = (
        out["unique_shared_SGBs"] / args.full_catalogue_size
    )
    out.to_csv(args.out_dir / "sup_table_R2Q9_scope_and_denominators.csv", index=False)

    event = (
        x["percent_compared"].ge(MIN_PERCENT_GENOME_COMPARED)
        & x["popANI"].ge(PRIMARY_POPANI)
    )
    distribution = (x.loc[event].groupby("genome", observed=True)
                    .agg(sharing_events=("sample_pair", "size"),
                         sample_pairs=("sample_pair", "nunique"),
                         source_pairs=("source_pair", "nunique"),
                         sampling_rounds=("time1", "nunique"))
                    .sort_values("sharing_events", ascending=False).reset_index())
    distribution.to_csv(args.out_dir / "sup_table_R2Q9_shared_SGB_distribution.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
