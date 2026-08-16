from __future__ import annotations

import ast
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
TEXT_SUFFIXES = {".py", ".r", ".md", ".txt", ".tsv", ".csv", ".yml", ".yaml"}
ALLOWED_DATA_FILES = {
    ROOT / "data" / "README.md",
    ROOT / "data" / "metadata_template.csv",
}
FORBIDDEN_PATTERNS = {
    "Windows user path": re.compile(r"[A-Za-z]:[\\/]Users[\\/]", re.IGNORECASE),
    "Unix home path": re.compile(r"/(?:home|Users)/[^/\s]+/"),
    "Deprecated non-farmer label": re.compile(
        r"\b" + "stu" + "dent" + r"\b", re.IGNORECASE
    ),
}
EXPECTED_THRESHOLD_SNIPPETS = {
    Path("code/analysis/02_batch_effects/batch_effects.R"): (
        "min_prevalence <- 0.10",
        "permutation_n <- 999L",
    ),
    Path("code/analysis/04_feature_associations/transformation_sensitivity.R"): (
        'genus_matrix, "Taxonomy", names(maps), min_prevalence = 0.10',
        'c("Human_nasal_vs_bird", "Human_gut_vs_bird"), min_prevalence = 0.05',
        "min_mean = 1e-5",
    ),
    Path("code/analysis/05_arg_associations/arg_class_associations.R"): (
        "min_primary_prevalence <- 0.10",
        'fit_matrix(copies_per_cell, "Copies per cell")',
        'fit_matrix(rpkm_tss, "TSS-scaled RPKM")',
    ),
    Path("code/analysis/06_strain_sharing/build_canonical_comparisons.py"): (
        "PRIMARY_POPANI = 0.99999",
        "MIN_SAMPLE_GENOME_BREADTH = 0.50",
        "MIN_PERCENT_GENOME_COMPARED = 0.50",
    ),
    Path("code/analysis/06_strain_sharing/intra_inter_baseline.py"): (
        "PRIMARY_POPANI = 0.99999",
        "MIN_PERCENT_GENOME_COMPARED = 0.50",
    ),
    Path("code/analysis/06_strain_sharing/lagged_sensitivity.py"): (
        "PRIMARY_POPANI = 0.99999",
        "MIN_PERCENT_GENOME_COMPARED = 0.50",
    ),
    Path("code/analysis/07_hgt_sensitivity/ribosomal_taxonomic_sensitivity.py"): (
        'default=99.5',
        'candidates.length.gt(500)',
        'MAX_MAG_PAIR_ANI = 95.0',
        '.lt(MAX_MAG_PAIR_ANI)',
    ),
    Path("code/analysis/08_population_genomics/01_mk_multiple_testing.R"): (
        "min_gene_coverage <- 5",
        "min_gene_breadth <- 0.5",
        "min_axis_counts <- 3",
    ),
    Path("code/analysis/08_population_genomics/02_asymptotic_mk.R"): (
        "min_site_coverage <- 5",
        "min_genome_coverage <- 5",
        "min_genome_breadth <- 0.5",
        "bootstrap_n <- 2000L",
    ),
    Path("code/analysis/08_population_genomics/03_recurrence_analysis.R"): (
        "min_site_coverage <- 5",
        "min_genome_coverage <- 5",
        "min_genome_breadth <- 0.5",
        "permutation_n <- 2000L",
    ),
}


def main() -> None:
    failures: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.is_relative_to(ROOT / "data") and path not in ALLOWED_DATA_FILES:
            failures.append(f"Unexpected file under data/: {path.relative_to(ROOT)}")
        if path.stat().st_size > 2_000_000:
            failures.append(f"Unexpected file larger than 2 MB: {path.relative_to(ROOT)}")
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for label, pattern in FORBIDDEN_PATTERNS.items():
            if pattern.search(text):
                failures.append(f"{label} in {path.relative_to(ROOT)}")

    for path in (ROOT / "code").rglob("*.py"):
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeDecodeError) as error:
            failures.append(
                f"Python syntax error in {path.relative_to(ROOT)}: {error}"
            )

    for relative_path, snippets in EXPECTED_THRESHOLD_SNIPPETS.items():
        path = ROOT / relative_path
        if not path.exists():
            failures.append(f"Missing threshold-bearing script: {relative_path}")
            continue
        text = path.read_text(encoding="utf-8")
        for snippet in snippets:
            if snippet not in text:
                failures.append(
                    f"Expected final threshold not found in {relative_path}: {snippet}"
                )

    if failures:
        print("Release audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        raise SystemExit(1)
    print("Release audit passed: no private data files, local user paths, or Python syntax errors detected.")


if __name__ == "__main__":
    main()
