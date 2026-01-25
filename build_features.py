#!/usr/bin/env python3
"""
build_features.py

Merge base test DB TSV with one or more feature TSVs (ESM, EVcouplings, etc.)
into a single feature table keyed by an ID column (default: test_id).

Backwards compatible:
  - Old workflow: requires --base --esm1v --esm2 --esm3
New workflow:
  - You can omit esm args and instead pass --feature_tsvs with any number of TSVs.

Examples:

# NEW: merge only EVcouplings features (0.3 + 0.7)
python build_features.py \
  --base data/testset_db_with_wtseq_wt123.tsv \
  --feature_tsvs evc_results_0.3/evc_features.tsv evc_results_0.7/evc_features.tsv \
  --out data/processed/test_features_evc_03_07.tsv \
  --id_col test_id

# OLD: original merge behavior
python build_features.py \
  --base data/testset_db_with_wtseq_wt123.tsv \
  --esm1v data/processed/esm1v_features.tsv \
  --esm2  data/processed/esm2_features.tsv \
  --esm3  data/processed/esm3_features.tsv \
  --out data/processed/test_features_all.tsv
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd


def _read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _read_tsv_infer(path: str) -> pd.DataFrame:
    """
    Read TSV but let pandas infer numeric types where possible.
    We still keep strings for IDs safely.
    """
    df = pd.read_csv(path, sep="\t")
    return df


def _ensure_id_col(df: pd.DataFrame, id_col: str, name: str) -> None:
    if id_col not in df.columns:
        raise ValueError(f"[build_features] {name} missing id_col='{id_col}'. Columns={list(df.columns)}")


def _safe_drop_dupe_cols(base: pd.DataFrame, feat: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """
    If a feature TSV repeats columns already in base (other than id_col),
    drop them from feat to avoid accidental overwrites.
    """
    overlap = [c for c in feat.columns if c in base.columns and c != id_col]
    if overlap:
        feat = feat.drop(columns=overlap)
    return feat


def merge_feature_tsvs(
    base_tsv: str,
    feature_tsvs: List[str],
    out_tsv: str,
    id_col: str = "test_id",
    verbose: bool = True,
) -> None:
    base = _read_tsv_infer(base_tsv)
    _ensure_id_col(base, id_col, "BASE")

    if verbose:
        print(f"[build_features] BASE rows={len(base)} cols={len(base.columns)} file={base_tsv}")

    merged = base.copy()

    for fp in feature_tsvs:
        feat = _read_tsv_infer(fp)
        _ensure_id_col(feat, id_col, f"FEATURE({fp})")

        feat = _safe_drop_dupe_cols(merged, feat, id_col)

        before = len(merged)
        merged = merged.merge(feat, on=id_col, how="left", validate="one_to_one")
        after = len(merged)

        if before != after:
            raise ValueError(
                f"[build_features] Merge changed rowcount for {fp}: before={before} after={after} "
                f"(likely duplicate IDs in feature TSV)."
            )

        if verbose:
            n_new_cols = len(feat.columns) - 1
            print(f"[build_features] + merged {fp}  (+{n_new_cols} cols) -> cols={len(merged.columns)}")

    out_path = Path(out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_tsv, sep="\t", index=False)

    if verbose:
        print(f"[build_features] WROTE: {out_tsv}")
        print(f"[build_features] FINAL rows={len(merged)} cols={len(merged.columns)}")


def main():
    p = argparse.ArgumentParser()

    # required base
    p.add_argument("--base", required=True, help="Base TSV (test DB) containing id_col.")
    p.add_argument("--out", required=True, help="Output TSV path.")
    p.add_argument("--id_col", default="test_id", help="Join key column name (default: test_id).")

    # new generic features
    p.add_argument(
        "--feature_tsvs",
        nargs="*",
        default=[],
        help="One or more feature TSVs to merge (must contain id_col).",
    )

    # legacy args (optional now)
    p.add_argument("--esm1v", default=None, help="Legacy ESM1v feature TSV (optional).")
    p.add_argument("--esm2", default=None, help="Legacy ESM2 feature TSV (optional).")
    p.add_argument("--esm3", default=None, help="Legacy ESM3 feature TSV (optional).")

    # keep old flags if your pipeline expects them later (no-op here unless you extend)
    p.add_argument("--llr-zero-as-missing", action="store_true", help="(legacy flag; ignored here)")
    p.add_argument("--auto-flip-pll", action="store_true", help="(legacy flag; ignored here)")
    p.add_argument("--zscore-pll-by-wtseq", action="store_true", help="(legacy flag; ignored here)")
    p.add_argument("--alpha-risk", type=float, default=0.0, help="(legacy flag; ignored here)")
    p.add_argument("--features", default=None, help="(legacy flag; ignored here)")
    p.add_argument("--disable-features", default=None, help="(legacy flag; ignored here)")

    args = p.parse_args()

    # Build the final list of feature TSVs to merge
    feature_tsvs: List[str] = []

    # Legacy ESM inputs: if provided, append them in order
    for legacy_fp in [args.esm1v, args.esm2, args.esm3]:
        if legacy_fp is not None:
            feature_tsvs.append(legacy_fp)

    # Generic feature TSVs appended after
    if args.feature_tsvs:
        feature_tsvs.extend(args.feature_tsvs)

    if len(feature_tsvs) == 0:
        raise SystemExit(
            "[build_features] ERROR: No feature TSVs provided.\n"
            "Provide at least one of:\n"
            "  --feature_tsvs <a.tsv> <b.tsv> ...\n"
            "or legacy:\n"
            "  --esm1v <...> --esm2 <...> --esm3 <...>\n"
        )

    merge_feature_tsvs(
        base_tsv=args.base,
        feature_tsvs=feature_tsvs,
        out_tsv=args.out,
        id_col=args.id_col,
        verbose=True,
    )


if __name__ == "__main__":
    main()