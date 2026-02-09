#!/usr/bin/env python
# scorer.py
# Usage: python scorer.py --in features_compiled.tsv --out scored.tsv --scores S1=llr_mean S2=delta_pll_mean S3="llr_mean-0.5*llr_std" S4="delta_pll_mean-0.5*delta_pll_std" [--drop-wt] [--id-col test_id] [--mutation-col mutation]

import argparse
import re
import sys
import numpy as np
import pandas as pd

SAFE_EXPR_RE = re.compile(r"^[A-Za-z0-9_\s\.\+\-\*\/\(\),]+$")

ID_CANDS = ["test_id", "_id", "name", "id", "variant_id", "seq_id"]
MUT_CANDS = ["_mutation_upper", "mutation", "mut_code", "mutation_code"]

def infer_col(df, cands):
    for c in cands:
        if c in df.columns:
            return c
    return None

def z(x):
    x = pd.to_numeric(x, errors="coerce")
    mu = np.nanmean(x.values)
    sd = np.nanstd(x.values)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - mu) / sd

def parse_scores(raw):
    out = []
    for s in raw:
        if "=" not in s:
            raise SystemExit(f"Bad --scores entry (need NAME=EXPR): {s}")
        name, expr = s.split("=", 1)
        name, expr = name.strip(), expr.strip()
        if not name or not expr:
            raise SystemExit(f"Bad --scores entry (empty name/expr): {s}")
        if not SAFE_EXPR_RE.match(expr):
            raise SystemExit(f"Unsafe characters in expr: {expr}")
        out.append((name, expr))
    return out

def eval_expr(df, expr):
    # Allow referencing columns by name + z(col) normalization.
    tokens = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expr))
    # tokens includes function name 'z' if used
    col_tokens = [t for t in tokens if t in df.columns]
    env = {c: pd.to_numeric(df[c], errors="coerce") for c in col_tokens}
    env["z"] = z
    try:
        return pd.eval(expr, local_dict=env, engine="python")
    except Exception as e:
        raise SystemExit(f"Failed to eval expr: {expr}\n{e}")

def rank_desc(x):
    x = pd.to_numeric(x, errors="coerce")
    return (-x).rank(method="average", na_option="bottom")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--scores", nargs="+", required=True)  # NAME=EXPR
    ap.add_argument("--drop-wt", action="store_true")
    ap.add_argument("--id-col", default="")
    ap.add_argument("--mutation-col", default="")
    args = ap.parse_args()

    df = pd.read_csv(args.in_path, sep="\t")

    id_col = args.id_col or infer_col(df, ID_CANDS)
    if not id_col:
        raise SystemExit("Could not infer id column; pass --id-col.")

    mut_col = args.mutation_col or infer_col(df, MUT_CANDS)

    if args.drop_wt:
        if "is_wt" in df.columns:
            df = df[pd.to_numeric(df["is_wt"], errors="coerce").fillna(0).astype(int) == 0].copy()
        else:
            if not mut_col:
                raise SystemExit("Requested --drop-wt but no mutation column and no is_wt column.")
            m = df[mut_col].astype(str).str.upper().str.strip()
            df = df[m.ne("WT")].copy()

    specs = parse_scores(args.scores)

    for name, expr in specs:
        df[name] = eval_expr(df, expr)
        df[f"rank_{name}"] = rank_desc(df[name])

    df.to_csv(args.out, sep="\t", index=False)
    print(f"[done] wrote {args.out} rows={len(df)}")

if __name__ == "__main__":
    main()