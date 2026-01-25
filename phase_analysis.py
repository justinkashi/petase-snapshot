#!/usr/bin/env python
# phase_analysis.py
# Usage: python phase_analysis.py --in INPUT.tsv --outdir phase_dashboard [--id-col _id] [--group-col _wt_seq] [--score-cols s1,s2,s3] [--label-cols y1,y2] [--label-higher-is-better|--label-lower-is-better] [--topk 200] [--feature-cols f1,f2,f3] [--max-hist 80] [--max-scatter-per-score 12] [--max-scatter-per-label 12]                          

#MODE A: Dashboard
#   --in                                    #for eg. features_compiled.tsv
#   --outdir phase_dashboard                # Output folder (writes: tables/, plots/, summary.json)
#   [--id-col _id]                          # Row identifier column (auto-inferred if omitted)
#   [--group-col _wt_seq]                   # Optional grouping (e.g., backbone/WT) for per-group counts/summaries
#   [--feature-cols f1,f2,f3]               # Optional: restrict features analyzed; default = auto-detect numeric cols excluding id/group/scores/labels
#   [--max-hist 80]                         # Max feature histograms to render

# MODE B: Score/rank diagnostics. 1) stability (topk overlap of scores). 2) candidate consistency top_candidates
#   --score-cols s1,s2,s3                   # Comma list of score columns to treat as ranking signals (enables topK overlap + top candidates + score-vs-feature plots)
#                                           # this will output topk * n_scoringlogics so lets say theres s1,s2 with topk=200, then youll see 400 rows in the top_candidates_by_score.tsv 
#   [--topk 200]                            # TopK used for overlap tables + top-candidate tables
#   [--max-scatter-per-score 12]            # Max scatter plots per score (score vs strongest-correlated features)

# MODE C: Evaluation (requires labels in the TSV)
#   --label-cols y1,y2                      # Comma list of ground-truth/proxy label columns (enables Pearson/Spearman + topK overlap vs labels)
#   [--label-higher-is-better]              # Default: higher label = better
#   [--label-lower-is-better]               # Use if label is a loss/error where lower = better
#   [--max-scatter-per-label 12]            # Max scatter plots per label (label vs best-aligned scores)

import argparse
import json
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def mkdirp(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def safe_cols(df: pd.DataFrame, cols: List[str]) -> List[str]:
    return [c for c in cols if c in df.columns]


def infer_id_col(df: pd.DataFrame, provided: Optional[str]) -> str:
    if provided and provided in df.columns:
        return provided
    for c in ["_id", "test_id", "name", "id", "variant_id", "seq_id"]:
        if c in df.columns:
            return c
    raise SystemExit("ERROR: could not infer id column; pass --id-col")


def infer_group_col(df: pd.DataFrame, provided: Optional[str]) -> Optional[str]:
    if provided and provided in df.columns:
        return provided
    for c in ["_wt_seq", "wt_seq", "backbone", "parent", "wt_name"]:
        if c in df.columns:
            return c
    return None


def infer_numeric_feature_cols(
    df: pd.DataFrame,
    id_col: str,
    group_col: Optional[str],
    score_cols: List[str],
    label_cols: List[str],
    user_feature_cols: Optional[List[str]] = None,
) -> List[str]:
    drop = {id_col}
    if group_col:
        drop.add(group_col)
    drop |= set(score_cols)
    drop |= set(label_cols)

    if user_feature_cols:
        cols = [c for c in user_feature_cols if c in df.columns]
    else:
        cols = [c for c in df.columns if c not in drop]

    numeric = []
    for c in cols:
        if pd.api.types.is_numeric_dtype(df[c]):
            numeric.append(c)
        else:
            x = pd.to_numeric(df[c], errors="coerce")
            if x.notna().sum() > 0:
                df[c] = x
                numeric.append(c)
    return numeric


def summarize_missingness(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = []
    n = len(df)
    for c in cols:
        miss = int(df[c].isna().sum())
        out.append((c, miss, miss / max(n, 1)))
    return (
        pd.DataFrame(out, columns=["col", "missing_n", "missing_frac"])
        .sort_values(["missing_frac", "col"], ascending=[False, True])
        .reset_index(drop=True)
    )


def corr_matrix(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    X = df[cols].apply(pd.to_numeric, errors="coerce")
    return X.corr(method="pearson", min_periods=50)


def save_corr_heatmap(corr: pd.DataFrame, outpath: str, title: str) -> None:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    im = ax.imshow(corr.values, aspect="auto")
    ax.set_title(title)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
    ax.set_yticklabels(corr.index, fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

def save_hist_grid(df: pd.DataFrame, cols: List[str], outpath: str, title: str, bins: int = 60) -> None:
    cols = [c for c in cols if c in df.columns]
    if not cols:
        return

    xs = []
    for c in cols:
        x = to_num(df[c]).dropna()
        xs.append(x if len(x) >= 50 else None)

    # drop columns with too few points
    cols2, xs2 = [], []
    for c, x in zip(cols, xs):
        if x is not None:
            cols2.append(c)
            xs2.append(x)
    cols, xs = cols2, xs2
    if not cols:
        return

    n = len(cols)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5.2 * ncols, 3.6 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    for i, (c, x) in enumerate(zip(cols, xs)):
        r, k = divmod(i, ncols)
        ax = axes[r, k]
        ax.hist(x.values, bins=bins)
        ax.set_title(c, fontsize=9)
        ax.set_xlabel(c, fontsize=8)
        ax.set_ylabel("count", fontsize=8)

    # hide unused subplots
    total = nrows * ncols
    for j in range(n, total):
        r, k = divmod(j, ncols)
        axes[r, k].axis("off")

    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def grouped_hist_specs() -> Dict[str, List[str]]:
    # Only includes columns if they exist in the dataframe at runtime.
    return {
        "esm_delta_pll__esm1v_esm2_esm3.png": [
            "esm1v_delta_pll", "esm2_delta_pll", "esm3_delta_pll"
        ],
        "esm_mutation_llr__esm1v_esm2_esm3.png": [
            "esm1v_mutation_llr", "esm2_mutation_llr", "esm3_mutation_llr"
        ],
        "esm_pseudo_likelihood__esm1v_esm2_esm3.png": [
            "esm1v_pseudo_likelihood", "esm2_pseudo_likelihood", "esm3_pseudo_likelihood"
        ],
        "esm_pll_wt__esm1v_esm2_esm3.png": [
            "esm1v_pll_wt", "esm2_pll_wt", "esm3_pll_wt"
        ],
        "llr_bundle.png": [
            # per-model + consensus summaries (keep short; remove/add as you prefer)
            "esm1v_mutation_llr", "esm2_mutation_llr", "esm3_mutation_llr",
            "llr_mean", "llr_std", "llr_range",
        ],
        "rank_bundle.png": [
            # ranks produced by ranks_basic and/or scorer outputs
            "rank_llr_mean", "rank_delta_pll_mean",
            "rank_llr_median", "rank_delta_pll_median",
            "rank_S1", "rank_S2",  # if you name score cols S1/S2 later
        ],
    }


def grouped_hist_exclusion_set() -> set:
    ex = set()
    for _, cols in grouped_hist_specs().items():
        ex |= set(cols)
    return ex

def save_histograms(
    df: pd.DataFrame,
    cols: List[str],
    outdir: str,
    max_cols: int = 80,
    exclude: Optional[set] = None,
) -> None:
    mkdirp(outdir)
    use = cols[:max_cols]
    exclude = exclude or set()

    for c in use:
        if c in exclude:
            continue
        x = to_num(df[c]).dropna()
        if len(x) < 50:
            continue
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.hist(x.values, bins=60)
        ax.set_title(c)
        ax.set_xlabel(c)
        ax.set_ylabel("count")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{c}.png"), dpi=200)
        plt.close(fig)


def top_corrs_with_score(df: pd.DataFrame, feature_cols: List[str], score_col: str, topn: int = 25) -> pd.DataFrame:
    s = to_num(df[score_col])
    out = []
    for c in feature_cols:
        x = to_num(df[c])
        m = x.notna() & s.notna()
        if int(m.sum()) < 200:
            continue
        r = float(np.corrcoef(x[m].values, s[m].values)[0, 1])
        if np.isfinite(r):
            out.append((c, r, int(m.sum())))
    res = pd.DataFrame(out, columns=["feature", "pearson_r", "n"]).sort_values("pearson_r", ascending=False)
    if len(res) == 0:
        return res
    top_pos = res.head(topn)
    top_neg = res.tail(topn).sort_values("pearson_r", ascending=True)
    return pd.concat([top_pos, top_neg], axis=0).drop_duplicates("feature").reset_index(drop=True)


def save_scatter(df: pd.DataFrame, xcol: str, ycol: str, outpath: str, nmax: int = 15000) -> None:
    x = to_num(df[xcol])
    y = to_num(df[ycol])
    m = x.notna() & y.notna()
    if int(m.sum()) < 200:
        return
    idx = np.where(m.values)[0]
    if len(idx) > nmax:
        idx = np.random.default_rng(0).choice(idx, size=nmax, replace=False)
    fig = plt.figure(figsize=(5.5, 4.5))
    ax = fig.add_subplot(111)
    ax.scatter(x.iloc[idx].values, y.iloc[idx].values, s=6, alpha=0.35)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(f"{ycol} vs {xcol}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def compute_topk_overlap(df: pd.DataFrame, id_col: str, score_cols: List[str], topk: int) -> pd.DataFrame:
    sets: Dict[str, set] = {}
    for sc in score_cols:
        s = to_num(df[sc])
        tmp = df[[id_col]].copy()
        tmp["_score"] = s
        tmp = tmp.dropna().sort_values("_score", ascending=False).head(topk)
        sets[sc] = set(tmp[id_col].astype(str).values)

    rows = []
    for i, a in enumerate(score_cols):
        for b in score_cols[i:]:
            A, B = sets.get(a, set()), sets.get(b, set())
            inter = len(A & B)
            union = len(A | B) if (A or B) else 0
            jac = (inter / union) if union else np.nan
            rows.append((a, b, inter, union, jac))
    return pd.DataFrame(rows, columns=["score_a", "score_b", "intersection", "union", "jaccard"])

def compute_topk_support_table(df: pd.DataFrame, id_col: str, score_cols: List[str], topk: int) -> pd.DataFrame:
    """
    For each score_col, take topK IDs (descending), then compute per-variant support_count
    = how many of the rankings include that variant in their topK.
    Returns a table of all variants that appear in >=1 topK list.
    """
    top_sets: Dict[str, set] = {}
    for sc in score_cols:
        s = to_num(df[sc])
        tmp = df[[id_col]].copy()
        tmp["_score"] = s
        tmp = tmp.dropna().sort_values("_score", ascending=False).head(topk)
        top_sets[sc] = set(tmp[id_col].astype(str).values)

    universe = sorted(set().union(*top_sets.values())) if top_sets else []
    rows = []
    for vid in universe:
        support = sum(1 for sc in score_cols if vid in top_sets[sc])
        rows.append((vid, support))

    out = pd.DataFrame(rows, columns=[id_col, "support_count"]).sort_values(
        ["support_count", id_col], ascending=[False, True]
    )

    # Add boolean membership columns (useful for inspection)
    for sc in score_cols:
        out[f"in_top{topk}__{sc}"] = out[id_col].astype(str).isin(top_sets[sc])

    return out


def save_support_hist(support_df: pd.DataFrame, outpath: str, topk: int, n_scores: int) -> None:
    """
    Histogram of support_count (1..n_scores).
    """
    if support_df is None or len(support_df) == 0:
        return
    vc = support_df["support_count"].value_counts().reindex(range(1, n_scores + 1), fill_value=0)

    fig = plt.figure(figsize=(7, 4))
    ax = fig.add_subplot(111)
    ax.bar(vc.index.astype(int), vc.values)
    ax.set_xlabel(f"support_count = # rankings where variant is in top {topk}")
    ax.set_ylabel("number of variants")
    ax.set_title(f"Top-{topk} consensus across {n_scores} rankings")
    ax.set_xticks(list(range(1, n_scores + 1)))
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def save_support_upset_like(support_df: pd.DataFrame, score_cols: List[str], outpath: str, topk: int, max_sets: int = 40) -> None:
    """
    Simple 'upset-like' bar chart: show the largest intersections (exact membership patterns).
    This avoids external libraries. If too many patterns, show top max_sets.
    """
    if support_df is None or len(support_df) == 0:
        return

    mem_cols = [f"in_top{topk}__{sc}" for sc in score_cols]
    for c in mem_cols:
        if c not in support_df.columns:
            return

    # Pattern key like "101010"
    pat = support_df[mem_cols].astype(int).astype(str).agg("".join, axis=1)
    vc = pat.value_counts().head(max_sets)

    fig = plt.figure(figsize=(10, 4.5))
    ax = fig.add_subplot(111)
    ax.bar(range(len(vc)), vc.values)
    ax.set_ylabel("number of variants")
    ax.set_xlabel("intersection pattern (in topK across score_cols)")
    ax.set_title(f"Top-{topk} intersections (largest {min(max_sets, len(vc))})")
    ax.set_xticks(range(len(vc)))
    ax.set_xticklabels(vc.index.tolist(), rotation=90, fontsize=7)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    
def write_top_candidates(df: pd.DataFrame, id_col: str, group_col: Optional[str], score_cols: List[str], topk: int, outpath: str) -> None:
    # choose mutation column if present
    mut_col = None
    for c in ["_mutation_upper", "mutation", "mut_code", "mutation_code"]:
        if c in df.columns:
            mut_col = c
            break

    keep = [id_col]
    if mut_col:
        keep.append(mut_col)
    keep += score_cols

    keep = safe_cols(df, keep)
    out = df[keep].copy()
    for sc in score_cols:
        out[sc] = to_num(out[sc])
    parts = []
    for sc in score_cols:
        tmp = out.dropna(subset=[sc]).sort_values(sc, ascending=False).head(topk).copy()
        tmp["score_col"] = sc
        parts.append(tmp)
    if parts:
        pd.concat(parts, axis=0).to_csv(outpath, sep="\t", index=False)


def score_vs_label_metrics(df: pd.DataFrame, score_col: str, label_col: str, topk: int, label_higher_is_better: bool) -> Dict[str, float]:
    s = to_num(df[score_col])
    y = to_num(df[label_col])
    m = s.notna() & y.notna()
    n = int(m.sum())
    if n < 50:
        return {"n": float(n), "pearson": np.nan, "spearman": np.nan, "topk_jaccard": np.nan}

    pearson = float(np.corrcoef(s[m].values, y[m].values)[0, 1])
    spearman = float(pd.Series(s[m].values).corr(pd.Series(y[m].values), method="spearman"))

    # TopK jaccard: top by score vs top by label
    tmp = df.loc[m, :].copy()
    tmp["_s"] = s[m].values
    tmp["_y"] = y[m].values
    top_s = set(tmp.sort_values("_s", ascending=False).head(topk).index.values.tolist())
    top_y = set(tmp.sort_values("_y", ascending=not label_higher_is_better).head(topk).index.values.tolist())
    inter = len(top_s & top_y)
    union = len(top_s | top_y) if (top_s or top_y) else 0
    jac = (inter / union) if union else np.nan

    return {"n": float(n), "pearson": pearson, "spearman": spearman, "topk_jaccard": float(jac)}


def per_group_metrics(df: pd.DataFrame, group_col: str, score_cols: List[str], label_cols: List[str], topk: int, label_higher_is_better: bool) -> pd.DataFrame:
    rows = []
    for g, sub in df.groupby(group_col):
        for sc in score_cols:
            for lc in label_cols:
                met = score_vs_label_metrics(sub, sc, lc, topk=topk, label_higher_is_better=label_higher_is_better)
                rows.append((str(g), sc, lc, int(met["n"]), met["pearson"], met["spearman"], met["topk_jaccard"]))
    return pd.DataFrame(rows, columns=["group", "score_col", "label_col", "n", "pearson", "spearman", "topk_jaccard"])


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (features or ranks table)")
    ap.add_argument("--outdir", default="phase_dashboard", help="Output directory")
    ap.add_argument("--id-col", default=None, help="ID column (default: infer)")
    ap.add_argument("--group-col", default=None, help="Group column for per-backbone summaries (default: infer if present)")
    ap.add_argument("--score-cols", default="", help="Comma list of score columns to analyze (optional)")
    ap.add_argument("--label-cols", default="", help="Comma list of label/ground-truth columns (optional)")
    ap.add_argument("--label-higher-is-better", action="store_true", help="If set, treat labels as higher=better (default: higher=better).")
    ap.add_argument("--label-lower-is-better", action="store_true", help="If set, treat labels as lower=better.")
    ap.add_argument("--feature-cols", default="", help="Comma list of feature cols (optional; default infer numeric)")
    ap.add_argument("--topk", type=int, default=200, help="TopK for overlap / candidate tables / evaluation overlap")
    ap.add_argument("--max-hist", type=int, default=80, help="Max histograms to write")
    ap.add_argument("--max-scatter-per-score", type=int, default=12, help="Max scatterplots per score")
    ap.add_argument("--max-scatter-per-label", type=int, default=12, help="Max score-vs-label scatters per label")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    mkdirp(args.outdir)
    mkdirp(os.path.join(args.outdir, "plots"))
    mkdirp(os.path.join(args.outdir, "plots", "hist"))
    mkdirp(os.path.join(args.outdir, "plots", "scatter"))
    mkdirp(os.path.join(args.outdir, "tables"))

    df = pd.read_csv(args.inp, sep="\t")
    id_col = infer_id_col(df, args.id_col)
    group_col = infer_group_col(df, args.group_col)

    score_cols = [c.strip() for c in args.score_cols.split(",") if c.strip()]
    score_cols = safe_cols(df, score_cols)

    label_cols = [c.strip() for c in args.label_cols.split(",") if c.strip()]
    label_cols = safe_cols(df, label_cols)

    label_higher_is_better = True
    if args.label_lower_is_better:
        label_higher_is_better = False
    if args.label_higher_is_better:
        label_higher_is_better = True

    user_feature_cols = [c.strip() for c in args.feature_cols.split(",") if c.strip()] or None
    feature_cols = infer_numeric_feature_cols(
        df,
        id_col=id_col,
        group_col=group_col,
        score_cols=score_cols,
        label_cols=label_cols,
        user_feature_cols=user_feature_cols,
    )

    summary = {
        "input": args.inp,
        "n_rows": int(len(df)),
        "id_col": id_col,
        "group_col": group_col,
        "score_cols": score_cols,
        "label_cols": label_cols,
        "label_higher_is_better": bool(label_higher_is_better),
        "n_feature_cols": int(len(feature_cols)),
    }
    if group_col:
        summary["n_groups"] = int(df[group_col].astype(str).nunique())

    # Missingness (features + scores + labels)
    miss_cols = feature_cols + score_cols + label_cols
    miss = summarize_missingness(df, miss_cols)
    miss.to_csv(os.path.join(args.outdir, "tables", "missingness.tsv"), sep="\t", index=False)
    summary["top_missing_cols"] = miss.head(15).to_dict(orient="records")

    # Feature correlation heatmap (subset)
    corr_cols = feature_cols[:80]
    if len(corr_cols) >= 10:
        corr = corr_matrix(df, corr_cols)
        save_corr_heatmap(
            corr,
            os.path.join(args.outdir, "plots", "corr_heatmap_features.png"),
            title="Feature correlation heatmap (subset)",
        )

    # Grouped histograms (reduce clutter)
    hist_dir = os.path.join(args.outdir, "plots", "hist")
    for fname, cols in grouped_hist_specs().items():
        outp = os.path.join(hist_dir, fname)
        save_hist_grid(df, cols, outp, title=fname.replace(".png", ""), bins=60)
    # Individual histograms (features), excluding the grouped columns
    exclude = grouped_hist_exclusion_set()
    save_histograms(df, feature_cols, hist_dir, max_cols=args.max_hist, exclude=exclude)

    # Score-centric diagnostics (no labels needed)
    if score_cols:
        if len(score_cols) >= 2:
            ov = compute_topk_overlap(df, id_col=id_col, score_cols=score_cols, topk=args.topk)
            ov.to_csv(os.path.join(args.outdir, "tables", "topk_overlap__scores.tsv"), sep="\t", index=False)
            summary["topk_overlap_scores"] = ov.to_dict(orient="records")

        write_top_candidates(
            df,
            id_col=id_col,
            group_col=group_col,
            score_cols=score_cols,
            topk=args.topk,
            outpath=os.path.join(args.outdir, "tables", "top_candidates__by_score.tsv"),
        )
        # --- NEW: consensus across topK for the provided score_cols ---
        support_df = compute_topk_support_table(df, id_col=id_col, score_cols=score_cols, topk=args.topk)
        support_df.to_csv(os.path.join(args.outdir, "tables", f"topk_support_counts__top{args.topk}.tsv"), sep="\t", index=False)

        save_support_hist(
            support_df,
            outpath=os.path.join(args.outdir, "plots", f"topk_support_hist__top{args.topk}.png"),
            topk=args.topk,
            n_scores=len(score_cols),
        )

        save_support_upset_like(
            support_df,
            score_cols=score_cols,
            outpath=os.path.join(args.outdir, "plots", f"topk_upset__top{args.topk}.png"),
            topk=args.topk,
            max_sets=40,
        )
        # Score-feature correlations + scatters
        for sc in score_cols:
            tc = top_corrs_with_score(df, feature_cols=feature_cols, score_col=sc, topn=25)
            tc.to_csv(os.path.join(args.outdir, "tables", f"score_corrs__{sc}.tsv"), sep="\t", index=False)
            if len(tc) > 0:
                tc2 = tc.copy()
                tc2["abs_r"] = tc2["pearson_r"].abs()
                tc2 = tc2.sort_values("abs_r", ascending=False).head(args.max_scatter_per_score)
                for _, row in tc2.iterrows():
                    fx = row["feature"]
                    outp = os.path.join(args.outdir, "plots", "scatter", f"{sc}__vs__{fx}.png")
                    save_scatter(df, xcol=fx, ycol=sc, outpath=outp)

    # Label evaluation (only if labels exist)
    if score_cols and label_cols:
        rows = []
        for sc in score_cols:
            for lc in label_cols:
                met = score_vs_label_metrics(
                    df,
                    score_col=sc,
                    label_col=lc,
                    topk=args.topk,
                    label_higher_is_better=label_higher_is_better,
                )
                rows.append((sc, lc, int(met["n"]), met["pearson"], met["spearman"], met["topk_jaccard"]))
        eval_df = pd.DataFrame(rows, columns=["score_col", "label_col", "n", "pearson", "spearman", "topk_jaccard"])
        eval_df.to_csv(os.path.join(args.outdir, "tables", "eval__score_vs_label.tsv"), sep="\t", index=False)
        summary["eval_score_vs_label"] = eval_df.to_dict(orient="records")

        # Score vs label scatters for the strongest |spearman|
        for lc in label_cols:
            sub = eval_df[eval_df["label_col"] == lc].copy()
            sub["abs_spear"] = sub["spearman"].abs()
            sub = sub.sort_values("abs_spear", ascending=False).head(args.max_scatter_per_label)
            for _, r in sub.iterrows():
                sc = r["score_col"]
                outp = os.path.join(args.outdir, "plots", "scatter", f"{lc}__vs__{sc}.png")
                save_scatter(df, xcol=sc, ycol=lc, outpath=outp)

        # Per-group evaluation (if group_col present)
        if group_col:
            gdf = per_group_metrics(
                df,
                group_col=group_col,
                score_cols=score_cols,
                label_cols=label_cols,
                topk=args.topk,
                label_higher_is_better=label_higher_is_better,
            )
            gdf.to_csv(os.path.join(args.outdir, "tables", "eval__score_vs_label__by_group.tsv"), sep="\t", index=False)

    # Group counts
    if group_col:
        g = df[group_col].astype(str)
        gtab = pd.DataFrame({"group": g.value_counts().index, "n": g.value_counts().values})
        gtab.to_csv(os.path.join(args.outdir, "tables", "group_counts.tsv"), sep="\t", index=False)

    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote dashboard to {args.outdir}")
    print(f"[rows] {len(df)}")
    print(f"[features] {len(feature_cols)}")
    print(f"[scores] {len(score_cols)}")
    print(f"[labels] {len(label_cols)}")


if __name__ == "__main__":
    main()