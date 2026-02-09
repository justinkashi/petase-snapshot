#!/usr/bin/env python3
# evc_writer.py: write EVcouplings features (sum/mean + coverage) for test variants from single_mutant_matrix.csv

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd


AA20 = set(list("ACDEFGHIKLMNPQRSTVWY"))


def eprint(*a):
    print(*a, file=sys.stderr)


def parse_mutation_string(mut_str: str):
    """
    Parse mutation strings like:
      - "V26A"
      - "V26A;Q29G"
      - "V26A,Q29G"
      - "V26A Q29G"
    Returns list of tokens like ["V26A", "Q29G"].
    """
    if mut_str is None:
        return []
    s = str(mut_str).strip()
    if s == "" or s.lower() in {"nan", "none"}:
        return []
    if s.upper() == "WT":
        return ["WT"]

    toks = re.split(r"[;, \t]+", s)
    toks = [t.strip() for t in toks if t.strip() != ""]
    return toks


def parse_mut_token(tok: str):
    """
    Parse single mutation token like "V26A" -> (pos:int, wt:str, sub:str).
    """
    m = re.fullmatch(r"([A-Za-z])(\d+)([A-Za-z])", tok.strip())
    if not m:
        raise ValueError(f"Bad mutation token: '{tok}'")
    wt = m.group(1).upper()
    pos = int(m.group(2))
    sub = m.group(3).upper()
    if wt not in AA20 or sub not in AA20:
        raise ValueError(f"Non-AA20 mutation token: '{tok}'")
    return pos, wt, sub


def find_run_dir_for_wt(evc_root: str, wt_id: str, target_glob: str):
    """
    Default mapping for your setup:
      tournament_wt_1 -> evc_root/evc_wt1/<TARGET>
      tournament_wt_2 -> evc_root/evc_wt2/<TARGET>
      tournament_wt_3 -> evc_root/evc_wt3/<TARGET>
    """
    m = re.fullmatch(r"tournament_wt_(\d+)", str(wt_id).strip())
    if not m:
        return None
    idx = int(m.group(1))
    if idx not in {1, 2, 3}:
        return None

    run_parent = Path(evc_root) / f"evc_wt{idx}"
    if not run_parent.exists():
        return None

    # target_glob like "TARGET_b0.3" or "TARGET_b*"
    hits = sorted(run_parent.glob(target_glob))
    if len(hits) == 0:
        return None
    if len(hits) > 1:
        # prefer exact match if exists
        exact = [h for h in hits if h.name == target_glob]
        if len(exact) == 1:
            return str(exact[0])
        # else pick first deterministically
    return str(hits[0])


def load_single_mutant_lut(run_dir: str):
    """
    Loads mutate/<target>_single_mutant_matrix.csv and builds lookup:
      key = (pos, sub) -> dict with epistatic/independent/frequency/colcon
    """
    run_dir = str(run_dir)
    mutate_dir = Path(run_dir) / "mutate"
    if not mutate_dir.exists():
        raise FileNotFoundError(f"Missing mutate dir: {mutate_dir}")

    csvs = sorted(mutate_dir.glob("*_single_mutant_matrix.csv"))
    if len(csvs) == 0:
        raise FileNotFoundError(f"No *_single_mutant_matrix.csv in {mutate_dir}")

    csv_path = str(csvs[0])
    df = pd.read_csv(csv_path)

    required = [
        "pos",
        "subs",
        "prediction_epistatic",
        "prediction_independent",
        "frequency",
        "column_conservation",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {csv_path}: {missing}")

    lut = {}
    for _, r in df.iterrows():
        pos = int(r["pos"])
        sub = str(r["subs"]).strip().upper()
        if sub not in AA20:
            continue
        lut[(pos, sub)] = {
            "epistatic": float(r["prediction_epistatic"]),
            "independent": float(r["prediction_independent"]),
            "freq": float(r["frequency"]),
            "colcon": float(r["column_conservation"]),
        }

    return lut, csv_path


def safe_mean(xsum, n):
    return (xsum / n) if (n is not None and n > 0) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test_db_tsv", required=True, help="TSV with test variants and mutation strings")
    ap.add_argument("--evc_root", required=True, help="Root folder containing evc_wt1/ evc_wt2/ evc_wt3")
    ap.add_argument("--target_glob", default="TARGET_b0.3", help='Run dir glob under evc_wt*/ (e.g. "TARGET_b0.3")')
    ap.add_argument("--out_tsv", required=True, help="Output TSV of EVC features")

    ap.add_argument("--wt_id_col", default="wt_id", help="Column containing WT id (e.g. tournament_wt_1)")
    ap.add_argument("--test_id_col", default="test_id", help="Column containing test variant id")
    ap.add_argument("--mutation_col", default="mutation", help="Column containing mutation string (e.g. V26A;Q29G)")

    ap.add_argument("--agg", default="sum", choices=["sum"], help="Aggregation mode (currently only sum)")

    args = ap.parse_args()

    eprint(f"[EVC_WRITER] Reading test DB: {args.test_db_tsv}")
    df = pd.read_csv(args.test_db_tsv, sep="\t", dtype=str)
    for c in [args.test_id_col, args.wt_id_col, args.mutation_col]:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in {args.test_db_tsv}. Have: {list(df.columns)}")

    evc_root = str(args.evc_root)

    # Determine which WT ids we need, but only keep wt_1..wt_3
    wt_ids_all = sorted(set(df[args.wt_id_col].astype(str).tolist()))
    wt_ids = [w for w in wt_ids_all if re.fullmatch(r"tournament_wt_[123]", str(w).strip())]
    if len(wt_ids) == 0:
        raise ValueError(f"No WT ids tournament_wt_1/2/3 found in {args.wt_id_col}")

    # Load LUTs for each WT
    wt_to_run = {}
    wt_to_lut = {}

    for wt_id in wt_ids:
        run_dir = find_run_dir_for_wt(evc_root, wt_id, args.target_glob)
        if run_dir is None or not Path(run_dir).exists():
            eprint(f"[EVC_WRITER][ERROR] Could not find EVcouplings run dir for {wt_id}: expected under {evc_root}/evc_wt*/{args.target_glob}")
            sys.exit(2)

        lut, csv_path = load_single_mutant_lut(run_dir)
        wt_to_run[wt_id] = run_dir
        wt_to_lut[wt_id] = lut
        eprint(f"[EVC_WRITER] Loading LUT for {wt_id} from {csv_path}")

    # Compute features
    rows = []
    n_unknown_wt = 0
    n_empty_mut = 0
    n_bad_mut = 0

    for _, r in df.iterrows():
        test_id = str(r[args.test_id_col])
        wt_id = str(r[args.wt_id_col])
        mut_str = r[args.mutation_col]

        muts = parse_mutation_string(mut_str)

        if len(muts) == 0:
            n_empty_mut += 1
            rows.append(
                {
                    args.test_id_col: test_id,
                    args.wt_id_col: wt_id,
                    args.mutation_col: mut_str if mut_str is not None else "",
                    "evc_n_muts": 0,
                    "evc_n_found": 0,
                    "evc_frac_found": 0.0,
                    "evc_missing_frac": 1.0,
                    "evc_epistatic_sum": None,
                    "evc_independent_sum": None,
                    "evc_freq_sum": None,
                    "evc_colcon_sum": None,
                    "evc_epistatic_mean": None,
                    "evc_independent_mean": None,
                }
            )
            continue

        if muts == ["WT"]:
            n_empty_mut += 1
            rows.append(
                {
                    args.test_id_col: test_id,
                    args.wt_id_col: wt_id,
                    args.mutation_col: "WT",
                    "evc_n_muts": 0,
                    "evc_n_found": 0,
                    "evc_frac_found": 0.0,
                    "evc_missing_frac": 1.0,
                    "evc_epistatic_sum": None,
                    "evc_independent_sum": None,
                    "evc_freq_sum": None,
                    "evc_colcon_sum": None,
                    "evc_epistatic_mean": None,
                    "evc_independent_mean": None,
                }
            )
            continue

        if wt_id not in wt_to_lut:
            n_unknown_wt += 1
            rows.append(
                {
                    args.test_id_col: test_id,
                    args.wt_id_col: wt_id,
                    args.mutation_col: mut_str if mut_str is not None else "",
                    "evc_n_muts": len(muts),
                    "evc_n_found": 0,
                    "evc_frac_found": 0.0,
                    "evc_missing_frac": 1.0,
                    "evc_epistatic_sum": None,
                    "evc_independent_sum": None,
                    "evc_freq_sum": None,
                    "evc_colcon_sum": None,
                    "evc_epistatic_mean": None,
                    "evc_independent_mean": None,
                }
            )
            continue

        lut = wt_to_lut[wt_id]

        ep_sum = 0.0
        indep_sum = 0.0
        freq_sum = 0.0
        colcon_sum = 0.0
        n_found = 0

        bad = False
        for tok in muts:
            try:
                pos, _wt, sub = parse_mut_token(tok)
            except Exception as ex:
                eprint(f"[EVC_WRITER][WARN] Bad mutation string for {test_id}: '{mut_str}' ({ex})")
                n_bad_mut += 1
                bad = True
                break

            key = (pos, sub)
            if key in lut:
                n_found += 1
                ep_sum += lut[key]["epistatic"]
                indep_sum += lut[key]["independent"]
                freq_sum += lut[key]["freq"]
                colcon_sum += lut[key]["colcon"]

        if bad:
            rows.append(
                {
                    args.test_id_col: test_id,
                    args.wt_id_col: wt_id,
                    args.mutation_col: mut_str if mut_str is not None else "",
                    "evc_n_muts": len(muts),
                    "evc_n_found": 0,
                    "evc_frac_found": 0.0,
                    "evc_missing_frac": 1.0,
                    "evc_epistatic_sum": None,
                    "evc_independent_sum": None,
                    "evc_freq_sum": None,
                    "evc_colcon_sum": None,
                    "evc_epistatic_mean": None,
                    "evc_independent_mean": None,
                }
            )
            continue

        n_muts = len(muts)
        frac = (n_found / n_muts) if n_muts > 0 else 0.0
        missing_frac = 1.0 - frac

        if n_found == 0:
            rows.append(
                {
                    args.test_id_col: test_id,
                    args.wt_id_col: wt_id,
                    args.mutation_col: mut_str if mut_str is not None else "",
                    "evc_n_muts": n_muts,
                    "evc_n_found": 0,
                    "evc_frac_found": 0.0,
                    "evc_missing_frac": 1.0,
                    "evc_epistatic_sum": None,
                    "evc_independent_sum": None,
                    "evc_freq_sum": None,
                    "evc_colcon_sum": None,
                    "evc_epistatic_mean": None,
                    "evc_independent_mean": None,
                }
            )
            continue

        rows.append(
            {
                args.test_id_col: test_id,
                args.wt_id_col: wt_id,
                args.mutation_col: mut_str if mut_str is not None else "",
                "evc_n_muts": n_muts,
                "evc_n_found": n_found,
                "evc_frac_found": frac,
                "evc_missing_frac": missing_frac,
                "evc_epistatic_sum": ep_sum,
                "evc_independent_sum": indep_sum,
                "evc_freq_sum": freq_sum,
                "evc_colcon_sum": colcon_sum,
                "evc_epistatic_mean": safe_mean(ep_sum, n_muts),
                "evc_independent_mean": safe_mean(indep_sum, n_muts),
            }
        )

    out_df = pd.DataFrame(rows)

    # Ensure stable column order
    out_cols = [
        args.test_id_col,
        args.wt_id_col,
        args.mutation_col,
        "evc_n_muts",
        "evc_n_found",
        "evc_frac_found",
        "evc_missing_frac",
        "evc_epistatic_sum",
        "evc_independent_sum",
        "evc_freq_sum",
        "evc_colcon_sum",
        "evc_epistatic_mean",
        "evc_independent_mean",
    ]
    out_df = out_df[out_cols]

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out_tsv, sep="\t", index=False)

    eprint(f"[EVC_WRITER] WROTE: {args.out_tsv}")
    eprint(f"[EVC_WRITER] Rows: {len(out_df)}")
    eprint(f"[EVC_WRITER] Unknown WT rows: {n_unknown_wt}")
    eprint(f"[EVC_WRITER] Empty mutation rows: {n_empty_mut}")
    eprint(f"[EVC_WRITER] Bad mutation strings: {n_bad_mut}")

    # Quick coverage stats
    cov = out_df["evc_frac_found"].dropna()
    if len(cov) > 0:
        eprint(f"[EVC_WRITER] evc_frac_found: mean={cov.mean():.4f} median={cov.median():.4f} n={len(cov)}")

    eprint("[EVC_WRITER] WT run dirs:")
    for k in sorted(wt_to_run.keys()):
        eprint(f"[EVC_WRITER]   {k} -> {wt_to_run[k]}")


if __name__ == "__main__":
    main()