#!/usr/bin/env python3
# add_wt_seq_to_testset_db.py
#
# Adds a wt_seq column to testset_db.tsv by mapping wt_id -> WT sequence
# from tournament_wt.fasta.
#
# Usage:
#   python add_wt.py --db testset_db.tsv --wt tournament_wt.fasta --out testset_db_with_wtseq.tsv
#
# Assumptions:
#   - testset_db.tsv is tab-separated and has a column named: wt_id
#   - tournament_wt.fasta headers contain the wt_id as the first token after '>'

import argparse
import pandas as pd

def load_fasta_dict(fasta_path: str) -> dict:
    d = {}
    with open(fasta_path, "r", encoding="utf-8") as f:
        cur_id = None
        buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    d[cur_id] = "".join(buf)
                cur_id = line[1:].split()[0]  # first token = ID
                buf = []
            else:
                buf.append(line)
        if cur_id is not None:
            d[cur_id] = "".join(buf)
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="Input TSV (testset_db.tsv)")
    ap.add_argument("--wt", required=True, help="WT FASTA (tournament_wt.fasta)")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    df = pd.read_csv(args.db, sep="\t", dtype=str)

    if "wt_id" not in df.columns:
        raise SystemExit(f"ERROR: input TSV missing required column 'wt_id'. Columns={list(df.columns)}")

    wt_map = load_fasta_dict(args.wt)

    df["wt_seq"] = df["wt_id"].map(wt_map)

    missing = df["wt_seq"].isna().sum()
    if missing:
        bad_ids = df.loc[df["wt_seq"].isna(), "wt_id"].unique().tolist()[:20]
        raise SystemExit(
            f"ERROR: {missing} rows could not map wt_id -> wt_seq from FASTA.\n"
            f"Example missing wt_id values: {bad_ids}\n"
            f"Fix: ensure tournament_wt.fasta headers match the wt_id strings exactly."
        )

    df.to_csv(args.out, sep="\t", index=False)
    print(f"OK: wrote {len(df)} rows with wt_seq -> {args.out}")

if __name__ == "__main__":
    main()