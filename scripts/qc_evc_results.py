#!/usr/bin/env python3
# qc_evc_results.py: QC EVcouplings result folders under evc_results/*/TARGET_b*/ and write a TSV summary.

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def count_globs(base: Path, patterns: List[str]) -> int:
    n = 0
    for pat in patterns:
        n += len(list(base.rglob(pat)))
    return n


def pick_first(base: Path, patterns: List[str]) -> str:
    for pat in patterns:
        hits = sorted(base.rglob(pat))
        if hits:
            return str(hits[0])
    return ""


def qc_one_run(run_dir: Path) -> Dict:
    align_dir = run_dir / "align"
    couplings_dir = run_dir / "couplings"
    mutate_dir = run_dir / "mutate"
    compare_dir = run_dir / "compare"

    align_file_count = 0
    if align_dir.exists():
        align_file_count = count_globs(
            align_dir,
            patterns=[
                "*.sto",
                "*.a2m",
                "*.fa",
                "*.fasta",
                "*.aln",
            ],
        )

    couplings_file_count = 0
    if couplings_dir.exists():
        couplings_file_count = count_globs(
            couplings_dir,
            patterns=[
                "*.pkl",
                "*.npz",
                "*.csv",
                "*.json",
                "*.txt",
            ],
        )

    mutate_file_count = 0
    if mutate_dir.exists():
        mutate_file_count = count_globs(
            mutate_dir,
            patterns=[
                "*.csv",
                "*.tsv",
                "*.txt",
                "*.npz",
            ],
        )

    # Try to find representative files (helpful for debugging)
    align_example = pick_first(align_dir, ["*.sto", "*.a2m", "*.fa", "*.fasta"]) if align_dir.exists() else ""
    couplings_example = pick_first(couplings_dir, ["*.pkl", "*.npz", "*.csv"]) if couplings_dir.exists() else ""
    mutate_example = pick_first(mutate_dir, ["*.csv", "*.tsv", "*.txt"]) if mutate_dir.exists() else ""

    ok_align = align_dir.exists() and align_file_count > 0
    ok_couplings = couplings_dir.exists() and couplings_file_count > 0
    ok_mutate = mutate_dir.exists() and mutate_file_count > 0

    notes = []
    if not run_dir.exists():
        notes.append("MISSING run_dir")
    if not align_dir.exists():
        notes.append("MISSING align/")
    if align_dir.exists() and align_file_count == 0:
        notes.append("EMPTY align/")
    if not couplings_dir.exists():
        notes.append("MISSING couplings/")
    if couplings_dir.exists() and couplings_file_count == 0:
        notes.append("EMPTY couplings/")
    if not mutate_dir.exists():
        notes.append("MISSING mutate/")
    if mutate_dir.exists() and mutate_file_count == 0:
        notes.append("EMPTY mutate/")

    return {
        "run_dir": str(run_dir),
        "ok": bool(ok_align and ok_couplings and ok_mutate),
        "ok_align": bool(ok_align),
        "ok_couplings": bool(ok_couplings),
        "ok_mutate": bool(ok_mutate),
        "has_compare": bool(compare_dir.exists()),
        "align_files": int(align_file_count),
        "couplings_files": int(couplings_file_count),
        "mutate_files": int(mutate_file_count),
        "align_example": align_example,
        "couplings_example": couplings_example,
        "mutate_example": mutate_example,
        "notes": "; ".join(notes) if notes else "",
    }


def find_runs(root: Path, target_glob: str) -> List[Path]:
    # Your structure: evc_results/evc_wt1/TARGET_b0.3/...
    # So we search for folders named TARGET_b* (or whatever you pass).
    runs = []
    for p in sorted(root.rglob(target_glob)):
        if p.is_dir():
            # must contain at least one expected stage folder to be considered a run
            if (p / "align").exists() or (p / "couplings").exists() or (p / "mutate").exists():
                runs.append(p)
    return runs


def main():
    ap = argparse.ArgumentParser(description="QC EVcouplings results folder tree.")
    ap.add_argument("--root", default="evc_results", help="Root folder containing evcouplings runs (default: evc_results)")
    ap.add_argument("--target_glob", default="TARGET_b*", help="Run folder name glob (default: TARGET_b*)")
    ap.add_argument("--out_tsv", default="", help="Output TSV path (default: <root>/evc_qc.tsv)")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"ERROR: root not found: {root}")

    runs = find_runs(root, args.target_glob)
    if not runs:
        raise SystemExit(f"ERROR: no run directories found under {root} matching {args.target_glob}")

    rows = [qc_one_run(r) for r in runs]
    df = pd.DataFrame(rows)

    # sort: failures first
    df = df.sort_values(["ok", "run_dir"], ascending=[True, True])

    out_tsv = Path(args.out_tsv) if args.out_tsv else (root / "evc_qc.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    # print compact summary
    n_ok = int(df["ok"].sum())
    n_total = df.shape[0]
    print(f"WROTE: {out_tsv}")
    print(f"OK: {n_ok}/{n_total}")
    print(df[["ok", "run_dir", "align_files", "couplings_files", "mutate_files", "has_compare", "notes"]].to_string(index=False))


if __name__ == "__main__":
    main()