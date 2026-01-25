#!/usr/bin/env python3
# Usage: python esm3_writer_mac.py --in_fasta data/test.fasta --out_tsv out/esm3_mac.tsv
# Requires your ESM3 install + weights. This is a wrapper skeleton.

import argparse
import sys
from pathlib import Path

import pandas as pd
import torch


def pick_device(device: str) -> torch.device:
    device = device.lower()
    if device == "auto":
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    if device == "mps":
        if not torch.backends.mps.is_available():
            print("[esm3_writer_mac][WARN] MPS not available; falling back to CPU", file=sys.stderr)
            return torch.device("cpu")
        return torch.device("mps")
    if device == "cpu":
        return torch.device("cpu")
    raise ValueError(f"Invalid --device {device}. Use auto|mps|cpu")


def read_fasta(path: str):
    records = []
    cur_id = None
    cur_seq = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            records.append((cur_id, "".join(cur_seq)))
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_fasta", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--id_col", default="seq_id")
    ap.add_argument("--device", default="auto", choices=["auto", "mps", "cpu"])
    ap.add_argument("--batch_size", type=int, default=1)
    ap.add_argument("--model", default="esm3", help="placeholder name for bookkeeping")
    args = ap.parse_args()

    device = pick_device(args.device)
    print(f"[esm3_writer_mac] device={device}", file=sys.stderr)

    # You MUST replace this with your real ESM3 load/inference code.
    # I am not guessing your ESM3 API because there are multiple incompatible setups.
    try:
        import esm  # noqa: F401
    except Exception:
        raise SystemExit(
            "[esm3_writer_mac] ESM3 dependencies not found in this env. "
            "If you're using Meta ESM3, install its package + weights, then adapt this wrapper."
        )

    records = read_fasta(args.in_fasta)
    if len(records) == 0:
        raise SystemExit(f"[esm3_writer_mac] No sequences found in {args.in_fasta}")

    # Placeholder output with NaNs so pipeline doesn't silently lie.
    rows = []
    for seq_id, seq in records:
        rows.append(
            {
                args.id_col: seq_id,
                "esm3_model": args.model,
                "seq_len": len(seq),
                "esm3_pll_sum": float("nan"),
                "esm3_pll_mean": float("nan"),
                "esm3_note": "TODO: implement ESM3 scoring for your install",
            }
        )

    df = pd.DataFrame(rows)
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[esm3_writer_mac] WROTE: {args.out_tsv} rows={len(df)}", file=sys.stderr)


if __name__ == "__main__":
    main()