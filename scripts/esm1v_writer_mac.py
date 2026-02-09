#!/usr/bin/env python3
# Usage: python esm1v_writer_mac.py --in_fasta data/test.fasta --out_tsv out/esm1v_mac.tsv

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

try:
    import esm
except ImportError:
    raise SystemExit(
        "[esm1v_writer_mac] Missing dependency: fair-esm. "
        "Install with: pip install fair-esm"
    )


def pick_device(device: str) -> torch.device:
    device = device.lower()
    if device == "auto":
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    if device == "mps":
        if not torch.backends.mps.is_available():
            print("[esm1v_writer_mac][WARN] MPS not available; falling back to CPU", file=sys.stderr)
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


@torch.no_grad()
def score_batch(model, alphabet, batch_converter, batch_records, device: torch.device):
    labels, seqs = zip(*batch_records)
    batch = list(zip(labels, seqs))
    batch_labels, batch_strs, batch_tokens = batch_converter(batch)
    batch_tokens = batch_tokens.to(device)

    out = model(batch_tokens, repr_layers=[], return_contacts=False)
    logits = out["logits"]  # [B, T, vocab]

    # token indices
    pad_idx = alphabet.padding_idx
    bos_idx = alphabet.cls_idx
    eos_idx = alphabet.eos_idx

    results = []
    for i in range(len(batch_records)):
        toks = batch_tokens[i]  # [T]
        logit = logits[i]       # [T, vocab]

        # valid positions = between BOS and EOS (exclude special/pad)
        valid = (toks != pad_idx) & (toks != bos_idx) & (toks != eos_idx)
        idx = torch.where(valid)[0]

        # log prob of true token at each position
        log_probs = torch.log_softmax(logit[idx], dim=-1)  # [L, vocab]
        true_tok = toks[idx]                               # [L]
        ll = log_probs.gather(1, true_tok[:, None]).squeeze(1)  # [L]

        pll_sum = float(ll.sum().cpu().item())
        pll_mean = float(ll.mean().cpu().item())
        L = int(len(idx))

        results.append((batch_labels[i], L, pll_sum, pll_mean))

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_fasta", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--id_col", default="seq_id")
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--device", default="auto", choices=["auto", "mps", "cpu"])
    ap.add_argument(
        "--model",
        default="esm1v_t33_650M_UR90S_1",
        help="One of: esm1v_t33_650M_UR90S_1..5 (default 1)",
    )
    args = ap.parse_args()

    device = pick_device(args.device)
    print(f"[esm1v_writer_mac] device={device}", file=sys.stderr)

    # load model
    model_name = args.model.strip()
    if not hasattr(esm.pretrained, model_name):
        raise SystemExit(
            f"[esm1v_writer_mac] Unknown model {model_name}. "
            f"Valid are esm1v_t33_650M_UR90S_1..5"
        )

    model_loader = getattr(esm.pretrained, model_name)
    model, alphabet = model_loader()
    model.eval()
    model = model.to(device)

    batch_converter = alphabet.get_batch_converter()

    records = read_fasta(args.in_fasta)
    if len(records) == 0:
        raise SystemExit(f"[esm1v_writer_mac] No sequences found in {args.in_fasta}")

    rows = []
    bs = max(1, args.batch_size)
    for i in range(0, len(records), bs):
        batch_records = records[i : i + bs]
        scored = score_batch(model, alphabet, batch_converter, batch_records, device)
        for seq_id, L, pll_sum, pll_mean in scored:
            rows.append(
                {
                    args.id_col: seq_id,
                    "esm1v_model": model_name,
                    "seq_len": L,
                    "esm1v_pll_sum": pll_sum,
                    "esm1v_pll_mean": pll_mean,
                }
            )
        if (i // bs) % 10 == 0:
            print(f"[esm1v_writer_mac] scored {min(i+bs, len(records))}/{len(records)}", file=sys.stderr)

    df = pd.DataFrame(rows)
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[esm1v_writer_mac] WROTE: {args.out_tsv} rows={len(df)}", file=sys.stderr)


if __name__ == "__main__":
    main()