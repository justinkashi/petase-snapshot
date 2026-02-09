#!/usr/bin/env python3
"""
Extract residues (and/or ungapped residue indices) from an MSA at reference (anchor) positions.

Key behavior:
- Anchor numbering is the anchor's UNGAPPED 1-based residue positions (gaps ignored).
- Validates that the anchor has the expected residue at each position (e.g., S160).
- For each other sequence, reports the residue aligned to that anchor position, plus that residue's
  UNGAPPED 1-based index in the target sequence (if not a gap).

"compact" mode matches your example:
- if target residue matches expected: output index only (e.g., 129)
- if mismatch: output *<AA><index> (e.g., *F61)
- if gap: output empty
"""

import argparse
import re
import sys
from collections import OrderedDict, defaultdict

GAPS = set("-.")

def read_fasta_alignment(path):
    seqs = OrderedDict()
    name = None
    buf = []
    opener = open(path, "r") if path != "-" else sys.stdin
    with opener as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                header = line[1:].strip()
                name = header.split()[0]  # first token
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    if not seqs:
        raise SystemExit("ERROR: no sequences found (FASTA alignment expected).")
    # sanity: all same length
    lens = {len(s) for s in seqs.values()}
    if len(lens) != 1:
        raise SystemExit(f"ERROR: alignment sequences have different lengths: {sorted(lens)}")
    return seqs

_pos_re = re.compile(r"^\s*([A-Za-z])?\s*(\d+)\s*$")

def parse_positions(lines):
    # returns list of (label, expected_aa_or_None, pos_int)
    out = []
    for raw in lines:
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        m = _pos_re.match(raw)
        if not m:
            raise SystemExit(f"ERROR: cannot parse position line: {raw!r} (expected like S160 or 160)")
        aa = m.group(1).upper() if m.group(1) else None
        pos = int(m.group(2))
        label = f"{aa}{pos}" if aa else str(pos)
        out.append((label, aa, pos))
    if not out:
        raise SystemExit("ERROR: no positions provided.")
    return out

def unique_headers(labels):
    counts = defaultdict(int)
    out = []
    for lab in labels:
        counts[lab] += 1
        out.append(lab if counts[lab] == 1 else f"{lab}_{counts[lab]}")
    return out

def find_anchor_name(seqs, anchor_query):
    if anchor_query in seqs:
        return anchor_query
    hits = [k for k in seqs.keys() if anchor_query in k]
    if len(hits) == 1:
        return hits[0]
    if len(hits) == 0:
        raise SystemExit(f"ERROR: anchor {anchor_query!r} not found (exact) and no substring matches.")
    raise SystemExit(f"ERROR: anchor query {anchor_query!r} is ambiguous; matches: {hits}")

def build_refpos_to_col(anchor_seq):
    refpos_to_col = {}
    refpos = 0
    for col, ch in enumerate(anchor_seq):
        if ch not in GAPS:
            refpos += 1
            refpos_to_col[refpos] = col
    return refpos_to_col, refpos

def cum_ungapped_counts(seq):
    # cum[i] = number of non-gap residues in seq[:i]
    cum = [0] * (len(seq) + 1)
    c = 0
    for i, ch in enumerate(seq):
        if ch not in GAPS:
            c += 1
        cum[i + 1] = c
    return cum

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--aln", required=True, help="Alignment FASTA file (use - for stdin).")
    ap.add_argument("--anchor", required=True, help="Anchor sequence name (exact or unique substring).")
    ap.add_argument("-p", "--positions", required=True,
                    help="Positions file (one per line like S160) or '-' to read from stdin.")
    ap.add_argument("--mode", choices=["compact", "residue", "pos", "both"], default="compact",
                    help="Output mode. Default: compact (matches your example).")
    ap.add_argument("--delim", default="\t", help="Output delimiter (default: tab).")
    ap.add_argument("--include", default=None,
                    help="Optional regex: only output sequence names matching this (anchor always included).")
    args = ap.parse_args()

    seqs = read_fasta_alignment(args.aln)
    anchor_name = find_anchor_name(seqs, args.anchor)
    anchor_seq = seqs[anchor_name]

    # positions input
    if args.positions == "-":
        pos_lines = sys.stdin.read().splitlines()
    else:
        with open(args.positions, "r") as fh:
            pos_lines = fh.read().splitlines()
    positions = parse_positions(pos_lines)  # list of (label, expected, pos)

    refpos_to_col, ref_len = build_refpos_to_col(anchor_seq)

    # validate anchor expected residues + map positions -> columns
    pos_to_col = []
    for label, expected, pos in positions:
        if pos not in refpos_to_col:
            raise SystemExit(f"ERROR: anchor position {pos} not present (anchor ungapped length={ref_len}).")
        col = refpos_to_col[pos]
        anchor_aa = anchor_seq[col].upper()
        if expected and anchor_aa != expected:
            raise SystemExit(
                f"ERROR: anchor mismatch at {label}: expected {expected} at ungapped pos {pos}, "
                f"but alignment has {anchor_aa} (col {col+1})."
            )
        pos_to_col.append((label, expected, pos, col))

    # optional include filter
    inc_re = re.compile(args.include) if args.include else None

    # header (unique if duplicates like I218)
    header_labels = unique_headers([x[0] for x in positions])
    print(args.delim.join(["name"] + header_labels))

    # precompute target cumcounts
    cumcounts = {name: cum_ungapped_counts(seq) for name, seq in seqs.items()}

    # output rows
    for name, seq in seqs.items():
        if name != anchor_name and inc_re and not inc_re.search(name):
            continue
        row = [name]
        cum = cumcounts[name]
        for (_, expected, _, col), _hdr in zip(pos_to_col, header_labels):
            aa = seq[col].upper()
            if aa in GAPS:
                # aligned gap at that reference column
                if args.mode == "residue":
                    row.append("-")
                else:
                    row.append("")
                continue
            idx = str(cum[col + 1])  # ungapped 1-based index
            if args.mode == "pos":
                row.append(idx)
            elif args.mode == "residue":
                row.append(aa)
            elif args.mode == "both":
                row.append(f"{aa}{idx}")
            else:  # compact
                if expected and aa != expected:
                    row.append(f"*{aa}{idx}")
                else:
                    row.append(idx)
        print(args.delim.join(row))

if __name__ == "__main__":
    main()