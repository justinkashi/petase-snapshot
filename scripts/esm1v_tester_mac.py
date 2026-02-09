#!/usr/bin/env python
# esm1v_tester_mac.py
# Test ESM-1v throughput on Mac CPU or MPS without writing outputs
# -------------------------------------------------------------
# Difference vs Windows+CUDA setup:
#   This version runs on Mac using CPU or MPS (Metal Performance Shaders). 
#   Unlike Windows+CUDA, where each process can spawn isolated GPU contexts, 
#   Mac MPS shares the GPU differently and does not allow true multiprocessing GPU 
#   acceleration. Batch size is smaller to avoid memory issues; performance is limited 
#   by single-device throughput rather than full GPU saturation.
# -------------------------------------------------------------

# Usage:
#   python esm1v_tester_mac.py --fasta data/esm1v_input.fasta  --batch 1 --device mps
#   Options:
#     --fasta    Path to input FASTA
#     --batch    Batch size (1–4 recommended for MPS)
#     --device   'cpu' or 'mps'

import os, time, psutil
import torch, esm

# ------------------ CPU/MEMORY ------------------
proc = psutil.Process(os.getpid())
def mem_cpu():
    mem_gb = proc.memory_info().rss / (1024**3)
    cpu = proc.cpu_percent(interval=None)
    return mem_gb, cpu

# ------------------ CONFIG ------------------
fasta_path = "data/esm1v_input.fasta"
device = "mps" if torch.backends.mps.is_available() else "cpu"
BATCH_SIZES = [1, 2, 4]           # try multiple batch sizes
N_REF = 50                        # test on first N sequences for speed

# ------------------ LOAD MODEL ------------------
model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
model = model.to(device).eval()
batch_converter = alphabet.get_batch_converter()

# ------------------ LOAD FASTA ------------------
seqs = []
with open(fasta_path) as f:
    name, buf = None, []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs.append((name, "".join(buf)))
            name = line[1:]
            buf = []
        else:
            buf.append(line)
    if name:
        seqs.append((name, "".join(buf)))

seqs = seqs[:N_REF]   # subset for testing
N = len(seqs)

print(f"Device: {device} | sequences: {N}")

# ------------------ TEST DIFFERENT BATCH SIZES ------------------
for BATCH_SIZE in BATCH_SIZES:
    print(f"\nTesting batch size {BATCH_SIZE} ...")
    t0 = time.time()

    with torch.no_grad():
        for i in range(0, N, BATCH_SIZE):
            batch = seqs[i:i+BATCH_SIZE]

            mem_before, _ = mem_cpu()
            _, _, tokens = batch_converter(batch)
            tokens = tokens.to(device)
            reps = model(tokens, repr_layers=[33])["representations"][33]
            emb = reps[:, 1:-1].mean(dim=1).cpu().numpy()  # just compute, don't save

            del tokens, reps, emb
            if device == "mps":
                torch.mps.empty_cache()
            mem_after, cpu = mem_cpu()
            delta = mem_after - mem_before

            print(f"[{i+1}/{N}] RAM {mem_after:.2f} GB (Δ {delta:+.2f}) | CPU {cpu:.1f}%")

    elapsed = time.time() - t0
    rate = N / elapsed
    print(f"Batch size {BATCH_SIZE} DONE | rate ~ {rate:.2f} seq/s | elapsed {elapsed:.1f}s")
