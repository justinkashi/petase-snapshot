#!/usr/bin/env python3
#usage: python3 /Users/bustin/petorg/scripts/run_foldx_parallel.py \
#  --foldx-bin /Users/bustin/petorg/bioapps/foldx5_1Mac/foldx_20251231 \
#  --wt-pdb /Users/bustin/petorg/data/wt1_af.pdb \
#  --mut-list /Users/bustin/petorg/data/wt1_individual_list.txt \
#  --out-dir /Users/bustin/petorg/results/foldx_parallel/wt1 \
#  --jobs 2 --chunk-size 20 &
import argparse
import os
import re
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

TOKEN_RE = re.compile(r"^[A-Z][A-Z][0-9]+[A-Z];$")

def load_mutations(mut_file: Path):
    muts = []
    for line in mut_file.read_text().splitlines():
        t = line.strip()
        if not t or t.lower() == "mutation":
            continue
        if not t.endswith(";"):
            t += ";"
        if not TOKEN_RE.match(t):
            raise ValueError(f"Bad FoldX token: {t}")
        muts.append(t)
    return muts

def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def run_cmd(cmd, cwd: Path, log_name: str):
    p = subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (cwd / log_name).write_text(p.stdout or "")
    if p.returncode != 0:
        raise RuntimeError(f"Command failed in {cwd}: {' '.join(cmd)}")

def repair_wt(foldx_bin: Path, wt_pdb: Path, repair_dir: Path):
    repair_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(wt_pdb, repair_dir / "input.pdb")
    cmd = [str(foldx_bin), "--command=RepairPDB", "--pdb=input.pdb", "--output-dir=."]
    run_cmd(cmd, repair_dir, "repair.log")
    repaired = repair_dir / "input_Repair.pdb"
    if not repaired.exists():
        raise RuntimeError("RepairPDB did not produce input_Repair.pdb")
    return repaired

def run_chunk(idx, chunk_muts, foldx_bin: Path, repaired_pdb: Path, jobs_dir: Path, out_pdb_dir: Path):
    job_dir = jobs_dir / f"chunk_{idx:04d}"
    if job_dir.exists():
        shutil.rmtree(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(repaired_pdb, job_dir / "input_Repair.pdb")
    (job_dir / "individual_list.txt").write_text("\n".join(chunk_muts) + "\n")

    cmd = [
        str(foldx_bin),
        "--command=BuildModel",
        "--pdb=input_Repair.pdb",
        "--mutant-file=individual_list.txt",
        "--numberOfRuns=1",
        "--output-dir=."
    ]
    run_cmd(cmd, job_dir, "buildmodel.log")

    produced = sorted(job_dir.glob("input_Repair_*.pdb"), key=lambda p: int(p.stem.split("_")[-1]))
    if len(produced) < len(chunk_muts):
        raise RuntimeError(f"{job_dir}: expected {len(chunk_muts)} mutant PDBs, got {len(produced)}")

    for i, mut in enumerate(chunk_muts):
        out_name = mut[:-1] + ".pdb"
        shutil.move(str(produced[i]), str(out_pdb_dir / out_name))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--foldx-bin", required=True)
    ap.add_argument("--wt-pdb", required=True)
    ap.add_argument("--mut-list", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--chunk-size", type=int, default=20)
    args = ap.parse_args()

    foldx_bin = Path(args.foldx_bin).resolve()
    wt_pdb = Path(args.wt_pdb).resolve()
    mut_list = Path(args.mut_list).resolve()
    out_dir = Path(args.out_dir).resolve()

    out_dir.mkdir(parents=True, exist_ok=True)
    repair_dir = out_dir / "00_repair"
    jobs_dir = out_dir / "01_jobs"
    out_pdb_dir = out_dir / "mutant_pdbs"
    out_pdb_dir.mkdir(parents=True, exist_ok=True)

    muts = load_mutations(mut_list)
    repaired = repair_wt(foldx_bin, wt_pdb, repair_dir)

    chunked = list(chunks(muts, args.chunk_size))
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [
            ex.submit(run_chunk, i, c, foldx_bin, repaired, jobs_dir, out_pdb_dir)
            for i, c in enumerate(chunked, start=1)
        ]
        for f in as_completed(futs):
            f.result()

    (out_dir / "done.txt").write_text(f"OK: {len(muts)} mutations\n")

if __name__ == "__main__":
    main()
