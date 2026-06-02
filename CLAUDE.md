# petorg — CLAUDE.md

Read this fully before doing anything. This project has a lot of context.

---

## What this project is

**Align 2025 competition** (currently paused, no labeled data returned yet) — zero-shot ranking of ~4,674 single-mutant PETase variants across 3 wildtype backbones.

**SciPy 2026 talk + proceedings manuscript** (1 month deadline) — "Simulation-Informed Machine Learning Workflows for PETase Engineering." This is the active goal.

**PETase** = plastic-degrading enzyme that hydrolyzes PET polymer into MHET/TPA.

---

## The 3 wildtype backbones

| ID | Enzyme | Structure | Notes |
|---|---|---|---|
| tournament_wt_1 | CaPETase | AlphaFold (X-ray 7YM9 exists but use AF for consistency) | Catalytic: Ser129, Asp175, His206 (estimated — verify) |
| tournament_wt_2 | Allorhizocola PETase | AlphaFold only | Catalytic: Ser130, Asp176, His207 (estimated) |
| tournament_wt_3 | PETaseSM14 | AlphaFold (X-ray 9HYD exists but use AF for consistency) | Catalytic: Ser131, Asp177, His209 (estimated) |

1,558 single-mutant variants per WT = 4,674 total tournament variants.

**Competition targets:** activity_1 (pH 5.5), activity_2 (pH 9.0), expression (E. coli BL21)
**Metric:** NDCG@10

---

## What is already computed and where

### The main feature table
**`.claude/worktrees/thirsty-davinci-951487/may11/results/mega_features/mega_features.tsv`**
- 4,674 rows × 76 features. This is the canonical feature table.
- Do NOT recompute features that already exist here.

**All 76 feature columns:**
- PETase flags (18): `flag_triad_broken`, `flag_oxyanion_disrupted`, `flag_cleft_landmark_mutated`, `flag_lipase_box_disrupted`, `flag_cys_introduced/removed`, `flag_disulfide_risk`, `flag_loop_rigidity_shock`, `flag_electrostatic_rewiring`, `flag_charge_flip_near_active`, `flag_product_release_sticky`, `flag_aromatic_clamp_change`, `flag_hydrophobic_core_change`, `flag_proline_introduced`, `flag_glycine_removed`, `flag_n_severe`, `flag_n_moderate`, `flag_penalty_score`, `flag_vector`
- ESM1v distances (6): `esm1v_l2`, `esm1v_l1`, `esm1v_cosine`, `esm1v_max_abs`, `esm1v_mean_abs`, `esm1v_rel_l2`
- EVcouplings (12): `evc_*` at b03 and b07 thresholds
- Stability: `dynamut2_ddg`, `foldx_ddg`
- Solubility: `netsolp_solubility`
- Fitness: `poet2_score`, `poet2_sub_score`, `mutation_llr`, `pseudo_likelihood`
- Aggregation: `a3d_mut_site`
- Structural: `dist_to_active_site`, `structure_risk_score`, `burial_status`, `near_active_site`, `secondary_structure`
- pH: `charge_change_pH5.5`, `charge_change_pH9.0`, `pH_differential`
- Expression: `cai`, `rare_codon_freq`, `gc_content`
- GNN: `gnn_activity_1`, `gnn_activity_2`, `gnn_expression`
- FireProt-derived: `fp_conservation`, `fp_correlation`, `fp_bfactor`, `fp_is_btc`
- MutCompute: `mc_wt_prob`, `mc_pred_prob`, `mc_avg_log_ratio`

### Rankings
**`.claude/worktrees/thirsty-davinci-951487/may11/results/mega_scored/scored.tsv`**
- 4,674 variants ranked by 16+ hand-written scoring formulas (weighted z-score combinations)
- **NOT ML predictions** — these are manually designed formulas
- Top consensus variant: tournament_test_1434 (A54G on wt1)

### ML benchmark validation (done, never applied to tournament variants)
**`.claude/worktrees/thirsty-davinci-951487/may11/results/`**
- `eval_fireprot_ml/eval_results.json` — Ridge/GBR on FireProtDB: GBR Spearman=0.638, NDCG@10=0.767
- `eval_stability_ml/eval_results.json` — expanded benchmark evaluation
- LOPO (leave-one-protein-out) Spearman=0.375 = honest generalization estimate
- **The ML model was trained on FireProtDB and evaluated there — it was NEVER applied to produce tournament variant scores. That is the critical missing piece.**

### Other computed results in `results/`
- `dynamut2/wt{1,2,3}_dynamut2.csv` — 1,559 rows each, complete
- `netsolp_test*.csv` — 8 batches, 4,988 complete solubility predictions
- `tournament_test_wt*_poet2_ranksequences.csv` — POET2 scores
- `esm1v_delta_features.tsv` — 4,674 rows
- `wt{1,2,3}_af_A3D.csv` — Aggrescan3D
- `evcouplings_results_wildtypeset/` — EVcouplings (sparse for many variants)
- `foldx_parallel/wt{1,2,3}/` — FoldX mutant PDB structures + ΔΔG

---

## Benchmark datasets for ML training

These are for training models to then apply to the 4,674 tournament variants.

### Stability (mutation-level ΔΔG/ΔTm → predicts activity_1, activity_2)
| Dataset | Size | Label | Local path |
|---|---|---|---|
| FireProtDB | 53K | ddG, dTm | `data/benchmark/stab_benchmark/fireprotdb1/fireprotdb_results_stability.csv` |
| ThermomutDB | 12K | ddG, dTm | `data/benchmark/stab_benchmark/thermomutdb.json` |
| Novozymes | 31K | Tm (whole protein) | `data/benchmark/novozymes-enzyme-stability-prediction/` |
| Meltome | 1M | Tm (noisy) | `data/benchmark/stab_benchmark/meltome/meltome_cross-species.csv` |
| ProteinGym | large | DMS fitness | `data/benchmark/proteingym/` (may need download) |
| MaveDB | varies | fitness | `data/benchmark/mavedb/` (may need download) |

### Solubility/expression (→ predicts expression target)
| Dataset | Size | Label | Local path |
|---|---|---|---|
| ProtSolM | 71K | binary solubility | `data/benchmark/protsolm_data/` |
| SoluProt | 11K | binary solubility | `data/benchmark/sol_benchmark/soluprot_data/` |
| NESG | 10K | exp + solubility | `data/benchmark/sol_benchmark/nesg/` |
| Price | 7K | usability (binary) | `data/benchmark/sol_benchmark/Price_usability_trainset.csv` |
| PSI Biology | 11K | binary solubility | `data/benchmark/sol_benchmark/PSI/` |

### Competition-format (MOST IMPORTANT — actual Align activity labels)
| Dataset | Notes |
|---|---|
| **Align2023** | 4 enzymes, focus on beta-glucosidaseB + alpha-amylase only, single mutants only. Path: `data/benchmark/align2023/**/train.csv`. Same experimental format as Align2025 PETase. |

**The benchmark.ipynb notebook** (`notebooks/benchmark.ipynb`) has complete loading functions for all datasets, CD-HIT clustering for deduplication, and WT-mutant pair building. Read it before writing any benchmark loading code.

---

## Key scripts

| Script | Purpose |
|---|---|
| `notebooks/benchmark.ipynb` | Full benchmark pipeline — loading, normalization, CD-HIT, pairs |
| `.claude/worktrees/.../benchmarks/eval_fireprot_ml.py` | ML training + eval on FireProtDB |
| `.claude/worktrees/.../benchmarks/score_mega.py` | Applies hand-written scoring formulas to mega_features.tsv |
| `.claude/worktrees/.../benchmarks/compile_all_features.py` | Builds mega_features.tsv from individual result files |
| `.claude/worktrees/.../benchmarks/build_consensus.py` | Consensus ranking across 16 methods |
| `scripts/run_foldx_parallel.py` | FoldX ΔΔG computation |
| `scripts/esm1v_writer_mac.py` / `esm2_writer_mac.py` | ESM embedding computation |
| `scripts/evc_writer.py` | EVcouplings feature extraction |
| `scripts/docking_protocol_code.txt` | AutoDock Vina batch docking scripts (shell + Python) |
| `scripts/docking protocol_v2.pdf` | Full 11-page docking protocol (Snakemake pipeline) |

---

## Docking protocol (not yet run)

**What it is:** Dock PET substrate into FoldX mutant structures to get catalytic geometry features. NOT FoldX — that's stability. Docking = ligand binding.

**Primary ligand:** PET trimer (2-HE(MHET)₃)
SMILES: `OCCOC(=O)c1ccc(cc1)C(=O)OCCOC(=O)c1ccc(cc1)C(=O)OCCO`
Secondary ligands: MHET, TPA (for completeness)

**Setup:** AutoDock Vina, 24Å box centered on catalytic Ser Oγ, 200 poses per run, filter poses by Ser–ester distance < 4Å. Run twice per variant: pH 5.5 and pH 9.0 (use PROPKA/pdb2pqr to protonate structures).

**What's blocking it:** Catalytic residue numbers (Ser129 etc.) are estimated — need verification from actual PDB files before running. Ligand PDBQT needs to be built. Batch runner script needs to be written.

**Full protocol:** See `scripts/docking protocol_v2.pdf` — 11 pages, covers Step 0 (filtering) through Step 7 (scoring aggregation) including Snakemake pipeline.

---

## What remains to do (priority order)

1. **Multi-dataset ML training experiment** — Train models on all benchmark datasets (FireProtDB, ThermomutDB, Align2023, ProteinGym, etc.) across stability and expression tracks. Multiple model types (Ridge, GBR, XGB, RF). Evaluate in parallel on GPU PC. Apply best model to 4,674 tournament variants → add `ml_stability_score` and `ml_expression_score` columns to mega_features.tsv. **This is the central ML contribution of the paper.**

2. **Docking pipeline** — Verify active site residue numbers, build PET trimer PDBQT, write batch runner, run on GPU PC overnight. Start with top 100–200 variants from scored.tsv. Produces catalytic geometry features to add to mega_features.tsv.

3. **Manuscript** — SciPy 2026 proceedings paper. Story: zero-shot pipeline → ML transfer learning (FireProtDB + Align2023) → docking validation → ranked variant list. Sections: intro, methods (feature pipeline), results (benchmark validation, ranking), discussion.

---

## Infrastructure

- **Mac (here):** Development, Claude Code, analysis, writing
- **GPU PC (remote):** Long compute jobs — ML training, docking batch runs, ESM inference. Access via AnyDesk or Tailscale SSH.
- **venv:** Always activate `/Users/bustin/petorg/venv/bin/activate` before any Python
- **Papers:** 46 PDFs in `petase_paperdb/` — use NotebookLM for paper Q&A
- **MCP server:** `scripts/paperqa_mcp.py` registered in `.mcp.json` (simple PDF text search, no API key)

---

## Things NOT to do

- Do not recompute features that exist in mega_features.tsv
- Do not use the worktree scripts in isolation — understand them first, they were written by a previous agent session
- Do not use X-ray structures — use AlphaFold for all 3 WTs for consistency
- Do not use BHET as the primary docking ligand — use PET trimer per the protocol
- Do not apply ML scores without first validating on a held-out benchmark split
- PROVEAN: abandoned — PSI-BLAST NaN bug, 200–300GB database requirement, not worth fixing
- Aggrescan3D CLI: requires Python 2.7, not usable — webserver results already in results/
