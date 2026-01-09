# Align 2025 PETase Engineering

## Current Status: Zero-Shot Phase 

**Deadline:** January 16, 2026  
**Goal:** Rank 4,988 single-variant PETase sequences for:
- **Expression / recoverable soluble titer (mg/mL)**  
- **Activity at pH 5.5**  
- **Activity at pH 9.0**
---

## Abstract (need to finish)
We are building a **hybrid ranking engine** that combines four types of features:

1. **General, family-agnostic features**  

2. **PETase-specific and evolutionary features**  

3. **3D modeling, docking/MD-sim features:**  


---

## Scoring System

### Stage 0 — Mutation Flagging
Before any ML, evolutionary scoring, or structure scoring, we annotate each variant with binary flags and severity penalties for failure modes that should dominate ranking regardless of downstream predictors.

#### 0.1 PETase catalytic architecture
- Catalytic triad integrity  
  Ser Asp His equivalents present at mapped landmark positions
- Oxyanion support  
  Y87-region equivalent present and M161-equivalent backbone-NH context not disrupted
- Cleft landmarks  
  W159 W185 S238 N241 equivalents flagged by presence and mutation severity

#### 0.2 PETase motif and structural integrity
- M-class M1 to M5 classification  
  block-based detection with required spacing and positional coupling
- Lipase box context  
  GxSxG neighborhood and local embedding context
- Disulfide liability  
  loss of native Cys gain of new Cys and Cys placement patterns consistent with mispairing risk

#### 0.3 Mechanistic mutation-type flags
- Loop rigidity shocks  
  Pro Gly swaps and bulky small swaps at loop and hinge sites near the gate and cleft
- Electrostatic rewiring  
  charge flips near catalytic residues salt-bridge networks or the cleft electrostatic field
- Product-release stickiness  
  new strong donor acceptor patterns near the exit path likely to trap products

**Outputs:** a per-variant flag vector and penalty score used as hard disqualifiers when severe and additive penalties in later stages.

## Feature Engineering

We compute three property tracks in parallel:

### A) Activity @ pH 5.5 and pH 9.0 
Activity is treated as a function of:
- catalytic geometry + cleft access (mechanistic constraints)
- electrostatics/protonation sensitivity (pH dependence)
- stability (active enzymes must stay folded under assay conditions)

### B) Expression 
Expression is treated as a competition between:
- translation + folding burden
- aggregation risk
- survivability post protein purification 

### C) Shared priors and constraints
These include PLM likelihoods, evolutionary couplings, and motif integrity.

---

## Tools Employed for Zero Shot Phase 

### 1) Evolutionary / Sequence-Based 
These generate a robust zero-shot prior for folding/function:
- **EVcouplings** likelihood / Δlog-likelihood (per lineage)
- **MutPSSM (Lu)** conservation-weighted mutation score
- **PLM log-likelihood signals** (ESM-1v LLR or equivalent)
- **Alam motif score** (M-class integrity + motif break penalties)

### 2) Biophysical stability 
We treat stability as a shared limiting factor for activity/expression:
- **ThermoProt**
- **ESM-1v stability proxy / LLR**
- **ddG predictors:** ddGemb, DeepDDG, RosettaΔΔG, FoldX (as available)
- Optional: Prostab / TemStaPro / other type-2 predictors

### 3) Solubility / aggregation
We explicitly penalize aggregation or insolubility failure modes:
- **NetSolP / ProteinSol / SoluProt**
- **Aggrescan3D** (if structure available; otherwise sequence proxy)
- Optional: SAP-like surface hydrophobicity proxy if we can afford it

### 4) pH / electrostatics 
We model pH dependence primarily through ionization changes near the cleft:
- **PROPKA** on WT structures at pH 5.5 and pH 9.0 (baseline)
- Variant pKa runs only for mutations near active site / charged networks

### 5) Structure/physics refinement (small subset only)
Run only after triage (top-K per lineage):
- Structure generation: AF2/ESMFold/ColabFold (WT + selected variants)
- Substrate docking: standardized protocol (e.g., DiffDock or internal)
- Constrained minimization: PyRosetta / Rosetta relax into a consistent pose family
- Metrics panel: interface energy proxies, clashes/packing, catalytic distances, cleft width proxies

---

## Property Equations (Work In Progress)

### ActivityScore(pH)
A weighted combination of:
- evolutionary priors (EVcouplings/PLM/PSSM)
- mechanistic penalties (motif/architecture gates)
- stability terms (ΔΔG ensemble)
- pH terms (ΔpKa / NetQ proxies)
- (optional) docking refinement terms for the top tranche

*We maintain separate weight for pH 5.5 vs pH 9.0 to reflect electrostatic shifts.

### ExpressionScore
A weighted combination of:
- solubility/aggregation predictors
- stability burden predictors (ΔΔG)
- sequence features tied to expression (codon/GC/CAI if relevant)
- purification survivability proxies (surface patches / tag exposure if modeled)

*Tournament protocols were fixed (expression, purification, assay), so the objective is not “universal expression” 

---

## 📋 Task Board (Updated: Jan 8, 2026)

1. WT and test-set mutation type flagging / annotation  
   - MSA zoom-ins for motif table  
   - Per-lineage M1–M5 placement + phylogenetic placement  

2. ThermoProt (Lu) → SANJU  
3. GRAPE (Lu) → SANJU  
4. MutCompute → SANJU  
5. Alam mutation score / motif integrity → JUSTIN  
6. MutPSSM (Lu) → JUSTIN  
7. EVcouplings → CHARLIE  
8. SoluProt + solubility consensus → JUSTIN / CHARLIE  
9. Rosetta/PyRosetta/FoldX pipeline (relax + ΔΔG) → SANJU  
10. Consolidate stability & solubility tool outputs into master feature table

---

## Project Timeline

| Date | Milestone | Focus |
| --- | --- | --- |
| **Jan 8** | First Rank Submission | Initial ranking from zero-shot priors + gates |
| **Jan 16** | Zero-Shot Deadline | Final explainable ranking + abstract |
| **Mar 9** | Predictive Deadline | Supervised refinement using training labels |
