# Align 2025 PETase Engineering

## Current Status: Zero-Shot Phase (No Tournament Labels)

**Deadline:** January 16, 2026  
**Goal:** Rank 4,988 single-variant PETase sequences for:
- **Expression / recoverable soluble titer (mg/mL)**  
- **Activity at pH 5.5**  
- **Activity at pH 9.0**
---

## Abstract: What We Are Doing (and Why It Works Zero-Shot)

We are building a **hybrid ranking engine** that combines four complementary pillars:

1. **General, model-agnostic biophysical features:**  
   Universal predictors of foldability and expressibility (stability/ΔΔG proxies, solubility/aggregation risk, sequence “naturalness”/complexity), used as a baseline prior across all variants.

2. **PETase-specific mechanistic constraints:**  
   Some mutations destroy PETase catalytic architecture; we explicitly detect and penalize these failures using PETase landmark and motif logic (triad integrity, oxyanion-hole geometry proxies, aromatic clamp/gate residues, disulfide liabilities, cleft-shaping residues, and motif/architecture breaks with positional coupling).

3. **Evolutionary priors (PETase family models):**  
   Variants that remain statistically consistent with functional PETase evolution are more likely to fold, express, and retain activity. We capture this with **co-evolutionary likelihood** and **PLM log-likelihood / pseudo-perplexity / Δlogprob** style scores, evaluated within the appropriate WT lineage.

4. **Physics/structure refinement on a small subset (Nimbus-style):**  
   We do not attempt docking/MD on all 4,988. Instead, we run **expensive structure modeling only after triage** to refine top candidates per lineage using standardized substrate placement, constrained minimization, and a compact set of geometric/energetic metrics.


---

## Dataset Structure: Backbone Lineages and Coordinate System

### Backbone grouping (lineage assignment)
We do not treat all 4,988 variants as one homogeneous space. Each variant is assigned to a parent backbone cluster (WT lineage) so that scoring is **relative to the correct reference** rather than mixing incompatible architectures.`

---

## Scoring System 

### Stage 0 — Mutation Flagging (Deterministic Gate Layer)
Before any ML, evolutionary, or structure scoring, we annotate each variant with **binary flags + severity penalties** for failure modes that are expected to dominate ranking regardless of downstream predictors.

#### 0.1 PETase Catalytic Architecture (landmark-mapped)
- **Catalytic triad integrity:** Ser–Asp–His equivalents present at mapped landmark positions  
- **Oxyanion support:** Y87-region equivalent present **and** M161-equivalent backbone-NH context preserved (no local breaker mutations)  
- **Cleft landmarks:** W159 / W185 / S238 / N241 equivalents (presence + mutation severity)

#### 0.2 PETase Motif / Structural Integrity (positional coupling enforced)
- **M-class (M1–M5) classification:** block-based detection with required spacing (not “motif exists somewhere”)  
- **Lipase box context:** GxSxG neighborhood (and embedding context where relevant)  
- **Disulfide liability:** loss of native Cys, gain of new Cys, or Cys count/placement changes consistent with mispairing risk

#### 0.3 Mechanistic Mutation-Type Flags (interpretable effect channels)
- **Loop rigidity shocks:** Pro↔Gly or bulky↔small swaps at loop/hinge sites (gate/W-loop/cleft-adjacent)  
- **Electrostatic rewiring:** charge flips near catalytic residues, salt-bridge networks, or cleft electrostatic field  
- **Product-release stickiness:** creation of strong donor/acceptor patterns near the exit path likely to trap TPA/MHET-like products

**Outputs:** a per-variant flag vector + penalty score used as (i) hard disqualifiers when severe, and (ii) additive penalties in later ranking stages.

## Feature Engineering Pipeline (Zero-Shot)

We compute three property tracks in parallel:

### A) Activity @ pH 5.5 and pH 9.0 (same features, different weights)
Activity is treated as a function of:
- catalytic geometry + cleft access (mechanistic constraints)
- electrostatics/protonation sensitivity (pH dependence)
- stability (active enzymes must stay folded under assay conditions)

### B) Expression / recoverable soluble titer
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

## 📋 Task Board (Updated: Jan 7, 2026)

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
