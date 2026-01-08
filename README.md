# 🧬 Align 2025: Mechanistic PETase Engineering

**Engineering the next generation of plastic-degrading enzymes through evolution-aware and biophysically-grounded Machine Learning.**

**Team:** Justin (Bioinformatics), Charlie (MLOps), Sanju (Docking/MD Dynamics), Aaisha (MLOps).

---

## 🚀 Current Status: The Zero-Shot Sprint

**Deadline:** January 16, 2026

**Objective:** Rank the Align Tournament Dataset (4,988 variants) for **Expression**, **Activity (pH 5.5)**, and **Activity (pH 9.0)** without training on tournament labels.

### Our Strategy: The Hybrid Heuristic-MLP Engine

We are not "blindly" predicting. We are building a **mechanistic scoring engine** that translates biological laws into ranking signals.

1. **Backbone Grouping:** Mapping all 4,988 variants to three core "ancestor" clusters identified in the test set: **CaPETase**, **WP_162908185.1**, and **WP_374935857.1**.
2. **Feature Enrichment:** Moving beyond simple sequence stats to include hard physics (MD/Docking), local chemistry (pKa/NetQ), and genetic bottlenecks (mRNA/Codon Bias).
3. **Hidden killers penalization:** Explicitly penalizing variants based on product inhibition, cysteine mismatches, and purification "ghosts" (IMAC-capture likelihood).

---

## 🛠️ The Feature Engineering Pipeline
For the expression note that they dont tell us a certain pH, stability and gene expression and gene purification are whats important here.  Whereas for activity, the effect of the point mutation on activity will depend on the PETase mutant's stability, catalytic mechanism (cleft, gate, bridges). 
We're looking for the mutation's effect on the activity at pH 5.5 and pH 9, effect on the expression. Note that all enzymes were expressed, purified, and assayed on PET the SAME way. So whats important is how the mutation changes key factors in gene expression, purification, and activity. We want to try to find this sort of golden equation that perfectly weights each feature to output activity (μmol [TPA]/min·mg [E]) and expression (mg/mL). The trick here is how do we find the function/weights of each feature to get to that equation? 
From all the information we have, we can either learn what features and weights to use to learn what best predicts the effect of a PETase's single point mutation. Or we can make a heuristic model. This requires analytical work for each feature and understanding its relationship with activity and expression. 

### **1. Activity Prediction (Unified Mechanistic Equation)**

Used for both Condition 1 (pH 5.5) and Condition 2 (pH 9.0) by varying the weights () based on the pH-dependent electrostatics of the active site.

**ActivityScore = (w1 * ESM-1v_LLR) + (w2 * MutMatrix) + (w3 * KLDiv) + (w4 * Vina) + (w5 * RMSF) - (w6 * pKaPen) + (w7 * NetQ) + (w8 * VolDelta) + (w9 * ddG) - (w10 * CoEvo) - (w11 * ProdRel) + (w12 * HydroDelta) + (w13 * AromaticAnchor)**

* **Aromatic Anchors (w13):** Rewarding Trp/Phe/Tyr in the cleft for surface adsorption to PET plastic.
* **Product Release (w11):** Penalizing new H-bond donors that "stick" to the TPA product, preventing turnover.
* **pH Discrimination:** Using **PROPKA** to calculate the charge (`NetQ`) and Histidine protonation (`pKaPen`) shifts between 5.5 and 9.0.

### **2. Expression Prediction (Total Titer Equation)**

Predicting **Soluble Titer (mg/mL)** by modeling the competition between translation speed, folding energy, and purification efficiency.

**ExpressionScore = (wA * RNA_DeltaG) + (wB * CAI) - (wC * Stall) + (wD * ddG) - (wE * SAP) + (wF * NetSolP) - (wG * SASA) - (wH * Protease) + (wI * NQ) + (wJ * Tag) + (wK * Surf) - (wL * Metal) - (wM * Cys) - (wN * Burden)**

* **The DNA Gatekeeper:** Using **RNAfold** to model the mRNA stability around the RBS/Start codon.
* **The Folding Burden:** Combining `ddG` (stability) and `SAP` (spatial aggregation) to predict inclusion body formation.
* **The Purification Logic:** Weighting `Tag_Exposure` and `Metal_Binding_Patches` to predict how much protein actually reaches the final tube after IMAC.

---

## 📋 Task Board (Updated: Jan 7, 2026)

- UPDATED PLAN FOR THE DAY 
	1. WT and Test set mutation type flagging/annotation 
	2. Thermoprot (lu) -> SANJU
	3. GRAPE (lu) -> SANJU 
	3. Mutcompute -> SANJU 
	4. Mut score (alam) -> JUSTIN 
	5. mutPSSM (Lu) -> JUSTIN 
	6. evocouplings -> CHARLIE 
	7. Soluprot -> JUSTIN CHARLIE 
	7. rosetta/pyrosetta/fastrelax/foldx method -> SANJU   
	7. STABILITY: esm-1v, ddgemb, rosettaddg, deepddg, mutcompute, rnafold,thermoprot, prostab, temstapro (type2) temberture (type2)
	8. SOLUBILITY: procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION 


---

## 🗓️ Project Timeline

| Date | Milestone | Focus |
| --- | --- | --- |
| **Jan 8** | **First Rank Submission** | Generate initial test set ranking based on MLP-derived weights and heuristic equation. |
| **Jan 16** | **Zero-Shot Deadline** | Finalized explainable ranking and biological abstract submission. |
| **Mar 9** | **Predictive Deadline** | Refine models using Supervised Track training data. |

---

## 🔬 Core Insights: What makes a PETase "Special"?

Our model prioritizes three critical structural motifs that distinguish true PETases from generic cutinases:

1. **The Aromatic Gate (W185):** Must be flexible enough to stack PET but stable enough to maintain the cleft.
2. **The  Loop:** Extended region (~238–260) that shapes the unique "open" cleft architecture.
3. **The DS1 Disulfide:** A PETase-specific bridge (C233-C282) that provides the structural integrity needed for high-turnover activity.

--- 