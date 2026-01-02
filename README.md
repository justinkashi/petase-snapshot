# 🧬 Align 2025: Mechanistic PETase Engineering

**Engineering the next generation of plastic-degrading enzymes through evolution-aware and biophysically-grounded Machine Learning.**

---

## 🚀 Current Status: The Zero-Shot Sprint

**Deadline:** January 16, 2026

**Objective:** Rank the Align Tournament Dataset (4,988 variants) for **Expression**, **Activity (pH 5.5)**, and **Activity (pH 9.0)** without training on tournament labels.

### Our Strategy: The Hybrid Heuristic-MLP Engine

We are not "blindly" predicting. We are building a **mechanistic scoring engine** that translates biological laws into ranking signals.

1. **Backbone Grouping:** Mapping all 4,988 variants to three core "ancestor" clusters identified in the test set: **CaPETase**, **WP_162908185.1**, and **WP_374935857.1**.
2. **Feature Enrichment:** Moving beyond simple sequence stats to include hard physics (MD/Docking), local chemistry (pKa/NetQ), and genetic bottlenecks (mRNA/Codon Bias).
3. **The "Hidden Killer" Audit:** Explicitly penalizing variants based on product inhibition, cysteine mismatches, and purification "ghosts" (IMAC-capture likelihood).

---

## 🛠️ The Feature Engineering Pipeline

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

## 📋 Task Board (Updated: Jan 2, 2026)
1. run all tools 
	* fetch remaining structures of testset
	* annotation: phobius, tmhmm 
	* STABILITY: esm-1v, ddgemb, rosettaddg, deepddg, mutcompute, rnafold,thermoprot, prostab, temstapro (type2) temberture (type2)
	* SOLUBILITY: procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION 
	* pH/pka: propka, 
	* compute mutation score suppinfo 
	* MDsim, docking biophysical features 
2. fetch activity/expression/ph/temp studies 
3. Code the features (transfer annotations)
4. MLP weights on masterdb 
OTHER 
5. Finish dataset statistics and annotation graphs 
6. ESM-2/3 fine-tuned: 
	* Trained on benchmark solubility and stability datasets 
	* Trained on PETase datasets 

---

## 🗓️ Project Timeline

| Date | Milestone | Focus |
| --- | --- | --- |
| **Jan 2 (Today)** | **Feature Freeze** | Finalize all 26 feature columns and start batch inference. |
| **Jan 8** | **First Rank Submission** | Generate initial CSV ranking based on MLP-derived weights. |
| **Jan 16** | **Zero-Shot Deadline** | Finalized explainable ranking and biological abstract submission. |
| **Mar 9** | **Predictive Deadline** | Refine models using Supervised Track training data. |

---

## 🔬 Core Insights: What makes a PETase "Special"?

Our model prioritizes three critical structural motifs that distinguish true PETases from generic cutinases:

1. **The Aromatic Gate (W185):** Must be flexible enough to stack PET but stable enough to maintain the cleft.
2. **The  Loop:** Extended region (~238–260) that shapes the unique "open" cleft architecture.
3. **The DS1 Disulfide:** A PETase-specific bridge (C233-C282) that provides the structural integrity needed for high-turnover activity.

---

**Team:** Justin (Bio-ML Architecture), Charlie (SQL/MLP Training), Sanju (Docking/MD Dynamics), Aaisha (Graph/Sequence Embeddings).