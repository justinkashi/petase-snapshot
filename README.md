# Welcome to the Align 2025 PETase Tournament Project

Welcome! We’re really happy you’re joining 😊  
---

## What is this project?

We’re taking part in the **Align Bio 2025 PETase Tournament**, a challenge where teams are asked to **rank PETase enzyme sequences** based on how well they might perform.

# 🗓️ Project timeline

**Align tournament phases**
- **Jan 16, 2025 — Zero-shot phase deadline**  
  Focus on understanding the tournament dataset, feature engineering, and explainable ranking without training on labeled tournament outcomes.
- **Mar 9, 2025 — Predictive phase deadline**  
  Incorporate predictive models and learned signals while maintaining interpretability and robustness.
- **Jul 20, 2025 — Generate phase (Round 1) deadline**  
  Propose and rank generated or modified sequences based on insights from earlier phases.
- **Nov 20, 2026 — Winners announced**  
  Final results, retrospectives, and broader dissemination.

**Beyond the tournament**
- **Manuscript preparation**  
  Formal write-up of the methodology, biological insights, and lessons learned from the tournament phases.
- **Talks & seminars**  
  Internal and external presentations (e.g. MILA, ML for Protein Engineering series, Formal Languages & LLM seminar, Valence Labs).
- **Conferences & workshops**  
  Potential submissions or presentations (e.g. BioML, ICML workshops, TechBioTransformers, related venues).
- **Community & outreach**  
  Sharing progress and insights via LinkedIn posts, Slack channels, and informal write-ups.

Work done early in the tournament directly feeds later phases, publications, and presentations—nothing is wasted.

---

## What kind of problem is this?

This is best thought of as:
- a **ranking problem**, not a perfect prediction problem
- an **understanding-first** problem, not a brute-force modeling problem
- a place where **clear reasoning beats fancy methods**

We are not trying to predict exact activity numbers.  
We are trying to decide *which sequences look better than others* and justify that decision.

---

## The two datasets you’ll hear about

### TOURNAMENTDB (the main thing we work on)
- Provided by Align
- This is **the dataset we analyze, rank, and submit**
- Everything we do should eventually connect back to this dataset

### MASTERDB (reference only)
- A curated collection of:
  - known PETases (IsPETase, CaPETase, BhrPETase)
  - LCCs and cutinase-like enzymes
  - mutants, structures, and literature annotations
- Used to:
  - understand what “real” PETases look like
  - define motifs and structural features
  - sanity-check our annotations
- This dataset is **frozen** and not the main modeling target

A helpful way to think about it:
> TournamentDB = what we’re judging  
> MasterDB = examples and background

---

## What makes a PETase “special”?

Not every enzyme that touches PET is the same.

### True PETases  
(e.g. IsPETase, CaPETase from *Cryptosporangium aurantiacum*, BhrPETase)

They usually share:
- An intact catalytic triad (IsPETase reference S160–D206–H237)
- **Trp185**, an aromatic residue that helps bind PET
- An **extended β8–α6 loop (~238–260)** that keeps the active site open
- PETase-specific disulfide bonds
- Activity at relatively mild temperatures

### LCC (leaf-branch compost cutinase)
- Much more thermostable
- Active site is more closed
- Works best near PET’s glass transition temperature

### Ancestral / cutinase-like hydrolases
- Have a catalytic triad
- Lack PETase-specific loops and disulfides
- Generally poorer or less specific PET activity

A lot of our work is about **telling these groups apart inside the tournament dataset**.

---

## What the main notebook actually does

Most work happens in `notebooks/main.ipynb`.

### 1. Dataset overview
We start by getting a feel for the data:
- sequence length distributions
- redundancy and clustering
- how many sequences fall in a PETase-like length range (~280–320 aa)
- ambiguous or unusual amino acids

### 2. PETase-specific features
For each tournament sequence, we annotate things like:
- catalytic triad integrity
- Trp185 identity (W / F / Y)
- presence and length of the β8–α6 loop
- disulfide architecture
- local sequence context around important residues:
  - flexibility (Gly/Pro)
  - charge
  - hydrogen-bond potential
  - steric bulk
- mutation-effect scores from supplement tables

### 3. Broader annotation
We also add:
- phylogenetic context
- signal peptide and transmembrane helix predictions
- taxonomy
- BLAST results (UniProt / InterPro / MasterDB)
- InterPro / Pfam domains
- AlphaFold2 structures (for selected sequences)
- ESM2 / ESM3 embeddings
- stability and fitness predictions (FoldX, MutCompute, etc.)

---

## How the GitHub repo is organized

```text
.
├── data
│   ├── master_db        # Curated reference data (frozen)
│   └── tournament_db    # Align tournament data + outputs
├── misc
│   └── suppinfo         # Supplementary tables, mutation scores, notes
├── notebooks
│   ├── main.ipynb       # Core analysis, annotation, ranking
│   └── tutorials        # Fine-tuning & past enzyme tournament examples
├── scripts
│   ├── bash             # FASTA handling, clustering, file ops
│   └── python           # Annotation and automation helpers