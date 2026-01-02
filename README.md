# Welcome to the Align 2025 PETase Tournament Project😊

# TASKS BOARD DO AS OF: JAN 1 2026
- NONCODING
	* Fetching activity data 
	* Define pH dependance of PETase 
	* Define temperature dependance of PETase 
	* Slides
	* Write manuscript, make lit review  
- CODING
	* Annotation 
		• Interproscan, MSA graphs, MEME motifs, PET-ase specific regions, 
		• BLASTP
		• Phobius
		• IQTREE
		• TMHMM
		• SignalP 

		• mutcompute
		• thermostability tools 
		• solubility prediction tools
		• Expression vector annotation (orf-level, dna-level)
	
		• Fetch pdb of tournament_wt + generate AF2 structures if no pdb + map tournament_test mutation codes + generate FoldX structures 
		
		• Try FoldX alternatives for tournament_wt/test, master_db
			•	Rosetta ΔΔG (ddg_monomer / cartesian_ddg) (slower)
			•	DeepDDG 
			•	ThermoNet (GNN)  — predict stability changes from structure graphs
			•	MAESTRO(ML)
			•	Sequence-based - EVEscape / ESM-Mut / ESM1b/GEMME/EVE PROVEAN / SIFT / PolyPhen

		• Docking, MD 


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
-- 

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