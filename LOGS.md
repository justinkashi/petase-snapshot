# **feb13** 
- BOLTZ2: 
- POET2: (1) substitution analysis of all 3 WT (2) rerun rank seq analysis but unmasking the remaining critical PETase motifs  
- TUNING POET2 RESULTS: try using strong activity+expression PETase variants for the CONTEXT. Also try fixing the very very important residues of wt1/2/3 in the query instead of just S D H triad. 
- POET2: Substitution analysis + Rank sequences analysis 
- context of POET2: In the context of the AlignBio 2025 PETase tournament, choosing the right context for PoET-2 is critical because the model uses "retrieval-augmentation." It doesn't just look at the query sequence; it uses the context to "prime" itself on the evolutionary constraints of that specific enzyme family.

Since you are targeting the zero-shot task for PETase variants, here is the recommended hierarchy for your context:

1. The "Gold Standard" Context: Curated PETase MSA
For the best performance, you should use a Multiple Sequence Alignment (MSA) specifically containing known PET-degrading enzymes and their close homologs.

Why: PoET-2 is designed to learn from "sequences-of-sequences." By providing an MSA of PETases, you are showing the model exactly which residues are conserved across plastic-degrading specialists versus general esterases.

What to include:

IsPETase (the wild type from Ideonella sakaiensis).

Known high-performers: FAST-PETase, LCC (Leaf-branch Compost Cutinase), and DuraPETase.

Pfam sequences: Use the Dienelactone hydrolase family (PF01738) or sequences from the PETase-like cluster.
- kept 1541/1557 sequences with S129/D175/H206 intact
# **feb12**
- CURRENTLY: NEED TO filter out the tournament_test from wt1/2/3 that dont have S D or H (the unmasked/fixed residues in poet2)
- running poet2 openprotein, with whole masterdb as context, query being wt1 and input the testset of wt1, with SDH triad fixed, but hten said seq 645 must have D206 meaning we must filter it out remove it 
# **feb12** 
- 
- What about fireprot, hotspot wizard? FireProt requires a protein structure (PDB) as input. It combines energy-based (FoldX + Rosetta) and evolution-based filtering to identify potential stabilizing mutations.Hotspot: It integrates structural, evolutionary, and functional information to identify positions likely to tolerate mutation without catastrophic loss of function. Structural hotspots (HotSpot Wizard) Use it to identify mutable positions and conserved positions to avoid before scoring with ESM or ΔΔG. Structural stability filters (FireProt) Focus on stabilizing mutations predicted by FireProt as likely beneficial for expression/stability, then combine with language model and ΔΔG scoring for activity.Use them only if you want structural stability guidance, not as your main zero-shot ranking engine. For both FireProt and HotSpot Wizard, you input the wild-type PETase structure (PDB or high-quality AlphaFold model), not your 5000 mutants; they analyze the WT and output (i) predicted stabilizing single mutations (FireProt) and/or (ii) mutable “hotspot” positions and mutability/stability assessments (HotSpot Wizard).Verify catalytic residues are correct for PETase. The list you show (129, 130, 154, 175, 206) does not match canonical IsPETase catalytic triad numbering (typically Ser160–His237–Asp206 in IsPETase numbering). If your PDB is renumbered, confirm positions correspond to the true catalytic residues. Do not proceed if these are wrong.

Keep all catalytic residues marked as essential. These should be fixed (non-mutable). That is correct.

Use a high-quality structure. Prefer X-ray PDB if available; otherwise a high-confidence AlphaFold model. Make sure chain and residue numbering match your intended WT.

Do not modify pocket/tunnel parameters unless you have a reason. Default probe radius (2.8 Å) and tunnel settings are fine for standard stability design.

Understand what FireProt will output. It will propose predicted stabilizing single mutations (ΔΔG-based + evolutionary filtering). It will not rank all 5000 mutants.
- Soluprot is not designed for single point mutation effect prediction -> why? 
- What about enzyminer -> no. design makes it structurally a homolog discovery and enzyme selection framework, not a single-point mutation effect predictor.
# **feb11**
- What produces NaN in PROVEAN’s code path

The only path to NaN is: SetSubjectSequences() returns early with no supporting sequences, so sum_weights_ stays 0 and delta_score_ /= 0 becomes NaN.
That happens if SetSequenceInfoFromBlastOut() returns -1 (parse error). In that case PROVEAN just returns early and keeps going, which yields NaNs without a crash. See SequenceDB.cpp (lines 388-394).

- provean returned nan for all 3 runs all 3wt mutations. what this means -> no alignment made by psi-blast? debug. 
- even with 2 provean running oin 2 terminal, the psiblast cpu usage is halved and total cpu iusgae on mac is same -> how can we make it use all cpu cores on mac ? 
- Your prompt sees C build artifacts and activates its C module, which displays:
- It did not run. The script is failing because macOS ships a BSD getopt, and PROVEAN requires GNU “enhanced” getopt.
- dynamut2 is directly relevant for us but (1) no open github available while the web server is broken/slow (2) specializing for certain species (bacteria/plants/etc.): What DynaMut2 Uses (from the paper)

Data: ProTherm single‑point mutations plus multiple‑point mutations. They add hypothetical reverse mutations to balance classes and exclude reverse entries with |ΔΔG| > 2.0 kcal/mol.
Structures: Prefer biological assemblies; if not available, use asymmetric unit (for NMR, use asymmetric unit).
Features: Protein dynamics from Normal Mode Analysis (NMA) using bio3D, wild‑type residue environment, substitution propensities, contact potential scores (AAINDEX), interatomic interactions (via Arpeggio), plus their graph‑based signatures (atoms as nodes, interactions as edges, pharmacophore classes and distance patterns summarized as cumulative distributions).
Model: Random Forest (scikit‑learn), with incremental stepwise greedy feature selection.
Outputs: ΔΔG for single and multiple mutations (up to three), plus NMA‑based dynamics visualizations in the web UI.
Replication Plan (local code, no web server)

Data curation
Use ProTherm for ΔΔG labels. Implement reverse mutations for balancing and filter out reverse entries with |ΔΔG| > 2.0 kcal/mol. Reproduce their train/test splits if possible.
Structure prep
Parse PDB/mmCIF, choose biological assembly, map chain and residue numbering, remove alternate locations, handle missing residues.
Feature extraction
Compute these per mutation using the wild‑type structure.
NMA dynamics: run bio3D (R) to extract normal modes and per‑residue dynamics features.
Graph‑based signatures: build atom‑level graph around the wild‑type residue, assign pharmacophore classes, summarize distance patterns as cumulative distribution functions.
Interatomic contacts: run Arpeggio and count hydrophobic and other contact types involving the wild‑type residue.
Contact potentials: compute AAINDEX‑based contact scores.
Substitution propensities: encode residue substitution scores (e.g., BLOSUM or other propensity measures if you follow their prior mCSM work).
Model training
Train RandomForestRegressor on feature vectors, apply stepwise greedy feature selection, and validate on non‑redundant and blind test sets.
Prediction CLI
Given a PDB + mutation list, compute features and output ΔΔG.
Gaps You’ll Need to Decide

The exact graph‑based signature recipe isn’t fully specified in this paper. It references their mCSM graph‑signature framework, which is not detailed here.
Exact feature list and RF hyperparameters aren’t enumerated; the paper notes stepwise greedy selection but doesn’t list final features.
Multiple‑mutation feature aggregation isn’t described in detail (sum/mean/pairwise distances), so you’ll need to choose a strategy or infer from earlier mCSM/DynaMut work.
Pragmatic path forward
If you want a working local replacement quickly, we can:

Implement a “DynaMut2‑inspired” pipeline that follows the same feature families (NMA + Arpeggio + AAINDEX + substitution propensities) and train an RF model on your curated dataset.
Add graph‑based signatures later if we can find or re‑create the mCSM signature algorithm.
If you want me to start coding this, tell me:

Which structures and mutation sets you want to run first (e.g., PETase variants).
Whether you prefer an approximate DynaMut2‑style model now, or you want to invest in matching the mCSM graph‑signature method more closely.
# **feb10** 
- why does provean need 200-300Gb of ncbi-nr data: 
PROVEAN’s score is based on evolutionary evidence: it runs BLAST/PSI‑BLAST to find homologs, clusters them, then computes how each mutation affects alignment scores. To get a reliable supporting set, it needs a very large protein database (nr) so it can find enough diverse homologs for many proteins. That’s why the dataset is huge.

If you use a smaller database (e.g., Swiss‑Prot), it can run faster, but you may get fewer homologs and less reliable or unstable scores for some proteins.------> the issue with all these tools is the yrely on ncbi-nr withotu any annotation to guide/tune to the specialization of the tool. Also they need better coverage of species in which they find homologs for eg. plants we need much more sequencing and marine too. 

- Short answer: with PROVEAN v1.1.5 (which uses BLAST+ 2.4.0 per the README you pasted), you should use dbV4. Support for dbV5 starts in BLAST+ 2.8.1+, so BLAST 2.4.0 won’t read dbV5; therefore PROVEAN 1.1.5 needs dbV4. This is an inference from BLAST’s version‑5 support timeline. (ncbi.nlm.nih.gov)

What dbV4/dbV5 are (and why it matters):

dbV4 = the legacy BLAST database format (GI‑based, no built‑in taxonomy). It’s no longer updated and lives under the /blast/db/v4/ directory on NCBI’s BLAST db FTP site; archives are named with _v4 (e.g., nr_v4.*.tar.gz). (ncbiinsights.ncbi.nlm.nih.gov)
dbV5 = the newer BLAST database format (taxonomy‑aware, accession‑based/GI‑less, supports multi‑character PDB chain IDs). It requires newer BLAST+ (2.8.1+/2.9.0+) and is now the default at the root /blast/db/ directory. (ncbi.nlm.nih.gov)

- Dynamut2 -> for MISSENSE mutations 
- When checking for chains in a x-ray structure why theres 2 chains: X‑ray structures, the asymmetric unit often has two copies of the protein. It’s different because you’re counting different things now. 507/479 was from len(residues) earlier, which counts everything in the chain (waters, ligands, modified residues, etc.). 257 is from one_letter_code() which only counts standard amino acids, and then you removed X (unknown/modified residues), which shortens it further. X‑ray structures often omit disordered residues, so the resolved chain is shorter than the full sequence (e.g., 299 aa).
- GEMMI: Red flags for “weird PDB”
Any chain IDs longer than 1 character (will be remapped/collide in PDB).
>99,999 atoms (PDB atom serial limit).
Residue numbers > 9999 (PDB residue limit).
Multiple models (some tools only read the first).
Huge number of chains (>62 single‑char IDs).

- .CIF -> .PDB: Most structures convert 1‑to‑1 with geometry preserved.
Losses/changes can happen because PDB has tighter limits (chain IDs, atom/residue numbering, alternate locations, long residue/ligand names, multiple models, metadata).
If the mmCIF uses features beyond PDB limits, the converter will truncate or remap those fields.
- dynamut 2, xray vs AF compare results 
- XRAY VS AF: Use the experimental X-ray structure if (a) it’s the correct wild-type sequence/backbone you’re mutating, (b) it has good local geometry around the mutated site(s), and (c) it’s reasonably complete (few missing residues near the region of interest). Otherwise use an AlphaFold model for the exact parent sequence you’re using in the tournament/backbone.An X-ray PDB is an experimentally determined structure; an AlphaFold (AF) model is a computational prediction. The main practical differences that matter for mutation-stability tools like DynaMut2:

What “reality” they reflect

X-ray: one (or a few) conformations consistent with crystallography data, influenced by crystal packing, constructs, ligands, temperature, etc.

AF: a predicted fold for the input sequence, closer to an “in-solution plausible” conformation, but not experimentally constrained.

Confidence / uncertainty

X-ray: quality is captured by resolution, R-factors, B-factors, electron density; local regions can still be ambiguous or missing.

AF: gives per-residue confidence (pLDDT) and sometimes predicted alignment error; low-confidence loops/termini can be wrong even if the core is right.

Completeness and modeling artifacts

X-ray: can have missing residues, alternate side-chain conformations, bound waters/ions/ligands; may have engineered mutations or truncations.

AF: usually full-length, but side-chain rotamers, loop placements, and local packing can be imperfect; no crystal waters unless you add them.

Local geometry at the mutation site (critical for ΔΔG predictions)

X-ray often has more realistic local packing if the region is well-resolved and matches your exact sequence.

AF can be better if the X-ray structure is missing that region, has different sequence, or has unusual strain from crystal contacts.

Comparability across backbones

Mixing X-ray and AF can shift baselines because they differ systematically in local strain, loop placement, and side-chain packing. That’s why standardizing the source (all AF or all X-ray) improves cross-parent comparability.
- docker run --rm -p 8070:8070 lfoppiano/grobid:0.7.2
- curl -s http://localhost:8070/api/isalive
- In that case you’ll need a true long‑context model (LED/LongT5/BigBird‑Pegasus), or a Pegasus‑X checkpoint that really exposes 16k in HF.
- BigBird‑Pegasus dates to 2020. (huggingface.co)
PEGASUS‑X was released in 2022 and supports up to 16k tokens when using a true PEGASUS‑X checkpoint. (huggingface.co)
LED‑16384 and LongT5‑16384 are long‑document summarizers that are still widely used, but they’re from the earlier long‑context era. (huggingface.co)
What’s newer (long‑context, open‑weights)

Llama 3.1 (128k context). (github.com)
Qwen2.5 (128k context; many sizes, including instruction‑tuned variants). (qwen2.org)
Qwen2.5‑Turbo (up to 1M context, hosted). (qwen2.org)
Llama 4 (1M–10M context according to Meta’s model table, but huge and not practical for most local runs). (github.com)

# **feb9**
- theres either webserver platform or githubs for compbio and some bioinformatic tools (theres a difference). 
- Working on codex to batch run all the github tools on our dataset without thinking too much 
	- Paper side: writing on codex the code needed to send the 
	- Github side: need to clone all tools 
- vector institute, mcgill bioinformatics hub, apple, PA elna medical 
- mcgill ventures, bioit, scipy, stanford rnafold 
- cleaning up repo (mac petorg)

 **feb8** 
- annoying to play with ESM and E1 models wiht python. openprotein (strong platform, no github, or its the old one) and esm (yes github but minimal platform) are using marketing/branding for competition and story in market they like 

# **feb3** 
- simulate discovery of CYP450 candidature in biosynthetic pathway of a plant by batch docking all CYP450s against all metabolites found via LCMS ?
# **feb2** 
- make a sandbox platform where can upload many types of bioinformatics files like boltzlab sandbox
# **feb1** 
- 2. GROUND TRUTH REFERENCE BENCHMARK
	align2023 amylase: wt var col 
	align2023 glucosidase: 
	ProteinGym: 
	MaveDB: 
	NESG Sol: 
	Soluprot Sol: 
	Price Sol: 
	PSI Sol:
	Meltome Stab:
	FireprotDB Stab:
	ThermomutDB Stab 
	CAFA-5 Kaggle: 142k var
	Novozymes kaggle: 31k var
	Protsol Sol: 71k 
	ProthermDB
test
1. RUNNING TOOLS
		Use structural encoding (AS) 
		Blosum matrix (AS)
		CataPro 
		Diffdock, rosettafold all atom, equibind, tankbind, GINNA, UniMol, flow dock, posebench/busters, gromacs, perturbation free energy toolkit, scheodinger FEP, Apodock, karmadock. 
		[EVcouplings]
		[ESM1v, ESM2, ESM3/c]
		(profluent E5, POET2, GeorgieV)
		(MutCompute)
		(Escalante/Mosaic)
		Envision
		DeepSequence / EVE / EVEscape
		SIFT / PolyPhen-2 / SNAP / SuSPect 
		ECNet, SCANEER
		UniKP, DLTKcat, CmpdEnzymPred 
		ProPRIME
		Provean 
		deeprankgnn
		geoevobuilder
		PETase specific: MutPSSM other paper 
		PETase specific: Alam mutation score /motif integrity
	Structural/docking/mdsim 		
		aqaffinity sandbox
		Boltz2  
		FoldX generate + ddg 
		RFjoint/RFold/Rosetta generate + 
		other sructure based: PoPMuSiC, mCSM, SDM2, SAAFEC 
		other sequence-based: I-Mutant 2.0, STRUM, SAAFEC-SEQ, BoostDDG,mGPfusion

	2. Stability 
		ThermoProt 
		ddgemb(https://ddgemb.biocomp.unibo.it/), deepddg (https://protein.org.cn/ddg.html), prostab, temstapro (type2) temberture (type2) 
		dynamut2 (https://biosig.lab.uq.edu.au/dynamut2/), enzyact ()
		
	3. Expression
		RiboDecode
		TIsigner, OSTIR, Salis RBS, DeepCodon
		Soluprot 
		Procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION
		
- run tools, build benchmarkdb+features - RUNNING ... benchmarkdb extracting wt-var pairs, to run all tools on it, to have ground truth reference to evaluate features on 
# **jan28** 
- (side) plan to learn rust: 
How to learn Rust efficiently

Weeks 1–2: fundamentals with immediate payoff: ownership/borrowing, iterators, error handling (Result), traits, lifetimes only as needed. Build 2–3 tiny CLIs that read/write FASTA/FASTQ/BAM-like data and emit TSV/Parquet.

Weeks 3–5: “bioinformatics-grade” engineering: streaming IO, zero-copy parsing, compression (BGZF), concurrency (rayon), and profiling (cargo flamegraph, criterion). Rust becomes valuable when you can process files without blowing RAM and can safely parallelize.

Weeks 6+: interoperability-first: wrap Rust as Python modules with pyo3/maturin so you can drop Rust into existing notebooks/pipelines without rewriting your stack. This is usually the highest ROI pattern in computational biology.

Where Rust is genuinely useful in compbio/bioinformatics/BioML

High-throughput file formats & streaming: spec-compliant readers/writers and transformations for genomics formats (SAM/BAM/CRAM/VCF/BED/GFF/GTF, BGZF, tabix, etc.)—noodles is a major ecosystem piece here.

Fast sequence parsing and k-mer/minimizer work: FASTA/FASTQ parsing and k-mer oriented processing with minimal copying—needletail is explicitly designed for this.

Embedding/search infrastructure: take PLM embeddings (ESM/ProtT5/etc.) and build a robust, fast ANN index + query service (HNSW, IVF-PQ) or offline batch join against metadata; Rust shines for latency, memory control, and safe concurrency.

Pipeline kernels: the “hot loops” in your feature extraction (k-mer stats, alignment summaries, pileup-like transforms, SIMD-friendly numeric transforms) as Rust extensions called from Python.

Production-ish tooling: reproducible single-binary CLIs, deterministic deployments, strong typing around schemas/metadata, and safer concurrency than ad-hoc Python multiprocessing.

Rust in ML (what’s real today)

Rust is not the center of model research/training in BioML, but it’s increasingly usable for inference and tooling. Options include Candle (Hugging Face) for performant inference (including GPU support) and Burn as a Rust-native DL framework. In practice: you’ll still train in PyTorch/JAX, then use Rust for serving, batch inference plumbing, or custom kernels.

A valuable/interesting Rust project in 2026 (practical + relevant)
Pick one that matches what you already do (PLM + large datasets):

“Genomics/Proteomics streaming feature engine”: a Rust CLI that ingests FASTA/FASTQ + optional BAM/CRAM + annotations, computes configurable features (k-mer/minimizers, QC stats, motif scans), and outputs Arrow/Parquet for downstream pandas/Polars; use noodles for formats and needletail for FASTX.

“Embedding index + retrieval service for proteins”: Rust service that stores ESM embeddings, supports nearest-neighbor queries, filters by metadata, and returns candidates + explanations; pair it with your existing Python that generates embeddings.

“Python→Rust acceleration package for your pipeline”: identify 1–2 bottlenecks (e.g., FASTA parsing + mutation featurization + big joins) and ship them as a pip install-able Rust extension. This is the most directly career-relevant because it demonstrates you can harden and speed up real BioML workflows.

What is ESM written in, and would it “need Rust”?

The mainstream ESM codebases are Python + PyTorch (e.g., Meta’s facebookresearch/esm for ESM-1/2/ESMFold) and the newer ESM3 codebase is also published as a Python repo.

ESM doesn’t need Rust; research iteration speed and CUDA kernel ecosystems keep it Python-centric. Rust is most sensible around ESM: fast preprocessing, dataset plumbing, embedding indexing, inference orchestration, or a high-performance service layer. Re-implementing ESM training in Rust is possible but usually negative ROI unless you’re specifically building infrastructure or kernels.
- what to do about expression from CDS-vector for align: If your PETase “expression” label is recombinant yield in E. coli, the realistic hackathon move
	1.	Fix the construct assumptions (at minimum: do you know the RBS/5′ UTR and start context, or do you only have CDS?).
	2.	Score each CDS variant with 5′ opening energy/accessibility and/or OSTIR translation initiation rate (this is usually the strongest signal).  ￼
	3.	Add lightweight CDS covariates: CAI/tAI, GC%, rare-codon runs/cluster metrics; optionally a solubility proxy at the AA level (orthogonal to “expression”).
	4.	If you need “one backbone,” treat OSTIR/TIsigner as the backbone score and the others as small additive features.
- Tool idea: takes input list of CDS variants with one vector backbone and predicts expression in different parameters like ecoli/etc. 
- PLM-gpu Project ideas:
	Local PLM embedding factory (MPS-accelerated)
	Build a CLI + cache that takes FASTA(s) and emits standardized embeddings (per-residue + pooled), with batching, resume, and metadata; supports ESM-family checkpoints you use; writes Parquet/Arrow; includes deterministic hashing so reruns are free. Value: turns “run ESM on 5k–500k sequences” into an industrial primitive you can reuse across every project.

	Mutation-effect benchmark harness for protein tasks
	A clean framework that ingests multiple deep mutational scanning datasets / assay tables and evaluates any scorer (PLL/LLR, embeddings + linear model, stability predictors) with consistent splits, metrics (Spearman, AUROC, NDCG), and leakage checks. Value: makes your modeling claims defensible and speeds iteration; immediately useful to teams.

	Structural surrogate pipeline without full AlphaFold
	Instead of full AF-at-scale, use “cheap structure signals”: ESM-based contact/structure heads (if available), secondary-structure/disorder predictors, coevolution proxies, and simple graph features; optionally add AF only for top-K. Value: 80/20 structure-informed ranking with 10–100× less compute.

	End-to-end “protein scoring service” (MLOps)
	A small service (FastAPI) that exposes endpoints: /embed, /score, /explain, with a job queue, artifact store, and versioned models; includes a web dashboard that shows distributions, outliers, and provenance. Value: converts notebooks into a shareable product; strong for real-world engineering credibility.

	Data-quality + provenance layer for bio pipelines
	A library that enforces schemas for FASTA/TSV/assay tables, validates sequence lengths/alphabet, detects duplicates/leakage, tracks transforms, and emits a “data card” report (hashes, counts, missingness, mutation density). Value: prevents silent errors; teams will reuse it.

	Active learning loop for variant selection
	Implement “pick next variants” with uncertainty + diversity (e.g., ensemble variance, embeddings clustering, constrained design regions), and simulate on historical assays to show reduced lab budget for same performance. Value: connects ML to experimental ROI.

	Apple-silicon-optimized inference cookbook (real, not hype)
	A reproducible guide + scripts showing when MPS helps (batch size curves, dtype choices, throughput), with profiling and fallbacks; includes a “device sanity check” to prove you’re on MPS. Value: saves you and others days; practical and shareable.
- Q: would generating alphafoldstructures (at scale, optimized pipeline) able to use metal mac
	A: AlphaFold2 is JAX-based; JAX can target Apple GPUs via Apple’s Metal plug-in (jax-metal), but that path is explicitly experimental and not all JAX functionality is guaranteed to work. In practice, most people doing serious batch structure generation still run on CUDA GPUs; on macOS you’re usually limited to CPU runs or experimental/community ports/forks (which can work for demos but are not a stable “optimized pipeline” for scale).
- Q: does esm models all that able to train on mac metal ? i know it does its thing on python but what about bare metal mac and would htere be improvements and is that useful? 
	- A: Yes—ESM (PyTorch) can run and be fine-tuned on Apple Silicon using Metal via PyTorch’s mps backend, i.e., you move the model/tensors to device="mps" and it executes on the M1 Max GPU. In practice, it’s most useful for inference/embedding generation and small fine-tunes; it’s usually not competitive with CUDA for large-scale training, and you can hit backend limitations (e.g., no float64, incomplete operator coverage, and single-device training constraints). For “bare metal” improvements: you’re basically choosing between (a) PyTorch MPS (easy drop-in, good for prototyping) and (b) rewriting/porting to Apple-first stacks (Core ML / MLX) which can be faster for some workloads but requires more engineering and won’t be a drop-in for ESM code. If your goal is your PETase workflow: MPS can be a real win for batch embedding generation locally; for heavy training or big sweeps, you’ll still typically want CUDA machines.
- run tools, build benchmarkdb+features 
# **jan27** 
- Pandas 3.0 just got released, what does this mean for computational systems biology: 
Core improvements in Pandas 3.0
	•	Copy-on-Write (CoW) becomes the default and only mode, making DataFrame modifications predictable and reducing unnecessary memory duplication. This eliminates the old SettingWithCopyWarning and can reduce memory footprint in large pipelines.  ￼
	•	A dedicated string data type (str) by default replaces the older NumPy-object storage for text columns. This improves type safety and performance on string operations.  ￼
	•	Better integration with Apache Arrow via PyArrow can accelerate I/O and in-memory operations, particularly on columnar data, yielding speed and memory improvements.  ￼
	•	Performance optimizations across group-by, window functions, and other core operations (30–40 % speedups in many cases).  ￼
	•	Removal of deprecated APIs and behaviors, requiring some code updates (e.g., append() removed, stricter indexing rules).  ￼

Practical impacts on AlignBio2025 and systems biology projects
	•	Large dataset handling: Many bioinformatics pipelines involve multi-GB tables (variants, annotations, expression matrices). Reduced memory overhead and CoW semantics help when you transform and subset these tables repeatedly. Faster group-by and joins can speed QC and feature engineering phases.  ￼
	•	Feature engineering pipelines: Nullable types make handling missing values more consistent. This reduces workarounds for NaN/None mixing, which matter when preparing features for ML models (e.g., annotation matrices or phenotype tables).  ￼
	•	IO performance: Tighter PyArrow integration benefits reading/writing Parquet or Arrow formats at scale (common in biobank-scale or multi-omics workflows).  ￼
	•	Reproducibility and debugging: More predictable indexing and copies mean fewer silent bugs in data transformations, improving reliability of preprocessing and analysis scripts.  ￼

Considerations before upgrading
	•	Code changes may be required to handle API removals and stricter behaviors.  ￼
	•	Test performance and results against your current pipelines to confirm gains outweigh migration effort.

# **jan26** 
- looking into Mosaic/Escalante
- RUNNING ... benchmarkdb extracting wt-var pairs, to run all tools on it, to have ground truth reference to evaluate features on 
- ran evcouplings on masterdb1, waiting for results, 
- After Amber MD (as a late-stage, top-N feature), the highest-leverage next 3–4 tools for zero-shot ranking are:
	1.	FoldX ΔΔG (stability/packing, fast, scalable to thousands)

	•	Output features: ΔΔG_fold (mutant−WT), local interaction energy around active-site/loop residues, clash flags.
	•	Why priority: cheapest structure-based signal that correlates with foldability and often with expression/activity viability when you’re exploring many mutants.

	2.	Solubility/expression predictor (SoluProt + one of Protsol / NetSolP / Protein-Sol-type)

	•	Output features: predicted soluble expression probability/score; aggregation propensity proxies.
	•	Why priority: Align’s targets include expression; this is usually orthogonal to PLM “fitness” scores and catches obvious liabilities.

	3.	Thermostability predictor (ThermoProt + one other stability model like TemStaPro/DynaMut2 class)

	•	Output features: ΔTm / stability score; consensus/mean across methods; disagreement as uncertainty.
	•	Why priority: gives a second, non-FoldX stability channel; combining a physics-ish ΔΔG with learned stability often improves rank robustness.

	4.	Docking / binding proxy to PET (quick docking for top-N, not all 5000)

	•	Output features: docking score, contact counts to catalytic pocket, substrate pose stability over short restrained MD, MM/GBSA ΔG_bind on a small snapshot set.
	•	Why priority: it’s your most direct “activity-ish” signal tied to PET interactions; do it only after triaging candidates with PLMs + stability + solubility.

Practical sequencing for the pipeline: (A) run FoldX + solubility + learned stability on all 5000, (B) take top 200–500 (plus a diversity slice), (C) do docking + short Amber MD/MMGBSA at pH 5.5 and 9, and use pH-differences as features.
- running masterdb1 into evcouplings 
- meeting with sanju 
	- went over interpreting evc results
		- run on masterdb1, benchmarkdb 
	- priority need to run 
	(1) amber MD-SIM 
	(2) FoldX DDG 
	(3) Soluprot, Protsol, NetsolP (4) Thermoprot, Dynamut2, TemStaPro 
	(5) Docking + MDSim  
# **jan25** 
2. GROUND TRUTH REFERENCE BENCHMARK
	align2023 amylase: wt var col 
	align2023 glucosidase: 
	ProteinGym: 
	MaveDB: 
	NESG Sol: 
	Soluprot Sol: 
	Price Sol: 
	PSI Sol:
	Meltome Stab:
	FireprotDB Stab
	ThermomutDB Stab 
	CAFA-5 Kaggle: 142k var
	Novozymes kaggle: 31k var
	Protsol Sol: 71k 
	ProthermDB
test
1. RUNNING TOOLS
		[EVcouplings]
		[ESM1v, ESM2, ESM3/c]
		(profluent E5, POET2, GeorgieV)
		(MutCompute)
		(Escalante/Mosaic)
		Hotspot Wizard
		Fireprot 
		Envision
		DeepSequence / EVE / EVEscape
		SIFT / PolyPhen-2 / SNAP / SuSPect 
		ECNet, SCANEER
		UniKP, DLTKcat, CmpdEnzymPred 
		ProPRIME
		Provean 
		deeprankgnn
		geoevotranfsormer
		PETase specific: MutPSSM other paper 
		PETase specific: Alam mutation score /motif integrity
	Structural/docking/mdsim 		
		aqaffinity sandbox
		Boltz2  
		FoldX generate + ddg 
		RFjoint/RFold/Rosetta generate + 
		other sructure based: PoPMuSiC, mCSM, SDM2, SAAFEC 
		other sequence-based: I-Mutant 2.0, STRUM, SAAFEC-SEQ, BoostDDG,mGPfusion

	2. Stability 
		ThermoProt 
		ddgemb(https://ddgemb.biocomp.unibo.it/), deepddg (https://protein.org.cn/ddg.html), prostab, temstapro (type2) temberture (type2) 
		dynamut2 (https://biosig.lab.uq.edu.au/dynamut2/), enzyact ()
		
	3. Expression
		RiboDecode
		TIsigner, OSTIR, Salis RBS, DeepCodon
		Soluprot 
		Procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION
		
- for the evcouplings align part: Use MMseqs2 instead of jackhmmer (often 10–100× faster)
- A evc_epistatic_sum
	•	What it is: predicted ΔE under the pairwise Potts model (PLMC) including couplings (J terms).
	•	Intuition: “Does this mutation fit the coevolution constraints?”
	•	Strength: captures context-dependent compatibility (good for function if function is constrained by networks of residues).
	•	Weakness: needs a decent MSA; if MSA is weak, scores are noisy or collapse.

B evc_independent_sum
	•	What it is: predicted ΔE under an independent model (fields only, no pair couplings).
	•	Intuition: “Is this mutation generally allowed at this site in evolution?”
	•	Strength: strong for “don’t break the protein” / general viability.
	•	Weakness: less sensitive to epistasis, can look like a conservation proxy.

C evc_freq_sum
	•	What it is: MSA observed frequency of the substituted AA at that position.
	•	Intuition: “Has evolution already tried this exact substitution?”
	•	Strength: very interpretable, strong “sanity check” feature.
	•	Weakness: sparse; many substitutions will be 0.

D evc_colcon_sum (column conservation)
	•	What it is: position-level conservation metric (same for all substitutions at same position).
	•	Intuition: “How constrained is this site overall?”
	•	Strength: good for “mutations in conserved sites are risky.”
	•	Weakness: does not distinguish V26A vs V26C (as you noticed).

- 	1.	Substitution never observed in the MSA at that position
Your frequency=0.0 for many of those rows confirms this.
When frequency=0, EVcouplings often assigns a similar “very bad” penalty for many rare AAs → same prediction_*.
	2.	Column is highly conserved
Your evc_colcon_sum is identical for all mutants at pos 26 because it’s a position-level value, not mutation-level.
So it will always be the same for all 19 substitutions at that position.
	3.	Regularization / pseudocount effects (plmc)
When the data doesn’t constrain a mutation, the model defaults to similar energies for multiple substitutions.
- Option A (recommended): treat “missing” as a feature, don’t force-fill. Option B: rebuild EVcouplings so it scores more positions (proper “fix”).  Biological interpretation: the MSA doesn’t support learning constraints there (not enough homologous evidence at that residue).. Biologically it means: EVcouplings has no reliable evolutionary constraint signal for that specific substitution at that position (in that WT background / alignment), so it cannot score it. It’s not “neutral” or “bad” by itself — it’s “out of model / unsupported”. EVcouplings only produces scores for positions that survive:
	•	focus region selection (region, first_index)
	•	gap/coverage filtering (minimum_sequence_coverage, minimum_column_coverage)
	•	redundancy filtering (hhfilter, seqid_filter)
	•	sometimes additional internal trimming / lowercase-column logic

So “not scanned” = not enough usable evolutionary information at that position, or it was excluded by preprocessing.

- You can now score any test variant by summing its mutations’ prediction_epistatic (or prediction_independent) from the correct WT table.

Lookup key

Use (pos, wt, subs) or just mutant string like G33A.

Example row:
	•	pos=33
	•	wt=G
	•	subs=A
	•	prediction_epistatic=1.206…

Multi-mutant scoring (simple baseline)

For a variant with mutations ["G33A","S160T","W185F"]:
	•	evc_score_epi_sum = sum(prediction_epistatic for each mutation)
	•	evc_score_ind_sum = sum(prediction_independent for each mutation)
Also store:
	•	min / max / mean
	•	count of mutations found vs missing
- using the reuslts of evcuoplings and its MSA and mutate matrix of all single point mutations instead of running locally but will try again locally once downloaded uniref90 instead of using swissprot for the evcouplings pipeline locally 
- notes on evcouplings:
 on engine: local these fields are mostly used to size internal parallelism and satisfy required config keys; they don’t “reserve” RAM/cores like a scheduler would, and (2) whether you actually get multi-core speedups depends on the underlying tools (notably whether plmc was compiled with OpenMP; your earlier Mac compile issues suggest it may be single-threaded there). Runtime varies widely with alignment size and the number of sequences; with align: existing the “align” stage is mostly preprocessing and should be relatively quick, while PLMC inference can dominate if the MSA is large.
Workflow (what each stage does): align = produce/ingest an MSA and derive a “focus” alignment around your target (region/first_index mapping, gap filtering, optional redundancy filtering); couplings = fit a Potts model (pairwise Markov random field) to the MSA; in practice EVcouplings calls plmc which does pseudo-likelihood maximization and writes a model (fields h and couplings J); mutate = score mutation(s) by computing ΔE (change in statistical energy under the fitted Potts model) for each mutation string you provide. The environment.time/memory/queue entries are “resource request” fields for cluster submitters; for engine: local they are not enforced reservations, but EVcouplings’ code path requires the keys to exist. You can boost cores/cpu in the YAML for the PC, but again: you only get real speedups if the underlying executable is threaded (e.g., OpenMP-enabled plmc).
Yes—you can provide your own MSA: set align.protocol: existing and point raw_alignment_file to your alignment (Stockholm/A2M/FASTA depending on what EVcouplings expects for that path). Compatibility requirements: the MSA must include the target sequence, the target should correspond to the region/first_index you specify, and it should be a clean amino-acid alignment (consistent length, gaps as -, no weird characters). EVcouplings may still build an internal “focus” alignment and apply some bookkeeping, but it won’t do the database search if you use existing.
frequencies_file is derived from the (focus) MSA: it’s the per-position residue frequency table used downstream by the couplings/contact machinery; it’s not “extra data” beyond the MSA. plmc is the core inference engine: it estimates the Potts parameters. lambda_J and lambda_h are regularization strengths (L2 penalties) on couplings J and fields h; higher values shrink parameters more (reducing overfitting but potentially underfitting).
- debugging evcouplings ...
- evcouplings needing jckhmmer, need to find whats at the root of jckhmmer and how many tools depend on it + for other dependancies than jckhmmer 
# **jan24**
- meeting with sanju
	- RUNNING TOOLS - EVCOUPLINGS 
		- look into the plmc requirement to oisntall evcouplings 
- note on results from the accuracy of esm1v/2/3 features on the labels of the train.csv 
	•	Activity signal: LLR-heavy scores (S_llr, S_risk) can be very strong on at least one enzyme (glucosidase) and are not purely one-hotspot (31 unique positions in top50).
	•	Stability/thermal signal: ΔPLL (S_dpll) gives a measurable melting-point lift and shows plausible stabilizing hotspots.
	•	Expression: mixed/LLR-leaning scores help on amylase; glucosidase expression labels aren’t informative.

# **jan23** 
- need to look more into the align2023 github for scoring approach and zero shot approaches  
- now also but not alot doing lit review, slides, github organization ....................  
- now also adding more PLM features (Profluent-E1, POET2, T5, GeorgieV) as writer.py files and test in the buildfeature phaseanalysis pipeline 
- now trying the label verification feature from phaseanalysis modeC, using features from esm1v/2/3 labels from integrated_benchmarkdb to fine tune the weights of the scorer.py module on the benchmarkdb and apply it to the align tournament test set. Also implement comparison statistical of the masterdb1 to the benchmarkdb to understand where the petase data lies in the distribution of benchmarkdb. 
- meeting tomorrow to implement a tool fully 
- note on esm1v/2/3 vs msa/alignents: For contrast with PLMs, classic MSA/covariation fitness models (e.g., EVmutation/EVcouplings-style) that also estimate evolutionary constraints but via alignments/statistical physics rather than neural LMs.
- phaseanalysismodeb + esm1v/2/3 results:
	Right now, you are ranking variants by ESM-derived signals, for example:
		•	llr_mean
	→ “On average, how much do ESM1v/2/3 prefer this mutation over WT?”
		•	delta_pll_mean
	→ “How much more likely does ESM think this mutant sequence is than the WT sequence?”
	
	Because you don’t know which signal corresponds to real activity or expression yet.
	So you ask:
		•	“What if I rank by LLR?”
		•	“What if I rank by ΔPLL?”
	Each choice produces a different ordering of the same variants.

	Mode B asks if I rank variants using signal A vs signal B, do I mostly get the same top variants?
	Ranking by llr_mean gives Top-200 set A, Ranking by delta_pll_mean gives Top-200 set B, Overlap = 34 variants. 

	Why this matters because you will submit a ranked list/ Even without knowing activity:
	•	Rank #1 is “best guess”
	•	Rank#200 is “still worth testing”
	•	Rank#3000 is “probably bad”

	The Jaccard ≈ 0.09 for delta PLL vs LLR -> Ranking by LLR and ranking by ΔPLL choose almost completely different mutants.so ranking defines which mutants get experimental attention.

- analyzing the results with phaseanalysis:
	- Raw esm*_pseudo_likelihood should look backbone-shifted (3 peaks in esm1v/esm3 and 2 in esm2)
	- esm*_delta_pll should be much more centered/overlapping.
- made the phase_analysis give grouped histograms now 
- 	•	LLR is usually more “local” (site-specific preference / mutation effect style).
	•	delta_pll is a “global” log-likelihood change relative to WT, and can reflect distributed effects (including context around the site) and how the model scores the whole sequence.
- organizing this 3 module workflow (1) build_features.py (2) phase_analysis.py (3) scorer.py
- workflow right now: 
	build_features.py → build_esm_consensus_features.tsv (features) 
	•	build_features.py stays a deterministic “feature assembly + QC normalization” step.

	rank.py → esm_only_prelim_ranks.tsv (define scoring logic)
	•	rank.py becomes the experimentation surface: you can change weights/penalties/objectives quickly without rebuilding ESM features every time. (same rows + added score columns like S1/S2/S3/S4). The way it scores is based on how you define it. It doesn’t automatically use all features. You decide which features matter by defining score formulas (initially simple; later weighted combos). The rest of the feature columns can remain in the output for debugging and analysis, but they don’t affect ranking unless referenced in a score. By the way you could also use just phasexanalysis to obtain scores for example python phasex_analysis.py --in esm_only_prelim_ranks.tsv --score-cols S1,S2,S3,S4 --topk 200. 
	
	phasex_analysis.py --score-cols ... → analyzes those score columns (topK overlap, top candidates per score, etc.)

- can we ask actual protein engineering ML teams and bioML teams how they organize their summary/reports/folders/meetings ?? -> ask charlie and sanju 
- can we look deeper into why esm1v works so well for variant calling and is different from esm2/3 ? 
- results on the buildfeatures(esm1v/2/3)+phase dashboard -> The blocks behave as expected:
	•	llr_mean/median/min/max are all highly correlated with the per-model LLR columns → they’re mostly redundant summaries of the same signal.
	•	delta_pll_mean/median are highly correlated with the per-model esm*_delta_pll → good, the WT mapping/delta construction is working.
	•	The “rank_*” columns are strongly anti-correlated with their underlying values because ranks are in the opposite direction (better value → smaller rank). That’s normal and not a bug.
- module 1 buildfeatures note: do not accidentally making one feature dominate by scale (handled later via normalization/weights).
- note on what the autoflip PLL thing we ran into was: auto_flip_pll_consensus exists because PLL-like outputs are not guaranteed to have the same sign convention across different writers/models/pipelines. If one model’s “pseudo_likelihood” is effectively the negative of another model’s (or you accidentally exported “loss” for one and “log-likelihood” for another), then averaging them is meaningless: they will anti-correlate and cancel.
# **jan22**
- today will work on 
	1. lit review
	2. play with features in petase_pipeline 
	3. play with weights in rank_from_esm_only 
- pipeline workflow is (1) build_esm_consensus -> data integration + feature engineering (2) phase_analysis -> what features we should keep/remove/weight differently 
- fixing build_esm_consensus and phase_analysis .py files to interpret the results / analyze the LLR and PLL features from esm1v/2/3
- Q: did align do a DMS to generate their test set of 5000 variants/ Most likely yes: it’s “DMS-like” (single mutants), but not a full DMS.
•	DMS-like = lots of single substitutions across many positions.
•	Not full DMS = they didn’t include every position × 19 AA (that would be far >5,000 unless they restricted positions heavily).
And no, it’s very unlikely they did iterative rounds (like evolve → pick winners → mutate again) for the test set generation. That kind of multi-round design usually produces:
	•	multi-mutants,
	•	clustered lineages,
	•	clear “families” of variants,
	•	and a non-uniform mutation-count distribution.
- (1) interpret rank_from_esm_only.py vs generate_rank.py 
- (2) work on features, ranking/weights, 
- (3) lit review
- (4) run tools 
# **jan21**
- confusing the submission requirement and the realtive ordering ranks 
	- *_rank_esm_only = true ranks (1 = best) from sorting the proxy score.
	- *_pred_esm_only in esm_only_prelim_submission_like.csv = submission-like numbers made by mapping those ranks/percentiles into the required output ranges: activity_1 → 0–10 activity_2 → 0–10 expression → 0–2 
	- They’re monotonic with the score (better score ⇒ higher “pred”), but not calibrated to real units.
- we can start getting distribution of mutations that are ranked top from different model iterations with diff features and weights to see if they land in hot spots? 
- starting preliminary ranking using esm1/2/3 llr pll and consensus stats as features for 3 ranking logics (act1 act2 expression). We can add extra features to feed into these 3 ranking logics, and play with the weights/logic  
-  was doing phasex analysis of esm1/2/3/ consensus rankings and i need to look more into	Most single mutants introduce at least one residue that the model thinks is less likely in that context, which drops the average PLL.And the PLL is averaged over the whole sequence, so even one bad position can noticeably reduce the mean.
- finished esm1v/2/3 consensus LLR and PLL, now looking at petase_pipeline 31 features. 
- put local esm3 usage in petase_pipeline.py
- nvidia-smi: is a client tool that talks to the NVIDIA driver through NVML (NVIDIA Management Library). It only shows real numbers if the process you’re running it from can “see” the physical GPU + driver stack.
- wrapping up esm1/2/3 writer and tester .py files. 
- CHARLIE SANJU MEETING NOTES 
	- for next few days, work more on features and tools 
	- MIT BIO IT stuff: 
		- feb 2 submission abstract, send an email about our idea cfde-trainingcenter@orau.org -> we did we asked the questions sent email 
	- 2. ENSEMBLE RANKING: need to tune the weights in charlie4/generate_submission.py  -> debug ranking 
	- 1. FEATURES: petase_pipeline.py -> gives features, we add features on that, everyone can contribute (eg. i need to try esm1v in there)
- now continuing with (1) lit review (2) tool notebook analysis runs (3) charlie4 and charlie2 analysis codes 
- wrote esmc_tester.py it works now, ok so now we have esm1v/2/3 testers that work. 
- huggingface is to give weights easily to donwload on ur machien isntead of using lets say forge api for esm3 
- need to finish this quick so i can get to looking at charlie2 and charlie4 folders 
- wrote esm3_concurrent.py, but it keeps failing with 401 error, the insight was So what happens is: thread A sends request with auth header thread B overwrites / clears session state requests start going out without auth server returns 401 for everything This is a classic “shared client across threads” failure.
	- CONCLUSION FINAL: JUST DO ASYNC THREADED LIKE THEY SAY IN THE GITHUB https://github.com/evolutionaryscale/esm  
	- CONCLUSION FOR NOW JAN21: Your account is rejecting parallel sessions (server-side), and you should stay with the non-concurrent esm3_tester.py (or reduce to --workers 1
	- esm.sdk.client("esm3") is NOT thread-safe and/or it mutates shared auth headers/state internall. That’s why it can fail even with workers=1 sometimes if the code path differs (warmup vs main loop timing), but with workers>1 it becomes consistent. When you add concurrency, requests race and the Authorization header/cookie/session gets clobbered → server sees “missing/invalid token” → 401.
	- The only difference is how many requests are in-flight at once: esm3_tester.py: 1 request at a time (serial) esm3_concurrent.py: N requests at a time (--workers N) (parallel)The difference between esm3_concurrent and esm3_tester is how the client is authenticated + how concurrency interacts with the ESM SDK. because m = client("esm3", token="DUMMY") meanas the function accepts the token but doesnt mean the client is authenticated with the key. So concurrency can increase throughput only if the server + your account allow parallel requests without rate-limiting/401/429.
	- so how to incerase throughput on esm3 remote? (1) esm 3.2.3 (latest ver) doesnt support batching in encode() (2) parallel requests using from concurrent.futures import ThreadPoolExecutor. 
	- esm3 remote throughput -> 21 seq/s To get real scaling, we need an API call that accepts multiple sequences in one request (true server batching). Your current SDK path is not doing that.
	- esm3_tester takes time to debug... 
	- esm2_tester works 
	- now writing esm2 and 3 writers and testers 
	- wrote esm1v_tester_mac 
	- ok so definitive on windwpws+ cuda the gpu throughput is 18 seq/s (of seq len around 250)
	- making esm1v_tester and esm1v_writer for WINDOWS + CUDA and then for MAC, and then for esm2 and esm3 
(1)charlie sbumision rankings 
(2) test esm3 batch mode for charlie code
(3) check charlie2 and charlie4 notebooks 
- annoying issue with connecting via ssh to tailscalepc from mac now its not working? also had issue with the right password to connect via termius via tailscalepc via iphone.. --> for the mac pc_tailscale thing, its cuz tailscale was not up on the mac-end. so do sudo tailscale up. then itll work. 
- workflow i like, testing the outputs from others, continuing lit review, running tools/deeper notebook analyses, cuda gpu stuff, slides 
- still got 1. lit review 2. run tools as my priorities 
- can use HPC from lab 
- first iteration of submission ranks analysis gave me all 5 submission types give the same top500 ranks ?? these plots/consensus outputs are telling you “the 5 submissions are basically the same ranking (at least at the top)”. That’s why it feels useless — because there’s no disagreement signal to visualize. --> treat the 5 submissions as one ranking, not five. They don’t add information right now.
- 
# **jan20**
- making biochatgpt but with CUDA-integrated workflows and code the whole way? -> more notes in cudaproject.txt 
- *logging proper chatgpt prompts and answers - to learn for autonomous systems that will be coming out *
- charlie's submission rankings 
- continuing lit review on manuscript_v1 
# **jan19** 
- lit review on manuscript_v1 on drive, then notebook analyses for tools below 
- leaderboard format possibility 
-  add rosettafold/rfjoint/rosetta to use generative modelling and generate the structures of the variants and analyze their ddG 
-  tools summary:
	0. Mutatagenesis/variant effect predictor specific 
		Envision
		DeepSequence / EVE / EVEscape
		SIFT / PolyPhen-2 / SNAP / SuSPect 
		MutCompute
		ECNet, SCANEER
		UniKP, DLTKcat, CmpdEnzymPred 
		ProPRIME
		MutCompute	
		MutComputeX
		EVcouplings
		Provean 
		deeprankgnn
		geoevotranfsormer
		PETase specific: MutPSSM (from that one paper i forgot author name)
		PETase specific: Alam mutation score / motif integrity
		aqaffinity sandbox
		proteingym
		Boltz2 
	1. DDG 
		FoldX generate + ddg 
		RFjoint/RFold/Rosetta generate + 
		other sructure based: PoPMuSiC, mCSM, SDM2, SAAFEC 
		other sequence-based: I-Mutant 2.0, STRUM, SAAFEC-SEQ, BoostDDG,mGPfusion

	2. Stability 
		ThermoProt 
		ddgemb(https://ddgemb.biocomp.unibo.it/), deepddg (https://protein.org.cn/ddg.html), prostab, temstapro (type2) temberture (type2) 
		dynamut2 (https://biosig.lab.uq.edu.au/dynamut2/), enzyact ()
	3. Expression
		Soluprot 
		Procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION
	4. PLMs
		* E5, Poet2, ESM1b, ESM1v, ESM2, ESM3, T5, georgiev
	GRAPE (?)
	consolidated results from all tools
	* docking output results from vina


# **jan18**
- set up the GPUs (I want to learn a lot more about low-level optimization of deep learning models to make the most of the available hardware) not necessarily making the most of whatever hardware I had at my disposal. As such, I would often wonder whether I should focus my attention on squeezing out better model performance or optimizing the runtime and resource usage. 
- Single point mutation summary by chatgpt
	1. understanding lit review so far on single point mutations on enzymes
	What the literature converges on is that single amino-acid substitutions usually act through a small number of physical “levers” (folding stability, conformational dynamics, active-site geometry/electrostatics, and cellular handling), and the observable phenotype (activity/expression/solubility) depends strongly on where the mutation is and what environment/host you measure in.
	Most single mutations are destabilizing (thermodynamically): across many proteins, the distribution of ΔΔG for random substitutions is skewed toward positive (destabilizing) values, and this statistical shape looks broadly similar across proteins (“universal” distribution idea).  ￼ Related, mutations that change or create new function are often destabilizing on average, because they perturb packing and active-site neighborhoods and/or increase dynamic heterogeneity; proteins can tolerate this only until a stability “buffer” is exhausted.  ￼

	Stability, expression, solubility, and activity are coupled but not identical traits. A destabilizing point mutation can reduce the folded fraction at equilibrium and/or slow folding/increase misfolding, which tends to lower soluble expression and increase aggregation/proteolysis in cells; however, stabilizing mutations do not guarantee better solubility or expression, and several analyses note that many stability-prediction workflows over-select designs that can have solubility liabilities (a practical “stability–solubility” tension).  ￼ Deep mutational scanning (DMS) studies show these couplings directly at scale: many mutations decrease measured “fitness” (which often conflates expression and activity), while some improve one axis while harming another, revealing tradeoffs and multiple mechanisms.  ￼

	For catalytic activity specifically, a single mutation can have (i) direct active-site effects (changing catalytic residues, substrate positioning, hydrogen-bond networks, electrostatics, water structure), or (ii) indirect allosteric/dynamical effects (shifting conformational ensembles, loop motions, access tunnels). This is why activity effects are harder to predict than stability: you can lose activity without large ΔΔG (pure mechanistic/dynamical disruption), or gain activity at a stability cost. Large experimental maps (thousands of missense variants) explicitly quantify these activity–stability landscapes and show that the tradeoff is common but not absolute: some regions/mutations improve both, many harm both, and the joint distribution depends on structural context.  ￼

	Context dependence is a first-order theme: the same single substitution can look neutral, beneficial, or deleterious depending on temperature, pH, substrate, cofactors, crowding, and expression level/host (proteostasis). Multi-condition mutational scans and “expression × coding” epistasis experiments show that environment and expression regime can reshape the mutation-effect landscape—important if you care about both activity and expression under specific assay conditions.  ￼

	2.  Understanding of effect of single point mutation on PETase enzymes 
	there are many PETase mutation-effect studies, but they are mostly (i) targeted single-point “rational” mutations and (ii) directed-evolution/ML-derived variants that are then dissected with a small set of single/revertant mutants. I’m not aware of a public, true deep-mutational-scan (all 19 substitutions at most positions with activity + expression + stability readouts) for IsPETase in the way DMS is done for many soluble enzymes; the PETase landscape is instead inferred from repeated convergence of mutations across engineering campaigns and mechanistic follow-ups.
	What we understand about the PETase landscape (IsPETase and engineered derivatives):
		1.	Stability is a gating constraint, and “apparent activity” often rises mainly because stability enables hotter assays where PET is more degradable. Many of the highest-performing PETases are built on a small set of thermostabilizing substitutions (classically S121E, D186H, R280A; and in some lineages disulfide-style N233C/S282C) that raise Tm and allow operation at higher temperatures, which boosts PET chain mobility and hydrolysis rates.  ￼
		2.	A small number of hotspots recur across independent lineages, implying strong structural “levers.” Reviews and primary engineering reports repeatedly mention mutations around (a) loop/active-site entrance regions, (b) positions that change local H-bonding/water-mediated networks, and (c) positions near the binding cleft and “subsite” shaping residues. For ThermoPETase (S121E/D186H/R280A), structural interpretations attribute gains largely to added hydrogen-bonding interactions and altered local hydration networks rather than dramatic active-site rewiring.  ￼
		3.	Activity changes can come from conformational-ensemble shifts, not just “better binding” in a static structure. A mechanistic ACS Catalysis study showed a single-point mutant can increase the population of productive ground-state conformations (a conformational-selection framing), which is consistent with PETase being sensitive to loop dynamics and binding-site conformational heterogeneity.  ￼
		4.	FAST-PETase and related “high-performance” variants illustrate a common pattern: combine a stability backbone with additional mutations that tune function at ~50–55 °C. Secondary sources summarizing FAST-PETase list the five substitutions relative to IsPETase (N233K, R224Q, S121E, D186H, R280A), and multiple later studies focus on why temperature and dynamics matter for this variant’s behavior.  ￼
		5.	Expression/solubility is context-dependent and can be engineered orthogonally (sometimes via PTMs). A notable thread is yeast expression + glycosylation engineering of (thermostable) PETase variants, where glycosylation state measurably affects thermostability and PET-degradation performance—i.e., “expression system biology” is part of the landscape, not just the amino-acid sequence.  ￼
		6.	High-concentration/crowded conditions reveal aggregation and non-ideal behavior that can dominate kinetics and “expression-like” outcomes in practice. A JBC 2024 kinetics/aggregation analysis reports aggregation tendencies under crowding for WT and a thermostable TS-PETase variant, which is a reminder that some mutations improve intrinsic stability but do not automatically solve formulation/aggregation constraints.  ￼

	If you want a practical “landscape prior” for PETase single mutants: (i) expect most random substitutions to reduce function; (ii) stability-raising mutations are enriched on surface/loop regions that improve packing, salt-bridges, and hydration networks without breaking catalysis; (iii) mutations near the cleft/loops often show large, condition-sensitive effects (beneficial or harmful) via dynamics and substrate engagement; and (iv) expression/solubility can be driven by the host system and post-translational effects at least as much as modest ΔΔG shifts.  ￼	

# **jan17**
- state machines kris
- randomized tests tiger beetlew
# **jan16**
- in the lit review, need to include how different classes of bioinformatic tools give what information in respect to single point mutations in enzymes (eg. what are the assumptions/limits/space of the algorithms behind each bioinformatic tools)
- 

# **jan15**
1. charlie sends ensemble ranking based on his pipeline -> will then analyze the top candidates from them according to literature 
2. repeat analysis script for EVcouplings ... remaining tools 
3. finishing esm1v.ipynb notebook 
	- comments and results 
	- compare esm1b/2/3, T5, georgiev, poet2, E5
- continue abstract on cudaproject.txt
- start literature review (witman han etc.)

- ADIL MEETING:
	Use fewer tools not more, start with simple, then can layer,
	You do need reference to effectively test out hypothesis on what makes a good model, small datasets have limited information
	Theres other datasets, if youre gonna make a new complex technique you need a good benchmarking way  
	What is the expected limitations and assumptions of all these tools, lot of models have not been trained to do single point mutation well or struggle eg. alphafold, number of the PLMs, docking
	Point mutations is not a trivial thing
	Maybe its not 
	Rationally layer on top
	The newer models, E1 

	Look at the datasets behind E1 POET start there 

	Do 2 models only, then think of maybe epxnading it  

	Dont think about all the infrastructure stuff, do simple notebook things get results we understand

	Then iterate iterate 


# **jan13** 
- need to make a detailed MSA of the wt and test set vs masterdb as well so that we have a good reference going forward 
- need the test set aligned residues positions mutations like we did for wt set found in masterdb 
- Remember to finish the mutations lists in masterdb2 to include Han cui erickson lu etc. 
- TO DO
	* E5, Poet2, 
	* esm1v 
	* esm3 
	* georgiev embedding 
	* docking output results from vina
	* evcouplings
	* deeprankgnn, geoevotranfsormer, 
	* Wittman2021 & Han & Lu & Alam & Cui & Erickson & Garcia & Hong & Norton 
	* mutcompute, mutcomputeX, provean, and protein-sol
	
	* docking, MDSim 

- test of sweep #batch size and #workers on TIME -> 
	- In a multi-process esm1v_mp.py batch-size sweep you should see the same qualitative behavior as single-process—runtime improves as --batch increases (better GPU occupancy, fewer kernel launches) and then plateaus—but you should not expect the same absolute timings or the same “best” batch size because 4 workers are sharing one GPU and will contend for SM time, memory bandwidth, and CUDA context scheduling; this typically makes the optimal per-worker batch smaller than in single-process. VRAM per worker will rise with batch due to larger activations, but weights dominate so it won’t scale linearly; total VRAM pressure is roughly procs × (model residency + activation(batch)), so you can hit OOM sooner than in single-process if you push batch too high. You should also expect more variance run-to-run from OS scheduling/context switching. Practically, you’ll usually see clear gains from batch 1→2→4, then diminishing returns somewhere around ~4–16 depending on your system, with multi-proc often preferring ~2–8 even if single-proc preferred ~8–20; the only reliable way to pick is wall-clock time for a representative subset or the full FASTA with the same --procs.
- now running multiple process esm1v script 
	- Total wall-clock runtime is governed by the slowest worker plus the one-time model load/warmup cost; batch size 10 from 5 improves throughput but does not significantly change VRAM because ESM-1v is weight-dominated, not activation-dominated
- the batch size sweep test show (single process)
	- increasing batch size doesnt speed up TIME that much. 
	- But the RAM does decrease which apparently is an artifact from RSS. RSS here would be that during the sweep at each batch we run the full loop, the loop allocates/frees/so on, and over time python sees that those arenas arent being reused, and pytorch sees the cached buffers arent needed, and OS sees the pages are cold. So we notice a drop in RAM mid-sweep. 
		- Batch size only affects temporary tensors during forward pass, GPU activations (not CPU RAM).CPU RAM dominated by Python interpreter + model metadata. 
	- And GPU alloc doesnt change. 
- so CUDA is 2x faster than MPS which is 5x faster than mac CPU (not multprocess,not fp16)
- wrote notes in the main.ipynb esm1v block 
- want to run esm1v using all C1-9 threads  
- TASK FOR JAN 14 MEETING: go through the higlighted tools and check their input/output/way its calculated 
	* ESM-1V: 
- notes on running ESM-1v is trying to allocate ~24 GB for attention tensors. On Apple Silicon (MPS), PyTorch uses unified memory, but there is still a hard per-allocation limit and fragmentation constraints. You exceeded it.Even with BATCH_SIZE = 1, you still get massive parallelism: Inside each transformer layer, Inside attention matrix multiplies, Inside feedforward blocks
That is handled by: Metal Performance Shaders (MPS) on Mac, and cuBLAS / cuDNN on NVIDIA GPUs

The two main drivers of memory here are:
	1.	Batch size (number of sequences at once)
	2.	Max sequence length in the batch

- Meeting with sanju charlie minutes same minutes and jan12_work compilation doc on the drive 
# **jan 12** 
- meeting with sanju minutes in https://docs.google.com/document/d/1rlT5fv3d3MEIsmlzb1Pukds4s6VxFLNMAoC4Ci3uVpo/edit?tab=t.0 
# **jan 11** 
- On the expression and the plasmid expression system of hte tournament: Same plasmid means the vector parts (promoter/RBS/origin/antibiotic/tag sequence if fixed) are constants, so they don’t help ranking. But you still care about features at the interface between the variant sequence and the expression system, because the mutations change the insert, and that changes how the host expresses/folds it even in the same plasmid. Concretely, the only plasmid-context features that matter are those whose value changes when the coding/protein sequence changes:

Translation / codon-level features (if CDS differs across variants)
Even with the same plasmid, different CDS → different CAI/rare codons/mRNA structure near start. If your test set is amino-acid mutants but the CDS is not specified (or all CDS are re-synthesized with standardized codons), then codon features are irrelevant; otherwise they can matter.

Signal peptide / secretion compatibility (if the construct uses secretion/periplasm)
Same signal peptide region in the plasmid doesn’t mean the processed protein is the same: mutations near the N-terminus (or cleavage site / early residues) can change secretion efficiency and processing. If mutations are only in the mature enzyme region and the signal peptide is fixed and not mutated, then this mostly becomes a constant and can be ignored.

Protein folding/solubility burden in the host
Even with identical promoter/RBS, mutations can:

	•	increase aggregation propensity
	•	destabilize the fold (lower soluble expression)
	•	expose hydrophobics, change pI, etc.
These are “expression” effects that do depend on sequence and are valid ranking features.

So the decision hinges on what exactly varies in your dataset:
	•	If variants are provided only as AA sequences and all are synthesized with the same codon optimization strategy → drop CAI/rare-codon/mRNA-folding features; keep protein-level solubility/stability/PLL features.
	•	If variants have real CDS that differs (your pipeline suggests you have CDS per variant) → codon and 5′ mRNA features can matter even with the same plasmid.


- Implementation pattern per tool (what you actually save): FoldX/Rosetta: ddg_kcal then define stability_score = -ddg_kcal; mutPSSM/LuPSSM/EvoCouplings: one scalar per variant like delta_pssm or delta_E (higher better); MutCompute: llr = log P(mutAA) − log P(wtAA) at mutated site(s), aggregate by sum for multi-mutants; ThermoProt/GRAPE/ProStab/etc.: whatever their output “Tm/ΔTm/stability probability” is, convert to “higher better” and keep raw value; Solubility tools: keep raw probability/score (higher better). Then merge everything into master_features.csv and compute the three submission scores as weighted sums with guardrail penalties: score_act1 = w*ESM_mutLLR + w*mutcompute_llr + w*stability_score + w*pssm + (optional pH5) − big_penalties(flags); score_act2 same but pH9; score_expr = w*PLL + w*CAI − w*rare_codons + w*solubility + (optional stability) − penalties. Don’t tune 20 weights; do 2–3 ablations only: baseline, +stronger stability/struct-risk suppression, +/− pssm/evo; pick the one whose top-200 doesn’t contain flag failures and looks stable. Operationally assign owners exactly as you listed, but enforce one delivery rule for everyone: “deliver features_<tool>.csv with columns: variant_id, backbone_id(optional), tool_score(s), and a README line stating whether higher is better.


- Priorities: yes to running/deriving features now; no to “filling masterdb with more labels” as a blocking task. With 6 days left, labels are useful only if they already exist and can be added quickly; otherwise they’ll eat your time and still be mismatched to the tournament assay.

- question?: turning masterdb into a ranked set, running charlie's code + tools/mutcomputeLupssm/mutationflagging on the masterdb to see if the ranked set makes sense compared to our ground truth? 
- need to support current ipynb progress via converting the LuPSSM + mutcompute + run tools + PETase specific mutation flagging into --> scoring system ran on the test set. then we integrate it. but how 
- does the flag of disulfidge bridge need 3d modeling 
- NDCG ranking means the ranking is what matters so we need to filter/triage the dataset alot 
- need to focus on one end-end pipeline that outputs an initial ranking.tsv after which we can 
- now prioritize how we can model all types of mutations such as those found in the test set instead of focusing on mutations covered by literature. Will implement the literature mutations later. 
- the mutations covered by the theoretical efficiency scoring system by alam et al. doesnt cover all the types of mutations we will find in the test set meaning we need that to be a feature with a weak weight unless one of the alam's mutation is found in the test set variant 
- Mapped the mutations from Alam ST1 to the wt set, now need to do it for test set depending on the test set variant's mutation, lets say eg.V26A, need to annotate if that falls in one of Alam's mutation. Specifically the extract_alignment_position.py script checks for an input MSA and input alam mutation list, and outputs the aligned position at those alam positions for all the sequences in the MSA, in this case being the wt set. Need to continue this to the test set. 
- need to do lit review it would have helped with focus from the beggining 
* Alam et al. theoretical efficiency scoring system from ST1 
* Running thermoprot, mutcompute, GRAPE, EVcouplings
* WT set annotation/alignment to literature 
* Go through Charlie's code
* Go through Sanju's docking 

# **jan9** 
- Run alam et al. Table S1 theoretical scoring system on the wt and test set and masterdb, masterdb is what we can use to make sure our methods make sense to interpret results from the  test ste 
- Want to also integrate in the project FAIR Data principles on top of open source principles 
- Now continuing off the wt set annotation, meaning characterizing as deeply as possible the 3 WT backbones compared to the literature, which means to IsPETase, CaPETase, and perhaps others that are well-characterized 
	- Need to set the official closest most-characterized PETase closest to wt_2 and wt_3
		wt2: 
		wt3: 
	- Need to list the mutation codes publicly studies from CaPETase IsPETase so on and convert them to mutation codes that would be for the 3WT 
	- wt annotation produce the graphs/results/outputs 


# **jan8** 
- break 
# **jan 7**
- Minutes of the  meeting: 
	- saw difference between the RSMD LB/UB between the top 2 conformations tehy are very different, do we tune the FE, UB 
	- (6) energy search is optimistic so we might remove (LESS LONG)
	- confused on (4) unbound systems energy 
	- she set the coordinates based on catalytic triad --> need to confirm cleft all that 
	- try different exhaustiveness levels, need to do it for top candidates of variants 
- what are other papers/attempts to predict enzyme activity/exprssion/stability change of large scale variants (>1000/wt) and in the 2023 align tournament with amylase: 

	• Team TUM Rostlab utilized ProtT5, an encoder-decoder PLM. Instead of using frozen embeddings, they fine-tuned the model using Low-Rank Adaptation (LoRA) on the provided experimental data (activity, expression, stability) to create specific regression predictors for each property. For the generation phase, they used these fine-tuned "oracles" to score variants suggested by Efficient Evolution (using ESM models) and EvoPlay (a reinforcement learning agent).
	• Team AI4PD employed ZymCTRL, a conditional language model. They fine-tuned ZymCTRL on sequences with high expression and stability values to generate new variants. They found that the model's perplexity (a measure of how "natural" or expected a sequence is to the model) correlated with high activity and stability, using it as a primary filter. They further ranked sequences using ESM-1v embeddings and ProteinMPNN log-likelihoods.
	• Team Medium Bio benchmarked various embeddings, including ESM-1v, ESM-1b, and Georgiev Embeddings, training machine learning models (like Random Forests) on these features to predict properties. They selected the model with the lowest mean squared error to guide their sequence selection.

2. Structure-Based & Physics-Based Modeling
Teams in this category modeled the 3D structure of the enzyme and its interaction with the substrate to predict performance, often focusing on the transition state of the reaction.
• Team Nimbus (the overall winner of the in silico round) used PyRosetta to model the pre-transition state of the enzyme-substrate complex. They predicted structures using ColabFold (AlphaFold2) or ESMFold and docked the substrate using DiffDock.
    ◦ Strategy: They extracted various energetic and geometric metrics from these models (e.g., interface energy, packing).
    ◦ Prediction: They trained linear regression and spline interpolation models to map these structural metrics to the experimental labels (activity, expression, stability) provided in the training set.
    ◦ Design: For the in vitro round, they used RFDiffusion and ProteinMPNN to design de novo structured insertions to interact with the substrate.
• Team SergiR1996 used a Genetic Algorithm for multi-objective optimization. They employed Rosetta to calculate substrate binding energy and total system energy as fitness functions.
    ◦ Activity: They additionally used DLKcat, a deep-learning tool designed to predict enzyme turnover numbers (k cat), to estimate activity against specific substrates.
    ◦ Expression: They used a consensus of SoluProt and MPEPE (machine learning predictors) to classify sequences as expressible or not.
3. Evolutionary & Sequence-Based Approaches
These strategies relied on the statistical patterns found in natural sequence alignments to infer functional constraints.
• Team Marks Lab used EVcouplings, a model that captures evolutionary couplings (co-variation) between residues in natural sequence alignments.
    ◦ Strategy: They scored variants based on their probability under the EVcouplings model (assuming natural-like sequences are functional).
    ◦ Bias: To improve performance, they "biased" the model using the provided experimental data—masking mutations known to be deleterious and boosting those known to be beneficial in the assay conditions.
• Team Medium Bio also employed a Greedy Recombination strategy. They programmatically identified mutations that improved all three properties (activity, stability, expression) based on the training data and exhaustively recombined them. They also used Markov Chain Monte Carlo (MCMC) sampling of the sequence space using their trained Random Forest evaluators.

- concept very important to understand and carry: what do motifs mean biologically is not that simple (eg. M4 is not: “there exists a [PG]G[YF] somewhere upstream and M3 somewhere downstream” M4 is: “the specific upstream oxyanion/clamp motif [PG]G[YF] that corresponds to the Tyr87 region occurs at the correct relative distance from the catalytic serine that defines the same catalytic architecture”) That implies positional coupling, not independent existence.
- need to annotate the phylogeneti time/divergence of each little motif/blocks (aa and dna level) to have the evolutionary model just how Nelson did for CYP450s 
- other useful tool is to improve ESM3 and PLMs at the attention and layer mechanism level
- one useful tool to build as part of this project is finding a way to scale MDsim on 5000 variants easily cuz cost explodes and we wish we could do it for more accurate results 
- (1) finish the MSA zoom in of the motif table in the slides, then (2) compare with papers/suppinfo (3) meet 
- easy workflow to quickly find the PEtase region/motifs in the 3wt backbone by doing an alignment against wtIsPETase on ncbi msa viewer you can just hover above and itll tell you the original position of that residue for each different query (3 wt)
- note that wt_6 iw wt_IsPETase but has LED added at the end. And the signal peptide trimmed is actually from 1-27 + another 18 trimmed 
- blitzing the Is vs 3wt alignment for mutation flag annotation (check supp too)'

- talk bit about the comps model 
# **jan 6** 
- 
- UPDATED PLAN FOR THE DAY
	- 
	- Han flags 
	1. WT and Test set mutation type flagging/annotation 
	2. Thermoprot (lu) -> SANJU
	3. GRAPE (lu) -> SANJU 
	3. Mutcompute -> SANJU 
	4. Mut score (alam) -> JUSTIN 
	5. mutPSSM (Lu) -> JUSTIN 
	6. evocouplings -> CHARLIE 
	7. Soluprot -> JUSTIN CHARLIE 
	7. rosetta/pyrosetta/fastrelax/foldx DDG method -> SANJU   
	7. STABILITY: esm-1v, ddgemb, rosettaddg, deepddg, mutcompute, rnafold,thermoprot, prostab, temstapro (type2) temberture (type2)
	8. SOLUBILITY: procesa/netsolp, protsol (ecoli), progsol/gatsol (type2), aggrescan3D, VECTOR ANNOTATION 

- try the flag --> score method from chatgpt 
- continue off with annotating the test set
- need to review the docking stuff from sanju 
- once we have reference numbers for the 3wt, we can compare the wt to the test set by converting the mutation# to the code below that outputs/distinguishes the M class 
- how can we distinguish the M class between the 3 wt? We can compare them structurally to the 5 classes by COMPARING:
	* MSA alignment key regions (alam) (chen) (cui)
	* central b-sheet region (bell)
	* electrostatics of binding cleft (bell) (erickson) (hong) (ma)
	* ensemble refinements (bell) 
	* surface patch docking icmpro (bell)
	* global superposition with zoom-ins (bell) (bugin) (chen) (erickson)
	* pH effect on structure (norton)

	* equilibriated michaelis complex (bugin)
	* acylation and deacylation free energy profile (bugin) (garcia )
	* Trp185 dihedral free energy profile (bugin) (Cui) (hong)
	* Trp185-Tyr87 Aromatic interaction (transition state snapshots)  (bugin)

	* RMSF of all atoms from 2 structures highlighting BhrPETase's H218/F222, F243 W104 (Cui) (garcia)
	* Time evolution dihedral angles between  CA-CB-CG-CG1 atoms of W190 
	* BhrPETase's Stacking interactions between W104-P258 (Cui)
	* Distances between S165 OG and substrate's labile carbonyl atom 20-100ns (Cui) 
	* Comprehensive one (Erickson)

- need to establish strong/detailed reference of the 3wt in 3D modeling compared to IsPETase to determine a reference activity for the 3. what is the model learning/referencing/attending to, can attend not only to its own set of mutations from Ca but can try to backtrack from other types of mutations if they converge structurally 

- Phylogeny
* M1: cutinase-like, low PET activity 
* M2: transitional class
* M3: PETase-like 
* M4: enhanced PETase
* M5: Hyperactive PETase

- if we can get an initial scoring using the mutationtype below, the next priorities will be: (1) run thermoprot (2) run mutcompute (lu) (3) calculate mutation score (alam) (4) calculate mutPSSM (Lu) (5) other stab/sol tools 

- on the question of is it regression or classification, thats a very good point. It's a triage workflow we need to do so it will be clustering first (classification) with regression on very close/good scoring variants based on a learned score on anything labeled we have  (masterdb)

- mutation types that single point mutations influence:

A. Cys / Disulfide Liability

A single mutation can introduce or remove a cysteine, creating risk of mis-paired disulfides, oxidation sensitivity, or disrupting an existing disulfide if that cysteine is part of one.

B. Active Site & Substrate Binding

A single mutation can change side-chain chemistry, size, or aromaticity of a residue that directly contacts PET or lines the binding cleft, altering binding affinity or orientation.

C. Surface Charge / Solubility Remodeling

A single mutation can add or remove a charged or hydrophobic residue on the surface, changing solubility, aggregation propensity, or long-range electrostatic attraction/repulsion.

D. pH Robustness

A single mutation can alter the local environment of ionizable residues (especially His/Asp/Glu), shifting their effective pKa and stabilizing or destabilizing catalysis across pH.

E. Loop Rigidification / Flexibility

A single mutation (especially Gly↔Pro or bulky↔small) can change backbone entropy or steric freedom of a loop, affecting dynamics and access to the active site.


F. Aromatic Clamp / W-loop “Wobbling”

A single mutation can add, remove, or reposition an aromatic side chain (e.g., Phe/Tyr/Trp), altering π-stacking and gating motions that control substrate positioning.

G. Local Salt-Bridge / H-bond Rewiring

A single mutation can create or break one salt bridge or hydrogen bond, locally stabilizing or destabilizing secondary structure or side-chain orientation.

H. Hydrophobic Core Packing

A single mutation can fill a cavity or create a void in the hydrophobic core, changing packing efficiency, stability, and thermal tolerance.

I. Oxyanion Hole Geometry Change

A single mutation can shift side-chain orientation or backbone positioning of residues contributing to oxyanion stabilization, altering transition-state stabilization.

J. M5 + Family Motif Integrity

A single mutation within a conserved motif can break a conserved interaction or geometry, even if chemically “small,” reducing catalytic efficiency or folding reliability.

K. Electrostatic Field Tuning

A single mutation can change the local or long-range electrostatic potential by introducing/removing a charged or polar side chain, thereby altering:
	•	orientation of dipoles,
	•	stabilization of charged intermediates,
	•	guidance of substrate/product along the active-site path.

This is not network rewiring, just field reshaping.

L. Metal / Ca Coordination Weakening

A single mutation can remove or alter a coordinating residue (Asp/Glu/Asn), weakening metal binding affinity or geometry without fully eliminating the site.

- methods thinking vs project/hollistic thinking 
- once finished with the test set, need to execute code on the masterdb as input, and give it to charlie so he sees if the features we give him correlate well with activity/stab/expression from public studies 
- first of all what are the types of PETases evolutionarily: 
	* 
- tackling the test set now and focusing on coding the docking protocol




# **jan 5** 
- must change approach to focus ONLY on test dataset 
- thermoprot, mutcompute, 
- almost done aggregating all docking methods from the papers
- add the vina script to scripting.nd
- Minutes:
	* meet end of week, we'll talk about bit about GNNs 
	* need to know what features   
	* 
- ideas developed throughout: 
	* tool to analyze/extract specific info from 30+ papers at scale
	* protein structure visualization and attention maps according to different labels/learnings 
	* a public PETase specific database of all enzyme labels studies with our way of integrating/denoising data from different experiment conditions 
- Meeting at 1pm, going over questions for the zero shot deadline 
# **jan 4**
- taking too long to fetch activity with docking info, just gonna fetch docking methods and graphs cuz docking is most critical at the moment with 12 days left.  
- lots of activity measurements with different conditions we dont need that level of details its hard to implement specific conditions in the model like that.
- wt_5 is LCC_wt 
- wt_1 is a standard type IIb petase like IsPETase (Tm 48), but CaPETase is not as fragile (Tm 67) because it relies on its rigid core, while Is has a flexible loop that melts easily. I really need to annotate each masterdb and tournamentdb PETases for each different regions. 
- found that tournament_wt_2 has an extra CBM35 domain fused to it naturally through evolution.
- 4pm update, working on 2. now 
	1. Going through papers (labels, methods)
	2. Generate remaining structure for wt, test and masterdb (foldx)
- going through each paper now, extracting docking methods, figures/tables for labels, writing notes on how to integrate different experiments, etc.  At bell rn. 
- will need to specify from each activity study what kind of PET substrate was used, because different (amorph, cyrst,bottle, etc.) have different properties as %crystallinity, Mw, Mn, Dispersity (Bell et al.)
- summarizing methods and figures from all papers to better guide docking 
- 1pm update
	1. Run all tools against masterdb and tournament_wt/test 
	2. Fetch remaining structure/seq/labels for mastedb and tournament_wt/test 

- Currently
	- Annotate mutant set by type of mutation (structurally/biologically beyond just mutation code)
	- OBTAIN REMAINING LABELS MASTERDB: go through 27 papers
	- OBTAIN REMAINING STRUCTURES: (ncbi ID) --> (taxa) + closest hit's pdb + fetch foldx/alphafold
- all 3 terrestrial bacteria no marine 
- wt_1 CaPETase
- wt_2 is WP_162908185: 
- wt_3 is WP_374935857: 
- figuring out the 3 + 313 wt backbone information, then figure out docking stuff for meeting with sanju today 
- Will need help coding features 	
# **jan 3** 
- working with gemini to figure out automating the docking for all 5000 test sets and what that looks like in code and in autodock and in pymol is it first pose then docking ? 
- downloading papers + suppinfo 
- Now working on 1. fetching remaining data for masterdb, aggregating all papers and studies possible. found 27 additional papers i did not look at, will need to quickly fetch any data and biological insights on engineering/mutations, but cant waste time, need to do gemini summary of biological insights, and quickly screening graphs/tables/suppinfo for labels, and will also need a note on how to consolidate data from different studies.
- Today will be the docking worksession w/ Sanju 
- Need to give charlie the masterdb with the columns filled up (1) seq wo sigpep (2) activity/expression/stab/ph labels (3) structures 
# **jan 2**
- set up team meetings and catching up 
- Activity Score = (w1 * ESM-1v_LLR) + (w2 * PETase_Specific_Mutation_Matrix_Score) + (w3 * KL_Divergence_PETase_vs_Cutinase) + (w4 * Docking_Vina_Affinity) + (w5 * MD_RMSF_W185_Gate) - (w6 * Active_Site_pKa_Penalty_at_pH) + (w7 * Cleft_Net_Charge_at_pH) + (w8 * Cleft_Volume_Delta) + (w9 * Global_Consensus_ddG) - (w10 * Coevolutionary_Disruption_Index) - (w11 * Product_Release_H_Bond_Penalty) + (w12 * Binding_Cleft_Hydrophobicity_Shift) + (w13 * Cleft_Aromatic_Anchor_Surface_Area)
- Expression Score = (wA * mRNA_Folding_DeltaG_5prime) + (wB * Codon_Adaptation_Index_EC) - (wC * Ribosomal_Stall_and_PolyProline_Motifs) + (wD * Global_Consensus_ddG) - (wE * Spatial_Aggregation_Propensity_SAP) + (wF * NetSolP_Consensus_Score) - (wG * Total_SASA_Hydrophobic_Residues) - (wH * Protease_Motif_Count) + (wI * N_Terminal_Charge_Density) + (wJ * Tag_Exposure_IUPred_Score) + (wK * Surface_Accessibility_Near_Tag) - (wL * Metal_Binding_Internal_Patch_Count) - (wM * Cysteine_Count_Mismatch_Penalty) - (wN * Predicted_Metabolic_Burden) + (wP * Dye_Assay_Composition_Bias)
- Caution and notes on manual equation versus learning the equation: 
	1. An MLP can  "figure out" the weights, but it needs labels. 
	Pros: It excels at capturing non-linear relationships.1 For example, it can learn the "Stability-Activity Trade-off" (where activity increases with stability up to a point, then crashes because the enzyme becomes too rigid). A linear equation struggles with this "sweet spot."
	Cons: It is a "black box." If your MasterDB is biased (e.g., mostly IsPETase data), the MLP might overfit to those specific patterns and fail on the more diverse WP-backbone variants in the tournament set.Requirement: You need at least ~100–500 high-quality labeled points in your MasterDB for a simple MLP to outperform a well-tuned manual equation.
	2. The "Heuristic" Approach has Total Interpretability to ensure that the "H237 pKa Penalty" is the dominant factor for pH 5.5 without worrying about the model getting distracted by "noise" features, but It assumes the relationship between features is linear and won't easily find the complex cross-talk between, say, "mRNA Stability" and "Metabolic Burden" unless you explicitly code a cross-term.

- Undersatnd the workflow for fast building of bioml tools, need to rapidly summarize aggregate all data and literature to draft an initial task score/equation/algorithm/rule and then rapidly test it with pipelines already made that visualize the progress in ML performance and undersatnding with interp and PCA t-sne and attention maps in structures and sequences alot of this automated (analytics side) for example i was supposed to bring in the most detailed/niche activity and expression algorithms from previous work. 

- Plan for today: 
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


- figuring out annoying xcode cant clone github organization's repo unless i fix the permissions.  Organizations themselves do not issue personal access tokens so need to use personal justinkashi one. forget it im not doing it with an org im gonna do a personal repo. 
- Want to meet with the team to speak about progress, fetching data (pH, temp, activity, expression), finalizing the features 
- interrupted IQTREE on tournament_test, taking too long (60h) but i got the most recent tree and checkpoint saved, anyway theyre all single point mutations so dont care about the perfect tree 
# **jan 1** 
- taking break will need to finalize the features + fetch labels data + ph dependance + temp dependance 
- Below is missing: (1) consensus score from stability/solubility prediction tools (2) mutation score suppinfo (3) coevolutionary disruptions (4) MSA and phylogeny feature (5) what kind of docking score exactly do we need (6) all expression features (7) pH dependance features 
- Tentative final feature list which minimally uses outputs from ML tools:
* Thermostability & Rigidity Balance
	* thermoprot/prostab type 1 of tools: These are excellent for rapid screening of single-point mutations. They are typically faster than Rosetta and provide
	* temstapro / temberture type 2 of tools: specifically designed to predict the Growth Temperature (OGT) or thermal niche of the protein to determine if a variant is "naturally" designed for the pH 9.0/Higher Temp environment. 
	* METL biophysical type: slower but much better at identifying the Activity-Stability Trade-off because it understands the energetic cost of moving loops.
	* Core_vs_Loop_Stability: Calculate ddg specifically for residues in the "Second Shell" (e.g., D186) versus the "Catalytic Loops."Logic: High Core Stability + Low Loop Stability = High Activity (The "HotPETase" signature).
	* W185_Rotamer_Freedom: Measure the distance between residues at 214 and 218. If these are mutated to bulky groups (His/Phe), the gate "locks."Output: A binary "Gate-Locked" penalty (0 or 1).
	
* Activity-Specific (pH & Evolution)
	* Binding_Face_Zeta_Potential: Sum the charges of all surface-exposed residues within 10Å of the cleft at pH 9.0.Output: Net charge (positive values = higher activity at pH 9.0).
	* Kullback-Leibler (KL) Divergence: Use your IQ-TREE MSA to see if a mutation is "Conserved in PETases" but "Rare in Cutinases."Output: A "PETase-Specialization Score."H237_pKa_Delta: The difference between the WT pKa and Mutant pKa (from PROPKA)
* Expression-Specific (The "Folding Burden")
	* TYPE 1: 
	NetSolP trained to predict whether a protein will be in the soluble fraction or inclusion bodies, deep learning on sequence features, it captures the "folding burden" that leads to aggregation.
	ProtSol specifically targets the E. coli expression system for identifying how surface-exposed hydrophobic residues act as "nucleation seeds" for aggregation.
	* TYPE 2:
	Prog-Sol at the "Progressive" or structural side of solubility for understanding if a mutation affects the stability of intermediate folding states.

	GatSol Graph Attention Networks (GAT) with structures for your test set, because it looks at the spatial arrangement of residues rather than just the sequence string to see if a mutation creates a hydrophobic patch on the surface even if the residues are far apart in the sequence.

	* custom Aggrescan3D_Solubility: Unlike a sequence-string search, this uses your AlphaFold structure to find spatial hydrophobic patches.Output: A "Solubility Index" (Higher = better expression)
	* N-Terminal_Charge_Density: The net charge of the first 15 residues, Highly charged N-termini often improve translocation and prevent early aggregation.

- Evolutionary Potential: A group of variants that share "Frustrated" catalytic loops but "Minimal Frustration" in the core is the "Elite Group"—they are stable enough to express but flexible enough to evolve high activity.
- "Highly Frustrated" Clusters: These are regions where the sequence is not optimized for stability. In PETases, these are often the active site loops.
- Need to define the way thermostability influences the equation in the equations predicting the 3 properties, higher thermostability = (1) +/- activity @ pH 5.5 (2) - or + activity @ pH 9 ? (3) + expression for sure
- How can we evaluate the activity potential and ceiling for a group, and how can we even define or categorize a group that has the same potential for evolution/optimization ? 
- Note: 
	* Q: How does thermostability affect activity and expression? 
	* A: To catalyze a reaction, an enzyme must be flexible enough to "breathe"—it needs to capture the substrate and release the product. If you make a PETase too stable (e.g., by adding too many salt bridges), you might lock the W185 Aromatic Gate in place. If the gate can't move, the PET polymer can't enter the active site. PETase is naturally a mesophilic enzyme (likes 30°C), if you give it the stability of a 70°C enzyme but test it at 30°C, it may be too rigid to function efficiently. 
	* Using thermostability as a feature: conditions are at 30˚C but PET 
	* In a zero-shot setting, a model can't "learn" that pH 9.0 is different unless you give it features that specifically change with pH.
	* Expression titer in E. coli BL21 (DE3) is often limited by how fast the ribosome moves versus how fast the protein folds.
	* Mutation Distance to WT: The total number of mutations. In zero-shot, activity usually decays exponentially as you move away from a functional WT.
- Need to finalize the remaining features Charlie is missing 
	* Thermostability features
		* TOO high stability = lower activity 
		* higher stability = longer kinetics, more proteins folded that are purified 
		
	* Activity-specific 
		* Consensus score of variant tools: esm-1v , ddgemb, rosettaddg, deepddg, mutcompute
		* ESM3 embeddings 
		* Phylogeny (IQTREE, MSA) relative to PETase specific fam
			* Relative Entropy (Kullback–Leibler Divergence) 
			* Co-evolutionary Disruptions 
		* Scoring according to key structural regions 
		* Mutation score & mutPSSM according to suppinfo  
		* Remaining activity labels of masterdb 
		* Docking scores with PET ligand 
		* MDsim
			* B-factor/RMSF of aromatic gate W158 
			* Cleft volume/occupancy delta, hydrophobic index shift 
		* Just obtain interpro & PDB for the wt set 
	* Expression-specific  
		* Screen spatial aggregation propensity (a single hydrophobic residue can act as nucleation seed leading to aggregation)
		* Screen aa for protease sites 
		* Screen CDS for ribosomal stall site 
		* Tournament_vector annotate and convert to features (like?)
			* 
		* RNAfold mRNA structure stability 
		* Total Solvent Accessible Surface Area of hydrophobic residues
	* Run thermostability tools, use Tm as a feature to help  
		* Active Site Net Charge: Calculate the net charge of the catalytic triad + the aromatic gate (W185) + residues within 5Å of the cleft at pH 5.5 vs 9.0 using PROPKA. 
		* Protonation State of H237: The catalytic Histidine's ability to act as a base is pH-dependent. Use the predicted $pKa$ of H237 as a feature.
		* Solvation Energy Change: How the solubility of the binding pocket changes at different pH levels.
	* 
		
- Seems that esm-1v is great for variant zero shot calling, but the model is not PETase specific it just thinks if a mutation is "unlikely" based on evolution that activity goes down. 

- Need to send signal peptide-trimmed masterdb dataset to charlie + fill up activity/expression labels 
- phylogenetic tree building has a certain scoring system, that in itself does it give an evolutionary ranking that is general, i can use the score of the test set as a feature, and if i use the score that is PETase specific you can infer better ? 
- switching to gemini, and waiting on IQTREE to finish running so i can finally push to git 
- IMPROVING CHARLIES MODEL
	- Separate PETase-specific functional priors from generic similarity
	- Penalize mutations in catalytic / cleft / disulfide regions
	- Reward PETase-specific motifs (Trp185, beta8-alpha6 loop, DS1)
	- Explicit clustering-aware normalization (3 WT backbones)
- Map the more variable regions in the seqlogo motifs found in the testset versus the distribution of regions of test set mutations --> make nice visual 
- I think my job is to translate all the biological insights into features to train 
- Can't finish alone too much coding work and testset analysis, best way forward is to delegate the tasks fast, which requires proper communication of this highly interdisciplinary field (biology - compsci, bioinfo-wetlab-softwareengineer-ML/AI scientist/engineer) 
- Need a realistic deadline and deliverable
	* BY JAN 8: FIRST TEST SET RANK PREDICTION  
	* BY JAN 26: MINIMAL ACTIVITY-EXPRESSION RANKING MODEL BUILT 
- how to build all this and future possible analysis question autonomously at scale 
- Now annotating tournament_wt vs PDB fastas (remember to turn the entire pipeline of annotating against PDB then generating alphafold of those missing)
- I need a tool/script that just gives me the hits of RCSB PDBs from an input fasta file, which would require a map PDB ID - PDB seq - input seq, because currently the reverse search is only one sequence at a time on the webapp, the forward serach from PDB ID to fasta seq does exist though. --> nvm the rcsb db fasta is not large at all 
- IQtree is almost done after what 2 days, now can do in parallel: ESM embedding inference on wt testset // activity data search // 
- Need to identify wt_1/2/3: wt_1 is CaPETase, wt_2 is WP_162908185.1, wt_3 is WP_374935857.1
- Because of nature of the test set blastp can map the wt and test easily. 4678 proteins are from 3 wt, 310 remaining variants are from the 310 remaining wt.   
- Will need to show that we can incentivise engineers to build and collaborate on one github which will have numbers to show as well
- Will need to later on go through the logs to see what step and tasks we can make an automated pipeline/tool to publish that will help the community 
- Scripting visualization of mutation hotspots of testset, using cd-hit and MSA here to map them 
- MEME apply PETase specific PSSM ?
	* To estimate a Markov background from sequences:  fasta-get-markov <sequences.faa> <bg.txt> and then run MEME with -bfile bg.txt
-  Tradeoffs to consider about the motif window search in MEME: Too wide → motifs smear together, overfit family composition. Too narrow → miss cooperative patterns that span multiple residues (e.g., helix–loop interfaces. 
- NONCODING
	* Fetching activity data 
	* Define pH dependance of PETase 
	* Define temperature dependance of PETase 
	* Slides
	* Write manuscript, make lit review  
- CODING
	* Annotation 
		• Interproscan, MSA graphs, MEME motifs,CDD, PET-ase specific regions, 
		• BLASTP
		• Phobius
		• IQTREE
		• TMHMM
		• SignalP 
		
		• ESM embeddings --> Attention maps, PCA clusters 
		• Structural superposition testset of key regions 


		• mutcompute, coevolution annotation 
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


- need motif structure file for itol tree

-  how to translate learned lessons across enzyme-specific studies into the model to learn (will need to implement things like signalP probability distribution as features for all enzymes) 
- Phobius is annoying stuck at line 408 problem, failed not working unless linux 

- SUMMARY 
	1. ACTIVITY ANNOTATIONS TO RUN tournamentdb 
		* interproscan, iqtree, TMHMM, phobius 
		* thermostability prediction tools 
		* activity mutation score, mutation matrix, annotate tournament_test against tournament_wt (eg. mutation code tournament_test = tournament_wt1_A234H/F111H etc.) 
		* mutcompute & other variant ML tools 
		* Structural features (docking, MD, PDB features) 

	2. MASTERDB
		* Fetch remaining activity studies 

	3. TOURNAMENTDB STATISTICS & ANNOTATION GRAPHS  
		* signalp graph, blastp, lendist, taxanomy, tree, mutcompute, mutation score, MSA graph with motifs shown, sequence logos of each region/motif. 

	4. EXPRESSION 
		- note that we dont want pure expression but rather things like soluble expression + correct folding (stays soluble), accessibility of His tag for binding (buried tags reduce capture), IMAC binding/elution efficiency (histidine-rich internal patches can compete; proteolysis removing tag lowers yield), dye assay sensitivity to contaminants (co-purified proteins, nucleic acids) despite IMAC/benzonas, so features that proxy “IMAC-capture likelihood”: N/C-terminus disorder (tags more exposed), predicted surface accessibility near tag, presence of metal-binding motifs, and predicted oligomerization/aggregation (reduces capture)

		* ANNOTATE tournament_wt & tournament_test with:
			* solubility prediction tools 
			* double check all vectors have the same sequence except for the cloned PETase sequence 
			* (ORF LEVEL)
				•	Length (aa), MW, pI; predicted extinction coefficient isn’t used (A660 is dye-based), but MW/length still correlates with folding burden and expression.
				•	Amino-acid composition features: % hydrophobics, aromatics, prolines, glycines, cysteines; low-complexity runs; net charge; predicted solubility proxies.
				•	Predicted folding stability / aggregation propensity: per-residue disorder (IUPred-like), aggregation hotspots (TANGO/AGGRESCAN-like), hydrophobic patch metrics, predicted oligomerization tendency.
				•	Predicted subcellular localization / secretion signals: signal peptide probability (SignalP/Phobius), transmembrane helices. Even “false” weak SP/TM calls often correlate with membrane targeting/mislocalization and low soluble yield.
				•	Cysteine/disulfide patterns: BL21 cytosol is reducing; extra cysteines can hurt soluble yield unless designed carefully.
				•	Protease susceptibility motifs (basic/hydrophobic cleavage-prone regions) as a proxy for degradation in lysate/overnight expression.
				•	Mutation features relative to reference PETase: count of substitutions, BLOSUM scores, “to/from” property changes (charge, polarity, size), proximity to catalytic residues (if known), and whether mutations cluster (local destabilization).
				•	Optional but strong: predicted thermostability / ΔΔG (FoldX/Rosetta) and predicted expression/solubility scores (Solubis/Protein-Sol–style). Even rough predictors help as features.
			* (DNA LEVEL)
				•	Codon usage indices for E. coli (CAI, tAI), rare codon counts, codon pair bias, GC%, mRNA folding near start codon (−30..+70 nt) and around RBS/start (translation initiation bottleneck), long homopolymers/repeats.
				•	Internal Shine–Dalgarno-like motifs, strong secondary structure in the 5’ region, and “ribosome stalling” motifs (polyproline, certain dipeptides) that affect elongation and folding.
				•	If every insert is in the same PET28a context, the RBS/promoter/UTRs are constant; the insert’s first ~50 codons matter disproportionately.



	
- While the annotations are running, will now work on fetching all possible activity studies of PETase --> priority is to train a minimal model ASAP to output activity predictions, urgent 
- Can do annotaions on usegalaxy why tree is running
- Need to make the Notion 
- Annotations done so far: (1) BLASTP (pfam,masterdb) (2) MSA (mafft, clustal) (3) Phylogenetic Tree (4) Clusters (5) interproscan (running) (6) Signalp 
- When building phylogenetic tree with IQTREE what does it compute:
	•	Dayhoff: very old, general protein matrix (PAM-based).
	•	mtREV: mitochondrial proteins (vertebrate mtDNA).
	•	mtART: mitochondrial arthropod proteins.
These define the relative exchange rates between amino acids.
Add-on modifiers (after the +):
	•	+F: empirical amino-acid frequencies estimated from your alignment (instead of fixed frequencies from the original matrix). Often improves fit.
	•	+I: proportion of invariant sites (a fraction of sites assumed to never change).
	•	+Rk: FreeRate model with k discrete rate categories (e.g. R5, R6, R7). This is IQ-TREE’s modern replacement for +G (gamma); more flexible, often better.

# **dec 31** 
- parsing the ncbip hit table awk -F'\t' '$0 !~ /^#/ && $3 == 100 && $5 == 0 {print $1 "\t" $2}' tournament_wt_ncbip.txt --> annoying there is (1) some hits 0 mismatches + 100% identity but not full alignment range, need to get % (col4/len(query)), (2) some queries have 2 hits () 
- Smarter not to sell the model, but to use the model to scale clients and functionalities rapidly 
- Need statistics on #mutations per MOTIF/REGION 
- Need to fill up masterdb with (1) more PETases sequences (pfam, ancestral, NCBI proteomes all kingdoms) (2) More activity and expression studies 
- Looks like masterdb is not comprehensive enough because the blastp of tournament_wt against masterdb shows alot of wt are from these WP proteins, major cluster is CaPETase, and theres a few LCC/Is/Thermo meaning that the set is mostly not IsPETase (which is to be expected we are trying to discover/optimize newer PETases) so it means our ML model would be biased towards Is, means maybe good generalization towards distant homologs, but nonetheless the model needs to stay robust when inferring on novel sequences (most data out there is on Is)
- Ran the MSA and tree for tournamentdb, waiting a while 
- Need to make the pipeline for (1) mafft-trimal-iqtree-weblogo (2) blastp, cd-hit, quick stats,  
- In tournament_test, that there are consistently 3 major clusters via cd-hit meaning that  the mutant dataset comes from mostly 3 species --> which need to be identified. In tournament_wt they are clean 
- Need nice visual of positional/residue alignment of tournament_test vs tournament_wt -> need to make the MSAs 
- The tournament_wt is clean - indeed mostly well-characterized ones (mostly 100% similarity to the PETase Interpro/PFAM dataset or IPR041127.fasta) 
- Need to refocus to getting results ASAP on the tournament_wt and tournament_test, at least final activities and expression  
- Need build some API 
- Need to at least build a library or a tool as a result of this project --> how to use a prior dataset (fragmented features, imbalance, batch effects, homologs, ancestral, etc.) to supplement a model to fine tune it to one enzyme family  
- Need to focus on snakemake/nextflow pipelines to prove scalability, but on what code/task --> customer-need, ML training, protein structure visualization, 
- Need to add a general intro to the current research and applications of AI 
- Need to add an interpretability section to the project 
- Industries: (1) medicine (2) energy/sustainability (3) new technology 
- Now need to watch compbio seminars and operational comps 
- Wrote an onboarding document 
- Finishing the annotation for the tournamentdb now 

# **dec 30** 
- Notes to keep in mind about 2. PETase-specific motifs: 

	* IsPETase/Ca/Bhr (reference PETase)
		•	Catalytic triad intact at ref S160–D206–H237
		•	Trp185 present and flexible (aromatic gate for PET π–π stacking)
		•	Extended β8–α6 loop (~238–260) shaping an open cleft
		•	Two disulfides present (near-active-site and distal; PETase-characteristic)
		•	Moderate thermostability; high activity at mild temperatures

	* LCC
		•	Closed or semi-closed cleft (shorter / absent β8–α6 loop)
		•	Often lacks PETase-specific disulfide architecture
		•	Catalytic triad present but PET binding geometry differs
		•	High thermostability, requires PET near Tg for high activity

	* Ancestral PETase / cutinase-like hydrolases
		•	Catalytic triad present
		•	Lack extended β8–α6 loop
		•	Often lack PETase-specific disulfide
		•	Cleft more closed; poorer PET specificity
		•	Lower or substrate-ambiguous PET activity

- Now working on main.ipynb for tournament db statistics, i want: 
	1. Statistics (len distribution, no unique, no clusters 100-90, % with len 280-320, presence of noncanonical/ambiguous/artifacts aa) 
	2. PETase-specific motifs 
		* Catalytic Triad (~S160–H237–D206)
		* Trp185 cleft gate: position and flexibility class W/F/Y (enable pi-pi stacking with PET's phenyl rings)
		* Cleft extended b8-a6 loop (238-260)  the opposing cleft gate around 245-255 --> shapes the cleft (of PETase, absent or small in cutinase)
			* open cleft = fewer bulky aromatics blocking entrance wh
		* Disulfides DS1 (unique to PETase) C233-C282 
		* Disulfides DS2 (ancestral, shared by cutinase)
		* Local context of key residues ^
			* Flexibility (Gly/Pro) (like H237 rotamer permissiveness), charge (D/E/K/R), H-bonds capacity (count donors H/K/R/N/Q/S/T acceptors D/E/N/Q penalize bulky hydrophobics), steric bulk 
		* Score according to PETase mutation table (suppinfo) 
		* Generate the mutation table for each tournament sequence (suppinfo)

	3. Annotation 
		* Phylogenetic tree evolutionary annotation  
		* TM Helix prediction 
		* Signal Peptide Prediction  
		* Taxanomy annotation 
		* BLASTP Annotation 
			* NCBI, Uniprot, masterdb, IPR041127 
		* Activity/stability/fitness tools
		* Mutcompute tool 
		
	* Generating Alphafold2 structures, esm2/esm3 embeddings 

	* Mutcompute Annotation 
- Looking at taxanomy distribution of IsPETase it is restricted heavily in bacteria --> hypothesis, evolution was rapid and adaptable in prokaryotes indeed but there is large potential for evolution in plants and marine bacteria --> can we simulate a richer evolutionary trajectory inspiring from plant evolution (eg. CYP450 or specialized metbaolism type of evolution) to generate a strong novel PETase ?? 
- Gathering the interpro pdb uniprot etc. annotations and sequence_db to blast against the tournamentdb (exhaustive search)
- **TOURNAMENTDB**: (1) Run analysis/statistics/annotation 
- **MASTERDB** : (1) find remaining signal peptides (2) remove/add signal peptide to foldx structures 
- To do Tournament Dataset Analysis, it's important to have (1) PETase biophysical property document with interpro/uniprot/pdb annotations explained
- Overinvesting in: 
	The MASTERDB expansion (more LCC variants, more studies, more historical data) is now diminishing returns for zero-shot. Keep it frozen as a reference, not an active axis of work.
	Deep docking/MD/QM/MM thinking is correct scientifically, but too slow and too fragile as a primary ranking engine for a 1-month hackathon unless it’s heavily down-scoped.
	Trying to reconcile exact activity numbers across papers (TPA vs BHET vs MHET vs conditions) will not pay off for zero-shot ranking. Relative ordering matters more. 
-  Q: Would the best method to start from analyzing/understanding the tournament dataset ranking it with all possible methods and then packaging that into a ML/AI model? Instead of starting from my curated database of public PETase Decided to shift/reset to (1) Tournament Dataset Analysis, understanding it, then this will guide the zero-shot phase. The end-goal here is simply to rank the tournament dataset based on 3 conditions no matter how. 
- Will need to shift focus to tournament database statistics soon 
- Can expand the database to include more studies for LCC variants, and more activity studies. 
- Can experiment different pdb templates for the P181A of IsPETase and others that have PDBs available for the mutant  
- Note: FoldX is not a structure prediction engine. It is a local side-chain + energy perturbation tool. For FASTPETase, TS_PETase, TM3, D1 which are multimutations --> pdb do not exist for them --> check alphafold structures --> 

- fetching pdb one-liner: 

for id in 6ILW 7OSB 6IJ6 6KY5 7SH6 7QVH 7YM9 5LUK 4EB0 6IJ5 7EOA; do
  curl -L -o ${id}.cif https://files.rcsb.org/download/${id}.cif
done
- around 4am now and I need to continue from where we left off on the masterdb 
# **dec 29** 
- It's now 6pm; need to finish, of masterdb, (1) structures (2) CDS (3) (4) interpro/uniprot (5) blast tournamentdb to masterdb (6) annotate tournamentdb  

- Finished the docking protocol by aggregating the methods from papers/PETASE folder 

- note on activity and finding enzyme families (hydrolase vs lipase vs ): TPA yield measures full PET mineralization, so success depends as much on MHET/BHET hydrolysis as on PET surface depolymerization. Because PETase --> BHET + MET while MHETase/esterase/BTA-hydrolase-like --> TPA + EG 
In other words try to:
Ensure PET attack is sufficient (baseline PETase activity)
	•	Strongly optimize MHETase-like activity or dual-function behavior
	•	Favor mutations that:
	•	reduce MHET binding without release
	•	increase catalytic turnover on monoesters
	•	improve thermostability so downstream steps proceed longer

- we can learn the rankings instead of Tm for each study meaning have ranks of data for each study. 

- trick on google sheets master_db joining signal peptides ="M"&A1

- some HIS tags are the front or end of certain pdb sequences of PETase mutants

- Found a weird error where the publication sequence has 2 mutation positions wrong; ./check_wt.sh "A214H/I168R/W159H/S188Q/R280A/A180I/G165A/Q119Y/L17F/T140D" 
"MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCS"
FAIL A214H → sequence[214]=S (expected WT=A)
OK  I168R → sequence[168]=I (WT matches)
OK  W159H → sequence[159]=W (WT matches)
OK  S188Q → sequence[188]=S (WT matches)
OK  R280A → sequence[280]=R (WT matches)
OK  A180I → sequence[180]=A (WT matches)
OK  G165A → sequence[165]=G (WT matches)
OK  Q119Y → sequence[119]=Q (WT matches)
FAIL L17F → sequence[17]=G (expected WT=L)
OK  T140D → sequence[140]=T (WT matches)

# **dec 28** 
- using pdb-tools to fetch fasta from a .pdbx/mmCIF file(removing the X padding): 

f=7osb.cif; name=${f%.cif}; { printf ">%s\n" "$name"; pdb_fromcif "$f" | pdb_selchain -A | pdb_tofasta | sed '/^[^>]/ s/X*$//' | grep -v '^>' | tr -d '\n' ; printf "\n"; } > "$name.fasta"

- finishing the masterdb columns: sequence with mutation. Adding an input mutation MUTS into an input SEQ, checking whether the mutation was succesfully added. 

mut_seq=$(./add_mut.sh "$MUTS" "$SEQ")
./check_mut.sh "$MUTS" "$mut_seq"

*note on behaiour of pdb-tools:

pdb_fromcif | pdb_selchain -A | pdb_tofasta:
	•	reads atomic coordinates only
	•	outputs only residues that have coordinates
	•	drops:
	•	signal peptide (not in structure)
	•	His-tag (not resolved)
	•	pads missing residues with X (which you removed)

- wrote add_mut.sh rm_subseq.sh mut_check.sh scripts for easy handling 

# **nov 24**  

* **Obtain structures of benchmark_db**
* **Obtain esm2/esm3 embeddings of benchmarkdb (precomputed or we run esm3)**

* **Get property tool results for PETase data**
* **fill activity ph labels** 
* **Finish the mutation table for each batch of petase** 
* **Check in on docking**

# **nov 13**  
Welcome Sanju ! :) 
Going through the google drive, tournament rules
Priorities for December 
Finish small MLP head for 3 properties on 
Finish docking all PET-PETase mutants, with ideally an automation step for when we receive the tournament dataset (which we will have to dock) 
Finalize the plan on the esm credits, finding precomputed embeddings and making a model on the big 150k proteins dataset 
NOTE: every Tm/activity/expression measurement is done at a specific pH ! 
We forgot to talk about fine-tuning ! 
Model learns the pH calibration curve according to the type of PETase (thermo/dura/etc.) liek in supp. fig 2 https://drive.google.com/file/d/1RkxzZKTijDyOmi3EHKebG7h3YG1Clnp_/view 
(Script) Mutation checker position 
(Csv) fi
(Csv) Mutation table integration 
(Csv) Column for aligned to tournament (diamond)
(Csv) Typically signal peptide where ? 
(Csv) AA sequence without peptide 
(Csv) AA sequence with with peptide 
(Csv) PDB of the tournament 
Feature engineering/enrichment (1) data imputation (2) normalization (3) PCA and t-SNE, feature selection/impact score 
Tool library that helps turn bi databases into modular elements for machine learning and visualization of stats 
Build tool library to quickly visualize protein structures 
(1) Integrating/aligning our master database embeddings and features into Align’s dataset 
PDB, length, signal peptides (eg. trace back the wt seqs to pdb using BLAST wtv) [JUSTIN] + I need to fill in the pH and activity + encode the vector plasmid map for expression features 
Verify or compare features from literature to predictive tools (e.g. Expasy, Justin’s folder of chemical tools) 
Feature normalization/Polish the masterdb so its ML-ready 
Min-max scaling 
Z score normalization 
Ensure consistent dtypes 
Log scaling 
Batch effects and reducing technical noise to max the biological signal 
Feature enrichment:
Master database labels (activity, stab, express) + PDB embeddings 
Expression: plasmid vector features (.gb), promoter sequences/SD, codon optimizations, protein solubility (global @ pH 5.5 and 9, etc…)
Activity 
Supervised phase labels @ pH 5.5 and 9
PETase-specific mutation table, substitution matrix (suppinfo) 
Docking Score affinity between PET and PETase 
51 features classical properties 
Amino acid properties wrt point mutation changes on features we care about = functional activity and stability of expression
Future model architectures/idaes:
Finetuning ESM3 for PETases
Check Justin’s ESM finetuning tutorial
Training benchmark (entire protein families) for expression in E. coli pre-specializing on PETases
Explore alternate embeddings (T5)
	•		•	Learn difference between alphafold/foldx/etc. And crystal structure to correct for batch effects in case 

Motivation, hypothesis 
Each graph global aspect is the enzyme. 
The GAN network part is to organize the graph of our dataset to organize it in a  graph, followed by RL 


# **nov 5** 

Charlie updates 
SQL database
Docking scores PET-PETase 
Aaisha updates on CNN/GNN
Reinforcing the model with graph and convolutional features (likely in the simulation features?) 
Justin
Talk about the remaining features we did not talk about (activity, mutation table, biophysical attributes) 
Still need to do: get function track labels for the PETase dataset 

Moving forward: 
Justin + Charlie: (1) wrap up docking scores and function labels, (2) work on initial MLP-head and fine-tuning 
Aaisha: Find graph features 
(To do) Embedding benchmark dataset, finding precomputed embeddings 
(To do) Embed our dataset with different PLMs (T5, BERT) 
(To do) Here we need to focus on activity not function pred, so look for mutant-activity datasets (like the kaggle tournament one) 
(To do) Write the readme github

# **OCT 27** 

Charlie: making SQL database of the “MASTER_DB” 150 PETase 
Foldx structure embeddings 
What part should be test set vs validation set vs training set 


What features to predict: 
Standardizing the database 
Function track: 
PETase specific function feature (docking)?
Objectives
(1) Linear regression + (2) Xgboost / bagging boosting 
(3) MLP/NN
(4) Fine-tuning 


Benchmark model vs PETase model --> P( | ) 

# **oct 20** 

Unifying the benchmark solubility and stability datasets 
Unifying it is basically done
Now for embedding them, we have 150k sequences/struct, too long, I asked ESM about our research for more credits
look for pretrained embeddings;  Need to look for embeddings from their training data they might have it for majority of our proteins 
Building the model
We can start with this until we hear back from ESM; start with 150 PETase dataset to predict 2 labels (solubility, stability)  --> read the 4-5 papers on ML for PETase 
Function track: 
Function key words (interpro, blast) of PETase db 
Docking/MD of 150 PET-PETase and embedding 
Presentation slides, lit review of manuscript 


Charlie embedding numpy files easier to train vs pickle 
Double check numpy vs pickle embeddings vs per residue vs per sequence vs mean embedding vs all mebeddinbg, Esm-3-open vs esmc-6b 
Make a standardized database with dataset description 

Important reminder that we need to include the pH in the labeled dataset because the same enzyme has different Tm under different pH (https://www.kaggle.com/competitions/novozymes-enzyme-stability-prediction/discussion/355209) 
Do t-SNE on current embeddings as result + simple data analysis of proteins  
Need help with building github later on 
Build model with 150 PETases 
Go through novozymes code
Go through esm-1v code 
Fine-tuning: go through nvidia, aws, and the petase.ipynb code
Embed novozyme dataset 
May need to embed locally if cant find precomputed embeddings / if API fetching esm2 embedding doesnt work 


Tasks 
Find public esm3 sequence and structure embeddings of proteins from our benchmark datasets (Aaisha) 
Finish function track annotations (interpro and blast)  (Justin) 
Work on docking of PET-PETase for function track (Zach?)
Build lightweight model with 150 PETase training data (Charlie)

tips logs only write 5-6 lines idk what the mental is but whatever it is for u make it 5-6 line long and understand deeply waht you wrote, dont do verbose lengthy texts in logs or u get lost, the point of logs is to not get lost as is the remaining resources or ai helpers or etc. 
also add actually good summaries from gemini/gpt like the single point mutation one i added on jan18