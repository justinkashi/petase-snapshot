# **dec 30** 
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




