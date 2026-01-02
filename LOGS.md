# **jan 2**
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



**other**
- looking at xcode can we use AR app, safari extension, sticker pack app or whatever else for something in bioinformatics/innovative ? 
- pivoting to instrumentation, dignaostics --> models for blood pressure/time-series big data, link meta, can we use visionOS 