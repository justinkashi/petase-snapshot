# **dec 29** 
- 

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


