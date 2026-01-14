# ￼BIOINFORMATICS SCRIPTING  
  
**PARSE CDHIT CLUSTERS BETWEEN WT VS MUTANTS**  
```


```
  
  
**SHOW SEQ ONLY PRESENT IN ONE OR IN BOTH FASTA FILES **  
```
pip install cdhit-reader

```
```
#or 
conda install -c bioconda -c conda-forge cdhit-reader

```
```
cdhit-compare data/input1.faa data/input2.faa  --id 0.99

```
```


```
  
  
**Tacca path**: desktop/tacca/mossbioinfo/transcriptomics/march25/deepec/venv/bin/activate   
  
**RUN LINUX TOOL VIA DOCKER**  
```
# unpack
tar -xzf phobius101_linux.tgz
cd phobius

# build minimal container
cat <<EOF > Dockerfile
FROM ubuntu:20.04
RUN apt-get update && apt-get install -y perl
COPY . /phobius
WORKDIR /phobius
ENV PHOBIUS_DIR=/phobius
EOF

docker build -t phobius .

# run on fasta
docker run --rm -v "$PWD:/data" phobius \
  ./phobius.pl /data/sequences.fasta > /data/phobius.out

##DEBUGGING Could not read provided fasta sequence at ./phobius.pl line 408.


```
```
--> The bundled decodeanhmm is a 32-bit Linux ELF hard-linked to

```
```
/lib/ld-linux.so.2.

```
```


```
  
Split fasta file into 2 or n files   
```
n=$(grep -c '^>' tournament_test.fasta); half=$(( (n+1)/2 )); awk -v half="$half" '/^>/{i++} {print > (i<=half ? "part1.fa" : "part2.fa")}' tournament_test.fasta 

## 

N=3; awk -v n="$N" '
/^>/{i++; f=int((i-1)*n/total)+1}
{print > ("part_" f ".fa")}
' total=$(grep -c '^>' tournament_wt.fasta) tournament_wt.fasta

```
  
  
  
**Fast visualization of an MSA**  
  
```
seqkit fx2tab aligned.fasta | less -S

jalview aligned.fasta


```
  
  
**Loop through cd-hit seq similarity **  
```
for c in {100..90}; do p=$(printf "%.2f" "$(echo "$c/100" | bc -l)"); cd-hit -i tournament_test.fasta -o tournament_test_cdhit${c}.fasta -c $p -n 5 -d 0; done


```
  
  
**Convert column of file into fasta file **  
```
name="myseq"; awk -F',' -v n="$name" 'NR>1 {print ">"n"_"(NR-1)"\n"$1}' input.csv > "$name.fasta"

#or if the fasta headers are known in a column

awk -F',' 'NR>1 {print ">"$1"\n"$6}' input.csv > output.fasta

```
  
**MAFFT MSA**  
  
```
mafft --auto --anysymbol {.fasta} > {nelson_cyp_aligned.fasta} 

```
  
  
**Trim post-msa **  
```
./../../bioapps/trimAl_MacOS_x86-64/trimal -in nelson_cyp_aligned.fasta -out nelson_cyp_trimAl.fasta -fasta -automated1    


```
  
  
**IQTREE**  
```
./../../bioapps/iqtree-2.3.6-macOS/bin/iqtree2 -s tctridadzp450_mafft_trim.fasta -B 1000


```
  
  
**weblogo**  
```
weblogo -f tc_cyp_hmmalign_trim.fasta -l 100 -u -o tc_cyp_hmmalign_weblogo.png -F png 


```
  
  
  
**CDHIT**  
```
cd-hit -i pooled_p450.fasta -o pooled_p450_99 -M 0 -T 0 -c 0.99 -n 5 -d 0

cd-hit -i pool_supertranscripts_round2.fasta -o pool_supertranscripts_round2_cdhit95 -M 0 -T 0 -c 0.95 -n 5 -d 0


```
  
  
**DIAMOND**   
```
diamond blastp -q tcsl1_zi1_unified_transdecoder_clean.fasta -d swissprot_db.dmnd -f 6 --more-sensitive -e 0.000001 -k 1 -o tcsl1_zi1_unified_transdecoder_swissprot.out

diamond blastp -q ../pooled_leaf_rhizome_analysis/tcsl1_zi1_unified_transdecoder_clean.fasta -d uniprot_db.dmnd -f 6 qseqid sseqid pident evalue length slen --more-sensitive -e 0.01 -k 2 -o tcsl1_zi1_uniprot.out

diamond makedb --in TCSZI2_trinity/TCSZI2_transdecoder.fasta  -d tcszi2_db 
diamond makedb --in TCSZO2_p450_cdhit.fasta  -d tcszo2_p450db


```
  
  
**PFAM**  
  
```
hmmsearch --tblout trinity19_hmm Pfam-A.hmm trinity19_withania_transdecoder.fasta
hmmsearch --tblout pooled_cdhit95_hmm_round2 Pfam-A.hmm pool_supertranscripts_round2_cdhit95.fasta

```
```
#better
hmmscan --domtblout tczl_pooled_pfam_hits.txt ../PfamScan/Pfam-A.hmm tcsl1_zi1_unified_transdecoder_clean.fasta


```
  
  
  
**PFAM output parsing different**  
```
grep "^p450" dalata_p450_pfam.txt > dalata_p450.txt
awk '{print $4}' dalata_p450.txt > dalata_p450_id.txt
sort dalata_p450_id.txt | uniq > dalata_p450_uniqueid.txt
./../../../../fetch_fasta.sh


```
  
**Parsing ncbi blastp hit table **  
```
awk -F'\t' '$0 !~ /^#/ && $3 == 100 && $5 == 0 {print $1 "\t" $2}' tournament_wt_ncbip.txt

######

awk -F'\t' '$0 !~ /^#/ && $3 == 100 && $5 == 0 {print $1 "\t" $2}' tournament_wt_ncbip.txt \
| awk '
{
  split($1,a,"_");
  id=a[3];
  hits[id] = (id in hits ? hits[id] "," $2 : $2);
}
END {
  for (i=11;i<=312;i++) {
    if (i in hits) print hits[i];
    else print "";
  }
}
' > ncbi_tournament_wt_column.txt


```
  
  
**Length distribution of fasta file **  
```
cat file.fa | awk '$0 ~ ">" {print c;c=0} $0 !~ ">" {c+=length($0);} END { print c; }' | sed '/^$/d' | sort | uniq -c | awk '{print $1,$2}'


```
  
  
**Average sequence length of fasta file (works even without strip, meaning there’s spaces or multiline fasta files) **  
```
awk '/^>/ {if (seqlen){sum+=seqlen; count++}; seqlen=0; next} {seqlen += length($0)} END {sum+=seqlen; count++; print sum/count}' test.fasta 


```
  
  
  
  
**Convert sto to fasta**  
[https://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php](https://sequenceconversion.bugaco.com/converter/biology/sequences/stockholm_to_fasta.php) or   
  
```
easel esl-reformat afa

```
  
  
**++Fetch pdb’s  RCSB API ++**  
```
for id in 6ILW 7OSB 6IJ6 6KY5 7SH6 7QVH 7YM9 5LUK 4EB0 6IJ5 7EOA; do
  curl -L -o ${id}.cif https://files.rcsb.org/download/${id}.cif
done

```
  
**Fetch fasta**  
```
protein_seq=$(grep -A 1 ">$gene_id" tcsl1_zi1_unified_transdecoder_clean.fasta | tail -1)

```
  
  
**Cat, touch, nano, grep, awk, module list, pip list, diff between 2 fasta files, **  
  
```
perl ../trinityrnaseq-v2.15.1/util/align_and_estimate_abundance.pl --transcripts ../pooled_rhizome_hits/tcsz_pooled_transdecoder_clean.fasta --est_method kallisto --trinity_mode --seqType fq --samples_file sample_expression.txt --output_dir pool_rhizome_abundance


```
  
     
**Create a folder in a certain path **  
mkdir -p   
  
**Add a column of values (file1 file2 etc.)**  
```
for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done


```
  
  
**Merge 2 columns tsv file **  
Paste a b (a to the left of b)  
cut -f1 inputs.tsv  
```
cut -f1 inputs.tsv
#Paste a b 


```
  
  
**Strip * end of fastas **  
```
sed 's/\*//g' input.fasta > output.fasta

```
  
**Number of newlines **  
  
```
tr -cd '\n' < data.txt | wc -c

```
  
  
**Cleanup multiline fasta of newlines:**  
```
awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' triterpene_actransf_pep.fasta > triterpene_actransf_pep_clean.fasta


```
  
  
**Fasta file writer **  
```
def newlinebreak_fasta -->  go on terminal and write awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' TCSZI1_transdecoder.fasta > stdout.txt


```
  
  
  
**Loop through id’s in a text file (gotten from cut -f1) to fetch pfam **  
  
```
cat test.txt| while read line;


```
  
**R packages**  
```
Press R, then install.packages(“ “)

```
  
  
**PFAM output PARSING **  
```
awk '{for (i=1; i<=NF; i++) if ($i ~ /TRINITY/) print $i}' tcszl_tps_pfam.txt > tcszl_tps_id.txt
awk '{for (i=NF-3; i<=NF; i++) printf "%s ", $i; print ""}' tcszl_tps_pfam.txt > tcszl_tps_description.txt


```
  
  
  
**Replace the IDs (first column) of TMM/TPM to another TPM/TMM tsv file**   
```
paste temp.tsv temp2.tsv > tcp450_kallistotpmclean.tsv
cut -f 2- tcp450_kallistotpm_reordered.tsv > temp2.tsv
awk 'NR==FNR {order[$1]=NR; next} $1 in order {print order[$1], $0}' sortedid.txt tcp450_kallistotpmcopy.txt | sort -k1,1n | cut -d' ' -f2- > tcp450_kallistotpm_reordered.tsv
cut -f1 tcp450_kallisto_TMM_reorderedcopy.tsv > temp.tsv


```
  
**Fetch fasta of a uniprot ID using a batch job with Rest api **  
```
IDS=$(paste -sd "," meltome_uniprot_ids.txt)

curl -X POST "https://rest.uniprot.org/idmapping/run" \
     -H "Content-Type: application/x-www-form-urlencoded" \
     -d "from=UniProtKB_AC-ID&to=UniProtKB&ids=$IDS" \
     -o job.json
JOB=$(jq -r '.jobId' job.json)
JOB="a2bWuH3MwR"

while true; do
    STATUS=$(curl -s "https://rest.uniprot.org/idmapping/status/$JOB" \
             | grep -o '"jobStatus":"[^"]*"' \
             | cut -d':' -f2 | tr -d '"')

    echo "Status: $STATUS"

    if [ "$STATUS" = "FINISHED" ]; then
        echo "Job complete — downloading FASTA..."
        curl -L "https://rest.uniprot.org/idmapping/stream/$JOB?format=fasta" \
            -o meltome_uniprot.fasta
        echo "Done."
        break
    fi

    if [ "$STATUS" = "FAILED" ]; then
        echo "Job failed."
        break
    fi

    sleep 60   # <— poll every 1 minute
done

```
  
***don’t forget to chmod +x new .sh files to make them executable **  
  
**Filter out rows who’s sum of columns2-7 (counts) less than 50 **  
```
wk -F'\t' '{sum=0; for (i=2; i<=7; i++) sum+=$i} sum >= 50' tc_genetmmmatrix_coexpression.tsv > tc_tmm_coexpression_filtered.tsv


```
  
**From uniprot id to fasta files **  
  
```
curl --request POST \
  --url "https://rest.uniprot.org/uniprotkb/stream?format=fasta" \
  --header "Content-Type: text/plain" \
  --data-binary @id_list.txt \
  -o all_uniprot_sequences.fasta

```
  
  
  
**REPLACE GENE ID OF FASTA WITH LIST’S**  
```
cat names.txt 
TR1|c0_g1_i1    scaf0432344_50037.734_wgs
TR6|c0_g1_i1    scaf0159424_10142.072_wgs
seqkit replace -p "(.+)" -r '{kv}|$1' -k names.txt seq.fa 


```
  
  
**Obtain dna from plantismash output convert to fasta**  
```
for i in $(seq 1 21); do
  # Adjust input filename format based on the value of i
  if [ "$i" -le 9 ]; then
    input_gbk="CM048818.1.cluster00${i}.gbk"
  else
    input_gbk="CM048818.1.cluster0${i}.gbk"
  fi

```
```


```
```
  tcbgc_dna="tcbgc${i}_dna.txt"
  tleontobgc_dna="tleontobgc${i}_dna.txt"

  awk '/ORIGIN/ {found=1} found' "$input_gbk" > temp1.txt
  sed 's/[0-9]//g' temp1.txt > temp2.txt
  tr -d " " < temp2.txt | tr -d "\n" | tr -d "ORIGIN" | tr -d '0-9' > "$tcbgc_dna"
  tr -d "/" < "$tcbgc_dna" > "$tleontobgc_dna"


```
  
  
###################  
  
```
awk '/ORIGIN/ {found=1} found' JBFCZA010000030.1.cluster026.gbk > temp1.txt; sed 's/[0-9]//g' temp1.txt; tr -d " " < temp1.txt > temp2.txt; tr -d "\n" < temp2.txt > temp3.txt; tr -d "ORIGIN" < temp3.txt > temp4.txt; tr -d '0-9' < temp4.txt > tcbgc26_dna.txt; tr -d "/" < tcbgc26_dna.txt > tleontobgc26_dna.txt; rm temp1.txt temp2.txt temp3.txt temp4.txt tcbgc26_dna.txt


```
  
  
**hmmer panther **  
download ascii database  
```
for hmm_file in $(find . -name "*.hmm"); do hmmpress $hmm_file done
find target/famlib/rel/PANTHER19.0_altVersion/ascii/PANTHER19.0/books/ -name "*.hmm" | xargs -n 10
00

```
  
  
**Annotation**  
- Use grep to filter pfam outputs for certain pfam   
  
**QC assembly**  
```
salloc --time=01:00:00 --cpus-per-task=4 --mem=8G

bowtie2-build [] []
pooled_leaf_rhizome_analysis/tcsl1_zi1_unified_trinity.fasta tcsl1_zi1_unified_trinity.fasta


```
  
  
**Expression abundance **  
```
abundance_estimates_to_matrix.pl --est_method kallisto \
    --gene_trans_map none \
    --out_prefix kallisto \
    --name_sample_by_basedir \
     TCSZI1_trinity/tcszi1_kallisto/abundance.tsv \
     TCSZI2_trinity/tcszi2_kallisto/abundance.tsv \
     TCSZO1_trinity/tcszo1_kallisto/abundance.tsv \
     TCSZO2_trinity/tcszo2_kallisto/abundance.tsv 

perl trinityrnaseq-v2.15.1/util/abundance_estimates_to_matrix.pl --est_method kallisto \
    --gene_trans_map none \
    --out_prefix tcszi2 \
    --name_sample_by_basedir \
     TCSZI2_trinity/tcszi2_kallisto/abundance.tsv 


```
  
  
**KALLISTO**   
```
kallisto quant -i tcszi1_kallisto_index -o tcszi1_kallisto_output_bootstrap -b 50 TCS-ZI1_2-3407463_S2_L001_R1_001.fastq TCS-ZI1_2-3407463_S2_L001_R2_001.fastq 


```
  
  
  
**TBLASTN**  
**BGC**  
```
for db in dzing_db.dmnd dalata_db.dmnd drotuntada_db.dmnd dcayensis_db.dmnd; do diamond blastp -k 1 -f 6 -e 1e-6 -d $db -q tleonto_bgc14.txt -o ${db}bgc14_tleonto_blastp.tsv
done

```
  
  
  
**HPC**  
- No conda no docker on hpc for security reasons …   
**RANDOM**   
- Cd-hit installation needs gcc openmp, but don’t forget to modify the makefile where it says CC=g++ to CC=g++-14 (the version on your system). For all the CC= lines.   
diskusage_report   
  
  
```
mkdir {gene_hits}

diamond blastp -d TCSZI1_transdecoder_db.dmnd -q qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

gene_id=$(awk -F'\t' '{print $2}' smt1_hits/tcszi1_hits_smt1.txt | head -n 1)


```
```
protein_seq = $(grep -A 1 ">$sequence_id" TCSZI1_transdecoder.fasta | tail -1)


```
  
  
**Average of TPM matrix first 3 columns and last 3 columns (after first sample name column) **  
```
awk '{ avg1=($2 + $3 + $4)/3; avg2=($5 + $6 + $7)/3; printf "%s\t%.3f\t%.3f\n", $1, avg1, avg2 }' pooled_leaf_rhizome_kallisto_tpm_phytosterolhits.txt


```
  
  
  
**EdgeR File sort by nth column of numbers **  
```
sort -nk4 tc_geneEdger.txt > sorted_edger.txt


```
  
  
******Grep file contents against another file*******  
```
while read -r pattern; 
	do grep -m 1 "$pattern" alignments/tcp450_algn/test2.txt 
done < sortedid.txt > search_results.txt


```
  
  
**Grep beginning of line**  
```
grep "^p450" filename

```
  
**PATH**  
```
PATH=$PATH:~/opt/bin


```
  
  
**Sort tsv file by value of fourth column**  
```
sort -k4,4 -n input.tsv > sorted_output.tsv


```
  
  
**UNIQUE TEXT FILE **  
  
```
sort gene_ids.txt | uniq > non_redundant_gene_ids.txt


```
  
  
  
Convert .txt to .tsv by replacing spaces with tab delimiter  
  
```
sed 's/ \+/\t/g' input.txt > output.tsv


```
  
  
  
**Reorder TMM.tsv file by order of a .txt file containing ordered list of gene id’s **  
```
awk 'NR==FNR {order[$1]=NR; next} $1 in order {print order[$1], $0}' sortedid.txt tcp450_kallisto_TMM.tsv | sort -k1,1n | cut -d' ' -f2- > reordered_TMM.tsv


```
  
  
**Extract gene ids from TMM or fasta file**  
```
cut -f1 kallisto.gene.TMM.EXPR.matrix.tsv > tmm_gene_ids.txt
grep ">" tcp450_cdhitclean.fasta | sed 's/>//' > tcp450_ids.txt

```
  
  
*****Trim the gene id of protein into just gene level *****  
```
awk -F'_i' '{print $1}' input > output 


```
  
**Convert pfam .txt file into .tsv to later use cut -f3 to get the gene IDs to then do pfam output parsing with fetch_fasta.sh script to get the fasta file**  
```
awk '{$1=$1; print}' OFS="\t" tl_p450_pfam.txt > tlp450pfam.tsv

```
  
  
  
**Sed tutorial guide  **  
sed -n -e 's/^.*stalled: //p'  
sed -n -e ’s/JBF //p’  
Detailed explanation:  
* -n means not to print anything by default.  
* -e is followed by a sed command.  
* s is the pattern replacement command.  
* The regular expression ^.*stalled:  matches the pattern you're looking for, plus any preceding text (.* meaning any text, with an initial ^ to say that the match begins at the beginning of the line). Note that if stalled:  occurs several times on the line, this will match the last occurrence.  
* The match, i.e. everything on the line up to stalled: , is replaced by the empty string (i.e. deleted).  
* The final p means to print the transformed line.  
  
**Remove nth columns of a .tsv file**  
cut -f1-4,6,8- tcp450_edger.txt > tcp450_edgerlogfc.txt  
  
****Filter a TMM gene matrix file for only the gene id’s of a fasta file:  ****  
```
awk 'NR==FNR {ids[$1]; next} $1 in ids' tcp450geneid.txt kallisto.gene.TMM.EXPR.matrix.tsv > filtered_kallisto_TMM.tsv


```
  
  
  
**Filter EdgeR File based on smaller than -2 as logFC threshold **  
```
awk -F'\t' -v margin=0.2 'FNR==2 {min=$4} FNR>1&&($4<=-2)' sorted_edger.txt > tcdgesorted.txt


```
  
**Hmmalign**  
```
hmmalign -o nelson_cyp_hmmalign.sto ../PF00067.hmm ../../../fasta_databases/cyp_db_full.fasta 


```
  
**Hmm fetch a pfam hmm model using a .txt containing your pfamID**  
```
hmmalign -o tc_cyp_hmmAl.sto PF00067.hmm tc_p450_cdhit.fasta


```
  
**Esl-reformat .sto to fasta (aligned) using git cloned easel **  
```
./miniapps/esl-reformat afa ../nelson_cyp_hmmalign.sto > ../nelson_cyp_hmmalign.fasta


```
  
**Filter non M-starting proteins fasta file**  
```
seqkit grep -r -s -p "^M" input -o output 


```
  
**Trim a fasta file by size/length **  
```
seqkit seq -m 380 -M 600 input.fasta -o filtered.fasta


```
  
**Trim a character from each line**  
```
sed 's/x//g' filename


```
  
**Curl vs wget**  
![Feature](Attachments/84D6FF2A-3671-403F-8514-5CE35A800ABF.png)  
  
Notes/things to study more officially   
- $PATH is just an environment variable containing directories/libraries you added. How do these variables differ from venv to venv?folder to folder, Jupyter kernel to another?   
- Rsync -avP, better than scp   
- **Using GPU, implementing it **  
  
  
—————   
Sep 22 2025 esm petase project  
  
  
**PARSING FOLDER OF RSCB PDB.GZ FILES TO FETCH TOTAL # SEQUENCES AND TOTAL LENGTH **  
```
!/bin/bash        
count=0
total=0
for f in *.pdb.gz; do
    len=$(zgrep "^SEQRES" "$f" | awk '{for(i=5;i<=NF;i++) c++} END{print c+0}')
    total=$((total + len))
    count=$((count + 1))
done
echo "Total proteins: $count"
echo "Sum of lengths: $total"
 

```
  
**OUTPUT TOTAL SEQUENCE LENGTH OF FASTA FILE ‘’ OF INTERPRO TSV FILE**  
```
 awk -F'\t' '       
NR>1 && $6 != "" {
  total += $6
  count++
}
END {
  print "Total sequences:", count
  print "Sum of lengths:", total
}' petase_db/IPR041127_sequence.tsv

awk '
  /^>/ { if (seq_len > 0) {
            total += seq_len; count++; 
            if (seq_len > 1024) long_count++; 
         }
         seq_len=0; next }
  { seq_len += length($0) }
  END {
    if (seq_len > 0) {
      total += seq_len; count++; 
      if (seq_len > 1024) long_count++;
    }
    print "Total sequences:", count;
    print "Sequences >1024:", long_count;
    print "Sum of lengths:", total;
  }
' input.fasta


```
  
  
**COMBINE FASTAS INTO ONE**  
```
for f in D1-PETase.fasta DuraN233K.fasta DuraPETase.fasta FAST-PETase.fasta HOT-PETase.fasta \
         IsPETase_S238FW159H.fasta IsPETase_W159H_F229Y.fasta IsPETasevariants_Lu.fasta \
         thermopetase.fasta TM3-PETase.fasta TS-PETase.fasta wtIsPETase_signaltrim.fasta \
         wtispetase.fasta; do
  cat "$f"
  echo
done > all_petase_variants.fasta


```
  
  
  
**FOLDX COMMANDS**  
```
. ../foldx5_1Mac/foldx_20251231 --command=RepairPDB --pdb=../petase_db/6ilw.pdb 

cat > mutations.txt <<EOF
A214H,A168R,A159W,A188Q,A280A,A180I,A165G,A119Y,A17F,A140D
EOF

foldx --command=BuildModel --pdb=6ilw_Repair.pdb --mutant-file=mutations.txt

```
  
**SPLIT MULTI FASTA INTO SINGLE FASTA FILES**  
```
awk '/^>/ {
    if (out) close(out);
    h = substr($0,2);             # remove ">"
    split(h, a, " ");             # split on space, take first word
    gsub(/[\/ ]/, "_", a[1]);     # replace "/" or space with "_"
    out = a[1] ".fasta";
}
{ print >> out }' LCC_variants_tournier.fasta 


```
  
[oct 17 2025]  
```
git init

```
```
git add . 
git commit -m "initial commit"

```
```
git branch -M main
git push -u origin main
git remote add origin https://github.com/justinkashi/petase.git


```
  
```

# 1. Remove any old container
docker rm -f fireprot

# 2. Start a new MySQL container with root password “root”
docker run -d --name fireprot -e MYSQL_ROOT_PASSWORD=root -v fireprot_data:/var/lib/mysql mysql:8.0

# 3. Wait 30–60s, then verify it’s ready
docker logs -f fireprot   # stop with Ctrl+C when you see “ready for connections”

# 4. Copy your SQL dump into the container
docker cp fireprotdb.sql fireprot:/fireprotdb.sql

# 5. Enter the MySQL shell
docker exec -it fireprot mysql -u root -p
# (enter password: root)

# 6. Inside MySQL, create and use the database, then import the dump
CREATE DATABASE fireprot;
USE fireprot;
SOURCE /fireprotdb.sql;
SHOW TABLES;


```

SOURCE JUSTIN KASHI NOTES 