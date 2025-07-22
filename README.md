# HMM-for-Kunitz-domains
## Prject Description
This project builds a Hidden Markov Model (HMM) to identify Kunitz-type protease inhibitor domains (Pfam: PF00014) in protein sequences. It integrates structural data from the Protein Data Bank (PDB) and curated sequences from UniProt to develop, train, and evaluate the model.

##### Key Features:
- Automated retrieval of Kunitz domain sequences from PDB and UniProt.
- Sequence clustering, cleaning, and alignment to construct a high-quality multiple sequence alignment (MSA).
- HMM profile generation using HMMER.
- Removal of high-similarity sequences to ensure dataset diversity.
- Construction of positive and negative datasets for performance evaluation.
- Scalable performance evaluation across multiple E-value thresholds.


###### This pipeline is useful for researchers aiming to detect Kunitz domains in novel sequences, benchmark domain prediction algorithms, or explore evolutionary patterns in protease inhibitors.


##  First Data Retrieval 
Download PDB entries matching the following criteria:
- Resolution ≤ 3.5 Å
- Sequence length between 45 and 80
- Pfam domain = PF00014
#### Output: `rcsb_pdb_custom_report_*.csv`
Note: the '*' in the output file stands for a specific number that you will get after downloading the file from PDB

## Convert CSV to FASTA
tr -d '"' < rcsb_pdb_custom_report_*.csv | \
awk -F ',' '{if (length($2)>0) {name=$2}; print name,$6,$8,$9}' | \
grep PF00014 | \
awk '{print ">"$1"_"$3; print $2}' > `pdb_kunitz.fasta`

## Clustering with CD-HID
- cd-hit -i pdb_kunitz.fasta -o pdb_kunitz_cluster.txt
#### output files:
- `pdb_kunitz_cluster.txt`
- `pdb_kunitz_cluster.txt.clstr`
- `pdb_kunitz_cluster.txt.clstr`

#### Converting pdb_kunitz_cluster.txt.clstr to tabular:
- clstr2txt.pl pdb_kunitz_cluster.txt.clstr > clusters_table.txt

## Multiple sequence alignment (MSA):
- grep '^>' pdb_kunitz_cluster.txt | sed 's/^>//' | sed 's/_/:/' > `pdb_kunitz_ids_25.txt`

- Manually remove sequences: 2ODY:E, 5JBT:Y, 4BQD:A
- Run MSA lignment again with now manually curated file  → *output: `pdbs_22_MSA.fasta`*

## Clean your MSA file:
Run in your terminal:
- python3 clean_fasta.py pdbs_22_MSA.fasta `pdbs_kunitz_MSA_clean.ali`

## Build a HMM profile
- hmmbuild    pdb_kunitz_HP_clean.hmm    pdbs_kunitz_MSA_clean.ali

output: `pdb_kunitz_HP_clean.hmm`

## Second Data retrieval:
Open Uniprot and download the results from the following queries:

•	`Human: (taxonomy_id:9606) AND (reviewed:true) AND (xref:pfam-PF00014)`

•	`Non-human: NOT (taxonomy_id:9606) AND (reviewed:true) AND (xref:pfam-PF00014)`

Download these two FASTA files from Uniprot

#### Merge the two Kunitz dataset files to form a unified collection of positive examples 

- cat uniprot_human_kunitz.fasta uniport_not_human_kunits.fasta > `kunitz_all.fasta`

## Remove high similarity Hits
#### Remove all sequences with a sequence identity ≥ 95% and a Nres ≥ 50 

- makeblastdb -in kunitz_all.fasta -input_type fasta -dbtype prot -out kunitz_all.fasta

- awk '/^>/ {print; next} {gsub("-", ""); print}' pdbs_22_cleaned.fasta > `pdbs_22_ungapped.fasta`

- blastp -query pdbs_22_ungapped.fasta -db kunitz_all.fasta -out kunitz_pdb22.blast -outfmt 7

- grep -v "^#" kunitz_pdb22.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > `high_match_22.txt`
- cut -d'|' -f2 high_match_22.txt > `to_remove.txt`

- grep ">" kunitz_all.fasta | cut -d"|" -f2 > kunitz_all.txt
- comm -23 <(sort kunitz_all.txt) <(sort to_remove.txt) > `kunitz_final.txt`

#### This is our positive evaluation set 

## Collect all SwissProt reviewed non-kunitz proteins (N= 573.832)
#### Download the Corresponding FASTA file:
 `all_swisprot.fasta`

#### Create a list of the IDs of the whole UniProtKB/SwissProt:
- grep ">" all_swiss.fasta | cut -d "|" -f2 > `all_swissprot.txt`
#### Remove the positive IDs from the list of the negatives:
- comm -23 <(sort all_swissprot.txt) <(sort kunitz_final.txt) > `negatives.txt`

## Subsetting Positive and Negative Datasets:
#### Random sorting IDs in the file:
- sort -R kunitz_final.txt > `kunitz_final_random.txt`
- sort -R negatives.txt > `negatives_random.txt`

#### Subsetting in half:
- head -n 183 kunitz_final_random.txt > `pos_1.txt`
- tail -n 182 kunitz_final_random.txt > `pos_2.txt`

- head -n 286416 negatives_random.txt > `neg_1.txt`
- tail -n 286416 negatives_random.txt > `neg_2.txt`

## Extract Sequences by ID:
#### downoad the `get_seq.py` file from attached to this repository and run the following commands
- python3 get_seq.py pos_1.txt all_swiss.fasta > `pos_1.fasta`
- python3 get_seq.py pos_2.txt all_swiss.fasta > `pos_2.fasta`
- python3 get_seq.py neg_1.txt all_swiss.fasta > `neg_1.fasta`
- python3 get_seq.py neg_2.txt all_swiss.fasta > `neg_2.fasta`
##### please make sure to change the file names according to you own files 
##### please make sure to have the file get_seq.py save din your working directory

## Run HMM search and Hit Mapping:
- hmmsearch -Z 1000 --max --tblout `pos_1.out` pdb_kunitz_HP_clean.hmm pos_1.fasta
- hmmsearch -Z 1000 --max --tblout `pos_2.out` pdb_kunitz_HP_clean.hmm pos_2.fasta

- hmmsearch -Z 1000 --max --tblout `neg_1.out` pdb_kunitz_HP_clean.hmm neg_1.fasta
- hmmsearch -Z 1000 --max --tblout `neg_2.out` pdb_kunitz_HP_clean.hmm neg_2.fasta

## Generate `.class` Files:
- grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > `pos_1.class`
- grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > `pos_2.class`

- grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > `neg_1.class`
- grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > `neg_2.class`

#### Add missing negatives (non-hits)
- comm -23 <(sort neg_1.txt) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class
- comm -23 <(sort neg_2.txt) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_2.class

## Combine Final Datasets
- cat pos_1.class neg_1.class > `set_1.class`
- cat pos_2.class neg_2.class > `set_2.class`

## Evaluate model Performance 
#### Download the file named `performance-2.py` and run the following commands in your terminal:

- for i in $(seq 1 12); do
    python3 performance-2.py set_1.class 1e-$i
done

- for i in $(seq 1 12); do
    python3 performance-2.py set_2.class 1e-$i
done


