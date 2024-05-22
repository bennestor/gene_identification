This script (find_genes_v2.sh) is used to identify homologous protein sequences from multiple gene families in multiple genomes using confirmation from BLAST statistics, HMMER, and checks to the genome sequence to ensure no incorrectly predicted proteins were missed.
The general pipeline for this search is shown below:

1. If using own query sequences, set up a directory named {QUERIES}/{TRANS}/ (e.g. blast_queries/nrt2) with the name '{REF}_{TRANS}_proteome_re.fas' (e.g. Arabidopsis_nrt2_proteome_re.fas). The ${QUERIES} varaiable, among other directory variables, is specified at the top of the script.
Otherwise, proteins can be extracted from the reference proteome based on annotations (regex will need to be specified inside the script under the first section '#Extract proteins from reference. Edit depending on annotations and proteins.').
2. Secondly, the extracted protein sequences are BLAST searched in each specified proteome using strict evalue (EVALUE1), coverage (COV1%), and length (LENGTH AA) cutoffs.
3. Thirdly, the full protein sequences from the initial BLAST are re-BLAST searched in the proteomes to retrieve more divergent sequences (similar to the PSI-BLAST tool). The filtering uses less strict evalue (EVALUE2) and coverage (COV2) cutoffs, but keeps the same lengths. This step can be skipped using the -R flag.
(Optional and CURRENTLY DISABLED - can be re-enabled by editing the procedure steps at the end of the script)
4. Fourthly, the retrieved proteins are searched for domains generated from the intial BLAST protein sequences using HMMsearch and filtered by evalue EVALUE3.                                                                                                         Domains from pfam can be used instead if placed as {HMM}/{TRANS}/{TRANS}.prof
5. Finally, the whole genome sequence is BLAST searched for the same protein family using EVALUE4 and any regions with significant hits that were not annotated are extracted within a specified flanking region (-f).
6. Later sections include automatic prediction of protein sequences using AUGUSTUS - these are currently commented out, but removing the hastags should make them active again.


Available flags:
-s <"species list"> #specify this as a quoted list separated by spaces (e.g. -s "Arabidopsis Gossypium Oryza")
-r <reference species> 
-q <"query gene family name"> #also specify as quoted list (e.g. -q "pht1 nrt2 npf")
-e <initial BLAST evalue>
-c <initial BLAST coverage> 
-E <re-BLAST evalue>
-C <re-BLAST coverage>
-l <min protein length> #anything below is filtered out
-v <HMM evalue>
-n <genome blast evalue>
-t <threads for each program>
-f <flanking range>
-i #Start from initial BLAST filtering to respecify evalue and coverage cutoff
-j #Start from re-BLAST filtering to respecify evalue and coverage cutoff
-m #Start from HMMsearch filtering to respecify evalue cutoff
-g #Start from genome BLAST search filtering

At the top of the script are lines for default variables to be set. Below that is a section for ordering the desired output file structure.
The script requires BLAST (blastp and blastn) and HMMER to be installed and available on $PATH

A phylogenetic tree making script is also available: make_trees_v2.sh, which requires RAXML.
SLURM scripts for running both scripts on HPCs are also available.
