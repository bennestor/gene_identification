#!/bin/bash
set -euxo pipefail

## DEFINE VARIABLES

#Set default variables. Use "" for lists
SPECIES_LIST="Amborella Eucalyptus Hakea Macadamia Medicago Nelumbo Oryza Protea Telopea Vitis Zea"
REF=Arabidopsis
TRANS_LIST=""
EVALUE1=1e-10
COV1=0
EVALUE2=1e-10
COV2=0
LENGTH=200
EVALUE3=1e-10
EVALUE4=1e-5
FLANK=900
THREADS=64

#Path variables
BLASTDB=$MYSCRATCH/2020_11_20_Transporters/1-data/blast_DB
PROTEOMES=$MYSCRATCH/2020_11_20_Transporters/1-data/proteomes
BLASTRES=$MYSCRATCH/2020_11_20_Transporters/2-blasts
TRANSFAS=$MYSCRATCH/2020_11_20_Transporters/3-trans_fas
HMM=${MYSCRATCH}/2020_11_20_Transporters/4-hmm
HMMFAS=${MYSCRATCH}/2020_11_20_Transporters/5-hmm_fas
FINALFAS=${MYSCRATCH}/2020_11_20_Transporters/6-final_fas
GENOMES=${MYSCRATCH}/2020_11_20_Transporters/1-data/genomes
QUERIES=${MYSCRATCH}/2020_11_20_Transporters/1-data/blast_queries
AUG=${MYSCRATCH}/2020_11_20_Transporters/0-scripts/Augustus/bin

#Steps
REBLAST=false
INITBLASTFILTER=false
REBLASTFILTER=false
HMMFILTER=false
GENOMEFILTER=false

#Usage command
print_usage () {
  echo -e "
Usage: $0 [-s <"species list">] [-r <reference species>] [-q <"query gene family name">] [-e <initial BLAST evalue>] \
[-c <initial BLAST coverage>] [-E <re-BLAST evalue>] [-C <re-BLAST coverage] [-l <min protein length>] [-v <HMM evalue>] \
[-n <genome blast evalue] [-t <threads for each program>] [-f <flanking range>] [-i] [-j] [-m] [-g] [-h] [-I]

This script identifies specified families of proteins in species using BLAST and HMMsearch tools. 

Firstly, proteins are extracted from the reference proteome based on annotations (regex will need to be specified inside the script). \
If using own query sequences, place in {QUERIES}/{TRANS}/ with the name '{REF}_{TRANS}_proteome_re.fas'.
Secondly, the extracted protein sequences are BLAST searched in each specified proteome using strict evalue (EVALUE1), coverage (COV1%),\
and length (LENGTH AA) cutoffs.
Thirdly, the full protein sequences from the initial BLAST are re-BLAST searched in the proteomes to retrieve more divergent \
sequences (similar to the PSI-BLAST tool). The filtering uses less strict evalue (EVALUE2) and coverage (COV2) cutoffs, but keeps the same lengths. This step can be skipped using the -R flag.
Fourthly, the retrieved proteins are searched for domains generated from the intial BLAST protein sequences using HMMsearch and filtered by evalue EVALUE3. \ 
Domains from pfam can be used instead if placed as {HMM}/{TRANS}/{TRANS}.prof
Finally, the whole genome sequence is BLAST searched for the same protein family using EVALUE4 and any regions with significant hits that \
were not annotated are extracted with a specified flanking region.

The following options are available to continue the analysis at certain steps:

  -i,   Start from initial BLAST filtering to respecify evalue and coverage cutoff
  -j,   Start from re-BLAST filtering to respecify evalue and coverage cutoff
  -m,   Start from HMMsearch filtering to respecify evalue cutoff
  -g,   Start from genome BLAST search filtering
  
  "
  exit 1
}


#Check for user defined variables
while getopts s:r:q:e:c:E:C:l:v:n:Rijmgf:t:h flag
do
  case "$flag" in
    s) SPECIES_LIST=${OPTARG} ;;
    r) REF=${OPTARG} ;;
    q) TRANS_LIST=${OPTARG} ;;
    e) EVALUE1=${OPTARG} ;;
    c) COV1=${OPTARG} ;;
    E) EVALUE2=${OPTARG} ;;
    C) COV2=${OPTARG} ;;
    l) LENGTH=${OPTARG} ;;
    v) EVALUE3=${OPTARG} ;;
    n) EVALUE4=${OPTARG} ;;
    R) REBLAST=true ;;
    i) INITBLASTFILTER=true ;;
    j) REBLASTFILTER=true ;;
    m) HMMFILTER=true ;;
    g) GENOMEFILTER=true ;;
    f) FLANK=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    h|*) print_usage ;;
  esac
done

## DEFINE FUNCTIONS

#Extract proteins from reference. Edit depending on annotations and proteins.
extract_reference () {
  echo -e "Extracting query sequences..."
  for TRANS in $TRANS_LIST ; do
    mkdir -p ${TRANSFAS}/${TRANS}
    if [ ! -f ${QUERIES}/${TRANS}/${REF}_${TRANS}_query.fa ] ; then
      seqkit grep --by-name --use-regexp --ignore-case --pattern "gn=${TRANS}" ${PROTEOMES}/${REF}_proteome.fa > ${QUERIES}/${TRANS}/${REF}_${TRANS}_query.fa  
      echo -e "\nExtracted $TRANS sequences from $REF"
    else
      echo -e "\nFound $TRANS query sequences in ${QUERIES}/${TRANS}/${REF}_${TRANS}_query.fa"
    fi
  done
}

#Template for iBLAST and reBLAST
BLAST () {  
  #Make directories
  mkdir -p ${BLASTDB}/${SPECIES}_proteome_DB ${BLASTRES}/${TRANS}/${SPECIES}

  #Make blast database if it doesn't exists. Assumes makeblastdb creates .psq files
  if [ ! -f ${BLASTDB}/${SPECIES}_proteome_DB/${SPECIES}_db.psq ] ; then
    makeblastdb -dbtype prot -in ${PROTEOMES}/${SPECIES}_proteome.fa -out ${BLASTDB}/${SPECIES}_proteome_DB/${SPECIES}_db
    echo -e "\nMade BLAST database for ${SPECIES}"
  fi
  
  #Blast proteome with query sequences extracted from reference species
  blastp -db ${BLASTDB}/${SPECIES}_proteome_DB/${SPECIES}_db -query ${QUERY} -out ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}.blast -evalue 1 -max_hsps 500 -outfmt "7 sseqid evalue qcovs pident length stitle" -num_threads ${THREADS}
}

iBLAST () {
  echo -e "\nStarting initial BLASTs..."
  for SPECIES in $SPECIES_LIST ; do
    for TRANS in $TRANS_LIST ; do  
      QUERY=${QUERIES}/${TRANS}/${REF}_${TRANS}_query.fa
      STAGE=init
      BLAST
      echo -e "\nFinished initial BLAST for $TRANS in $SPECIES"
    done
  done
  #for TRANS in $TRANS_LIST ; do
    #cd ${BLASTRES}/${TRANS}
    #echo -e "\nNumber of unfiltered $TRANS hits:"
    #grep -Ev "^#" -c */*_${TRANS}_proteome_${STAGE}.blast >&1
  #done
}

reBLAST () {
  if [ $REBLAST = true ] ; then  
    echo -e "\nStarting re-BLASTs..."
    for SPECIES in $SPECIES_LIST ; do
      for TRANS in $TRANS_LIST ; do  
        QUERY=${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_init.fas
        STAGE=re
        BLAST
        echo -e "\nFinished re-BLAST for $TRANS in $SPECIES"
      done
    done
  fi
  #for TRANS in $TRANS_LIST ; do
    #cd ${BLASTRES}/${TRANS}
    #echo -e "\nNumber of unfiltered $TRANS hits:"
    #grep -Ev "^#" -c */*_${TRANS}_proteome_${STAGE}.blast >&1
  #done
}


#Template for filtering iBLAST and reBLAST hits based on evalue and coverage then extract full sequences
BLASTfilter () {
  for SPECIES in $SPECIES_LIST ; do
    for TRANS in $TRANS_LIST ; do
      
      #Make list of sequences too short
      seqkit seq -g -M $LENGTH ${PROTEOMES}/${SPECIES}_proteome.fa | grep ">" | cut -f 1 -d ' ' | tr -d '>' > ${BLASTRES}/${TRANS}/${SPECIES}/remove_seq.txt

      #Filter BLAST output to sort and deduplicate hits into pass and not pass (based on best evalue and coverage)
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}.blast | sort -k1,1 -k2,2g | sort -k1,1 -u | grep -vFw -f ${BLASTRES}/${TRANS}/${SPECIES}/remove_seq.txt | awk -v a="$EVALUE" -v b="$COV" -v c="$STAGE" -v d="${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_" '{if ($1 !~ /^#/ && $2 <= a && $3 >= b) print > d"proteome_"c"_pass.txt"
      else if ($1 !~ /^#/) print > d"proteome_"c"_notpass.txt"}' 
     
      #Print out highest e-values for pass
      if [ -f ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass.txt ] ; then
        echo -e "\nE-values, coverage and annotation of passed hits for ${TRANS} with highest E-values in $SPECIES:"
        cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass.txt | grep -vFw -f ${BLASTRES}/${TRANS}/${SPECIES}/remove_seq.txt | sort -k2g | cut -f 1,2,3,6 | tail -n 5
	echo ${SPECIES} >> ${BLASTRES}/${TRANS}/${TRANS}_pass_statistics.txt
	cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass.txt | sort -k2g | cut -f 1,2,6 | tail -n 1 | tr ' ' '_' >> ${BLASTRES}/${TRANS}/${TRANS}_pass_statistics.txt
      fi

      #Print out lowest e-values for not_pass
      if [ -f ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_notpass.txt ] ; then
        echo -e "\nE-values, coverage and annotation of filtered out hits for ${TRANS} with lowest E-values in $SPECIES:"
        cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_notpass.txt | sort -k2g | cut -f 1,2,3,6 | sed -n 1,5p
        echo ${SPECIES} >> ${BLASTRES}/${TRANS}/${TRANS}_notpass_statistics.txt
        cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_notpass.txt | sort -k2g | cut -f 1,2,6 | sed -n 1p | tr ' ' '_' >> ${BLASTRES}/${TRANS}/${TRANS}_notpass_statistics.txt
      fi

      #Trim pass files into ids
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass.txt | cut -f 1 | sort -u > ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass_ids.txt

      #Extract full sequences of hits from the proteome
      seqkit grep -f ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass_ids.txt ${PROTEOMES}/${SPECIES}_proteome.fa | seqkit rmdup -s | seqkit seq -g -m $LENGTH > ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_${STAGE}.fas
    
      #Remove uneeded files
      rm ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_proteome_${STAGE}_pass_ids.txt
    done
  done
}  

iBLASTfilter () {
  echo -e "\nStarting initial BLAST filtering..."
  STAGE=init
  EVALUE=${EVALUE1}
  COV=${COV1}
  BLASTfilter
  echo -e "\nFiltered initial BLAST results with E-value = ${EVALUE1} and coverage = ${COV1}%"
  for TRANS in $TRANS_LIST ; do
    #cp ${TRANSFAS}/${TRANS}/${REF}_${TRANS}_proteome_re.fas ${TRANSFAS}/${TRANS}/${REF}_${TRANS}_proteome_init.fas
    echo -e "\nStatistics for $TRANS identified by initial BLAST:"
    cd ${TRANSFAS}/${TRANS}
    seqkit stats *${TRANS}_proteome_${STAGE}.fas
  done
}

reBLASTfilter () {
  if [ $REBLAST = true ] ; then
    echo -e "\nStarting re-BLAST filtering..."
    STAGE=re
    EVALUE=${EVALUE2}
    COV=${COV2}
    BLASTfilter
    echo -e "\nFiltered re-BLAST results with E-value = ${EVALUE2} and coverage = ${COV2}%"
    for TRANS in $TRANS_LIST ; do
      echo -e "\nStatistics for $TRANS identified by re-BLAST:"
      cd ${TRANSFAS}/${TRANS}
      seqkit stats *${TRANS}_proteome_${STAGE}.fas
    done
  else
    for SPECIES in $SPECIES_LIST ; do
      for TRANS in $TRANS_LIST ; do
        cp ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_init.fas ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_re.fas
      done
    done
  fi
}

#Rename fasta headers to keep just Species and Annotation
#rename_sequences () {  
  #for SPECIES in $SPECIES_LIST $REF ; do
    #for TRANS in $TRANS_LIST ; do
      #sed -r "s/>(\S+ )(.+)(OS.+|\(.+)/>${SPECIES^} \2/" ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_re.fas \
      #| tr -d '()' | sed "s/LOW QUALITY PROTEIN//" > ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_rnd.fas  
      #echo -e "Renamed fastas for $TRANS in $SPECIES"
    #done
    #echo -e "Renamed fastas for $SPECIES"  
  #done
  #echo -e "Renamed all fastas"
#}

#Search for HMM protein domains in proteomes using profiles generated from iBLAST sequences
HMM_search () {
  echo -e "\nStarting HMMER searches..."
  for TRANS in $TRANS_LIST ; do
    #Check if profile is already there
    if [ ! -f ${HMM}/${TRANS}/${TRANS}.prof ] ; then
    
      #Make directories
      mkdir -p ${HMM}/${TRANS}

      #Align all iBLAST transporters with MUSCLE
      cat ${TRANSFAS}/${TRANS}/*proteome_init.fas > ${TRANSFAS}/${TRANS}/all_proteome_init.fas 
      muscle -align ${TRANSFAS}/${TRANS}/all_proteome_init.fas -output ${HMM}/${TRANS}/${TRANS}.afa > /dev/null 2>&1
      rm ${TRANSFAS}/${TRANS}/all_proteome_init.fas

      #Convert to stockholm format
      esl-reformat stockholm ${HMM}/${TRANS}/${TRANS}.afa > ${HMM}/${TRANS}/${TRANS}.sto

      #Build a HMM profile
      hmmbuild --cpu ${THREADS} ${HMM}/${TRANS}/${TRANS}.prof ${HMM}/${TRANS}/${TRANS}.sto > /dev/null 2>&1
      echo -e "\nMade HMM profile for $TRANS"
    else echo -e "\nFound HMM profile for $TRANS in ${HMM}/${TRANS}/${TRANS}.prof"
    fi

    #Search for the profile in each proteome  
    for SPECIES in $SPECIES_LIST ; do
      mkdir -p ${HMM}/${TRANS}/${SPECIES}
      hmmsearch --cpu ${THREADS} --noali --notextw --tblout ${HMM}/${TRANS}/${SPECIES}/seq.tbl ${HMM}/${TRANS}/${TRANS}.prof ${PROTEOMES}/${SPECIES}_proteome.fa > ${HMM}/${TRANS}/${SPECIES}/output.hmm
      echo -e "\nFinished HMM search for $TRANS in $SPECIES"
    done
  done
}

#Filter HMM_search by evalue
HMM_search_filter () {
  echo -e "\nStarting HMMER result filtering..."
  for TRANS in $TRANS_LIST ; do

    #Make HMM fasta directory and copy reference sequences into there
    mkdir -p ${HMMFAS}/${TRANS}
    #cp ${TRANSFAS}/${TRANS}/${REF}_${TRANS}_proteome_re.fas ${HMMFAS}/${TRANS}
    
    #Filter and extract hmm hits for all species
    for SPECIES in $SPECIES_LIST ; do    
      
      #Trim seq.tbl to tab delimited, filter for short sequences and extract sequence ids based on evalues of all domain and best domain hits
      cat ${HMM}/${TRANS}/${SPECIES}/seq.tbl | sed 's/ \+/\t/g' | grep -vFw -f ${BLASTRES}/${TRANS}/${SPECIES}/remove_seq.txt | awk -v a="$EVALUE3" -v d="${HMM}/${TRANS}/${SPECIES}/" '{if ($1 !~ /^#/ && $5 <= a) print > d"pass.txt"
      else if ($1 !~ /^#/) print > d"not_pass.txt"}'

      #Print out highest e-values for pass
      if [ -f ${HMM}/${TRANS}/${SPECIES}/pass.txt ] ; then
        echo -e "\nHighest E-values of passed hits for ${TRANS} in $SPECIES:"
        cat ${HMM}/${TRANS}/${SPECIES}/pass.txt | sort -k5g | cut -f 1,5 | tail -n 5
	echo ${SPECIES} >> ${HMM}/${TRANS}/${TRANS}_pass_statistics.txt
        cat ${HMM}/${TRANS}/${SPECIES}/pass.txt | sort -k5g | cut -f 1,5 | tail -n 1 >> ${HMM}/${TRANS}/${TRANS}_pass_statistics.txt
      fi

      #Print out lowest e-values for not_pass
      if [ -f ${HMM}/${TRANS}/${SPECIES}/not_pass.txt ] ; then
        echo -e "\nLowest E-values of filtered out hits for ${TRANS} in $SPECIES:"
        cat ${HMM}/${TRANS}/${SPECIES}/not_pass.txt | sort -k5g | cut -f 1,5 | sed -n 1,5p
        echo ${SPECIES} >> ${HMM}/${TRANS}/${TRANS}_notpass_statistics.txt
        cat ${HMM}/${TRANS}/${SPECIES}/not_pass.txt | sort -k5g | cut -f 1,5 | sed -n 1p >> ${HMM}/${TRANS}/${TRANS}_notpass_statistics.txt	
      fi

      #Extract full sequences of hits from the proteome
      cut -f 1 ${HMM}/${TRANS}/${SPECIES}/pass.txt > ${HMM}/${TRANS}/${SPECIES}/pass_ids.txt
      seqkit grep -f ${HMM}/${TRANS}/${SPECIES}/pass_ids.txt ${PROTEOMES}/${SPECIES}_proteome.fa | seqkit rmdup -s | seqkit seq -g -m $LENGTH > ${HMMFAS}/${TRANS}/${SPECIES}_${TRANS}_hmm.fas

      #Remove uneeded files
      rm ${HMM}/${TRANS}/${SPECIES}/pass_ids.txt
    done
  done
  echo -e "\nFiltered HMM outputs with E-value $EVALUE3"
  
  #HMMER result statistics
  for TRANS in $TRANS_LIST ; do
    cd ${HMMFAS}/${TRANS}
    echo -e "\nStatistics for $TRANS proteins identified with HMMER:"
    seqkit stats *${TRANS}_hmm.fas
  done
}

#BLAST of genome sequence to find unannotated genes
check_genomes () { 
  echo -e "\nStarting genome BLASTs"
  for SPECIES in $SPECIES_LIST ; do
    for TRANS in $TRANS_LIST ; do
      
      #Make directories
      mkdir -p ${BLASTDB}/${SPECIES}_genome_DB ${BLASTRES}/${TRANS}/${SPECIES} ${TRANSFAS}/${TRANS}

      #Make blast db if it doesn't exists. Assumes makeblastdb creates .psq files
      if [ ! -f ${BLASTDB}/${SPECIES}_genome_DB/${SPECIES}_db.nsq ] ; then
        makeblastdb -dbtype nucl -in ${GENOMES}/${SPECIES}_genome.fa -out ${BLASTDB}/${SPECIES}_genome_DB/${SPECIES}_db
      fi

      #Blast genome with initial BLAST transporter query sequences
      if [ ! -f ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome.blast ] ; then
        tblastn -db ${BLASTDB}/${SPECIES}_genome_DB/${SPECIES}_db -query ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_init.fas -out $BLASTRES/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome.blast -evalue 1 -max_hsps 500 -outfmt 7 -num_threads ${THREADS}
      fi
      echo -e "\nFinished genome BLAST for $TRANS in $SPECIES"
    done
  done
}

#Filter genome BLASTs based on evalue and coverage and check if they overlap with annotations in the GFF3 file
check_genomes_filter () {
  echo -e "\nStarting genome BLAST filtering at E-value $EVALUE4 and protein prediction"
  for SPECIES in $SPECIES_LIST ; do
    for TRANS in $TRANS_LIST ; do

      #Filter BLAST output to sort and deduplicate hits into pass and not pass (based on best evalue and coverage) -u in sort no 2
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome.blast | awk -v a="$EVALUE4" -v d="${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_" '{if ($1 !~ /^#/ && $11 <= a) print > d"genome_pass.txt"
      else if ($1 !~ /^#/) print > d"genome_notpass.txt"}'

      #Print out highest e-values for pass
      echo -e "\nE-value and length of passed hits with highest E-values in $SPECIES:"
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass.txt | sort -k11g | tail -n 5 | cut -f 1,2,4,11

      #Switch coordinates of hits if inverted, print out chromosome coordinates, check if they overlap with the gff file, and merge coordinates
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass.txt | awk '{if ($10>$9) ($14=$9) && ($15=$10); else ($14=$10) && ($15=$9); print $2, $14, $15, OFS = "\t"}' | sed -r 's/\t+$//' | bedtools intersect -v -nonamecheck -a stdin -b ${GENOMES}/${SPECIES}.gff | sort -k1,1 -k2,2n | bedtools merge > ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated.txt

      #Add flanks and print out coordinate file for protein sequence prediction
      cat ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated.txt | awk -v a="$FLANK" 'FS="\t" {($4=$2-a) && ($5=$3+a); print $1, $4, $5, OFS = "\t"}' | sed -r 's/\t+$//' | sed 's/\t/:/g' | sed -r 's/:([^:]*)$/-\1/' > ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.txt

      if [ -s ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.txt ] ; then 
        #Extract regions using samtools faidx (skipping augustus by renaming to _genome.fas not _flanks.fa
        samtools faidx ${GENOMES}/${SPECIES}_genome.fa -r ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.txt > ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_genome.fas 
     
        #Run Augustus
        #${AUG}/augustus --protein=on --species=arabidopsis ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.fa > ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.gff 
    
        #Extract protein sequences
        #getAnnoFasta.pl ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.gff

        #Copy to transfas directory
        #cp ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.aa ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_genome.fas
      else
        cp ${BLASTRES}/${TRANS}/${SPECIES}/${SPECIES}_${TRANS}_genome_pass_unannotated_flanks.txt ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_genome.fas 
      fi
    echo -e "\nFiltered genome BLASTn results for $TRANS in $SPECIES"
    done
  done
  for TRANS in $TRANS_LIST ; do
    cd ${TRANSFAS}/${TRANS}
    echo -e "\nStatistics of unannotated regions which have BLASTn hits above ${EVALUE4}"
    seqkit stats *${TRANS}_genome.fas
  done
}

combine () {
  for SPECIES in $SPECIES_LIST ; do
    for TRANS in $TRANS_LIST ; do
    
    #Make directory and copy over reference proteins
    mkdir -p ${FINALFAS}/${TRANS} 
    #cp ${TRANSFAS}/${TRANS}/${REF}_${TRANS}_proteome_re.fas ${FINALFAS}/${TRANS}    

    #Cat all result files, remove duplicates, filter by length, and rename sequences 
    cat ${TRANSFAS}/${TRANS}/${SPECIES}_${TRANS}_proteome_re.fas ${HMMFAS}/${TRANS}/${SPECIES}_${TRANS}_hmm.fas | seqkit rmdup -s | seqkit seq -g -m $LENGTH | sed -r "s/>(\S+)(.+)/>${SPECIES^}_${TRANS^^}_ \1/" | seqkit replace -p "^(.+) " -r "\${1}{nr} " > ${FINALFAS}/${TRANS}/${SPECIES}_${TRANS}_final.fa
    done
  done
  
  #Statistics of final protein sequences
  for TRANS in $TRANS_LIST ; do
    cd ${FINALFAS}/${TRANS}
    echo -e "\nStatistics of final $TRANS protein sequences:"
    seqkit stats *_final.fa
  done
}

## BEGIN STEP CHECKING

if [ $INITBLASTFILTER = true ] ; then
  iBLASTfilter ; reBLAST ; reBLASTfilter ; HMM_search ; HMM_search_filter ; check_genomes ; check_genomes_filter ; combine
elif [ $REBLASTFILTER = true ] ; then
  reBLASTfilter ; HMM_search ; HMM_search_filter ; check_genomes ; check_genomes_filter ; combine
elif [ $HMMFILTER = true ] ; then
  HMM_search_filter ; check_genomes ; check_genomes_filter ; combine
elif [ $GENOMEFILTER = true ] ; then
  check_genomes_filter ; combine
else 
  extract_reference ; iBLAST ; iBLASTfilter ; reBLAST ; reBLASTfilter ; HMM_search ; HMM_search_filter ; combine
fi
#Removed check_genomes and check_genomes_filter 02/06/23
