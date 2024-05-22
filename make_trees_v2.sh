#!/bin/bash
set -euo pipefail

## DEFINE VARIABLES

#Set default variables. Use "" for lists
SPECIES_LIST="Amborella Arabidopsis Eucalyptus Hakea Macadamia Medicago Nelumbo Oryza Protea Telopea Vitis Zea"
NAME="all"
REF=Arabidopsis
TRANS=nrt2
BOOTSTRAPS=1000
OUTGROUP=AtrNRT2
THREADS=4

#Path variables
QUERIES=${MYSCRATCH}/2020_11_20_Transporters/1-data/blast_queries
FINALFAS=${MYSCRATCH}/2020_11_20_Transporters/6-final_fas

#Check for user defined variables
while getopts s:n:g:r:b:o:t:h flag
do
  case "$flag" in
    s) SPECIES_LIST=${OPTARG} ;;
    n) NAME=${OPTARG} ;;
    g) TRANS=${OPTARG} ;;
    r) REF=${OPTARG} ;;
    b) BOOTSTRAPS=${OPTARG} ;;
    o) OUTGROUP=${OPTARG} ;;
    t) THREADS=${OPTARG} ;;
    h|*) print_usage ;;
  esac
done

## DEFINE FUNCTIONS

#Make the tree
make_tree () {
	#Make folder for trees and copy all final.fa (also from queries) into it
	mkdir -p ${FINALFAS}/${TRANS}/${NAME}_tree_building
	cp ${QUERIES}/${TRANS}/${REF}_${TRANS}_query.fa ${FINALFAS}/${TRANS}/${REF}_${TRANS}_final.fa
	if [ -f ${QUERIES}/${TRANS}/*_final.fa ] ; then
	  cp ${QUERIES}/${TRANS}/*_final.fa ${FINALFAS}/${TRANS}/
	fi
	for SPECIES in $SPECIES_LIST ; do
		cat ${FINALFAS}/${TRANS}/${SPECIES}_${TRANS}_final*.fa >> ${FINALFAS}/${TRANS}/${NAME}_tree_building/${TRANS}_${NAME}.fa
	done
	cd ${FINALFAS}/${TRANS}/${NAME}_tree_building
	
        #Align sequences with MUSCLE
	echo -e "\nAligning $TRANS sequences for ${NAME}"
	mafft --thread ${THREADS} --maxiterate 1000 --localpair ${TRANS}_${NAME}.fa > ${TRANS}_${NAME}.afa
	echo -e "\nFinished aligning $TRANS sequences"

        #Make tree with RAxML
	echo -e "\nBuilding phylogenetic tree of $TRANS sequences for ${NAME}"
	raxmlHPC-AVX2 -T ${THREADS} -o ${OUTGROUP} -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -N ${BOOTSTRAPS} -s ${TRANS}_${NAME}.afa -n ${TRANS}_${NAME}
	echo -e "\nFinished phylogenetic tree of $TRANS sequences" 
}

make_tree
