#!/bin/bash -l

# Script for running SLURM jobs using multithreaded programs on Pawsey Setonix (2022)
# Copyright statement: Copyright (c) 2022 Applied Bioinformatics Group, UWA, Perth WA, Australia

# Available partitions:
##	Name	Time limit (h)	Cores/node	No. nodes	MEM/node (Gb)	Mem/core (Gb)	
## 	work	24		128		316		230		~2
## 	long	96		128		8		230		~2
##	highmem 24		128		8		980		~8
##	copy	24		64		8		118		~2
##	debug	1		128		8		230		~2

# User defined SLURM commands. #CPUS per task (no. threads program uses) should be a multiple of 8.
#SBATCH --job-name=find_genes
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=24
#SBATCH --partition=work
#SBATCH --output=1_find_genes_job
#SBATCH --error=1_find_genes_err
#SBATCH --mail-user=21979717@student.uwa.edu.au

# Generic SLURM commands
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --clusters=setonix
#SBATCH --account=pawsey0149
#SBATCH --mail-type=ALL
#SBATCH --export=NONE 
source /software/projects/pawsey0149/groupEnv/ivec/groupResource.cfg

# SLURM useful commands: sbatch, squeue, scancel
# Run this script using sbatch [slurm].sh

echo "========================================="
echo "SLURM_JOB_ID = $SLURM_JOB_ID"
echo "SLURM_NODELIST = $SLURM_NODELIST"
if [ ! -z $SLURM_ARRAY_TASK_ID ]; then
	echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
fi
echo "========================================="

# Script to run (srun -m command recommended by Pawsey to pack threads) 
cd /scratch/pawsey0149/bnestor/2020_11_20_Transporters/0-scripts
conda activate blast
module load samtools/1.15--h3843a85_0
module load blast/2.12.0--pl5262h3289130_0
module load bedtools/2.30.0--h468198e_3

#time srun -m block:block:block bash find_genes_v2.sh -s "Nelumbo Medicago Eucalyptus Vitis Zea" -q nrt2 > 1_find_genes_out/nrt2_other.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Nelumbo Medicago Eucalyptus Vitis Zea" -q nar2 -l 150 > 1_find_genes_out/nar2_other.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Hakea Telopea Macadamia" -q npf > results_npf_proteaceae
#time srun -m block:block:block bash find_genes_v2.sh -s "Nelumbo Medicago Eucalyptus Vitis Zea" -q npf -t 64 > 1_find_genes_out/npf_other.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Amborella Oryza" -q "pht2" -e 1e-25 -v 1e-40 -l 300 -t 24 > 1_find_genes_out/pht2.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Amborella Oryza" -q "pht3" -e 1e-25 -v 1e-40 -l 150 -t 24 > 1_find_genes_out/pht3.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Amborella Oryza" -q "pht4" -e 1e-25 -v 1e-40 -l 200 -t 24 > 1_find_genes_out/pht4.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Amborella Oryza" -q "pht5" -e 1e-25 -v 1e-40 -l 350 -t 24 > 1_find_genes_out/pht5.txt
#time srun -m block:block:block bash find_genes_v2.sh -s "Amborella" -q "pho1" -r Arabidopsis_Oryza -e 1e-40 -v 1e-50 -t 24 > 1_find_genes_out/pho1.txt
#time srun -m block:block:block bash find_genes_v2.sh -q "amt1 amt2" -r Arabidopsis_Oryza -e 1e-40 -v 1e-50 -t 24 > 1_find_genes_out/amt.txt
time srun -m block:block:block bash find_genes_v2.sh -q "clc slac" -r Arabidopsis -e 1e-40 -v 1e-50 -t 24 > 1_find_genes_out/clc_slac.txt
