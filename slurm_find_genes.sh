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

#Example
time srun -m block:block:block bash find_genes_v2.sh -q "pht1 nrt2" -r Arabidopsis -e 1e-40 -v 1e-50 -t 24 > 1_find_genes_out/nrt2_pht1.txt
