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
#SBATCH --job-name=make_tree
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=work
#SBATCH --output=2_make_trees_job
#SBATCH --error=2_make_trees_err
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
conda activate raxml

#time srun -m block:block:block make_trees.sh -n proteaceae -s "Hakea Telopea Macadamia Arabidopsis Amborella" -g nrt2 -o AtrNRT2
#time srun -m block:block:block make_trees.sh -n proteaceae -s "Hakea Telopea Macadamia Arabidopsis Amborella" -g nar2 -o AtrNAR2
#time srun -m block:block:block make_trees.sh -n all -s "Hakea Telopea Macadamia Nelumbo Medicago Eucalyptus Vitis Zea Arabidopsis Amborella" -g nrt2 -o AtrNRT2
#time srun -m block:block:block make_trees.sh -n all -s "Hakea Telopea Macadamia Nelumbo Medicago Eucalyptus Vitis Oryza Zea Arabidopsis Amborella" -g nrt2 -o AtrNRT2 -b 500
#time srun -m block:block:block make_trees.sh -n all -s "Hakea Telopea Macadamia Nelumbo Medicago Eucalyptus Vitis Oryza Zea Arabidopsis Amborella" -g nar2 -o AtrNAR2 -b 1000

#time srun -m block:block:block make_trees.sh -t 120 -n proteaceae -s "Hakea Telopea Macadamia Arabidopsis Amborella" -g npf -o AtrNPF6.1 -b 50
#time srun -m block:block:block make_trees.sh -t 120 -n all -s "Hakea Telopea Macadamia Nelumbo Medicago Eucalyptus Vitis Zea Arabidopsis Amborella" -g npf -o AtrNPF6.1 -b 10
#time srun -m block:block:block make_trees.sh -t 128 -n all -g pht2 -o AtrPHT2.1 -b 1000
#time srun -m block:block:block make_trees.sh -t 4 -n all -g pht3 -o AtrPHT3.1 -b 1000
#time srun -m block:block:block make_trees.sh -t 8 -n all -g pht4 -o AtrPHT4.1 -b 500
#time srun -m block:block:block make_trees.sh -t 4 -n all -g pht5 -o AtrPHT5.1 -b 1000
time srun -m block:block:block make_trees.sh -t 16 -n all -g pho1 -o AtrPHO1-1 -b 500
#time srun -m block:block:block make_trees.sh -t 60 -n all -g amt1 -o Amborella_AMT1_1 -b 10
#time srun -m block:block:block make_trees.sh -t 60 -n all -g amt2 -o Amborella_AMT2_1 -b 10
#time srun -m block:block:block make_trees.sh -t 60 -n all -g clc -o Amborella_CLC_1 -b 10
#time srun -m block:block:block make_trees.sh -t 60 -n all -g slac -o Amborella_SLAC_1 -b 10
