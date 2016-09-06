#!/bin/tcsh
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4

srun $Rscript simulate_sets.R
