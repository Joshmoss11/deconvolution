#!/bin/csh
#SBATCH --mem=8000
#SBATCH --time=00:20:00

if ($#argv != 1) then
    echo "Usage: $0 < sim id>"
        exit 0
endif
srun $Rscript  /cs/icore/joshua.moss/scripts/deconvolution/nnls_sim.R $1
