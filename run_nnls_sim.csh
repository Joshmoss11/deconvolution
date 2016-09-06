#!/bin/csh
mkdir nnls_sim_out
foreach s (`seq 1 100`)
	sbatch --output=nnls_sim_out/set_${s}.out /cs/icore/joshua.moss/scripts/deconvolution/nnls_sim.csh $s
end
