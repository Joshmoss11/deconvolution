#!/bin/csh
mkdir sim_out
foreach s (`seq 1 100`)
	sbatch --output=sim_out/set_${s}.out /cs/icore/joshua.moss/scripts/deconvolution/simulate_single_set.csh $s
end
