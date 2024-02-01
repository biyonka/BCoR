#!/bin/bash
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --time=
#SBATCH --partition=  # Partition to submit to
#SBATCH --mem= # Memory pool for all cores (see also --mem-per-cpu) can lower memory
#SBATCH -o output/hostname_%j.out # File to which STDOUT will be written
#SBATCH -e err/hostname_%j.err # File to which STDERR will be written

module load python/3.10.9-fasrc01
source activate gurobi_env

Rscript bayesrmab_greedy.R ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}
Rscript bayesrmab_whittle.R ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}
Rscript tsgreedy.R ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}
python ucwhittle/main_ucw.py -N ${N} -B ${budget} -H ${t} -T 1 --data ns -E 1 --iter ${E} --k ${k} --dim_time_spline ${d} --sigma_a ${sigma_a} --sigma_b ${sigma_b} --mu_beta ${mu_beta} --sigma_mu_beta ${sigma_mu_beta} --mu_B ${mu_B} --sigma_X ${sigma_X} --sigma_B ${sigma_B} --sigma ${sigma} --tau ${tau}
Rscript compile_csvs.R ${total_iterations} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}

