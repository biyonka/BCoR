#all effects
for total_iterations in 1000; do
mkdir -p data
for N in 400; do
for k in 8; do
for d in 3; do
for t in 50; do
for budget in 10; do
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
    mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
for sigma_a in 0.1; do #zero out for no info between arms
for sigma_b in 0.1; do #zero out for no info between arms
for mu_beta in 0; do
for sigma_mu_beta in 0.3; do #zero out for no covariate effects
for mu_B in 0; do
for sigma_X in 0.1; do #zero out for no covariate effects
for sigma_B in 0.3; do #zero out for no time effects
for sigma in 1; do
for tau in 100; do
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
    mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
    mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_data
        mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_plots
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_runs
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_params
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_rewards
    #mkdir data/N_${N}_T_${t}_B_${B}_k_${k}_d_${d}/posterior_samples/sigmaa_${sigma_a}_sigmab_${sigma_b}_muX_${mu_X}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/plots
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/posterior_samples
    mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/final_data
for E in $(seq 1 $total_iterations); do
#for r in $(seq 2 $m); do
#for k in $(seq 0 $m); do

 echo "${total_iterations} ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}"
export total_iterations E N t budget k d sigma_a sigma_b mu_beta sigma_mu_beta mu_B sigma_X sigma_B sigma tau


sbatch -o output/ts_out_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
-e err/ts_err_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
--job-name="ts_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}" \
batch.sh
#batch_final_br.sh
#
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done


#Zero out the appropriate pieces of the model to run the misspecified settings. For example:

# #no time effects
# for total_iterations in 1000; do
# mkdir -p data
# for N in 400; do
# for k in 8; do
# for d in 3; do
# for t in 50; do
# for budget in 10; do
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
# for sigma_a in 0.1; do #zero out for no info between arms
# for sigma_b in 0.1; do #zero out for no info between arms
# for mu_beta in 0; do
# for sigma_mu_beta in 0.3; do #zero out for no covariate effects
# for mu_B in 0; do
# for sigma_X in 0.1; do #zero out for no covariate effects
# for sigma_B in 0; do #zero out for no time effects
# for sigma in 1; do
# for tau in 100; do
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_data
#         mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_plots
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_runs
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_params
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_rewards
#     #mkdir data/N_${N}_T_${t}_B_${B}_k_${k}_d_${d}/posterior_samples/sigmaa_${sigma_a}_sigmab_${sigma_b}_muX_${mu_X}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/plots
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/posterior_samples
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/final_data
# for E in $(seq 1 $total_iterations); do
# #for r in $(seq 2 $m); do
# #for k in $(seq 0 $m); do

#  echo "${total_iterations} ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}"
# export total_iterations E N t budget k d sigma_a sigma_b mu_beta sigma_mu_beta mu_B sigma_X sigma_B sigma tau


# sbatch -o output/ts_out_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
# -e err/ts_err_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
# --job-name="ts_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}" \
# batch_final.sh
# #batch_final_br.sh
# #
# #
# sleep 1 # pause to be kind to the scheduler
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done


# #no covariate effects
# for total_iterations in 1000; do
# mkdir -p data
# for N in 400; do
# for k in 8; do
# for d in 3; do
# for t in 50; do
# for budget in 10; do
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}
# for sigma_a in 0.1; do
# for sigma_b in 0.1; do
# for mu_beta in 0; do
# for sigma_mu_beta in 0; do
# for mu_B in 0; do
# for sigma_X in 0; do
# for sigma_B in 0.3; do
# for sigma in 1; do
# for tau in 100; do
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_data
#         mkdir -p ANONYMIZED FILEPATH/data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/transition_plots
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_runs
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_params
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/iteration_rewards
#     #mkdir data/N_${N}_T_${t}_B_${B}_k_${k}_d_${d}/posterior_samples/sigmaa_${sigma_a}_sigmab_${sigma_b}_muX_${mu_X}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/plots
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/posterior_samples
#     mkdir -p data/N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}/sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}/final_data
# for E in $(seq 1 $total_iterations); do
# #for r in $(seq 2 $m); do
# #for k in $(seq 0 $m); do

#  echo "${total_iterations} ${E} ${N} ${t} ${budget} ${k} ${d} ${sigma_a} ${sigma_b} ${mu_beta} ${sigma_mu_beta} ${mu_B} ${sigma_X} ${sigma_B} ${sigma} ${tau}"
# export total_iterations E N t budget k d sigma_a sigma_b mu_beta sigma_mu_beta mu_B sigma_X sigma_B sigma tau


# sbatch -o output/out_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
# -e err/err_E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}.stdout.txt \
# --job-name="E_${E}_N_${N}_T_${t}_B_${budget}_k_${k}_d_${d}_sigmaa_${sigma_a}_sigmab_${sigma_b}_mubeta_${mu_beta}_sigmamubeta_${sigma_mu_beta}_muB_${mu_B}_sigmaX_${sigma_X}_sigmaB_${sigma_B}_sigma_${sigma}_tau_${tau}" \
# batch_final.sh
# #batch_final_br.sh
# #
# #
# sleep 1 # pause to be kind to the scheduler
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done
# done



