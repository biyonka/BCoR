# A Bayesian Approach to Online Learning for Contextual Restless Bandits with Applications to Public Health

We present Bayesian Learning for Contextual RMABs (BCoR), an online RL approach for restless bandits that novelly combines techniques in Bayesian modeling with Thompson sampling to flexibly model a wide range of complex restless bandit, such as contextual and non-stationary RMABs, which we design for a real-world public health intervention.

This repository contains the implementation of the BCoR algorithm and all files for reproducing the results of Section 5.2 in the paper. Though we cannot release the real data simulator or anonymized dataset used in Section 5.3, the implementation of all methods under comparison were the same for both Sections 5.2 and 5.3, and hence, the implementations provided here are equivalent to what was used in the real data example. All scripts were run on a high-performance computing cluster from 2023-2024. 

 This repository contains following directories and files:

1. `submit.sh`, `batch.sh`: These are the bash script specifying the N, T, B, and prior configurations used for submitting each of the four settings used explored in Section 5.2 to the Slurm Workload Manager on Harvard University's computing cluster.

2. `bcor_whittle.R`: This is the R script used to run BCoR with the Whittle index policy.

3. `bcor_greedy.R`: This is the R script used to run BCoR with the greedy policy.

4. `tsgreedy.R`: This is the R script used to run TS with the greedy policy.

5. `compile_csvs.R`: This is the R script for compiling all data files and generating the plots shown in Section 5.2.

6. `bcor_model.stan`: This is the stan file used for BCoR, for both its Whittle index policy and Greedy policy implementations.

7. `ucwhittle`: We cloned this directory from the public online repository for implementing UCWhittle (provided here: https://github.com/lily-x/online-rmab) in order to run UCWhittle for our experimental results. We added implementations of the Whittle oracles described in the paper, and TS with the Whittle index policy, and adapted their RMAB simulator to handle the non-stationary transition dynamics we considered in our paper. No modifications were made to the implementation of UCWhittle, which remains the original authors' implementation. 