library(ggplot2)
library(dplyr)
library(readr)
library(matrixStats)

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('total_iterations', 'N', 't', 'budget', 'k', 'd', 'sigma_a', 'sigma_b', 'mu_beta', 'sigma_mu_beta', 'mu_B', 'sigma_X', 'sigma_B', 'sigma', 'tau')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  total_iterations = 1000
  N = 400 
  t = 50 
  budget = 10
  k = 4
  d = 3
  sigma_a = 0.1
  sigma_b =0.1
  sigma=1
  tau=100 
  mu_beta = 0
  sigma_mu_beta = 0.3
  sigma_X = 0.1
  sigma_B = 0.3
  mu_B = 0
} else {
  for (i in 1:length(args)) {
    eval (parse (text = paste(arguments[i]," = ",args[[i]],sep="") ))
  }
}

base_path =  paste0("data/N_",N,'_T_',t,'_B_',budget,'_k_',k,'_d_',d,'/sigmaa_',sigma_a,'_sigmab_',sigma_b,'_mubeta_',mu_beta,'_sigmamubeta_',sigma_mu_beta,'_muB_',
                    mu_B,'_sigmaX_',sigma_X,'_sigmaB_',sigma_B,'_sigma_',sigma,'_tau_',tau)


files = list.files(path=paste0(base_path, '/iteration_runs'), full.names = TRUE) 

if (length(files) == 3*total_iterations){
  df <- files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 
  
 
  ts_whittle_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_ts_whittle'), full.names = TRUE)
  ts_whittle_files = ts_whittle_files[sort(as.numeric(gsub("\\D", "", ts_whittle_files )), index=TRUE)$ix]
  ucw_ts_whittle =  ts_whittle_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 

  qp_whittle_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_ucw_qp'), full.names = TRUE)
  qp_whittle_files = qp_whittle_files[sort(as.numeric(gsub("\\D", "", qp_whittle_files )), index=TRUE)$ix]
  ucw_qp  =  qp_whittle_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 
  
  random_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_random'), full.names = TRUE)
  random_files = random_files[sort(as.numeric(gsub("\\D", "", random_files )), index=TRUE)$ix]
  random =  random_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 
  
  optimal_all_time_avg_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_optimal_all_time_avg'), full.names = TRUE)
  optimal_all_time_avg_files =  optimal_all_time_avg_files[sort(as.numeric(gsub("\\D", "",  optimal_all_time_avg_files)), index=TRUE)$ix]
  optimal_all_time_avg =  optimal_all_time_avg_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 
  
  optimal_current_time_avg_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_optimal_current_time_avg'), full.names = TRUE)
  optimal_current_time_avg_files =  optimal_current_time_avg_files[sort(as.numeric(gsub("\\D", "",  optimal_current_time_avg_files)), index=TRUE)$ix]
  optimal_current_time_avg =  optimal_current_time_avg_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 
  
  optimal_current_time_files = list.files(path=paste0(base_path, '/final_data/whittle_methods_optimal_current_time'), full.names = TRUE)
  optimal_current_time_files =  optimal_current_time_files[sort(as.numeric(gsub("\\D", "",  optimal_current_time_files)), index=TRUE)$ix]
  optimal_current_time =  optimal_current_time_files %>% 
    lapply(read_csv, show_col_types=FALSE) %>% 
    bind_rows 

  #take avg of rewards across the seeds
  #excludes initial "reward" of starting state
  random_rewards = sapply(seq(14, 14+t-1), function(i){random[[i]]})
  ucw_qp_rewards = sapply(seq(14, 14+t-1), function(i){ucw_qp[[i]]})-random_rewards
  ts_whittle_rewards = sapply(seq(14, 14+t-1), function(i){ ucw_ts_whittle[[i]]})-random_rewards
  optimal_all_time_avg_rewards = sapply(seq(14, 14+t-1), function(i){optimal_all_time_avg[[i]]})-random_rewards
  optimal_current_time_avg_rewards = sapply(seq(14, 14+t-1), function(i){optimal_current_time_avg[[i]]})-random_rewards
  optimal_current_time_rewards = sapply(seq(14, 14+t-1), function(i){optimal_current_time[[i]]})-random_rewards
  
 
  #random_rewards is sum of rewards at each timestep, t x total_iterations matrix
  rand_reward_sub = rep(as.vector(random_rewards), 3)
  

  reward_df_test = df %>% group_by(Method, time, iteration) %>% summarize(reward_per_it = sum(y))
  reward_df_test['reward_minus_random'] = reward_df_test$reward_per_it - rand_reward_sub
  
  rewards_df = reward_df_test %>% group_by(Method, time) %>% summarize(avg_reward = sum(reward_minus_random)/total_iterations)%>% mutate(cum_avg_reward = cumsum(avg_reward)) %>%
    mutate(time_avg_reward = cumsum(avg_reward)/time)
  
  ses_reward_per_it = reward_df_test %>% group_by(Method, time) %>% summarize(se_reward = sd(reward_minus_random))
  
  ses_cumsum_per_it = reward_df_test %>% group_by(Method, iteration) %>%
    mutate(cumsum_reward_per_it = cumsum(reward_minus_random)) %>% group_by(Method, time) %>% summarize(se_cum_reward = sd(cumsum_reward_per_it))
  
  ses_time_avg_reward_per_it =  reward_df_test %>% group_by(Method, iteration) %>%
    mutate(time_avg_reward_per_it = cumsum(reward_minus_random)/time)  %>% group_by(Method, time) %>% summarize(se_time_avg_reward = sd(time_avg_reward_per_it))

  rewards_df['Method_Type'] =  c(rep('BCoR', t),rep('Greedy Oracle', t), rep('BCoR', t), rep('TS', t))
  #puts method type column first
  print('made it here 2')
  rewards_df <- rewards_df[, c(6, 1, 2, 3, 4, 5)]
  print('made it here 3')
  rewards_df['ses_reward'] = ses_reward_per_it$se_reward/sqrt(total_iterations)
  rewards_df['ses_cumsum'] = ses_cumsum_per_it$se_cum_reward/sqrt(total_iterations)
  rewards_df['ses_timeavg'] = ses_time_avg_reward_per_it$se_time_avg_reward/sqrt(total_iterations)
  
  final_path = paste0(base_path, '/final_data')
  write.csv(df, file = paste0(final_path, '/full_rmab_run.csv'), row.names = FALSE)
  print('download final data')

  random_rewards_sub = random_rewards - random_rewards
  whittle_rewards = data.frame(
    Method_Type = c(
      rep('UCW', t),
      rep('TS', t),
      rep('Random', t),
      rep('Whittle Oracles', t),
      rep('Whittle Oracles', t),
      rep('Whittle Oracles', t)
    ),
    Method = c(
      rep('qp', t),
      rep('ts_whittle', t),
      rep('random', t),
      rep('all_time_avg', t),
      rep('current_time_avg', t),
      rep('current_time', t)
    ),
    time = rep(seq(1, t, 1), 6),
    avg_reward = c(
      colSums(ucw_qp_rewards )/total_iterations,
      colSums(ts_whittle_rewards)/total_iterations,
      colSums(random_rewards_sub)/total_iterations,
      colSums(optimal_all_time_avg_rewards)/total_iterations,
      colSums(optimal_current_time_avg_rewards)/total_iterations,
      colSums(optimal_current_time_rewards)/total_iterations
    ),
    cum_avg_reward = c(
      cumsum(colSums(ucw_qp_rewards )/total_iterations),
      cumsum(colSums(ts_whittle_rewards)/total_iterations),
      cumsum(colSums(random_rewards_sub)/total_iterations),
      cumsum(colSums(optimal_all_time_avg_rewards)/total_iterations),
      cumsum(colSums(optimal_current_time_avg_rewards)/total_iterations),
      cumsum(colSums(optimal_current_time_rewards)/total_iterations)
    ),
    time_avg_reward =  c(
      cumsum(colSums(ucw_qp_rewards )/total_iterations)/seq(1, t),
      cumsum(colSums(ts_whittle_rewards)/total_iterations)/seq(1, t),
      cumsum(colSums(random_rewards_sub)/total_iterations)/seq(1, t),
      cumsum(colSums(optimal_all_time_avg_rewards)/total_iterations)/seq(1, t),
      cumsum(colSums(optimal_current_time_avg_rewards)/total_iterations)/seq(1, t),
      cumsum(colSums(optimal_current_time_rewards)/total_iterations)/seq(1, t)
    ),
    ses_reward = c(
      apply(ucw_qp_rewards, 2, sd)/sqrt(total_iterations),
      apply(ts_whittle_rewards, 2, sd)/sqrt(total_iterations), 
      apply(random_rewards_sub, 2, sd)/sqrt(total_iterations),
      apply(optimal_all_time_avg_rewards, 2, sd)/sqrt(total_iterations),
      apply(optimal_current_time_avg_rewards, 2, sd)/sqrt(total_iterations),
      apply(optimal_current_time_rewards, 2, sd)/sqrt(total_iterations)
    ),
    ses_cumsum = c(
      apply(rowCumsums(ucw_qp_rewards), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(ts_whittle_rewards), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(random_rewards_sub), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_all_time_avg_rewards), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_current_time_avg_rewards), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_current_time_rewards), 2, sd)/sqrt(total_iterations)
    ),
    ses_timeavg = c(
      apply(rowCumsums(ucw_qp_rewards)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(ts_whittle_rewards)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(random_rewards_sub)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_all_time_avg_rewards)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_current_time_avg_rewards)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations),
      apply(rowCumsums(optimal_current_time_rewards)/matrix(rep(seq(1, t), total_iterations), c(total_iterations, t), byrow = TRUE), 2, sd)/sqrt(total_iterations)
    )
  )
  
  print('whittle df made')
  final_rewards_df = rbind(rewards_df, whittle_rewards)
  print('final rewards df made')

  final_rewards_df['Method_Type'] =  as.factor(final_rewards_df$Method_Type)
  final_rewards_df$Method_Type = factor(final_rewards_df$Method_Type, levels = c("Greedy Oracle" ,"Random", "Whittle Oracles", "UCW", 'TS', "BCoR"))
  final_rewards_df$Method = factor(final_rewards_df$Method, levels = c("greedy_oracle" ,"all_time_avg", 'current_time_avg', "current_time",
                                                                       'ts_la', 'bcor_Whittle',
                                                                       "ts_whittle", 'tsgreedy',
                                                                       "qp",
                                                                       'random'))
  
  

  all_methods_error_bars = ggplot(data= final_rewards_df, aes(x=time, y = cum_avg_reward, linetype = Method, color = Method_Type, fill=Method_Type)) + geom_line(size = 0.75) +
    geom_ribbon(aes(ymin = cum_avg_reward - ses_cumsum , ymax = cum_avg_reward + ses_cumsum ),   alpha=0.2, colour = NA) + theme_minimal() +
    scale_linetype_manual(values=c(1, 1, 2, 3, 1,2, 1, 2, 1, 1), labels=c( 'Greedy Oracle',
                                                                           'Whittle Oracle: Time Avg','Whittle Oracle: Current Time Avg', 'Whittle Oracle: Current Time',
                                                                           'BCoR-Greedy', 'BCoR-Whittle',
                                                                           "TS-Whittle", 'TS-Greedy',
                                                                           'UCW-Penalty',
                                                                           'Random')) +
    xlab('Time') +
    ylab('Cumulative Average Reward')
  
  
 
  all_methods_time_avg_reward = ggplot(data=final_rewards_df, aes(x=time, y = time_avg_reward,  linetype = Method, color = Method_Type, fill=Method_Type)) + geom_line(size = 0.75) +
    geom_ribbon(aes( ymax=time_avg_reward + 2*ses_timeavg, ymin = time_avg_reward- 2*ses_timeavg),
                alpha=0.2, colour = NA)+
    theme_minimal() +
    scale_linetype_manual(values=c(1, 1, 2, 3, 1,2, 1, 2, 1, 2), labels=c( 'Greedy Oracle',
                                                                           'Whittle Oracle: Time Avg','Whittle Oracle: Current Time Avg', 'Whittle Oracle: Current Time',
                                                                           'BCoR-Greedy','BCoR-Whittle',
                                                                           "TS-Whittle", 'TS-Greedy',
                                                                           'UCW-Penalty',
                                                                           'Random')) +
    xlab('Time') +
    ylab('Time Average Reward') + theme(axis.title =element_text(size=20), axis.text = element_text(size = 17), legend.title = element_text(size = 17), legend.text =  element_text(size = 15)) 
  
  
  write.csv(final_rewards_df, file = paste0(final_path, '/final_rewards_df.csv'), row.names = FALSE)
  
  
  plot_path = paste0(base_path, '/plots/')
  
  plot_str = paste0('_N_',N,'_T_',t,'_B_',budget,'_k_',k,'_d_',d,'_sigmaa_',sigma_a,'_sigmab_',sigma_b,'_mubeta_',mu_beta,'_sigmamubeta_',sigma_mu_beta,'_muB_',
                    mu_B,'_sigmaX_',sigma_X,'_sigmaB_',sigma_B,'_sigma_',sigma,'_tau_',tau)
  
  ggsave(filename = paste0('random_sub_cumulative_reward_errors' ,plot_str,'.png'),
         plot = all_methods_error_bars, path = plot_path, bg = 'white')
  
  # 
  ggsave(filename = paste0('random_sub_time_avg_reward', plot_str,'.png'),
         plot =   all_methods_time_avg_reward , path = plot_path, bg = 'white')
 
} else {
  print('waiting for full iterations')
}

