library(invgamma)
library(dplyr)
library(dlnm)
#library(fitdistrplus)
# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('E', 'N', 't', 'budget', 'k', 'd', 'sigma_a', 'sigma_b', 'mu_beta', 'sigma_mu_beta', 'mu_B', 'sigma_X', 'sigma_B', 'sigma', 'tau')

if (length (args) == 0) {
  print("No arguments supplied.")
  ## supply default values
  E = 1
  N = 400 
  t = 50 
  budget = 3
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

#set a unique random seed
random_seed = 122+E
set.seed(random_seed)

#path to write data files, change to desired file path
csv = paste0("data/N_",N,'_T_',t,'_B_',budget,'_k_',k,'_d_',d,'/sigmaa_',sigma_a,'_sigmab_',sigma_b,'_mubeta_',mu_beta, '_sigmamubeta_',sigma_mu_beta,'_muB_',
             mu_B,'_sigmaX_',sigma_X,'_sigmaB_',sigma_B,'_sigma_',sigma,'_tau_',tau)



gen_X = function(N, k){
  #generate covariates
  age = round(rnorm(N, 22, 2))
  X = matrix(c(rnorm(N*(k-2), 0, 1), #some continuous predictors
               rbinom(N, size=1, prob = 0.5), #some binary predictor
               (age-mean(age))/sd(age) #mean center age variable
  ), ncol = k, byrow = FALSE)
  return(X)
}

#generate spline basis for time effect
gen_B = function(t, d){
  B = ps(seq(1, t, 1), df=d, knots=NULL, degree=3)
  return(B)
}

X = gen_X(N, k)
B = gen_B(t, d)

generate_transitions = function(X, B, N, k, d, t, sigma_b0, sigma_b1, sigma, tau, mu_beta, tau_mu_beta, tau_X, mu_B, tau_B){
  n=1
  b0 = rnorm(n, mean=0, sd = sigma_b0)
  b1 = rnorm(n, mean=0, sd = sigma_b1)
  mu_beta_1 = rnorm(k, mean=mu_beta, sd = tau_mu_beta)
  tau_alpha_2_00 = rinvgamma(n, shape = tau, rate = sigma) 
  tau_alpha_2_10 = rinvgamma(n, shape = tau, rate = sigma) 
  tau_alpha_2_01 = rinvgamma(n, shape = tau, rate = sigma) 
  tau_alpha_2_11 = rinvgamma(n, shape = tau, rate = sigma) 
  #coefficients corresponding to covariates
  beta_00 = mu_beta_1 + rnorm(k, 0, tau_X) 
  beta_10 = mu_beta_1 + rnorm(k, 0, tau_X) 
  beta_01 = mu_beta_1 + rnorm(k, 0, tau_X) 
  beta_11 = mu_beta_1 + rnorm(k, 0, tau_X)
  #time varying effects
  eta_00 = rnorm(d, mu_B, tau_B) 
  eta_10 = rnorm(d, mu_B, tau_B)
  eta_01 = rnorm(d, mu_B, tau_B)
  eta_11 = rnorm(d, mu_B, tau_B)
  #random effects when a=0
  alpha_00 = rnorm(N, mean=0, sd=sqrt(tau_alpha_2_00))
  alpha_10 =rnorm(N, mean=0, sd=sqrt(tau_alpha_2_10))
  #generate the N x t matrix where each ij-entry is the transition probability
  #of woman i to enter s'=1 from s=0, a=0 at time j
  transitions_s0_a0 = sapply(seq(1, t, 1), function(timestep){D_t = cbind(X, t(replicate(N, B[timestep, ])), diag(N))
  p_t_s0_a0 = pnorm(D_t %*% c(beta_00, eta_00, alpha_00))
  return(p_t_s0_a0)})
  transitions_s1_a0 = sapply(seq(1, t, 1), function(timestep){D_t = cbind(X, t(replicate(N, B[timestep, ])), diag(N))
  p_t_s1_a0 = pnorm(D_t %*% c(beta_10, eta_10, alpha_10))
  return(p_t_s1_a0)})
  #random effects for when a=1
  delta_0 = alpha_00
  delta_1 = alpha_10
  alpha_01 = rnorm(N, mean=b0*delta_0 + b1*delta_1, sd=sqrt(tau_alpha_2_01))
  alpha_11 =rnorm(N, mean=b0*delta_0 + b1*delta_1, sd=sqrt(tau_alpha_2_11))
  #transition probabilities for when a=1
  transitions_s0_a1 = sapply(seq(1, t, 1), function(timestep){D_t = cbind(X, t(replicate(N, B[timestep, ])), diag(N))
  p_t_s0_a1 = pnorm(D_t %*% c(beta_01, eta_01, alpha_01))
  return(p_t_s0_a1)})
  transitions_s1_a1 = sapply(seq(1, t, 1), function(timestep){D_t = cbind(X, t(replicate(N, B[timestep, ])), diag(N))
  p_t_s1_a1 = pnorm(D_t %*% c(beta_11, eta_11, alpha_11))
  return(p_t_s1_a1)})
  return (list(transitions_s0_a0, transitions_s1_a0, transitions_s0_a1, transitions_s1_a1, 
               list('b0'=b0, 'b1'=b1, 'alpha_00' = alpha_00, 'alpha_01' = alpha_01, 'alpha_10' = alpha_10, 'alpha_11' = alpha_11,
                    'tau_alpha_2_00'=tau_alpha_2_00, 'tau_alpha_2_01'=tau_alpha_2_01, 'tau_alpha_2_10'=tau_alpha_2_10, 'tau_alpha_2_11'=tau_alpha_2_11,
                    'eta_00'=eta_00, 'eta_01'=eta_01,'eta_10'=eta_10, 'eta_11'=eta_11, 'mu_beta' = mu_beta_1,
                    'beta_00'=beta_00, 'beta_01'=beta_01,'beta_10'=beta_10, 'beta_11'=beta_11
               ))#also output the true parameter values
  )
}

true_transition_list = generate_transitions(X, B, N, k, d, t, sigma_a, sigma_b, sigma, tau, mu_beta, sigma_mu_beta, sigma_X, mu_B, sigma_B)
transitions_s0_a0 = true_transition_list[[1]]
transitions_s1_a0 = true_transition_list[[2]]
transitions_s0_a1 = true_transition_list[[3]]
transitions_s1_a1 = true_transition_list[[4]]
true_params = true_transition_list[[5]]

#generate data frame of state, action, and transition data
true_transitions_df = data.frame(person = rep(seq(1, N), 4*t),
                                 time = rep(rep(seq(1, t, 1), each = N), 4),
                                 state = rep(c(0, 1, 0, 1), each = N*t),
                                 action = rep(c(0, 0, 1, 1), each = N*t),
                                 transitions = c(transitions_s0_a0[1:(N*t)], transitions_s1_a0[1:(N*t)],
                                                 transitions_s0_a1[1:(N*t)],  transitions_s1_a1[1:(N*t)])
)

run_ts_greedy = function(true_transitions_df, N, t, budget){
  #pre-allocate final data matrix
  df <- matrix(data=NA,nrow=(N*t),ncol=5)
  colnames(df) = c('person', 'time', 'state', 'action', 'y')
  #initialize df to keep track of ts draws, initialize draws using runif
  ts_draws = cbind(true_transitions_df[true_transitions_df$time==1,c('person', 'state', 'action')],rbeta(N*4,1,1))
  colnames(ts_draws)[4] = 'est_transitions'
  #wrap this in a loop over time until end of time horizon
  for(current_time in seq(1, t)){
    print(current_time)
    #get current state of woman (next state from previous time step if t>1, else, randomly initialize)
    if(current_time == 1){
      #initialize with a draw from the prior
      s = round(runif(N))
    }else{
      s = next_state
    }
    
    if(current_time == 1){
      #just select a random set of woman to give action
      pull_idx = sample(seq(1:N), budget, replace = FALSE)
    }else{
      #generate draw from the posterior given data seen up until this point
      current_df = df[1:(N*(current_time-1)), ]
      update_table = data.frame(current_df) %>% group_by(person, state, action) %>% summarise(count = n(), heads = sum(y))
      updates = sapply(seq(1, N, 1), function(i){
        #for each person, state, action, count up number of observed "heads" aka engaging transitions and number of total times in that state-action pair
        x = c(0, 0, 0, 0)
        j = 1
        for(actions in c(0, 1)){
          for(states in c(0, 1)){
            #if I have seen this person's state actio transitions before
            check = update_table[update_table$person == i & update_table$state == states & update_table$action==actions,]
            #print(c(i, states, actions))
            if(dim(check)[1] == 0){
              #if I haven't seen that person's state-action transitions before, just use unif[0,1]
              x[j] =  rbeta(1,1,1)
              j = j+1
            }else{
              #update according to beta binomial model
              x[j] = rbeta(1, 1+check$heads, check$count-check$heads+1)
              j = j+1
            }
          }
        }
        return(x)})
      
      ts_draws$est_transitions = as.vector(t(updates))
      
      te = sapply(seq(1, N, 1), function(i){
        transition_est_a1 = filter(ts_draws, person == i, state==s[i], action==1)$est_transitions
        transition_est_a0 = filter(ts_draws, person == i, state==s[i], action==0)$est_transitions
        return(transition_est_a1 - transition_est_a0)
      })
      #pull women with budget highest treatment effects
      #get indices of women to pull
      pull_idx = sort(te, index.return = TRUE, decreasing= TRUE)$ix[1:budget]
    }
    a = rep(0, N)
    a[pull_idx] = 1
    #pull out correct current_time, state, action for each women in true_transitions_df and generate bernoulli draw from true transitions
    next_state = sapply(seq(1, N), function(i){
      underlying_transition = filter(true_transitions_df, person == i, state==s[i], action==a[i], time == current_time)$transitions
      return(rbinom(n = 1, size = 1, underlying_transition))
    })
    #add new transition data to growing data table for new t
    df[((N*current_time)- (N-1)):(N*current_time),] = matrix(c(seq(1, N), rep(current_time, N), s, a, next_state), nrow = N, byrow = FALSE)
    #save intermediate data files
    write.csv(df, file = paste0(csv, '/iteration_rewards/ts_greedy_run_E_',E,'.csv'), row.names = FALSE)
  }
  return(list(df, update_table))
}


ts_greedy_run = run_ts_greedy(true_transitions_df, N, t, budget)


tsgreedy_run_df = data.frame(ts_greedy_run[[1]])
update_table = ts_greedy_run[[2]]
tsgreedy_run_df['Method'] = 'tsgreedy'
tsgreedy_run_df['iteration'] = E

# #generate rewards df for this epoch
rewards_df = tsgreedy_run_df %>% group_by(Method, time) %>% summarize(reward = sum(y))%>% mutate(cum_reward = cumsum(reward)) %>%
  mutate(time_avg_reward = cumsum(reward)/time)
rewards_df['iteration'] = E
#
# #save full run df and rewards
write.csv(tsgreedy_run_df, file = paste0(csv,"/iteration_runs/ts_greedy_run_E_",E,".csv"),  row.names = FALSE)
write.csv(rewards_df, file = paste0(csv,"/iteration_rewards/ts_greedy_rewards_E_",E,".csv"),  row.names = FALSE)

