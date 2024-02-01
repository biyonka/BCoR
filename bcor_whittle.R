library(reticulate)
library(invgamma)
library(rstan)
library(dplyr)
library(dlnm)

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))# args is now a list of character vectors
# First check to see if arguments are passed.
# Then cycle through each element of the list and evaluate the expressions.
arguments = c('E', 'N', 't', 'budget', 'k', 'd', 'sigma_a', 'sigma_b', 'mu_beta', 'sigma_mu_beta', 'mu_B', 'sigma_X', 'sigma_B', 'sigma', 'tau')

stan_filepath = '{INSERT PATH TO FILE}/bcor_model.stan'
source_python('{INSERT PATH TO FILE}/compute_whittle.py') 

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
  sigma_b = 0.1
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
burnin = 500

#path to write data files
csv = paste0("data/N_",N,'_T_',t,'_B_',budget,'_k_',k,'_d_',d,'/sigmaa_',sigma_a,'_sigmab_',sigma_b,'_mubeta_',mu_beta, '_sigmamubeta_',sigma_mu_beta,'_muB_',
             mu_B,'_sigmaX_',sigma_X,'_sigmaB_',sigma_B,'_sigma_',sigma,'_tau_',tau)

gen_X = function(N, k){
  #generate covariates
  age = round(rnorm(N, 22, 2))
  X = matrix(c(rnorm(N*(k-2), 0, 1), #some continuous predictors
               rbinom(N, size=1, prob = 0.5), #some binary predictor
               (age-mean(age))/sd(age)
               # round(rnorm(N, 22, 2))-mean(round(rnorm(N, 22, 2))) #mean-centered and standardized simulated age variable
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
  eta_00 = rnorm(d, mu_B, tau_B) #assume tau_B^2 * I_d covariance matrix
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


start = proc.time()
model = stan_model(file = stan_filepath)
end = proc.time()
print(end - start)

#takes in data up to the current time, runs Stan model and outputs a single posterior draw
get_posterior_draw = function(df, stanmodel, burnin, N_sampl = 1, num_chains = 1, alg = 'NUTS', seed = random_seed, verb = FALSE){
  start = proc.time()
  stan_samples <- sampling(object = stanmodel, data = df, iter = burnin+N_sampl,
                           warmup=burnin, chains = num_chains, cores = 4,
                           init="0", algorithm=alg ,seed=random_seed, verbose=verb, refresh=0)
  end = proc.time()
  print(end - start)
  return(stan_samples)
}

discount = 0.9
min_chosen_subsidy = 0
whittle_threshold = 1e-4

bcor_Whittle = function(E, true_transitions_df, model, X, B, N, t, k, d,  budget, burn, seed=123){
  #pre-allocate final data matrix
  df <- matrix(data=NA,nrow=(N*t),ncol=5)
  pulls = matrix(data=NA, ncol=budget, nrow=t)
  colnames(df) = c('person', 'time', 'state', 'action', 'y')
  #wrap this in a loop over time until end of time horizon
  for(current_time in seq(1, t)){
    print(current_time)
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
      data = list(N=N, K=k, d=d, time_horizon = t,
                  t_c=current_time-1, #-1 since we only observed until current_time-1
                  y=matrix(current_df[,'y'], nrow = N, byrow = FALSE),
                  S=matrix(current_df[,'state'] + 1, nrow = N, byrow = FALSE), #need to add 1 for stan indexing
                  A=matrix(current_df[, 'action'] + 1, nrow = N, byrow = FALSE),
                  X=X,
                  B=B
      )
      est_parameters = get_posterior_draw(df=data, stanmodel=model, burnin=burn, N_sampl = 1, num_chains = 1, alg = 'NUTS', seed = seed, verb = FALSE)
      post_params = extract(est_parameters)
      # just take one sample for the transition est
      alpha_HMC = post_params$alpha[1,,,]#
      beta_HMC = post_params$beta[1,,,]#
      eta_HMC = post_params$eta[1,,,]#
      b0_HMC = post_params$b0[1]# 
      b1_HMC = post_params$b1[1]# 
      mu_beta_HMC = post_params$mu_beta[1,]
      tau_HMC = post_params$tau[1,,]
      params_t = list('alpha' = post_params$alpha,
                      'beta' =post_params$beta,
                      'eta' = post_params$eta,
                      'b0' = post_params$b0,
                      'b1' = post_params$b1,
                      'mu_beta' = post_params$mu_beta,
                      'tau' =post_params$tau)
  
      #get estimated transitions for current_time from estimated parameters for each woman and her state in s
      #estimated transitions for state=0, action = 1 for all woman
      compute_whittle = sapply(seq(1, N, 1), function(i){
        current_state = as.integer(s[i])
        transition_est_s0a1 = pnorm(X[i,]%*%beta_HMC[1, 2, ] + B[current_time,]%*%eta_HMC[1, 2, ] + alpha_HMC[1, 2,i])
        transition_est_s0a0 = pnorm(X[i,]%*%beta_HMC[1, 1, ] + B[current_time,]%*%eta_HMC[1, 1, ] + alpha_HMC[1, 1,i])
        transition_est_s1a1 = pnorm(X[i,]%*%beta_HMC[2, 2, ] + B[current_time,]%*%eta_HMC[2, 2, ] + alpha_HMC[2, 2,i])
        transition_est_s1a0 = pnorm(X[i,]%*%beta_HMC[2, 1, ] + B[current_time,]%*%eta_HMC[2, 1, ] + alpha_HMC[2, 1,i])
        test_transitions = array(c(transition_est_s0a0, transition_est_s1a0, transition_est_s0a1, transition_est_s1a1), dim=c(2, 2))
        wi = arm_compute_whittle(test_transitions, current_state, discount, min_chosen_subsidy, eps=whittle_threshold)
        return(wi)
      })
      
      #get top B whittle indices
      pull_idx = sort(compute_whittle, index.return = TRUE, decreasing= TRUE)$ix[1:budget]
 
    }
    pulls[current_time,] = pull_idx
    #action vector
    a = rep(0, N)
    a[pull_idx] = 1
    #pull out correct current_time, state, action for each women in true_transitions_df and generate bernoulli draw from true transitions
    next_state = sapply(seq(1, N), function(i){
      underlying_transition = filter(true_transitions_df, person == i, state==s[i], action==a[i], time == current_time)$transitions
      return(rbinom(n = 1, size = 1, underlying_transition))
    })
    #add new transition data to growing data table for new t
    df[((N*current_time)- (N-1)):(N*current_time),] = matrix(c(seq(1, N), rep(current_time, N), s, a, next_state), nrow = N, byrow = FALSE)

    write.csv(df, file = paste0(csv, '/iteration_rewards/bw_run_E_',E,'.csv'), row.names = FALSE)
  }
  return(df)
}


bcor_Whittle_run = bcor_Whittle(E, true_transitions_df,
                                          model = model,
                                          X, B, N, t, k, d, budget, burn = burnin, seed = random_seed)

warnings()

rmab_run_df = data.frame(bcor_Whittle_run)

rmab_run_df['Method'] = 'bcor_Whittle'
rmab_run_df['iteration'] = E

# #generate rewards df for this epoch
rewards_df = rmab_run_df %>% group_by(Method, time) %>% summarize(reward = sum(y))%>% mutate(cum_reward = cumsum(reward)) %>%
  mutate(time_avg_reward = cumsum(reward)/time)

# #save full run df and rewards
write.csv(rmab_run_df, file = paste0(csv,"/iteration_runs/bw_run_E_",E,".csv"),  row.names = FALSE)
write.csv(rewards_df, file = paste0(csv,"/iteration_rewards/bw_rewards_E_",E,".csv"),  row.names = FALSE)

