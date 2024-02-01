import numpy as np
import pandas as pd
import scipy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.special import logit, expit
import math

class Memoizer:
    """ improve performance of memoizing solutions (to QP and WI value iteration) """
    def __init__(self, method):
        self.method = method
        self.solved_p_vals = {}

    def to_key(self, input1, input2):
        """ convert inputs to a key

        QP: inputs: LCB and UCB transition probabilities
        UCB and extreme: inputs - estimated transition probabilities and initial state s0 """
        if self.method in ['lcb_ucb', 'QP', 'QP-min']:
            lcb, ucb = input1, input2
            p_key = (np.round(lcb, 4).tobytes(), np.round(ucb, 4).tobytes())
        elif self.method in ['tswhittle','p_s', 'optimal', 'UCB', 'extreme', 'ucw_value']:
            transitions, state = input1, input2
            p_key = (np.round(transitions, 4).tobytes(), state)
        elif self.method in ['lcb_ucb_s_lamb']:
            lcb, ucb = input1
            s, lamb_val = input2
            p_key = (np.round(lcb, 4).tobytes(), np.round(ucb, 4).tobytes(), s, lamb_val)
        else:
            raise Exception(f'method {self.method} not implemented')

        return p_key

    def check_set(self, input1, input2):
        p_key = self.to_key(input1, input2)
        if p_key in self.solved_p_vals:
            return self.solved_p_vals[p_key]
        return -1

    def add_set(self, input1, input2, wi):
        p_key = self.to_key(input1, input2)
        self.solved_p_vals[p_key] = wi

def get_valid_lcb_ucb(arm_p_lcb, arm_p_ucb):
    n_states, n_actions = arm_p_lcb.shape

    # enforce validity constraints
    assert n_actions == 2  # these checks only valid for two-action
    for s in range(n_states):
        # always better to act
        if arm_p_ucb[s, 0] > arm_p_ucb[s, 1]:  # move passive UCB down
            arm_p_ucb[s, 0] = arm_p_ucb[s, 1]
        if arm_p_lcb[s, 1] < arm_p_lcb[s, 1]:  # move active LCB up
            arm_p_lcb[s, 1] = arm_p_lcb[s, 1]

    assert n_states == 2  # these checks only valid for two-state
    for a in range(n_actions):
        # always better to start in good state
        if arm_p_ucb[0, a] > arm_p_ucb[1, a]:  # move bad-state UCB down
            arm_p_ucb[0, a] = arm_p_ucb[1, a]
        if arm_p_lcb[1, a] < arm_p_lcb[0, a]:  # move good-state LCB up
            arm_p_lcb[1, a] = arm_p_lcb[0, a]

    # these above corrections may lead to LCB being higher than UCBs... so make the UCB the optimistic option
    if arm_p_ucb[0, 0] < arm_p_lcb[0, 0]:
        print(f'ISSUE 00!! lcb {arm_p_lcb[0, 0]:.4f} ucb {arm_p_ucb[0, 0]:.4f}')
        arm_p_ucb[0, 0] = arm_p_lcb[0, 0] # p_ucb[i, 0, 0]
    if arm_p_ucb[0, 1] < arm_p_lcb[0, 1]:
        print(f'ISSUE 01!! lcb {arm_p_lcb[0, 1]:.4f} ucb {arm_p_ucb[0, 1]:.4f}')
        arm_p_ucb[0, 1] = arm_p_lcb[0, 1] # p_ucb[i, 0, 1]
    if arm_p_ucb[1, 0] < arm_p_lcb[1, 0]:
        print(f'ISSUE 10!! lcb {arm_p_lcb[1, 0]:.4f} ucb {arm_p_ucb[1, 0]:.4f}')
        arm_p_ucb[1, 0] = arm_p_lcb[1, 0] # p_ucb[i, 1, 0]
    if arm_p_ucb[1, 1] < arm_p_lcb[1, 1]:
        print(f'ISSUE 11!! lcb {arm_p_lcb[1, 1]:.4f} ucb {arm_p_ucb[1, 1]:.4f}')
        arm_p_ucb[1, 1] = arm_p_lcb[1, 1] # p_ucb[i, 1, 1]

    return arm_p_lcb, arm_p_ucb


def get_ucb_conf(cum_prob, n_pulls, t, alpha, episode_count, delta=1e-3):
    """ calculate transition probability estimates """
    n_arms, n_states, n_actions = n_pulls.shape

    with np.errstate(divide='ignore'):
        n_pulls_at_least_1 = np.copy(n_pulls)
        n_pulls_at_least_1[n_pulls == 0] = 1
        est_p               = cum_prob / n_pulls_at_least_1
        est_p[n_pulls == 0] = 1 / n_states  # where division by 0

        conf_p = np.sqrt( 2 * n_states * np.log( 2 * n_states * n_actions * n_arms * ((episode_count+1)**4 / delta) ) / n_pulls_at_least_1 )
        conf_p[n_pulls == 0] = 1
        conf_p[conf_p > 1]   = 1  # keep within valid range

    # if VERBOSE: print('conf', np.round(conf_p.flatten(), 2))
    # if VERBOSE: print('est p', np.round(est_p.flatten(), 2))

    return est_p, conf_p


def get_glm_conf(t, X_sa, X, which_women_sa_t, next_state_sa, n_pulls, prev_est_p, prev_conf_p, alpha = 0.05):
    """ calculate new transition probability estimates, update est_p and conf_p estimates from previous iteration
    t: timestep
    X_sa: dictionary of context np.array matrices for each (s, a) pair up to time t
    X: N x k context matrix for the women (no repeated rows)
    next_states_sa: dictionary of 1-d np.array vectors of next states for each (s, a) pair up to time t
    which_women_sa_t: dictionary of 1-d np.arrays containing the indices of which women (1, ..., N) are associated with each (s, a) at time t
    alpha: significance level of CI
    n_pulls: number of pulls for each (s, a) across all N women
    """
    n_arms, n_states, n_actions = n_pulls.shape
    #initialize est_p and conf_p, set all est to 0.5 and all conf_p to [0, 1] for now
    est_p = prev_est_p.copy()#np.full((n_arms, n_states, n_actions),  1/n_states)
    conf_p = prev_conf_p.copy()# np.zeros((n_arms, n_states, n_actions, 2))
    if t==1:
        return est_p, conf_p
    #for each (s, a)
    for s in np.arange(n_states):
        for a in np.arange(n_actions):
            #get context matrix associated with (s, a)
            X_sa_train = X_sa[(s, a)]
            #check that it isn't empty
            if (len(X_sa_train) > 0)  & (X_sa_train.shape[0] > X_sa_train.shape[1]):
                X_sa_train_df = pd.DataFrame(X_sa_train)
                #add intercept because silly statsmodel doesn't do this automatically
                X_sa_train_df = sm.add_constant(X_sa_train_df)
                #get response vector of next states associated with each (s, a)
                y_sa_train = next_state_sa[(s, a)]
                #run logistic regression
                try:
                    log_reg = sm.GLM(y_sa_train, X_sa_train_df, family=sm.families.Binomial())
                    result = log_reg.fit()
                    contexts_sa_t = X[which_women_sa_t[(s, a)], :]
                    #get estimated p for each woman at (s, a) at timestep t-1
                    est_p_sa = result.predict(sm.add_constant(contexts_sa_t))
                    #get predictions on link scale
                    fit_link = logit(est_p_sa)
                    #estimate asymptotic var-cov matrix on link scale
                    cov_mat_est_link_scale = np.matmul(np.matmul(sm.add_constant(contexts_sa_t), np.array(result.cov_params())),  sm.add_constant(contexts_sa_t).transpose())
                    ses_link= np.sqrt(cov_mat_est_link_scale.diagonal())
                    #get CIs of predictions on response scale
                    conf_p_lower = expit(fit_link - scipy.stats.norm.ppf(1-(alpha/2))*ses_link)
                    conf_p_upper = expit(fit_link + scipy.stats.norm.ppf(1-(alpha/2))*ses_link)
                    #map back to est_p based on indices of women in X_sa_train
                    est_p[which_women_sa_t[(s, a)], s, a] = est_p_sa
                    conf_p[which_women_sa_t[(s, a)], s, a, 0] = conf_p_lower
                    conf_p[which_women_sa_t[(s, a)], s, a, 1] = conf_p_upper
                except:
                    #if an error occurs, just exit the loop and return the previous estimates
                    print('An error occured, likely perfect separation')
            else:
                print(str((s, a))+': not enough rows')
                #if no data in matrix, then no women were in state s and received action a
                #cannot estimate W_sa, so just make confident intervals [0, 1]
#                 est_p[which_women_sa_t[(s, a)], s, a] = 0.5
#                 conf_p[which_women_sa_t[(s, a)], s, a, :] = np.array([0, 1])
    #conf_p is length of interval; may need to adjust uc_whittle step functions to take in an interval instead
    #can check if conf_p is a number or a list, then adjust logic to accomodate in each step function
    return est_p, conf_p
