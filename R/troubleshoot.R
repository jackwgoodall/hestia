
################# FUNCTIONS ###################################

softmax <- function(x) {
  return(exp(x)/sum(exp(x)))
}

get_diagonal_element <- function(m, i){
  out = 1
  for(j in 1:nrow(m)) {
    if(j != i) {
      out = out - m[j, i]  
    }
  }
  return(out)
}

normalize_cols <- function(m){
  out <- matrix(nrow = nrow(m), ncol = ncol(m))
  
  for(i in 1:ncol(m)) {
    out[,i] = m[,i]/sum(m[,1])
  }
  return(out)
}

################# DATA ##########################################
{n_states <- dat_stan$n_states
trans <- dat_stan$trans
n_inf_states <- dat_stan$n_inf_states
inf_states <- dat_stan$inf_states
n_trans_fit <- dat_stan$n_trans_fit
param_index <- dat_stan$param_index
trans_index <- dat_stan$trans_index
source_states <- dat_stan$source_states
multiplier <- dat_stan$multiplier
n_params <- dat_stan$n_params
n_hh <- dat_stan$n_hh
hh_size <- dat_stan$hh_size
n_obs <- dat_stan$n_obs
n_obs_type <- dat_stan$n_obs_type
n_unique_obs <- dat_stan$n_unique_obs
y <- dat_stan$y
part_id <- dat_stan$part_id
t_day <- dat_stan$t_day
obs_per_hh <- dat_stan$obs_per_hh
hh_start_ind <- dat_stan$hh_start_ind
hh_end_ind <- dat_stan$hh_end_ind
hh_tmin <- dat_stan$hh_tmin
hh_tmax <- dat_stan$hh_tmax
obs_prob <- dat_stan$obs_prob
init_probs <- dat_stan$init_probs
epsilon <- dat_stan$epsilon}

################### PARAMETERS ################################
params <- c()
beta_eh <- logit(0.01)
beta_ih <- logit(0.05)

############## TRANSFORMED PARAMETERS ##############################
llik <- matrix(0, nrow = sum(hh_size), ncol = max(hh_tmax)-min(hh_tmin) + 1) # lik contribution for enrolled per participant and time
llik_final <- numeric(n_hh) # sum of logalpha for final timestep
logalpha <- matrix(nrow = sum(hh_size)*n_states, ncol = max(hh_tmax)-min(hh_tmin) + 1) # log forward probability
alpha <- matrix(nrow = sum(hh_size)*n_states, ncol = max(hh_tmax)-min(hh_tmin) + 1) # forward prob, normalized
obs <- array(dim = c(n_obs_type, sum(hh_size), n_states)) # observation component for enrolled memebrs, set to 1 if no observation for this time step
trans_temp <- matrix(nrow = n_states, ncol = n_states) 

ih_prob <- inv_logit(beta_ih)
eh_prob <- inv_logit(beta_eh)
trans_temp <- trans

for(h in 1:n_hh) { # Loop through households
  
  i_rows <- matrix(0,nrow = hh_size[h], ncol = n_states)
  
  if(h == 1) {
    last_lik = 0
  } else {
    last_lik = sum(hh_size[1:(h-1)])
  }
  
  # Subset to data only for the given HH
  y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),]
  t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])]
  part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])]
  
  index = 1
  
  for(i in 1:hh_size[h]) { # Loop through participants for t = 1
    ref <- numeric(n_states)
    
    for(k in 1:n_states) {ref[k] = n_states*last_lik+n_states*(i-1)+k} 
    
    obs_switch = 0 
    
    if(t_day_hh[index] == 1) {
      if(part_id_hh[index] == i) {
        obs_switch = 1
      }
    }
    
    if(obs_switch == 1) {
      for(k in 1:n_obs_type) {
        if(y_hh[index, k] != -1) {
          obs[k, last_lik+i, ] = obs_prob[k, y_hh[index, k],]
        } else {
          obs[k, last_lik+i, ] = rep_row_vector(1, n_states)
        }
      }
    } else {
      obs[,last_lik+i,] = 1
    }
    
    if(obs_switch == 1) {
      index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1)  
    }
    
    # Fill in starting probability for SIR states
    logalpha[ref, 1] = log(init_probs)
    for(k in 1:n_obs_type) {
      logalpha[ref, 1] = logalpha[ref, 1] + log(obs[k,last_lik+i,])
    }
    for(s in inf_states) {
      i_rows[i, s] = n_states*last_lik+n_states*(i-1)+s  
    }
    
    llik[last_lik + i, 1] = log(sum(exp(logalpha[ref,1])))
    
    # normalize and convert to the probability scale
    alpha[ref, 1] = softmax(logalpha[ref,1])
    
  } # End participant loop for t = 1
  
  for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {
    for(p in 1:hh_size[h]) {
      no_inf_prob <- numeric(n_states) # probability of avoiding all infections
      no_hh_inf_prob <- matrix(nrow = hh_size[h], ncol = n_states) # probability of avoiding infection from each HH member
      ref <- numeric(n_states)
      
      for(k in 1:n_states) {ref[k] = n_states*last_lik+n_states*(p-1)+k}
      logalpha_temp = logalpha[ref,tt-1]
      
      obs_switch = 0 
      
      if(t_day_hh[index] == tt) {
        if(part_id_hh[index] == p) {
          obs_switch = 1
        }
      }
      
      if(obs_switch == 1) {
        for(k in 1:n_obs_type) {
          if(y_hh[index, k] != -1) {
            obs[k,last_lik+p,] = obs_prob[k, y_hh[index, k], ]
          } else {
            obs[k,last_lik+p,] = 1
          }
        }
      } else {
        obs[,last_lik+p,] = 1
      }
      
      if(obs_switch == 1) {
        index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1)  
      }
      
      for(s in 1:n_states) {
        if(s %in% inf_states) {
          no_hh_inf_prob[,s] = alpha[i_rows[, s], tt-1]*(1-ih_prob) + (1 - alpha[i_rows[,s], tt-1]) # Pr of avoiding infection from each household member
          no_hh_inf_prob[p, s] = 1 # Particpant can't infect themselves  
          } else {
            no_hh_inf_prob[,s] = 1 
          }
          no_inf_prob[s] = prod(no_hh_inf_prob[,s]) # Probability of avoiding infection from all household members
        }
        
        # fill in tranistions that are being fit
        for(m in 1:n_trans_fit) {
          if(sum(source_states[m,]) == 0) {
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = params[param_index[m]]*multiplier[m]
          } else {
            no_inf = 1
            for(s in 1:n_states) {
              if(source_states[m,s] == 1) {
                no_inf = no_inf*no_inf_prob[s]
              }
            }
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = (1-no_inf)*(1-eh_prob)*multiplier[m]
          }
        }

        # fill in diagonals (columns must sum to one)
        for(i in 1:ncol(trans_temp)) {
          trans_temp[i,i] = get_diagonal_element(trans_temp, i)
        }
        
        # normalize
        trans_temp = normalize_cols(trans_temp)
        
        # Compute the probability of each epidemiological state
        logalpha[ref, tt] = log(trans_temp %*% exp(logalpha_temp))
        for(k in 1:n_obs_type) {
          logalpha[ref, tt] = logalpha[ref, tt] + log(obs[k,last_lik+p,])
        }
        
        # normalize and convert to probability scale
        alpha[ref, tt] = softmax(logalpha[ref,tt])
        
        llik[last_lik + p, tt] = log(sum(exp(logalpha[ref,tt])))
        
      } # end participant loop - update logalpha with observation probability
      
      if(tt == (hh_tmax[h] - hh_tmin[h] + 1)) {
        llik_final[h] = sum(llik[(last_lik+1):(last_lik+hh_size[h]),tt]) 
      }
        
    } # end time loop

} # End household loop



