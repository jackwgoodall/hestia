data {
  
  // Model information
  int  n_states; // number of latent states
  matrix[n_states, n_states] trans; // transition probabilities, cols = starting state, rows = ending states
  int n_inf_states; // number of infectious states
  int inf_states[n_inf_states]; // infectious states
  int n_trans_fit; // number of transitions to fit
  int trans_index[n_trans_fit, 2]; // row/col indices of transition matrix corresponding to each parameter to be fit
  int inf_mult[n_trans_fit]; // foi multiplier for infection transistions
  int n_params; // number of additional (non-infection) parameters to fit
  
  // Household information
  int n_hh; // number of households
  int hh_size[n_hh]; // household size
    
  // Data 
  int n_obs; // number of observations
  int n_obs_type; // number of observation types
  int n_unique_obs; // number of unique outcomes for enrolled individuals
  int y[n_obs, n_obs_type]; // outcome vector, ordered by 1) household, then 2) time, then 3) individual
  int part_id[n_obs]; // particpants associated with the observation
  int t_day[n_obs]; // observation times for enrolled (day)
  int obs_per_hh[n_hh]; // total number of observations per HH for enrolled
  int hh_start_ind[n_hh]; // starting index for the HH
  int hh_end_ind[n_hh]; // ending index for the HH
  int hh_tmin[n_hh]; // minimum day for which there is a HH observation
  int hh_tmax[n_hh]; // maximum day for which there is a HH observation

  // Initial state and observation probabilities
  matrix[n_unique_obs, n_states] obs_prob[n_obs_type]; // observation process for SIR states
  vector[n_states] init_probs; // starting state probabilities
  
  real epsilon; // small number to avoid log(0) issues
}

parameters {
  real params[n_params];
  real beta_eh; // monthly intercepts for extra-household probabilities
  real beta_ih; // intra-household probability

}

transformed parameters {
  matrix[sum(hh_size), max(hh_tmax)-min(hh_tmin) + 1] llik; // lik contribution for enrolled per participant and time
  vector[n_hh] llik_final; // sum of logalpha for final timestep
  real ih_prob;
  real eh_prob;
  matrix[sum(hh_size)*n_states, max(hh_tmax)-min(hh_tmin) + 1] logalpha; // log forward probability
  matrix[sum(hh_size)*n_states, max(hh_tmax)-min(hh_tmin) + 1] alpha; // forward prob, normalized
  matrix[sum(hh_size), n_states] obs[n_obs_type]; // observation component for enrolled memebrs, set to 1 if no observation for this time step
  matrix[n_states, n_states] trans_temp;
  
  llik = rep_matrix(0, sum(hh_size), max(hh_tmax)-min(hh_tmin) + 1);
  ih_prob = inv_logit(beta_ih);
  eh_prob = inv_logit(beta_eh);
  trans_temp = trans;
  
  for(h in 1:n_hh) { // loop through household
    
    int y_hh[obs_per_hh[h], n_obs_type];
    int part_id_hh[obs_per_hh[h]];
    int t_day_hh[obs_per_hh[h]];
    int t_week_hh[obs_per_hh[h]];
    int index; // index for next observation
    int i_rows[hh_size[h], n_inf_states];// rows in alpha corresponding to infectious states
    int last_lik;
    int obs_switch; // indicator for whether there is an observation corresponding to this time step
    
    if(h == 1) {
      last_lik = 0;
    } else {
      last_lik = sum(hh_size[1:(h-1)]);
    }
    
    // subset to data only for the given HH
    y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),];
    t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
    part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];
    
    index = 1;
    
    { // START FORWARD ALGORITHM
    
    // fill first column of alpha using starting probabilities
    for(i in 1:hh_size[h]) {
      int ref[n_states];
    
      for(k in 1:n_states) {ref[k] = n_states*last_lik+n_states*(i-1)+k;} 
      
      obs_switch = 0; 
      
      if(t_day_hh[index] == 1) {
        if(part_id_hh[index] == i) {
          obs_switch = 1;
        }
      }
      
      if(obs_switch == 1) {
        for(k in 1:n_obs_type) {
          if(y_hh[index, k] != -1) {
            obs[k, last_lik+i, ] = obs_prob[k, y_hh[index, k],];
          } else {
            obs[k, last_lik+i, ] = rep_row_vector(1, n_states);
          }
        }
      } else {
        obs[,last_lik+i,] = rep_array(rep_row_vector(1, n_states), n_obs_type);
      }
        
      if(obs_switch == 1) {
       index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1);  
      }
      
      // Fill in starting probability for SIR states
      logalpha[ref, 1] = log(init_probs);
      for(k in 1:n_obs_type) {
        logalpha[ref, 1] = logalpha[ref, 1] + to_vector(log(obs[k,last_lik+i,]));  
      }
      for(s in 1:n_inf_states) {
        i_rows[i, s] = n_states*last_lik+n_states*(i-1)+inf_states[s];  
      }
      
      llik[last_lik + i, 1] = log_sum_exp(logalpha[ref,1]);
      
      // normalize and convert to the probability scale
      alpha[ref, 1] = softmax(logalpha[ref,1]);
    
    } // end participant loop - t=1, update logalpha with observation probability

    for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {
      
      for(p in 1:hh_size[h]) {
        real no_inf_prob; // probability of avoiding all infections
        matrix[hh_size[h], n_inf_states] no_hh_inf_prob; // probability of avoiding infection from each HH member
        int ref[n_states];
        vector[n_states] logalpha_temp; // log forward probability
        int offset;
        
        for(k in 1:n_states) {ref[k] = n_states*last_lik+n_states*(p-1)+k;}
        logalpha_temp = logalpha[ref,tt-1];
        
        obs_switch = 0; 
      
        if(t_day_hh[index] == tt) {
          if(part_id_hh[index] == p) {
            obs_switch = 1;
          }
        }
      
        if(obs_switch == 1) {
          for(k in 1:n_obs_type) {
            if(y_hh[index, k] != -1) {
              obs[k,last_lik+p,] = obs_prob[k, y_hh[index, k], ];
            } else {
              obs[k,last_lik+p,] = rep_row_vector(1, n_states);
            }
          }
        } else {
          obs[,last_lik+p,] = rep_array(rep_row_vector(1, n_states), n_obs_type);
        }
        
        if(obs_switch == 1) {
          index = min(index + 1, hh_end_ind[h]-hh_start_ind[h]+1);  
        }
        
        for(s in 1:n_inf_states) {
          no_hh_inf_prob[,s] = to_vector(alpha[i_rows[, s], tt-1])*(1-ih_prob) + (1 - to_vector(alpha[i_rows[,s], tt-1])); // Pr of avoiding infection from within the home 
          no_hh_inf_prob[p, s] = 1; // Particpant can't infect themselves  
        }
        no_inf_prob = (1-eh_prob)*prod(no_hh_inf_prob); // Pr avoiding all infections
        
        offset = 0;
        
        // fill in tranistions that are being fit
        for(m in 1:n_trans_fit) {
          if(inf_mult[m] == 0) {
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = params[m-offset];
          } else {
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = (1-no_inf_prob)*inf_mult[m];
            offset = offset + 1;
          }
        }

        // fill in diagonals
        
        // normalize
        
        
        // Compute the probability of each epidemiological state
        logalpha[ref, tt] = log(trans_temp*exp(logalpha_temp));
        for(k in 1:n_obs_type) {
          logalpha[ref, tt] = logalpha[ref, tt] + to_vector(log(obs[k,last_lik+p,]));
        }
        
        // normalize and convert to probability scale
        alpha[ref, tt] = softmax(logalpha[ref,tt]);
        
        llik[last_lik + p, tt] = log_sum_exp(logalpha[ref,tt]);
        
      } // end participant loop - update logalpha with observation probability
      
      if(tt == (hh_tmax[h] - hh_tmin[h] + 1)) {
        llik_final[h] = sum(llik[(last_lik+1):(last_lik+hh_size[h]),tt]); 
      }
        
    } // end time loop
    } // END FORWARD ALGORITHM
    
  } // end household loop
  
}

model {
  
  beta_eh ~ normal(-3,3);
  beta_ih ~ normal(-3,3);
  
  // Only increment by final alpha
  target += sum(llik_final);

}

