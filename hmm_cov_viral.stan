functions {

  // Return the diagonal element for column i so that each column sums to 1.
  // In this model, transition matrix columns represent the "from" state,
  // and rows represent the "to" state, so for each column i:
  //   m[i,i] = 1 - sum_{j != i} m[j,i]
  real get_diagonal_element(matrix m, int i){
    real out;
    out = 1;
    for(j in 1:rows(m)) {
      if(j != i) {
        out = out - m[j, i];
      }
    }
    return out;
  }

  // Stan helper comparable to R's `%in%` for integer arrays.
  // Returns 1 if `pos` appears in `pos_var`, else 0.
  int is_in(int pos, array[] int pos_var) {
    int pos_match;
    array[size(pos_var)] int all_matches;

    for (p in 1:(size(pos_var))) {
      all_matches[p] = (pos_var[p] == pos);
    }

    if(sum(all_matches) > 0) {
      pos_match = 1;
      return pos_match;
    } else {
      pos_match = 0;
      return pos_match;
    }
  }

  // Normalize each column so the column sums to 1.
  // Used after transition edits to ensure valid probabilities.
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;

    for(i in 1:cols(m)) {
      out[,i] = m[,i] / sum(m[,i]);
    }
    return out;
  }

  // Replace exact zeros with epsilon to avoid log(0) and numerical issues.
  matrix replace_zeroes(matrix m, real epsilon) {
    matrix[rows(m), cols(m)] out;
    out = m;
    for(i in 1:rows(m)) {
      for(j in 1:cols(m)) {
        if(m[i,j] == 0) {
          out[i,j] = epsilon;
        }
      }
    }
    return(out);
  }

}


data {

  // -------------------------
  // Core HMM / transition data
  // -------------------------
  int n_states; // number of latent epidemiological states
  // Transition matrix template: columns = state at t-1, rows = state at t.
  matrix[n_states, n_states] trans;

  int n_inf_states; // number of states considered infectious
  array[n_inf_states] int inf_states; // indices of infectious states

  // Sparse description of transition entries that are estimated instead of fixed.
  int n_trans_fit;
  // Which parameter controls each fitted transition entry.
  // 0 indicates an infection transition computed from force-of-infection terms.
  array[n_trans_fit] int param_index;
  // [row, col] location in transition matrix for each fitted transition.
  array[n_trans_fit, 2] int trans_index;
  // For infection transitions, marks which source infectious states contribute.
  // All 0 means this is not infection-driven.
  array[n_trans_fit, n_states] int source_states;

  int n_params; // number of non-infection transition parameters

  // -------------------------
  // Transition multipliers
  // -------------------------
  // Multipliers allow splitting one baseline transition into multiple branches.
  matrix[n_states, n_states] transition_multiplier;
  int n_mult_fit; // number of multiplier entries that are estimated
  int n_mult_params; // number of unique multiplier parameters
  // Index map for multiplier parameters.
  // Positive value: use mult_params[index]
  // Negative value: use (1 - mult_params[abs(index)]) style complement.
  array[n_mult_fit] int mult_param_index;
  // [row, col] locations where multipliers are applied.
  array[n_mult_fit, 2] int mult_index;

  // -------------------------
  // Household structure
  // -------------------------
  int n_hh; // number of households
  array[n_hh] int hh_size; // members per household

  // -------------------------
  // Observation data
  // -------------------------
  int n_obs; // total observation records (all households/times/participants)
  int n_obs_type; // number of observation channels (e.g., symptoms, tests)
  int n_unique_obs; // number of unique outcomes per observation type

  // y is ordered by household, then time, then participant.
  // y[r, k] = observed category for row r and observation type k.
  // Convention in this model: y == -1 means missing for that channel.
  array[n_obs, n_obs_type] int y;

  array[n_obs] int part_id; // participant id within household for each row
  array[n_obs] int t_day; // day index for each row

  // Household-specific indexing into the long observation arrays.
  array[n_hh] int obs_per_hh; // number of rows belonging to each household
  array[n_hh] int hh_start_ind; // first row index for each household
  array[n_hh] int hh_end_ind; // last row index for each household
  array[n_hh] int hh_tmin; // first modeled day for each household
  array[n_hh] int hh_tmax; // last modeled day for each household

  // -------------------------
  // Covariates for infection pressure
  // -------------------------
  int k_ih; // number of intra-household covariates
  // One row per participant in the full dataset (across all households).
  matrix[sum(hh_size), k_ih] x_ih;

  int k_eh; // number of extra-household covariates
  matrix[sum(hh_size), k_eh] x_eh;

  // -------------------------
  // HMM initial/observation model
  // -------------------------
  // obs_prob[k] is a matrix with rows = observed category, cols = latent state,
  // so obs_prob[k][obs_value, state] = P(observed value | latent state).
  array[n_obs_type] matrix[n_unique_obs, n_states] obs_prob;
  vector[n_states] init_probs; // prior state probabilities at first modeled day

  real epsilon; // tiny floor for numerical stability
  // Number of distinct intra-household infection probability logits.
  // Either 1 shared across infectious states, or one per infectious state.
  int n_inf_prob;
}

parameters {
  // Unconstrained parameters transformed via inv_logit to (0,1).
  array[n_params] real logit_params;
  array[n_mult_params] real logit_mult_params;

  // Logistic regression coefficients for infection pressure.
  vector[k_eh] beta_eh; // extra-household effects
  vector[k_ih] beta_ih; // intra-household effects
  real beta0_eh; // extra-household intercept
  array[n_inf_prob] real beta0_ih; // intra-household intercept(s)

}

transformed parameters {
  // Household-level contribution to log-likelihood.
  // Each entry is sum over members of log-sum-exp(alpha) at final time.
  vector[n_hh] llik_final;

  // Participant-level probabilities from logistic models.
  // ih_prob[, c] = probability of transmission from infectious-state group c.
  matrix[sum(hh_size), n_inf_prob] ih_prob;
  vector[sum(hh_size)] eh_prob; // extra-household infection probability

  // Stores per-person log forward probabilities over time.
  // Rows are grouped by participant/state (household by household).
  matrix[sum(hh_size) * n_states, max(hh_tmax) - min(hh_tmin) + 1] logalpha;

  // Working transition matrix updated each step.
  matrix[n_states, n_states] trans_temp;

  // Probability-scale transition parameters.
  array[n_params] real params;
  array[n_mult_params] real mult_params;

  // Map unconstrained real values to (0,1).
  params = inv_logit(logit_params);
  mult_params = inv_logit(logit_mult_params);

  // Intra-household infection probabilities (possibly one column or many).
  for(i in 1:n_inf_prob) {
    ih_prob[,i] = inv_logit(beta0_ih[i] + x_ih * beta_ih);
  }

  // Extra-household infection probability per participant.
  eh_prob = inv_logit(beta0_eh + x_eh * beta_eh);

  // Start from baseline transition matrix template.
  trans_temp = trans;


  // Iterate households independently in the forward algorithm.
  for(h in 1:n_hh) {

    // alpha is the normalized forward probability on probability scale.
    matrix[hh_size[h] * n_states, max(hh_tmax) - min(hh_tmin) + 1] alpha;

    // llik[p, t] stores log normalizing constants for each person/time.
    matrix[hh_size[h], max(hh_tmax) - min(hh_tmin) + 1] llik;

    // Household-specific slices of observation arrays.
    array[obs_per_hh[h], n_obs_type] int y_hh;
    array[obs_per_hh[h]] int part_id_hh;
    array[obs_per_hh[h]] int t_day_hh;

    int index; // pointer to next observation row for this household

    // For each participant i and infectious state s,
    // i_rows[i,s] points to that participant/state row in alpha/logalpha.
    array[hh_size[h], n_states] int i_rows;

    // Offset in participant indexing across previous households.
    int last_lik;

    // 1 if current participant/time has an observation row, else 0.
    int obs_switch;

    llik = rep_matrix(0, hh_size[h], max(hh_tmax) - min(hh_tmin) + 1);

    if(h == 1) {
      last_lik = 0;
    } else {
      last_lik = sum(hh_size[1:(h-1)]);
    }

    // Subset long vectors/matrices to this household's observation rows.
    y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),];
    t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
    part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];

    index = 1;

    { // START FORWARD ALGORITHM

    // -------------------------
    // Initialization at first modeled day (t = 1)
    // -------------------------
    for(i in 1:hh_size[h]) {
      array[n_states] int ref;

      // Observation likelihood contribution for this person/time.
      // Default is 1 (no information) when there is no observation/missingness.
      matrix[n_obs_type, n_states] obs;

      // Ref gives the row block in global logalpha for participant i.
      ref = linspaced_int_array(n_states,
                                n_states * last_lik + n_states * (i-1) + 1,
                                n_states * last_lik + n_states * (i-1) + n_states);

      obs_switch = 0;

      // Detect whether the next observation row belongs to (day 1, participant i).
      if(t_day_hh[index] == 1) {
        if(part_id_hh[index] == i) {
          obs_switch = 1;
        }
      }

      // Build observation likelihood vectors by state.
      if(obs_switch == 1) {
        for(k in 1:n_obs_type) {
          if(y_hh[index, k] != -1) {
            obs[k, ] = obs_prob[k][y_hh[index, k],];
          } else {
            obs[k, ] = rep_row_vector(1, n_states);
          }
        }
      } else {
        obs = rep_matrix(1, n_obs_type, n_states);
      }

      // Advance observation pointer only when an obs row was consumed.
      if(obs_switch == 1) {
       index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1);
      }

      // Initialize forward log-probability with prior state probabilities.
      logalpha[ref, 1] = log(init_probs);

      // Add log-likelihood contributions from each observation type.
      for(k in 1:n_obs_type) {
        logalpha[ref, 1] = logalpha[ref, 1] + to_vector(log(obs[k,]));
      }

      // Precompute row indices of infectious states for participant i.
      for(s in inf_states) {
        i_rows[i, s] = n_states * (i-1) + s;
      }

      // log-sum-exp is the normalization constant for this participant/time.
      llik[i, 1] = log_sum_exp(logalpha[ref,1]);

      // Softmax gives normalized forward state probabilities.
      alpha[(n_states * (i-1) + 1):(n_states * (i-1) + n_states), 1] =
        softmax(logalpha[ref,1]);

    } // end participant loop (initialization)

    // -------------------------
    // Recursion for t = 2..T
    // -------------------------
    for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

      for(p in 1:hh_size[h]) {
        // no_inf_prob[s] = probability participant p avoids infection pressure
        // associated with infectious state s from all household members.
        array[n_states] real no_inf_prob;

        // no_hh_inf_prob[j,s] = probability p avoids infection from member j
        // via infectious state s.
        matrix[hh_size[h], n_states] no_hh_inf_prob;

        array[n_states] int ref;
        vector[n_states] logalpha_temp; // previous-time log forward state probs
        matrix[n_obs_type, n_states] obs;
        matrix[n_states, n_states] mult_temp;

        ref = linspaced_int_array(n_states,
                                  n_states * last_lik + n_states * (p-1) + 1,
                                  n_states * last_lik + n_states * (p-1) + n_states);

        logalpha_temp = logalpha[ref, tt-1];

        obs_switch = 0;

        // Detect whether next observation row belongs to (day tt, participant p).
        if(t_day_hh[index] == tt) {
          if(part_id_hh[index] == p) {
            obs_switch = 1;
          }
        }

        if(obs_switch == 1) {
          for(k in 1:n_obs_type) {
            if(y_hh[index, k] != -1) {
              obs[k,] = obs_prob[k][y_hh[index, k], ];
            } else {
              obs[k,] = rep_row_vector(1, n_states);
            }
          }
        } else {
          obs = rep_matrix(1, n_obs_type, n_states);
        }

        if(obs_switch == 1) {
          index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1);
        }

        // ct indexes which intra-household infection intercept column to use.
        int ct = 1;

        // Compute probability of avoiding infection contribution by state.
        for(s in 1:n_states) {
          if(is_in(s, inf_states)) {

            // For each potential source j:
            // P(avoid from j via state s) =
            //   P(j in infectious state s) * (1 - ih_prob[p,ct])
            //   + P(j not in state s)
            no_hh_inf_prob[,s] =
              to_vector(alpha[i_rows[, s], tt-1]) * (1 - ih_prob[last_lik + p, ct])
              + (1 - to_vector(alpha[i_rows[,s], tt-1]));

            ct += 1;

            // No self-infection contribution.
            no_hh_inf_prob[p, s] = 1;

          } else {
            // Non-infectious states contribute no infection pressure.
            no_hh_inf_prob[,s] = rep_vector(1, hh_size[h]);
          }

          // Combine independent source contributions multiplicatively.
          no_inf_prob[s] = prod(no_hh_inf_prob[,s]);
        }

        // Fill transition entries designated as estimated.
        for(m in 1:n_trans_fit) {
          if(sum(source_states[m,]) == 0) {
            // Non-infection transition: direct parameter.
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = params[param_index[m]];
          } else {
            // Infection transition: combine selected no-infection terms,
            // then include extra-household infection pressure.
            real no_inf;
            no_inf = 1;
            for(s in 1:n_states) {
              if(source_states[m,s] == 1) {
                no_inf = no_inf * no_inf_prob[s];
              }
            }
            // 1 - [prob no household infection * prob no external infection]
            trans_temp[trans_index[m, 1],trans_index[m, 2]] =
              1 - (no_inf * (1 - eh_prob[last_lik + p]));
          }
        }

        // Apply estimated transition multipliers.
        // Reset each time because updates can be self-referential.
        mult_temp = transition_multiplier;
        for(m in 1:n_mult_fit) {
          if(mult_param_index[m] > 0) {
            mult_temp[mult_index[m, 1],mult_index[m, 2]] =
              mult_params[mult_param_index[m]];
          } else {
            mult_temp[mult_index[m, 1], mult_index[m, 2]] +=
              -1 * mult_params[-mult_param_index[m]];
          }
        }

        // Transition split: element-wise scaling.
        trans_temp = trans_temp .* mult_temp;

        // Recompute diagonals so each column remains stochastic.
        for(i in 1:cols(trans_temp)) {
          trans_temp[i,i] = get_diagonal_element(trans_temp, i);
        }

        // Stabilize and renormalize after edits.
        trans_temp = replace_zeroes(trans_temp, epsilon);
        trans_temp = normalize_cols(trans_temp);

        // Forward recursion:
        // predicted state prob = trans * previous state prob,
        // then incorporate observation likelihood in log-space.
        logalpha[ref, tt] = log(trans_temp * exp(logalpha_temp));
        for(k in 1:n_obs_type) {
          logalpha[ref, tt] = logalpha[ref, tt] + to_vector(log(obs[k,]));
        }

        // Normalize for stable recursion and to obtain filtering probabilities.
        alpha[(n_states * (p-1) + 1):(n_states * (p-1) + n_states), tt] =
          softmax(logalpha[ref,tt]);

        // Save log normalization constant for likelihood.
        llik[p, tt] = log_sum_exp(logalpha[ref,tt]);

      } // end participant loop at time tt

      // This implementation uses only final-time contributions per household.
      if(tt == (hh_tmax[h] - hh_tmin[h] + 1)) {
        llik_final[h] = sum(llik[,tt]);
      }

    } // end time recursion
    } // END FORWARD ALGORITHM

  } // end household loop

}

model {

  // Weakly informative priors for infection covariate effects.
  beta_eh ~ normal(-3, 3);
  beta_ih ~ normal(-3, 3);

  // Add household log-likelihood contributions (computed in transformed parameters).
  target += sum(llik_final);

}
