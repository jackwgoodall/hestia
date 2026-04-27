functions {

  // -------------------------------------------------------------------------
  // Helper functions (unchanged from original)
  // -------------------------------------------------------------------------

  // Given an almost-complete transition matrix column, recover the diagonal
  // entry needed to make that column sum to 1.
  real get_diagonal_element(matrix m, int i) {
    real out;
    out = 1;
    for (j in 1:rows(m)) {
      if (j != i) {
        out = out - m[j, i];
      }
    }
    return out;
  }

  // Minimal helper equivalent to R's `%in%` for integer arrays.
  int is_in(int pos, array[] int pos_var) {
    int pos_match;
    array[size(pos_var)] int all_matches;
    for (p in 1:(size(pos_var))) {
      all_matches[p] = (pos_var[p] == pos);
    }
    if (sum(all_matches) > 0) {
      pos_match = 1;
      return pos_match;
    } else {
      pos_match = 0;
      return pos_match;
    }
  }

  // Force each column of a matrix to sum to 1.
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;
    for (i in 1:cols(m)) {
      out[, i] = m[, i] / sum(m[, i]);
    }
    return out;
  }

  // Replace exact zeros by a small positive value.
  matrix replace_zeroes(matrix m, real epsilon) {
    matrix[rows(m), cols(m)] out;
    out = m;
    for (i in 1:rows(m)) {
      for (j in 1:cols(m)) {
        if (m[i, j] == 0) {
          out[i, j] = epsilon;
        }
      }
    }
    return out;
  }

  // -------------------------------------------------------------------------
  // Partial log-likelihood function for reduce_sum
  //
  // Runs the HMM forward algorithm for a contiguous slice of households
  // (indexed start:end) and returns their combined log-likelihood.
  // reduce_sum calls this function in parallel across threads.
  // -------------------------------------------------------------------------
  real partial_log_lik(
    // Slice of household indices — reduce_sum passes these automatically.
    array[] int slice_hh,
    int start,
    int end,
    // ---- Latent-state transition structure ----
    int n_states,
    int n_inf_states,
    array[] int inf_states,
    int n_trans_fit,
    array[] int param_index,
    array[, ] int trans_index,
    array[, ] int source_states,
    matrix trans,
    matrix transition_multiplier,
    int n_mult_fit,
    array[] int mult_param_index,
    array[, ] int mult_index,
    // ---- Household layout ----
    array[] int hh_size,
    // ---- Observation data ----
    int n_obs_type,
    int n_unique_obs,
    array[, ] int y,
    array[] int part_id,
    array[] int t_day,
    array[] int obs_per_hh,
    array[] int hh_start_ind,
    array[] int hh_end_ind,
    array[] int hh_tmin,
    array[] int hh_tmax,
    // ---- Time-varying covariates (pre-computed probabilities) ----
    int T_global,
    array[] matrix ih_prob,   // [n_inf_prob] matrix[sum(hh_size), T_global]
    matrix eh_prob,           // matrix[sum(hh_size), T_global]
    // ---- Observation model and initialisation ----
    array[] matrix obs_prob,
    vector init_probs,
    real epsilon,
    int n_inf_prob,
    // ---- Fitted parameters (probability scale) ----
    array[] real params,
    array[] real mult_params
  ) {

    // Accumulator for the log-likelihood contributions of this slice.
    real llik_sum = 0;

    // Loop over only the households assigned to this thread.
    for (h in start:end) {

      // Normalized forward probabilities for this household.
      matrix[hh_size[h] * n_states, max(hh_tmax) - min(hh_tmin) + 1] alpha;

      // Per-participant, per-time log-normalising constants.
      matrix[hh_size[h], max(hh_tmax) - min(hh_tmin) + 1] llik;

      // Household-specific slices of the observation arrays.
      array[obs_per_hh[h], n_obs_type] int y_hh;
      array[obs_per_hh[h]] int part_id_hh;
      array[obs_per_hh[h]] int t_day_hh;

      // Log forward probabilities — local to this household.
      matrix[hh_size[h] * n_states, max(hh_tmax) - min(hh_tmin) + 1] logalpha_hh;

      // Working transition matrix (rebuilt each time step).
      matrix[n_states, n_states] trans_temp;

      // Pointer into the household observation arrays.
      int index;

      // Row-index lookup: i_rows[i, s] gives the row in alpha for
      // participant i in state s. Used to read other members' filtering
      // distributions when computing infection pressure.
      array[hh_size[h], n_states] int i_rows;

      // Number of participants in households before household h —
      // the person-level offset into ih_prob / eh_prob.
      int last_lik;

      int obs_switch;

      llik = rep_matrix(0, hh_size[h], max(hh_tmax) - min(hh_tmin) + 1);

      if (h == 1) {
        last_lik = 0;
      } else {
        last_lik = sum(hh_size[1:(h - 1)]);
      }

      y_hh        = y[(hh_start_ind[h]):(hh_end_ind[h]), ];
      t_day_hh    = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
      part_id_hh  = part_id[(hh_start_ind[h]):(hh_end_ind[h])];

      index = 1;

      { // START FORWARD ALGORITHM

        // ---- Initialisation at the first modelled day ----
        for (i in 1:hh_size[h]) {

          array[n_states] int ref;
          matrix[n_obs_type, n_states] obs;

          ref = linspaced_int_array(
            n_states,
            n_states * (i - 1) + 1,
            n_states * (i - 1) + n_states
          );

          obs_switch = 0;
          if (t_day_hh[index] == 1) {
            if (part_id_hh[index] == i) {
              obs_switch = 1;
            }
          }

          if (obs_switch == 1) {
            for (k in 1:n_obs_type) {
              if (y_hh[index, k] != -1) {
                obs[k, ] = obs_prob[k][y_hh[index, k], ];
              } else {
                obs[k, ] = rep_row_vector(1, n_states);
              }
            }
          } else {
            obs = rep_matrix(1, n_obs_type, n_states);
          }

          if (obs_switch == 1) {
            index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1);
          }

          logalpha_hh[ref, 1] = log(init_probs);
          for (k in 1:n_obs_type) {
            logalpha_hh[ref, 1] = logalpha_hh[ref, 1] + to_vector(log(obs[k, ]));
          }

          for (s in inf_states) {
            i_rows[i, s] = n_states * (i - 1) + s;
          }

          llik[i, 1] = log_sum_exp(logalpha_hh[ref, 1]);
          alpha[(n_states * (i - 1) + 1):(n_states * (i - 1) + n_states), 1] =
            softmax(logalpha_hh[ref, 1]);

        } // end initialisation participant loop

        // ---- Forward recursion for later days ----
        for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

          int actual_day;
          actual_day = hh_tmin[h] + tt - 1;

          for (p in 1:hh_size[h]) {

            array[n_states] real no_inf_prob;
            matrix[hh_size[h], n_states] no_hh_inf_prob;
            array[n_states] int ref;
            vector[n_states] logalpha_temp;
            matrix[n_obs_type, n_states] obs;
            matrix[n_states, n_states] mult_temp;

            ref = linspaced_int_array(
              n_states,
              n_states * (p - 1) + 1,
              n_states * (p - 1) + n_states
            );

            logalpha_temp = logalpha_hh[ref, tt - 1];

            obs_switch = 0;
            if (t_day_hh[index] == tt) {
              if (part_id_hh[index] == p) {
                obs_switch = 1;
              }
            }

            if (obs_switch == 1) {
              for (k in 1:n_obs_type) {
                if (y_hh[index, k] != -1) {
                  obs[k, ] = obs_prob[k][y_hh[index, k], ];
                } else {
                  obs[k, ] = rep_row_vector(1, n_states);
                }
              }
            } else {
              obs = rep_matrix(1, n_obs_type, n_states);
            }

            if (obs_switch == 1) {
              index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1);
            }

            int ct = 1;

            for (s in 1:n_states) {
              if (is_in(s, inf_states)) {
                no_hh_inf_prob[, s] =
                  to_vector(alpha[i_rows[, s], tt - 1])
                    * (1 - ih_prob[ct][last_lik + p, actual_day])
                  + (1 - to_vector(alpha[i_rows[, s], tt - 1]));
                ct += 1;
                no_hh_inf_prob[p, s] = 1;
              } else {
                no_hh_inf_prob[, s] = rep_vector(1, hh_size[h]);
              }
              no_inf_prob[s] = prod(no_hh_inf_prob[, s]);
            }

            // Update fitted transition entries.
            trans_temp = trans;
            for (m in 1:n_trans_fit) {
              if (sum(source_states[m, ]) == 0) {
                trans_temp[trans_index[m, 1], trans_index[m, 2]] =
                  params[param_index[m]];
              } else {
                real no_inf;
                no_inf = 1;
                for (s in 1:n_states) {
                  if (source_states[m, s] == 1) {
                    no_inf = no_inf * no_inf_prob[s];
                  }
                }
                trans_temp[trans_index[m, 1], trans_index[m, 2]] =
                  1 - (no_inf * (1 - eh_prob[last_lik + p, actual_day]));
              }
            }

            // Update splitting multipliers.
            mult_temp = transition_multiplier;
            for (m in 1:n_mult_fit) {
              if (mult_param_index[m] > 0) {
                mult_temp[mult_index[m, 1], mult_index[m, 2]] =
                  mult_params[mult_param_index[m]];
              } else {
                mult_temp[mult_index[m, 1], mult_index[m, 2]] +=
                  -1 * mult_params[-mult_param_index[m]];
              }
            }

            trans_temp = trans_temp .* mult_temp;

            // Fill diagonals.  Guard against negative diagonals (possible
            // during warmup if the off-diagonal sums temporarily exceed 1).
            for (i in 1:cols(trans_temp)) {
              real d = get_diagonal_element(trans_temp, i);
              trans_temp[i, i] = d > 0 ? d : epsilon;
            }

            trans_temp = replace_zeroes(trans_temp, epsilon);
            trans_temp = normalize_cols(trans_temp);

            logalpha_hh[ref, tt] = log(trans_temp * exp(logalpha_temp));
            for (k in 1:n_obs_type) {
              logalpha_hh[ref, tt] =
                logalpha_hh[ref, tt] + to_vector(log(obs[k, ]));
            }

            alpha[(n_states * (p - 1) + 1):(n_states * (p - 1) + n_states), tt] =
              softmax(logalpha_hh[ref, tt]);

            llik[p, tt] = log_sum_exp(logalpha_hh[ref, tt]);

          } // end participant loop at time tt

        } // end time recursion

      } // END FORWARD ALGORITHM

      // Accumulate the final-day log-likelihood for this household.
      llik_sum += sum(llik[, hh_tmax[h] - hh_tmin[h] + 1]);

    } // end household slice loop

    return llik_sum;
  }

} // end functions block


data {

  // ---- Latent-state transition structure ----
  int n_states;
  matrix[n_states, n_states] trans;
  int n_inf_states;
  array[n_inf_states] int inf_states;
  int n_trans_fit;
  array[n_trans_fit] int param_index;
  array[n_trans_fit, 2] int trans_index;
  array[n_trans_fit, n_states] int source_states;
  int n_params;

  // ---- Transition multipliers ----
  matrix[n_states, n_states] transition_multiplier;
  int n_mult_fit;
  int n_mult_params;
  array[n_mult_fit] int mult_param_index;
  array[n_mult_fit, 2] int mult_index;

  // ---- Household layout ----
  int n_hh;
  array[n_hh] int hh_size;

  // ---- Observation data ----
  int n_obs;
  int n_obs_type;
  int n_unique_obs;
  array[n_obs, n_obs_type] int y;
  array[n_obs] int part_id;
  array[n_obs] int t_day;
  array[n_hh] int obs_per_hh;
  array[n_hh] int hh_start_ind;
  array[n_hh] int hh_end_ind;
  array[n_hh] int hh_tmin;
  array[n_hh] int hh_tmax;

  // ---- Covariates ----
  int k_ih;
  int k_eh;
  int T_global;
  array[T_global] matrix[sum(hh_size), k_ih] x_ih;
  array[T_global] matrix[sum(hh_size), k_eh] x_eh;

  // ---- Observation model priors (Beta hyperparameters) ----
  // obs_prob_alpha[k, s] and obs_prob_beta[k, s] are the alpha and beta
  // parameters of the Beta prior on P(positive obs | test k, state s).
  array[n_obs_type, n_states] real<lower=0> obs_prob_alpha;
  array[n_obs_type, n_states] real<lower=0> obs_prob_beta;

  // ---- Bounds on obs_params, used to break label-switching symmetry ----
  // For "negative" states (test should be negative)  set [0,   0.5]
  // For "positive" states (test should be positive)  set [0.5, 1  ]
  // The R helper derives these automatically from the prior mean.
  array[n_obs_type, n_states] real<lower=0, upper=1> obs_lb;
  array[n_obs_type, n_states] real<lower=0, upper=1> obs_ub;

  // ---- Initialisation ----
  vector[n_states] init_probs;
  real epsilon;
  int n_inf_prob;

}

parameters {
  array[n_params] real logit_params;
  array[n_mult_params] real logit_mult_params;
  vector[k_eh] beta_eh;
  vector[k_ih] beta_ih;
  real beta0_eh;
  array[n_inf_prob] real beta0_ih;

  // Re-parameterised observation probabilities on [0,1].  The actual
  // probability obs_params is built in transformed parameters by linear
  // scaling to [obs_lb, obs_ub].  This breaks label-switching symmetry by
  // hard-bounding each (test, state) probability above or below 0.5
  // according to the prior mean.
  array[n_obs_type, n_states] real<lower=0, upper=1> obs_raw;
}

transformed parameters {

  // Pre-compute person- and time-varying infection probabilities.
  array[n_inf_prob] matrix[sum(hh_size), T_global] ih_prob;
  matrix[sum(hh_size), T_global] eh_prob;
  array[n_params] real params;
  array[n_mult_params] real mult_params;

  // Linearly rescale obs_raw ∈ [0,1] to obs_params ∈ [obs_lb, obs_ub].
  array[n_obs_type, n_states] real<lower=0, upper=1> obs_params;

  // Build the obs_prob array from estimated obs_params.
  // obs_prob[k][1, s] = P(negative | test k, state s) = 1 - obs_params[k, s]
  // obs_prob[k][2, s] = P(positive | test k, state s) =     obs_params[k, s]
  array[n_obs_type] matrix[n_unique_obs, n_states] obs_prob;

  params      = inv_logit(logit_params);
  mult_params = inv_logit(logit_mult_params);

  for (tt in 1:T_global) {
    for (i in 1:n_inf_prob) {
      ih_prob[i][, tt] = inv_logit(beta0_ih[i] + x_ih[tt] * beta_ih);
    }
    eh_prob[, tt] = inv_logit(beta0_eh + x_eh[tt] * beta_eh);
  }

  for (k in 1:n_obs_type) {
    for (s in 1:n_states) {
      obs_params[k, s] = obs_lb[k, s]
                       + (obs_ub[k, s] - obs_lb[k, s]) * obs_raw[k, s];
      obs_prob[k][1, s] = 1 - obs_params[k, s];
      obs_prob[k][2, s] =     obs_params[k, s];
    }
  }

}

model {

  // Weakly informative priors on covariate coefficients.
  beta_eh ~ normal(-3, 3);
  beta_ih ~ normal(-3, 3);

  // Beta priors on observation probabilities.  Priors are placed on
  // obs_params (the actual probability scale).  The change-of-variables
  // Jacobian from obs_raw → obs_params is constant in the parameters and
  // is therefore omitted; using target += beta_lpdf(...) suppresses Stan's
  // warning about sampling statements applied to transformed quantities.
  for (k in 1:n_obs_type) {
    for (s in 1:n_states) {
      target += beta_lpdf(obs_params[k, s] | obs_prob_alpha[k, s],
                                              obs_prob_beta[k, s]);
    }
  }

  // Integer array to slice over households.
  array[n_hh] int hh_indices;
  for (h in 1:n_hh) hh_indices[h] = h;

  // Parallelised forward algorithm.
  target += reduce_sum(
    partial_log_lik,
    hh_indices,
    1,
    // ---- shared data ----
    n_states, n_inf_states, inf_states,
    n_trans_fit, param_index, trans_index, source_states,
    trans, transition_multiplier,
    n_mult_fit, mult_param_index, mult_index,
    hh_size,
    n_obs_type, n_unique_obs,
    y, part_id, t_day,
    obs_per_hh, hh_start_ind, hh_end_ind,
    hh_tmin, hh_tmax,
    T_global,
    ih_prob, eh_prob,
    obs_prob, init_probs, epsilon, n_inf_prob,
    // ---- fitted parameters ----
    params, mult_params
  );

}

generated quantities {

  // Recompute per-household log-likelihoods for LOO.
  vector[n_hh] llik_final;

  array[n_hh] int hh_indices;
  for (h in 1:n_hh) hh_indices[h] = h;

  for (h in 1:n_hh) {
    llik_final[h] = partial_log_lik(
      hh_indices[h:h], h, h,
      n_states, n_inf_states, inf_states,
      n_trans_fit, param_index, trans_index, source_states,
      trans, transition_multiplier,
      n_mult_fit, mult_param_index, mult_index,
      hh_size,
      n_obs_type, n_unique_obs,
      y, part_id, t_day,
      obs_per_hh, hh_start_ind, hh_end_ind,
      hh_tmin, hh_tmax,
      T_global,
      ih_prob, eh_prob,
      obs_prob, init_probs, epsilon, n_inf_prob,
      params, mult_params
    );
  }

}
