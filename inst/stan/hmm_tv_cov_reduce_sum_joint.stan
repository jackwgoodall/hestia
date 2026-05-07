// =============================================================================
// Joint SIRS-SIS Hidden Markov Model (now with reduce_sum parallelisation.........)
//
// Joint state space (6 states):
//   1: (S_v, S_b)   Viral susceptible,  Bacterial susceptible
//   2: (S_v, I_b)   Viral susceptible,  Bacterial infected
//   3: (I_v, S_b)   Viral infected,     Bacterial susceptible
//   4: (I_v, I_b)   Viral infected,     Bacterial infected  [co-infected]
//   5: (R_v, S_b)   Viral recovered,    Bacterial susceptible
//   6: (R_v, I_b)   Viral recovered,    Bacterial infected
//
// Cross-immunity effects modelled:
//   cross_ih_trans  — co-infected transmitter (state 4) has enhanced bacterial
//                     transmissibility (logit-scale additive effect)
//   cross_ih_susc   — virally-infected recipient (states 3/4) has enhanced
//                     bacterial susceptibility, both IH and EH (logit-scale)
//
// Viral-dependent bacterial assay sensitivity:
//   obs_params[bac_test, 4] > obs_params[bac_test, 2] encodes that viral
//   co-infection boosts sub-detection bacterial load above assay threshold.
// =============================================================================

functions {

  // ---------------------
  // Helper: diagonal entry that makes column i sum to 1.
  // ---------------------
  real get_diagonal_element(matrix m, int i) {
    real out = 1;
    for (j in 1:rows(m)) {
      if (j != i) out -= m[j, i];
    }
    return out;
  }

  // ---------------------
  // Helper: normalise each column of a matrix to sum to 1.
  // ---------------------
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;
    for (i in 1:cols(m)) out[, i] = m[, i] / sum(m[, i]);
    return out;
  }

  // ---------------------
  // Helper: replace exact zeros with epsilon (avoids log(0)).
  // ---------------------
  matrix replace_zeroes(matrix m, real epsilon) {
    matrix[rows(m), cols(m)] out = m;
    for (i in 1:rows(m))
      for (j in 1:cols(m))
        if (m[i, j] == 0) out[i, j] = epsilon;
    return out;
  }

  // ---------------------
  // Partial log-likelihood for reduce_sum.
  //
  // Runs the HMM forward algorithm for households [start, end] and returns
  // their summed log-likelihood.  All prob matrices are pre-computed
  // in transformed params.
  //
  // `slice_hh` is required by reduce_sum's signature but is not used inside
  // the function body — the loop indexes households via `start:end` directly.
  // ---------------------
  
  real partial_log_lik(
    array[] int slice_hh,
    int start,
    int end,
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
    // ---- Time-varying infection probabilities (pre-computed) ----
    int T_global,
    matrix ih_prob_vir,         // viral IH:               [N, T]
    matrix ih_prob_bac,         // bacterial IH base:       [N, T]
    matrix ih_prob_bac_trans,   // bac IH, enhanced trans:  [N, T]
    matrix ih_prob_bac_susc,    // bac IH, enhanced susc:   [N, T]
    matrix ih_prob_bac_both,    // bac IH, both enhanced:   [N, T]
    matrix eh_prob_vir,         // viral EH:                [N, T]
    matrix eh_prob_bac,         // bacterial EH base:       [N, T]
    matrix eh_prob_bac_susc,    // bac EH, enhanced susc:   [N, T]
    // ---- Fixed rates (probability scale) ----
    real gamma_v,
    real gamma_b,
    real rho,
    // ---- Observation model and initialisation ----
    array[] matrix obs_prob,    // [n_obs_type] matrix[n_unique_obs, 6]
    vector init_probs,          // length 6
    real epsilon
  ) {

    real llik_sum = 0;

    for (h in start:end) {

      int T_hh = max(hh_tmax) - min(hh_tmin) + 1;

      // Forward probability matrices (rows = person*state, cols = time).
      matrix[hh_size[h] * 6, T_hh] alpha;
      matrix[hh_size[h] * 6, T_hh] logalpha_hh;

      // Per-participant, per-time log-normalising constants.
      matrix[hh_size[h], T_hh] llik;

      // Household slices of the observation arrays.
      array[obs_per_hh[h], n_obs_type] int y_hh;
      array[obs_per_hh[h]] int part_id_hh;
      array[obs_per_hh[h]] int t_day_hh;

      // Row-index lookup: i_rows[i, s] = row in alpha for person i, state s.
      array[hh_size[h], 6] int i_rows;

      // Offset into the global person-level probability matrices.
      int last_lik;

      int obs_switch;
      int index;

      llik = rep_matrix(0, hh_size[h], T_hh);

      last_lik = (h == 1) ? 0 : sum(hh_size[1:(h - 1)]);

      y_hh       = y[(hh_start_ind[h]):(hh_end_ind[h]), ];
      t_day_hh   = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
      part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];
      index = 1;

      { // START FORWARD ALGORITHM

        // ------------
        // Initialisation at the first modelled day
        // ------------
        for (i in 1:hh_size[h]) {

          array[6] int ref = linspaced_int_array(6, 6*(i-1)+1, 6*i);
          matrix[n_obs_type, 6] obs;

          obs_switch = 0;
          if (t_day_hh[index] == 1 && part_id_hh[index] == i) obs_switch = 1;

          if (obs_switch == 1) {
            for (k in 1:n_obs_type) {
              if (y_hh[index, k] != -1)
                obs[k, ] = obs_prob[k][y_hh[index, k], ];
              else
                obs[k, ] = rep_row_vector(1, 6);
            }
          } else {
            obs = rep_matrix(1, n_obs_type, 6);
          }

          if (obs_switch == 1)
            index = min(index + 1, obs_per_hh[h]);

          logalpha_hh[ref, 1] = log(init_probs);
          for (k in 1:n_obs_type)
            logalpha_hh[ref, 1] += to_vector(log(obs[k, ]));

          // Initialise i_rows for all 6 states.
          for (s in 1:6) i_rows[i, s] = 6*(i-1) + s;

          llik[i, 1] = log_sum_exp(logalpha_hh[ref, 1]);
          alpha[(6*(i-1)+1):(6*i), 1] = softmax(logalpha_hh[ref, 1]);

        } // end init participant loop

        // ------------
        // Forward recursion for t = 2, ..., this household's last day.
        // Note: T_hh is the *global* maximum window (used only for matrix
        // sizing); the loop must stop at this household's own length.
        // ------------
        
        for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

          int actual_day = hh_tmin[h] + tt - 1;

          for (p in 1:hh_size[h]) {

            array[6] int ref = linspaced_int_array(6, 6*(p-1)+1, 6*p);
            vector[6] logalpha_temp = logalpha_hh[ref, tt - 1];
            matrix[n_obs_type, 6] obs;
            matrix[6, 6] trans_temp;

            // Observation at (p, tt)
            obs_switch = 0;
            if (t_day_hh[index] == tt && part_id_hh[index] == p) obs_switch = 1;

            if (obs_switch == 1) {
              for (k in 1:n_obs_type) {
                if (y_hh[index, k] != -1)
                  obs[k, ] = obs_prob[k][y_hh[index, k], ];
                else
                  obs[k, ] = rep_row_vector(1, 6);
              }
            } else {
              obs = rep_matrix(1, n_obs_type, 6);
            }

            if (obs_switch == 1) index = min(index + 1, obs_per_hh[h]);

            // --------
            // Compute household force of infection for person p
            //
            //  no_vir_inf      — P(p not virally infected by household)
            //  no_bac_inf_base — P(p not bac-infected by household)
            //                    for base susceptibility (p in S_v or R_v)
            //  no_bac_inf_susc — "" but p is virally infected (state 3)
            // --------
            real no_vir_inf      = 1.0;
            real no_bac_inf_base = 1.0;
            real no_bac_inf_susc = 1.0;

            for (q in 1:hh_size[h]) {
              if (q != p) {
                // Viral source: states 3 (I_v,S_b) and 4 (I_v,I_b)
                real prob_Iv = alpha[i_rows[q, 3], tt-1]
                             + alpha[i_rows[q, 4], tt-1];
                no_vir_inf *= 1.0 - prob_Iv * ih_prob_vir[last_lik + p, actual_day];

                // Bacterial sources:
                //   base transmitters: states 2 (S_v,I_b) and 6 (R_v,I_b)
                //   enhanced transmitter: state 4 (I_v,I_b — co-infected)
                real prob_bac_base = alpha[i_rows[q, 2], tt-1]
                                   + alpha[i_rows[q, 6], tt-1];
                real prob_coinf    = alpha[i_rows[q, 4], tt-1];

                // Base susceptibility (p not virally infected)
                no_bac_inf_base *= 1.0
                  - prob_bac_base * ih_prob_bac[last_lik + p, actual_day]
                  - prob_coinf    * ih_prob_bac_trans[last_lik + p, actual_day];

                // Enhanced susceptibility (p is virally infected, state 3)
                no_bac_inf_susc *= 1.0
                  - prob_bac_base * ih_prob_bac_susc[last_lik + p, actual_day]
                  - prob_coinf    * ih_prob_bac_both[last_lik + p, actual_day];
              }
            }

            // --------
            // Build transition matrix (columns = source states)
            // Simultaneous transitions are approximated as zero (??valid for
            // small daily time steps).
            // ---------
            trans_temp = rep_matrix(0.0, 6, 6);

            // Viral infection probability (same for all S_v source states)
            real p_vir = 1.0 - no_vir_inf * (1.0 - eh_prob_vir[last_lik + p, actual_day]);

            // Bacterial infection probability: base susceptibility (states 1, 5)
            real p_bac_base = 1.0 - no_bac_inf_base * (1.0 - eh_prob_bac[last_lik + p, actual_day]);

            // Bacterial infection probability: enhanced susceptibility (state 3)
            real p_bac_susc = 1.0 - no_bac_inf_susc * (1.0 - eh_prob_bac_susc[last_lik + p, actual_day]);

            // Viral acquisition: S_v -> I_v
            trans_temp[3, 1] = p_vir;   // (S_v,S_b) -> (I_v,S_b)
            trans_temp[4, 2] = p_vir;   // (S_v,I_b) -> (I_v,I_b)

            // Bacterial acquisition: S_b -> I_b
            trans_temp[2, 1] = p_bac_base;  // (S_v,S_b) -> (S_v,I_b)
            trans_temp[4, 3] = p_bac_susc;  // (I_v,S_b) -> (I_v,I_b)  [i.e. enhanced/reduced]
            trans_temp[6, 5] = p_bac_base;  // (R_v,S_b) -> (R_v,I_b)

            // Bacterial recovery: I_b -> S_b
            trans_temp[1, 2] = gamma_b;   // (S_v,I_b) -> (S_v,S_b)
            trans_temp[3, 4] = gamma_b;   // (I_v,I_b) -> (I_v,S_b)
            trans_temp[5, 6] = gamma_b;   // (R_v,I_b) -> (R_v,S_b)

            // Viral recovery: I_v -> R_v
            trans_temp[5, 3] = gamma_v;   // (I_v,S_b) -> (R_v,S_b)
            trans_temp[6, 4] = gamma_v;   // (I_v,I_b) -> (R_v,I_b)

            // Viral waning immunity: R_v -> S_v
            trans_temp[1, 5] = rho;       // (R_v,S_b) -> (S_v,S_b)
            trans_temp[2, 6] = rho;       // (R_v,I_b) -> (S_v,I_b)

            // Fill diagonal so each column sums to 1.  Guard against
            // negative diagonals (possible during warmup if the sum of
            // off-diagonals exceeds 1).
            for (s in 1:6) {
              real d = get_diagonal_element(trans_temp, s);
              trans_temp[s, s] = d > 0 ? d : epsilon;
            }

            trans_temp = replace_zeroes(trans_temp, epsilon);
            trans_temp = normalize_cols(trans_temp);

            // Log forward step
            logalpha_hh[ref, tt] = log(trans_temp * exp(logalpha_temp));
            for (k in 1:n_obs_type)
              logalpha_hh[ref, tt] += to_vector(log(obs[k, ]));

            alpha[(6*(p-1)+1):(6*p), tt] = softmax(logalpha_hh[ref, tt]);
            llik[p, tt] = log_sum_exp(logalpha_hh[ref, tt]);

          } // end participant loop at tt
        } // end time recursion

      } // END FORWARD ALGORITHM

      llik_sum += sum(llik[, hh_tmax[h] - hh_tmin[h] + 1]);

    } // end household slice loop

    return llik_sum;
  }

} 


// =============================================================================
data {

  // ---- Household layout ----
  int n_hh;
  array[n_hh] int hh_size;

  // ---- Observation data ----
  int n_obs;
  int n_obs_type;     // total number of test types (viral + bacterial)
  int n_unique_obs;   // 2 for binary tests
  array[n_obs, n_obs_type] int y;   // 1 = negative, 2 = positive, -1 = missing
  array[n_obs] int part_id;
  array[n_obs] int t_day;
  array[n_hh] int obs_per_hh;
  array[n_hh] int hh_start_ind;
  array[n_hh] int hh_end_ind;
  array[n_hh] int hh_tmin;
  array[n_hh] int hh_tmax;

  // ---- Covariates (shared design matrix for both pathogens) ----
  // Separate coefficient vectors per pathogen are estimated.
  int k_ih;     // number of IH covariate columns
  int k_eh;     // number of EH covariate columns
  int T_global;
  array[T_global] matrix[sum(hh_size), k_ih] x_ih;
  array[T_global] matrix[sum(hh_size), k_eh] x_eh;

  // ---- Observation model Beta prior hyperparameters ----
  // obs_prob_alpha[k, s]: k = test type (1..n_obs_type), s = joint state (1..6)
  array[n_obs_type, 6] real<lower=0> obs_prob_alpha;
  array[n_obs_type, 6] real<lower=0> obs_prob_beta;

  // ---- Bounds on obs_params, used to break the pesky label-switching symmetry ----
  // For "negative" states (test should be negative)  set [0,   0.5]
  // For "positive" states (test should be positive)  set [0.5, 1  ]
  // (The R helper derives these automatically from the prior mean.)
  array[n_obs_type, 6] real<lower=0, upper=1> obs_lb;
  array[n_obs_type, 6] real<lower=0, upper=1> obs_ub;

  // ---- Initialisation ----
  simplex[6] init_probs;
  real<lower=0> epsilon;

}


// =============================================================================
parameters {

  // ---- Viral dynamics ----
  real logit_gamma_v;         // viral recovery rate
  real logit_rho;             // viral waning immunity rate
  real beta0_ih_vir;          // viral IH intercept
  real beta0_eh_vir;          // viral EH intercept
  vector[k_ih] beta_ih_vir;   // viral IH covariate effects
  vector[k_eh] beta_eh_vir;   // viral EH covariate effects

  // ---- Bacterial dynamics ----
  real logit_gamma_b;         // bacterial recovery rate
  real beta0_ih_bac;          // bacterial IH intercept
  real beta0_eh_bac;          // bacterial EH intercept
  vector[k_ih] beta_ih_bac;   // bacterial IH covariate effects
  vector[k_eh] beta_eh_bac;   // bacterial EH covariate effects

  // ---- Cross-immunity ----
  // Logit-scale additive effects on bacterial transmission/susceptibility.
  real cross_ih_trans;   // co-infected transmitter -> modifiy bacterial transmissibility
  real cross_ih_susc;    // I_v recipient -> modify. bacterial susceptibility (IH + EH)

  // ---- Observation model ----
  // Re-parameterised on [0,1]; the actual probability obs_params is
  // obtained in transformed parameters by linear scaling to [obs_lb, obs_ub].
  // This (I think) breaks label-switching symmetry by hard-bounding each (test, state)
  // probability above or below 0.5 according to the prior's mean.
  array[n_obs_type, 6] real<lower=0, upper=1> obs_raw;

}


// =============================================================================
transformed parameters {

  // Rates on probability scale.
  real<lower=0, upper=1> gamma_v = inv_logit(logit_gamma_v);
  real<lower=0, upper=1> gamma_b = inv_logit(logit_gamma_b);
  real<lower=0, upper=1> rho     = inv_logit(logit_rho);

  // Pre-computed person×time infection probability matrices.
  // Eight variants to cover all cross-immunity combinations.
  matrix[sum(hh_size), T_global] ih_prob_vir;
  matrix[sum(hh_size), T_global] ih_prob_bac;
  matrix[sum(hh_size), T_global] ih_prob_bac_trans;  // + cross_ih_trans
  matrix[sum(hh_size), T_global] ih_prob_bac_susc;   // + cross_ih_susc
  matrix[sum(hh_size), T_global] ih_prob_bac_both;   // + both
  matrix[sum(hh_size), T_global] eh_prob_vir;
  matrix[sum(hh_size), T_global] eh_prob_bac;
  matrix[sum(hh_size), T_global] eh_prob_bac_susc;   // + cross_ih_susc

  for (tt in 1:T_global) {
    ih_prob_vir[, tt]       = inv_logit(beta0_ih_vir + x_ih[tt] * beta_ih_vir);
    ih_prob_bac[, tt]       = inv_logit(beta0_ih_bac + x_ih[tt] * beta_ih_bac);
    ih_prob_bac_trans[, tt] = inv_logit(beta0_ih_bac + x_ih[tt] * beta_ih_bac + cross_ih_trans);
    ih_prob_bac_susc[, tt]  = inv_logit(beta0_ih_bac + x_ih[tt] * beta_ih_bac + cross_ih_susc);
    ih_prob_bac_both[, tt]  = inv_logit(beta0_ih_bac + x_ih[tt] * beta_ih_bac + cross_ih_trans + cross_ih_susc);
    eh_prob_vir[, tt]       = inv_logit(beta0_eh_vir + x_eh[tt] * beta_eh_vir);
    eh_prob_bac[, tt]       = inv_logit(beta0_eh_bac + x_eh[tt] * beta_eh_bac);
    eh_prob_bac_susc[, tt]  = inv_logit(beta0_eh_bac + x_eh[tt] * beta_eh_bac + cross_ih_susc);
  }

  // Linearly rescale obs_raw E [0,1] to obs_params E [obs_lb, obs_ub].
  array[n_obs_type, 6] real<lower=0, upper=1> obs_params;
  for (k in 1:n_obs_type) {
    for (s in 1:6) {
      obs_params[k, s] = obs_lb[k, s]
                       + (obs_ub[k, s] - obs_lb[k, s]) * obs_raw[k, s];
    }
  }

  // Observation probability matrices: obs_prob[k][row, state]
  // Row 1 = P(negative | state), row 2 = P(positive | state).
  array[n_obs_type] matrix[n_unique_obs, 6] obs_prob;
  for (k in 1:n_obs_type) {
    for (s in 1:6) {
      obs_prob[k][1, s] = 1 - obs_params[k, s];      // P(negative test k | state s)
      obs_prob[k][2, s] = obs_params[k, s];         // P(positive test k | state s)
    }
  }

}


// =============================================================================
model {

  // ---- Priors ----

  // Rates: centred on plausible daily probabilities
  logit_gamma_v ~ normal(-2, 1);   // prior mean gamma_v ~= 0.12
  logit_gamma_b ~ normal(-2, 1);   // prior mean gamma_b ~= 0.12
  logit_rho     ~ normal(-3, 1);   // prior mean rho     ~= 0.05

  // Baseline IH/EH probabilities: weakly informative ( logit scale)
  beta0_ih_vir ~ normal(-4, 2);
  beta0_ih_bac ~ normal(-4, 2);
  beta0_eh_vir ~ normal(-4, 2);
  beta0_eh_bac ~ normal(-4, 2);

  // Covariate effects
  beta_ih_vir ~ normal(0, 1);
  beta_ih_bac ~ normal(0, 1);
  beta_eh_vir ~ normal(0, 1);
  beta_eh_bac ~ normal(0, 1);

  // Cross-immunity
  cross_ih_trans ~ normal(0, 1);
  cross_ih_susc  ~ normal(0, 1);

  // Observation model : Beta priors are placed on the
  // probability scale (obs_params), not on obs_raw.
  for (k in 1:n_obs_type) {
    for (s in 1:6) {
      target += beta_lpdf(obs_params[k, s] | obs_prob_alpha[k, s],
                                              obs_prob_beta[k, s]);
    }
  }

  // ---- Parallelised forward algorithm ----
  array[n_hh] int hh_indices;
  for (h in 1:n_hh) hh_indices[h] = h;

  target += reduce_sum(
    partial_log_lik,
    hh_indices,
    1,
    // shared data
    hh_size,
    n_obs_type, n_unique_obs,
    y, part_id, t_day,
    obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
    T_global,
    ih_prob_vir, ih_prob_bac, ih_prob_bac_trans, ih_prob_bac_susc, ih_prob_bac_both,
    eh_prob_vir, eh_prob_bac, eh_prob_bac_susc,
    gamma_v, gamma_b, rho,
    obs_prob, init_probs, epsilon
  );

}


// =============================================================================
generated quantities {

  // Per-household log-likelihoods (for LOO-CV.)
  vector[n_hh] llik_final;

  {
    array[n_hh] int hh_indices;
    for (h in 1:n_hh) hh_indices[h] = h;

    for (h in 1:n_hh) {
      llik_final[h] = partial_log_lik(
        hh_indices[h:h], h, h,
        hh_size,
        n_obs_type, n_unique_obs,
        y, part_id, t_day,
        obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
        T_global,
        ih_prob_vir, ih_prob_bac, ih_prob_bac_trans, ih_prob_bac_susc, ih_prob_bac_both,
        eh_prob_vir, eh_prob_bac, eh_prob_bac_susc,
        gamma_v, gamma_b, rho,
        obs_prob, init_probs, epsilon
      );
    }
  }

}
