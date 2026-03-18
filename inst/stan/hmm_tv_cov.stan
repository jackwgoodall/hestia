functions {

  // Given an almost-complete transition matrix column, recover the diagonal entry
  // needed to make that column sum to 1. In this model:
  // - columns correspond to the state at time t-1 ("from" state)
  // - rows correspond to the state at time t ("to" state)
  // So for column i we set:
  //   trans[i, i] = 1 - sum_{j != i} trans[j, i]
  // This is used after we update off-diagonal transition probabilities and need
  // the self-transition to absorb whatever probability mass is left.
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

  // Minimal helper equivalent to R's `%in%` for integer arrays.
  // It is used to ask questions like "is latent state s one of the infectious
  // states?" and returns 1 for yes, 0 for no.
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

  // Force each column of a matrix to sum to 1.
  // Transition columns must represent proper probability distributions over the
  // next state, so this is a final cleanup step after edits/multipliers.
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;

    for(i in 1:cols(m)) {
      out[,i] = m[,i] / sum(m[,i]);
    }
    return out;
  }

  // Replace exact zeros by a small positive value.
  // This avoids numerical problems later when the code takes `log(...)`.
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
  // Latent-state transition structure
  // -------------------------

  // Total number of latent epidemiological states in the HMM.
  // Example interpretation could be susceptible / exposed / infectious / recovered,
  // but the model is written generically and only needs the state count and the
  // transition layout supplied from R.
  int n_states;

  // Baseline transition matrix template.
  // Columns are the state at time t-1 and rows are the state at time t.
  // Some entries are fixed by the data passed in, while others are overwritten
  // inside the model using fitted parameters.
  matrix[n_states, n_states] trans;

  // Number of latent states that count as infectious sources for transmission.
  int n_inf_states;

  // Integer indices of those infectious states.
  // These are used to decide which latent states can generate within-household
  // infection pressure on other members.
  array[n_inf_states] int inf_states;

  // Number of transition-matrix entries that are not fixed and will be updated
  // during model evaluation.
  int n_trans_fit;

  // For each fitted transition entry:
  // - a positive value means use `params[param_index[m]]`
  // - 0 means this is an infection-driven transition whose probability is built
  //   from household and extra-household infection pressure instead.
  array[n_trans_fit] int param_index;

  // Location of each fitted transition in the matrix.
  // `trans_index[m, 1]` is the destination row and `trans_index[m, 2]` is the
  // origin column for the m-th editable transition entry.
  array[n_trans_fit, 2] int trans_index;

  // For infection-driven transitions only, this row vector says which infectious
  // latent states contribute to the infection hazard for that transition.
  // A row of all zeros means the transition is not infection-driven.
  array[n_trans_fit, n_states] int source_states;

  // Number of non-infection transition probabilities that are estimated on the
  // logit scale in `logit_params`.
  int n_params;

  // -------------------------
  // Transition multipliers / splits
  // -------------------------

  // Baseline matrix of multiplicative modifiers applied element-wise to `trans`.
  // This lets one baseline transition be split across several destinations.
  matrix[n_states, n_states] transition_multiplier;

  // Number of multiplier entries that are themselves estimated.
  int n_mult_fit;

  // Number of unique free multiplier parameters.
  int n_mult_params;

  // Map from each editable multiplier entry to a free parameter.
  // Positive values mean "use mult_params[index]".
  // Negative values mean "use a complementary 1 - parameter style update" via
  // the additive expression later in the code.
  array[n_mult_fit] int mult_param_index;

  // Row/column locations of editable multiplier entries.
  array[n_mult_fit, 2] int mult_index;

  // -------------------------
  // Household layout
  // -------------------------

  // Number of households.
  int n_hh;

  // Household sizes. The sum over this array is the total number of people.
  array[n_hh] int hh_size;

  // -------------------------
  // Observation data
  // -------------------------

  // Total number of observation records across all households, all time points,
  // and all enrolled participants.
  int n_obs;

  // Number of distinct observation channels per record, for example symptoms,
  // PCR, serology, or other measurement types.
  int n_obs_type;

  // Number of possible categorical outcomes for each observation channel.
  int n_unique_obs;

  // Observed outcomes, ordered by:
  // 1. household
  // 2. time
  // 3. participant within household
  // `y[r, k]` is the observed category for row r and channel k.
  // A value of -1 means that channel is missing for that record.
  array[n_obs, n_obs_type] int y;

  // For each observation row, which participant within the household it belongs to.
  array[n_obs] int part_id;

  // For each observation row, the observed day index.
  array[n_obs] int t_day;

  // Household-specific indexing helpers so the long observation arrays can be
  // sliced into one household at a time.
  array[n_hh] int obs_per_hh;
  array[n_hh] int hh_start_ind;
  array[n_hh] int hh_end_ind;

  // Minimum and maximum modeled day for each household.
  // These determine how many time steps are run in the forward algorithm.
  array[n_hh] int hh_tmin;
  array[n_hh] int hh_tmax;

  // -------------------------
  // Covariates driving infection pressure
  // -------------------------

  // Number of covariates for intra-household infection probability.
  int k_ih;

  // Number of covariates for extra-household infection probability.
  int k_eh;

  // Total number of calendar days spanned by the study (= max(hh_tmax)).
  // Used to size the time-varying covariate arrays below.
  int T_global;

  // Intra-household covariate array.
  // `x_ih[t]` is a matrix of dimension sum(hh_size) x k_ih giving each
  // person's covariate values on calendar day t.
  array[T_global] matrix[sum(hh_size), k_ih] x_ih;

  // Extra-household covariate array with the same indexing as `x_ih`.
  array[T_global] matrix[sum(hh_size), k_eh] x_eh;

  // -------------------------
  // Observation model and initialization
  // -------------------------

  // Observation probability tables.
  // For observation type k:
  // `obs_prob[k][obs_value, state] = P(observed category obs_value | latent state)`
  // Rows therefore index the observed outcome, and columns index the hidden state.
  array[n_obs_type] matrix[n_unique_obs, n_states] obs_prob;

  // Prior probabilities over the latent state at the first modeled day.
  vector[n_states] init_probs;

  // Small positive constant used to avoid exact zeros before taking logarithms.
  real epsilon;

  // Number of distinct intra-household infection probabilities.
  // This is either:
  // - 1, meaning the same infection probability is used for all infectious states
  // - or one value per infectious state
  int n_inf_prob;
}

parameters {
  // Free transition probabilities on the unconstrained real line.
  // These are transformed with `inv_logit` into (0, 1) inside transformed parameters.
  array[n_params] real logit_params;

  // Same idea for transition multipliers that must live on the probability scale.
  array[n_mult_params] real logit_mult_params;

  // Regression coefficients for the extra-household infection probability.
  vector[k_eh] beta_eh;

  // Regression coefficients for the intra-household infection probability.
  vector[k_ih] beta_ih;

  // Intercept for extra-household infection probability.
  real beta0_eh;

  // Intercept(s) for intra-household infection probability.
  // There may be one shared intercept or one per infectious-state-specific
  // infection probability, depending on `n_inf_prob`.
  array[n_inf_prob] real beta0_ih;

}

transformed parameters {
  // Final log-likelihood contribution for each household.
  // The model block later adds `sum(llik_final)` directly to the target.
  vector[n_hh] llik_final;

  // Person- and time-varying within-household infection probabilities.
  // `ih_prob[c][person, t]` is the infection probability for person `person`
  // on day `t` via the c-th infectious-state probability slot.
  array[n_inf_prob] matrix[sum(hh_size), T_global] ih_prob;

  // Person- and time-varying extra-household infection probabilities.
  // `eh_prob[person, t]` is the extra-household infection probability for
  // person `person` on day `t`.
  matrix[sum(hh_size), T_global] eh_prob;

  // Global storage for log forward probabilities.
  // For each person we store `n_states` rows, and columns correspond to time steps.
  // Row blocks are stacked household by household, participant by participant.
  matrix[sum(hh_size) * n_states, max(hh_tmax) - min(hh_tmin) + 1] logalpha;

  // Working transition matrix that is repeatedly edited as the recursion proceeds.
  matrix[n_states, n_states] trans_temp;

  // Probability-scale versions of `logit_params` and `logit_mult_params`.
  array[n_params] real params;
  array[n_mult_params] real mult_params;

  // Map unconstrained reals to probabilities.
  params = inv_logit(logit_params);
  mult_params = inv_logit(logit_mult_params);

  // Convert the logistic regressions to person- and time-specific probabilities.
  // Outer loop runs over days; inner loop over infection-probability slots.
  for(tt in 1:T_global) {
    for(i in 1:n_inf_prob) {
      ih_prob[i][,tt] = inv_logit(beta0_ih[i] + x_ih[tt] * beta_ih);
    }
    eh_prob[,tt] = inv_logit(beta0_eh + x_eh[tt] * beta_eh);
  }

  // Start from the baseline transition template.
  trans_temp = trans;


  // Households are conditionally independent given the parameters, so the
  // forward algorithm can be run one household at a time.
  for(h in 1:n_hh) {

    // Normalized forward probabilities for this household only.
    // Each participant contributes a block of `n_states` rows.
    matrix[hh_size[h] * n_states, max(hh_tmax) - min(hh_tmin) + 1] alpha;

    // Per-participant, per-time log normalizing constants from the forward pass.
    // These are later used to build the household log-likelihood.
    matrix[hh_size[h], max(hh_tmax) - min(hh_tmin) + 1] llik;

    // Household-specific slices of the observation arrays.
    array[obs_per_hh[h], n_obs_type] int y_hh;
    array[obs_per_hh[h]] int part_id_hh;
    array[obs_per_hh[h]] int t_day_hh;

    // Pointer to the next observation row to be consumed for this household.
    int index;

    // For participant i and latent state s, `i_rows[i, s]` stores the row index
    // in `alpha` that corresponds to that participant-state combination.
    // This is mainly used to pull out the probability that another household
    // member is currently in infectious state s.
    array[hh_size[h], n_states] int i_rows;

    // Number of participants in all households before household h.
    // This provides the offset into person-level vectors like `ih_prob` and `eh_prob`.
    int last_lik;

    // Indicator for whether the current participant/day actually has an observation
    // row waiting in the household-specific observation arrays.
    int obs_switch;

    llik = rep_matrix(0, hh_size[h], max(hh_tmax) - min(hh_tmin) + 1);

    if(h == 1) {
      last_lik = 0;
    } else {
      last_lik = sum(hh_size[1:(h-1)]);
    }

    // Pull just this household's observations out of the long stacked arrays.
    y_hh = y[(hh_start_ind[h]):(hh_end_ind[h]),];
    t_day_hh = t_day[(hh_start_ind[h]):(hh_end_ind[h])];
    part_id_hh = part_id[(hh_start_ind[h]):(hh_end_ind[h])];

    index = 1;

    { // START FORWARD ALGORITHM

    // -------------------------
    // Initialization at the first modeled day
    // -------------------------
    // For each participant, start from the prior state distribution `init_probs`
    // and then multiply by any observation likelihood available at day 1.
    for(i in 1:hh_size[h]) {
      // `ref` gives the row block in the global `logalpha` matrix corresponding
      // to participant i in this household.
      array[n_states] int ref;

      // Observation likelihood contributions for this participant at this time.
      // `obs[k, s]` will hold P(observation type k | latent state s).
      // If there is no observation for a channel, the entry is set to 1 so it
      // does not change the latent-state probabilities.
      matrix[n_obs_type, n_states] obs;

      ref = linspaced_int_array(n_states,
                                n_states * last_lik + n_states * (i-1) + 1,
                                n_states * last_lik + n_states * (i-1) + n_states);

      obs_switch = 0;

      // Check whether the next stored observation row belongs to participant i
      // on the first modeled day. Because the input is ordered by household,
      // then time, then participant, a single forward-moving pointer is enough.
      if(t_day_hh[index] == 1) {
        if(part_id_hh[index] == i) {
          obs_switch = 1;
        }
      }

      // Build the observation-likelihood matrix.
      // For observed channels, pull the appropriate row from `obs_prob`.
      // For missing channels (coded -1), leave the contribution neutral at 1.
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

      // Consume the observation row only if it was actually matched here.
      if(obs_switch == 1) {
       index = min(index + 1, hh_end_ind[h] - hh_start_ind[h] + 1);
      }

      // Start the forward recursion from the initial state probabilities.
      logalpha[ref, 1] = log(init_probs);

      // Add observation information on the log scale, one observation type at a time.
      for(k in 1:n_obs_type) {
        logalpha[ref, 1] = logalpha[ref, 1] + to_vector(log(obs[k,]));
      }

      // Record the alpha-row locations for infectious states for this participant.
      // Later, when participant p's infection hazard is computed, these indices let
      // the model look up each household member's probability of occupying an
      // infectious latent state at the previous time.
      for(s in inf_states) {
        i_rows[i, s] = n_states * (i-1) + s;
      }

      // `log_sum_exp` is the normalization constant for this participant at time 1.
      llik[i, 1] = log_sum_exp(logalpha[ref,1]);

      // Convert the unnormalized log forward values into normalized probabilities.
      alpha[(n_states * (i-1) + 1):(n_states * (i-1) + n_states), 1] =
        softmax(logalpha[ref,1]);

    } // end participant loop for initialization

    // -------------------------
    // Forward recursion for later days
    // -------------------------
    for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

      // Convert the household-relative time index to the global day index
      // (1..T_global). Different households may start on different days, so
      // the offset hh_tmin[h] - 1 is needed to look up the correct covariate
      // values for this time step.
      int actual_day;
      actual_day = hh_tmin[h] + tt - 1;

      // Update one participant at a time, conditioning on the other household
      // members' filtering distributions from the previous day.
      for(p in 1:hh_size[h]) {
        // For each latent state s, `no_inf_prob[s]` will hold the probability that
        // participant p avoids all within-household infection pressure associated
        // with source state s during this time step.
        array[n_states] real no_inf_prob;

        // `no_hh_inf_prob[j, s]` is the probability that participant p avoids
        // infection from household member j through infectious state s.
        // The eventual `prod(...)` across j assumes independent avoidance across
        // household members conditional on the latent state probabilities.
        matrix[hh_size[h], n_states] no_hh_inf_prob;

        // Row block in `logalpha` for participant p.
        array[n_states] int ref;

        // Previous-time forward probabilities for participant p on the log scale.
        vector[n_states] logalpha_temp;

        // Observation likelihood matrix for participant p at time tt.
        matrix[n_obs_type, n_states] obs;

        // Working copy of the transition-multiplier matrix.
        matrix[n_states, n_states] mult_temp;

        ref = linspaced_int_array(n_states,
                                  n_states * last_lik + n_states * (p-1) + 1,
                                  n_states * last_lik + n_states * (p-1) + n_states);

        logalpha_temp = logalpha[ref, tt-1];

        obs_switch = 0;

        // Check whether the next stored observation belongs to participant p on day tt.
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

        // Counter over the within-household infection-probability columns.
        // It advances only when the current latent state is infectious.
        int ct = 1;

        // Build the probability of avoiding infection from the household.
        for(s in 1:n_states) {
          if(is_in(s, inf_states)) {

            // For a given potential source member j and infectious state s:
            // - with probability alpha[j,s] they are in infectious state s, so
            //   avoidance contributes (1 - ih_prob[p, ct])
            // - otherwise they are not in state s, which contributes 1
            // Summing those two cases gives the marginal avoidance probability
            // from member j via source state s.
            no_hh_inf_prob[,s] =
              to_vector(alpha[i_rows[, s], tt-1]) * (1 - ih_prob[ct][last_lik + p, actual_day])
              + (1 - to_vector(alpha[i_rows[,s], tt-1]));

            ct += 1;

            // A person cannot infect themselves.
            no_hh_inf_prob[p, s] = 1;

          } else {
            // Non-infectious latent states create no infection pressure.
            no_hh_inf_prob[,s] = rep_vector(1, hh_size[h]);
          }

          // Multiply across household members to get the total probability of
          // avoiding infection pressure associated with state s.
          no_inf_prob[s] = prod(no_hh_inf_prob[,s]);
        }

        // Update whichever transition probabilities are being estimated.
        for(m in 1:n_trans_fit) {
          if(sum(source_states[m,]) == 0) {
            // Ordinary non-infection transition:
            // pull a direct probability parameter from `params`.
            trans_temp[trans_index[m, 1],trans_index[m, 2]] = params[param_index[m]];
          } else {
            // Infection-driven transition:
            // combine the relevant no-infection terms across the chosen source
            // states, then combine that with the probability of avoiding infection
            // from outside the household.
            real no_inf;
            no_inf = 1;
            for(s in 1:n_states) {
              if(source_states[m,s] == 1) {
                no_inf = no_inf * no_inf_prob[s];
              }
            }

            // Probability of infection = 1 - probability of avoiding both:
            // - all relevant household infection routes
            // - extra-household infection
            trans_temp[trans_index[m, 1],trans_index[m, 2]] =
              1 - (no_inf * (1 - eh_prob[last_lik + p, actual_day]));
          }
        }

        // Update any transition-splitting multipliers that are being estimated.
        // `mult_temp` is reset from the baseline matrix each time step because
        // some updates depend on the current unmodified multiplier values.
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

        // Apply multiplier-based transition splits element-wise.
        trans_temp = trans_temp .* mult_temp;

        // Recompute the diagonal entries so each "from-state" column still sums to 1.
        for(i in 1:cols(trans_temp)) {
          trans_temp[i,i] = get_diagonal_element(trans_temp, i);
        }

        // Small numerical safeguards before moving to log space.
        trans_temp = replace_zeroes(trans_temp, epsilon);
        trans_temp = normalize_cols(trans_temp);

        // Standard HMM forward step:
        // 1. propagate yesterday's filtering distribution through the transition matrix
        // 2. multiply by today's observation likelihood
        // The multiplication is done on the probability scale and then logged.
        logalpha[ref, tt] = log(trans_temp * exp(logalpha_temp));
        for(k in 1:n_obs_type) {
          logalpha[ref, tt] = logalpha[ref, tt] + to_vector(log(obs[k,]));
        }

        // Normalize to recover the filtering distribution over latent states.
        alpha[(n_states * (p-1) + 1):(n_states * (p-1) + n_states), tt] =
          softmax(logalpha[ref,tt]);

        // Store the normalizing constant, which is the participant/time-point
        // log-likelihood contribution from the forward recursion.
        llik[p, tt] = log_sum_exp(logalpha[ref,tt]);

      } // end participant loop at time tt

      // The likelihood contribution used by this model is the sum of the final-day
      // participant log normalizing constants for the household.
      if(tt == (hh_tmax[h] - hh_tmin[h] + 1)) {
        llik_final[h] = sum(llik[,tt]);
      }

    } // end time recursion
    } // END FORWARD ALGORITHM

  } // end household loop

}

model {

  // Weakly informative priors for the infection-probability regressions.
  // Because these coefficients are on the logit scale, a normal(-3, 3) prior
  // places substantial mass on small probabilities while still allowing wide variation.
  beta_eh ~ normal(-3, 3);
  beta_ih ~ normal(-3, 3);

  // Add the household log-likelihoods that were accumulated during the forward pass.
  target += sum(llik_final);

}
