// =============================================================================
// Joint SIRS-SIS Hidden Markov Model — 8-state Viral Inflammatory Phase (VIP)
//
// State space (8 states):
//   1: (S_v, S_b, novip)   Viral susceptible,  Bacterial susceptible
//   2: (S_v, I_b, novip)   Viral susceptible,  Bacterial infected
//   3: (I_v, S_b, vip)     Viral infected,     Bacterial susceptible  [VIP always active]
//   4: (I_v, I_b, vip)     Viral infected,     Bacterial infected     [VIP always active]
//   5: (R_v, S_b, vip)     Viral recovered,    Bacterial susceptible  [VIP persists post-clearance]
//   6: (R_v, I_b, vip)     Viral recovered,    Bacterial infected     [VIP persists post-clearance]
//   7: (R_v, S_b, novip)   Viral recovered,    Bacterial susceptible
//   8: (R_v, I_b, novip)   Viral recovered,    Bacterial infected
//
// Viral Inflammatory Phase (VIP):
//   VIP is active throughout Iv infection (states 3-4) and persists into early
//   Rv (states 5-6), waning to novip at rate lambda.  Sv is always novip
//   (inflammation resolves long before immunity wanes fully).
//
//   VIP drives three cross-pathogen effects:
//     cross_ih_trans — VIP+Ib transmitters (states 4, 6) shed bacteria at enhanced rate
//     cross_ih_susc  — VIP+Sb recipients   (states 3, 5) have enhanced bac susceptibility
//     delta_bac      — VIP+Ib states       (4, 6) have enhanced bac PCR detectability
//
// Simplified emission model (5 parameters):
//   sens_vir        P(viral PCR+ | Iv)          -> states 3, 4
//   fpr_vir         P(viral PCR+ | not Iv)     -> states 1, 2, 5, 6, 7, 8
//   sens_bac        P(bac  PCR+ | Ib, novip)   -> states 2, 8
//   sens_bac_coinf  P(bac  PCR+ | Ib, VIP)     -> states 4, 6
//                   = inv_logit(logit(sens_bac) -> delta_bac)
//   fpr_bac         P(bac  PCR+ | Sb)            -> states 1, 3, 5, 7
//
//   Label-switching resolved by parameter constraints:
//     logit_sens_* >= 0  →  sens > 0.5
//     logit_fpr_*  <= 0  →  fpr  < 0.5
// =============================================================================

functions {

  // ---------------------------------------------------------------------------
  // Helper: diagonal entry so that column i sums to 1.
  // ---------------------------------------------------------------------------
  real get_diagonal_element(matrix m, int i) {
    real out = 1;
    for (j in 1:rows(m))
      if (j != i) out -= m[j, i];
    return out;
  }

  // ---------------------------------------------------------------------------
  // Helper: replace exact zeros with epsilon (avoids log(0)).
  // ---------------------------------------------------------------------------
  matrix replace_zeroes(matrix m, real epsilon) {
    matrix[rows(m), cols(m)] out = m;
    for (i in 1:rows(m))
      for (j in 1:cols(m))
        if (m[i, j] == 0) out[i, j] = epsilon;
    return out;
  }

  // ---------------------------------------------------------------------------
  // Helper: normalise each column of a matrix to sum to 1.
  // ---------------------------------------------------------------------------
  matrix normalize_cols(matrix m) {
    matrix[rows(m), cols(m)] out;
    for (i in 1:cols(m)) out[, i] = m[, i] / sum(m[, i]);
    return out;
  }

  // ---------------------------------------------------------------------------
  // Partial log-likelihood for reduce_sum.
  //
  // Runs the 8-state HMM forward algorithm for households [start, end] and
  // returns their summed log-likelihood.
  //
  // Force-of-infection terms:
  //   no_bac_inf_base  — escape probability for novip recipient (states 1, 7)
  //   no_bac_inf_susc  — escape probability for VIP   recipient (states 3, 5)
  //   VIP transmitters  — states 4 (Iv_Ib_vip) and 6 (Rv_Ib_vip)
  //   Base transmitters — states 2 (Sv_Ib)      and 8 (Rv_Ib_novip)
  // ---------------------------------------------------------------------------
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
    matrix ih_prob_vir,           // viral IH:                    [N, T]
    matrix ih_prob_bac,           // bacterial IH, base:          [N, T]
    matrix ih_prob_bac_trans,     // bacterial IH, VIP transmitter:[N, T]
    matrix ih_prob_bac_susc,      // bacterial IH, VIP recipient: [N, T]
    matrix ih_prob_bac_both,      // bacterial IH, VIP both:      [N, T]
    matrix eh_prob_vir,           // viral EH:                    [N, T]
    matrix eh_prob_bac,           // bacterial EH, base:          [N, T]
    matrix eh_prob_bac_susc,      // bacterial EH, VIP recipient: [N, T]
    // ---- Fixed daily transition rates ----
    real gamma_v,
    real gamma_b,
    real rho,
    real lambda,                  // VIP waning rate: Rv_vip -> Rv_novip
    // ---- Observation model and initialisation ----
    array[] matrix obs_prob,      // [n_obs_type] matrix[n_unique_obs, 8]
    vector init_probs,            // length 8
    real epsilon
  ) {

    real llik_sum = 0;

    for (h in start:end) {

      int T_hh = max(hh_tmax) - min(hh_tmin) + 1;

      // Forward probability matrices (rows = person * 8 states, cols = time).
      matrix[hh_size[h] * 8, T_hh] alpha;
      matrix[hh_size[h] * 8, T_hh] logalpha_hh;

      // Per-participant, per-time log-normalising constants.
      matrix[hh_size[h], T_hh] llik;

      // Household slices of observation arrays.
      array[obs_per_hh[h], n_obs_type] int y_hh;
      array[obs_per_hh[h]] int part_id_hh;
      array[obs_per_hh[h]] int t_day_hh;

      // Row-index lookup: i_rows[i, s] = row in alpha for person i, state s.
      array[hh_size[h], 8] int i_rows;

      // Offset into global person-level probability matrices.
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

        // ----------
        // Initialisation at t = 1
        // ----------
        for (i in 1:hh_size[h]) {

          array[8] int ref = linspaced_int_array(8, 8*(i-1)+1, 8*i);
          matrix[n_obs_type, 8] obs;

          obs_switch = 0;
          if (t_day_hh[index] == 1 && part_id_hh[index] == i) obs_switch = 1;

          if (obs_switch == 1) {
            for (k in 1:n_obs_type) {
              if (y_hh[index, k] != -1)
                obs[k, ] = obs_prob[k][y_hh[index, k], ];
              else
                obs[k, ] = rep_row_vector(1, 8);
            }
          } else {
            obs = rep_matrix(1, n_obs_type, 8);
          }

          if (obs_switch == 1) index = min(index + 1, obs_per_hh[h]);

          logalpha_hh[ref, 1] = log(init_probs);
          for (k in 1:n_obs_type)
            logalpha_hh[ref, 1] += to_vector(log(obs[k, ]));

          for (s in 1:8) i_rows[i, s] = 8*(i-1) + s;

          llik[i, 1] = log_sum_exp(logalpha_hh[ref, 1]);
          alpha[(8*(i-1)+1):(8*i), 1] = softmax(logalpha_hh[ref, 1]);

        } // end init loop

        // ----------
        // Forward recursion for t = 2, ..., T_hh
        // ----------
        for (tt in 2:(hh_tmax[h] - hh_tmin[h] + 1)) {

          int actual_day = hh_tmin[h] + tt - 1;

          for (p in 1:hh_size[h]) {

            array[8] int ref = linspaced_int_array(8, 8*(p-1)+1, 8*p);
            vector[8] logalpha_temp = logalpha_hh[ref, tt - 1];
            matrix[n_obs_type, 8] obs;
            matrix[8, 8] trans_temp;

            // Observation at (p, tt)
            obs_switch = 0;
            if (t_day_hh[index] == tt && part_id_hh[index] == p) obs_switch = 1;

            if (obs_switch == 1) {
              for (k in 1:n_obs_type) {
                if (y_hh[index, k] != -1)
                  obs[k, ] = obs_prob[k][y_hh[index, k], ];
                else
                  obs[k, ] = rep_row_vector(1, 8);
              }
            } else {
              obs = rep_matrix(1, n_obs_type, 8);
            }

            if (obs_switch == 1) index = min(index + 1, obs_per_hh[h]);

            // --------
            // Household force of infection for person p
            //
            //   no_vir_inf      P(p escapes viral infection from household)
            //   no_bac_inf_base P(p escapes bac infection), novip recipient [states 1, 7]
            //   no_bac_inf_susc P(p escapes bac infection), VIP  recipient [states 3, 5]
            //
            //   VIP transmitters : states 4 (Iv_Ib_vip) and 6 (Rv_Ib_vip)
            //   Base transmitters: states 2 (Sv_Ib)      and 8 (Rv_Ib_novip)
            // --------
            real no_vir_inf      = 1.0;
            real no_bac_inf_base = 1.0;
            real no_bac_inf_susc = 1.0;

            for (q in 1:hh_size[h]) {
              if (q != p) {

                // Viral source: states 3 (Iv_Sb_vip) and 4 (Iv_Ib_vip)
                real prob_Iv = alpha[i_rows[q, 3], tt-1]
                             + alpha[i_rows[q, 4], tt-1];
                no_vir_inf *= 1.0 - prob_Iv * ih_prob_vir[last_lik + p, actual_day];

                // Bacterial transmitters
                real prob_bac_base  = alpha[i_rows[q, 2], tt-1]   // Sv_Ib
                                    + alpha[i_rows[q, 8], tt-1];   // Rv_Ib_novip
                real prob_vip_trans = alpha[i_rows[q, 4], tt-1]   // Iv_Ib_vip
                                    + alpha[i_rows[q, 6], tt-1];   // Rv_Ib_vip

                // Novip recipient (states 1, 7): base susceptibility
                no_bac_inf_base *= 1.0
                  - prob_bac_base  * ih_prob_bac[last_lik + p, actual_day]
                  - prob_vip_trans * ih_prob_bac_trans[last_lik + p, actual_day];

                // VIP recipient (states 3, 5): enhanced susceptibility
                no_bac_inf_susc *= 1.0
                  - prob_bac_base  * ih_prob_bac_susc[last_lik + p, actual_day]
                  - prob_vip_trans * ih_prob_bac_both[last_lik + p, actual_day];
              }
            }

            // --------
            // Daily infection probabilities
            // --------
            // Viral: same for Sv_Sb (1) and Sv_Ib (2)
            real p_vir = 1.0 - no_vir_inf * (1.0 - eh_prob_vir[last_lik + p, actual_day]);

            // Bacterial, novip recipient: states 1 (Sv_Sb) and 7 (Rv_Sb_novip)
            real p_bac_base = 1.0 - no_bac_inf_base * (1.0 - eh_prob_bac[last_lik + p, actual_day]);

            // Bacterial, VIP recipient: states 3 (Iv_Sb_vip) and 5 (Rv_Sb_vip)
            real p_bac_susc = 1.0 - no_bac_inf_susc * (1.0 - eh_prob_bac_susc[last_lik + p, actual_day]);

            // --------
            // 8x8 transition matrix (columns = source state).
            // One-event-per-day approximation: simultaneous transitions = 0.
            // --------
            trans_temp = rep_matrix(0.0, 8, 8);

            // Viral acquisition: Sv -> Iv_vip
            trans_temp[3, 1] = p_vir;       // Sv_Sb_novip -> Iv_Sb_vip
            trans_temp[4, 2] = p_vir;       // Sv_Ib_novip -> Iv_Ib_vip

            // Bacterial acquisition
            trans_temp[2, 1] = p_bac_base;  // Sv_Sb_novip  -> Sv_Ib_novip [base]
            trans_temp[4, 3] = p_bac_susc;  // Iv_Sb_vip    -> Iv_Ib_vip   [VIP]
            trans_temp[6, 5] = p_bac_susc;  // Rv_Sb_vip    -> Rv_Ib_vip   [VIP]
            trans_temp[8, 7] = p_bac_base;  // Rv_Sb_novip  -> Rv_Ib_novip [base]

            // Bacterial recovery: Ib -> Sb
            trans_temp[1, 2] = gamma_b;     // Sv_Ib_novip  -> Sv_Sb_novip
            trans_temp[3, 4] = gamma_b;     // Iv_Ib_vip    -> Iv_Sb_vip
            trans_temp[5, 6] = gamma_b;     // Rv_Ib_vip    -> Rv_Sb_vip
            trans_temp[7, 8] = gamma_b;     // Rv_Ib_novip  -> Rv_Sb_novip

            // Viral recovery: Iv -> Rv_vip (VIP persists into recovery)
            trans_temp[5, 3] = gamma_v;     // Iv_Sb_vip    -> Rv_Sb_vip
            trans_temp[6, 4] = gamma_v;     // Iv_Ib_vip    -> Rv_Ib_vip

            // VIP waning: Rv_vip -> Rv_novip
            trans_temp[7, 5] = lambda;      // Rv_Sb_vip    -> Rv_Sb_novip
            trans_temp[8, 6] = lambda;      // Rv_Ib_vip    -> Rv_Ib_novip

            // Waning immunity: Rv -> Sv (both VIP and novip)
            trans_temp[1, 5] = rho;         // Rv_Sb_vip    -> Sv_Sb_novip
            trans_temp[2, 6] = rho;         // Rv_Ib_vip    -> Sv_Ib_novip
            trans_temp[1, 7] = rho;         // Rv_Sb_novip  -> Sv_Sb_novip
            trans_temp[2, 8] = rho;         // Rv_Ib_novip  -> Sv_Ib_novip

            // Diagonal: 1 - sum of off-diagonals; guard against warmup overflow.
            for (s in 1:8) {
              real d = get_diagonal_element(trans_temp, s);
              trans_temp[s, s] = d > 0 ? d : epsilon;
            }

            trans_temp = replace_zeroes(trans_temp, epsilon);
            trans_temp = normalize_cols(trans_temp);

            // Log forward step
            logalpha_hh[ref, tt] = log(trans_temp * exp(logalpha_temp));
            for (k in 1:n_obs_type)
              logalpha_hh[ref, tt] += to_vector(log(obs[k, ]));

            alpha[(8*(p-1)+1):(8*p), tt] = softmax(logalpha_hh[ref, tt]);
            llik[p, tt] = log_sum_exp(logalpha_hh[ref, tt]);

          } // end participant loop

        } // end time recursion

      } // END FORWARD ALGORITHM

      llik_sum += sum(llik[, hh_tmax[h] - hh_tmin[h] + 1]);

    } // end household loop

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
  int n_obs_type;       // number of test types (viral + bacterial = 2)
  int n_unique_obs;     // 2 for binary tests
  array[n_obs, n_obs_type] int y;   // 1 = negative, 2 = positive, -1 = missing
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

  // ---- Initialisation ----
  simplex[8] init_probs;
  real<lower=0> epsilon;

}


// =============================================================================
parameters {

  // ---- Viral dynamics ----
  real logit_gamma_v;
  real logit_rho;
  real beta0_ih_vir;
  real beta0_eh_vir;
  vector[k_ih] beta_ih_vir;
  vector[k_eh] beta_eh_vir;

  // ---- Bacterial dynamics ----
  real logit_gamma_b;
  real beta0_ih_bac;
  real beta0_eh_bac;
  vector[k_ih] beta_ih_bac;
  vector[k_eh] beta_eh_bac;

  // ---- VIP waning ----
  real logit_lambda;             // daily Rv_vip -> Rv_novip transition probability

  // ---- Cross-immunity (logit-scale additive effects) ----
  real cross_ih_trans;           // VIP+Ib transmitters: enhanced bac transmissibility
  real cross_ih_susc;            // VIP+Sb recipients:   enhanced bac susceptibility

  // ---- Emission model ----
  // Label-switching resolved by constraints:
  //   logit_sens >= 0  =>  sens > 0.5
  //   logit_fpr  <= 0  =>  fpr  < 0.5
  real<lower=0> logit_sens_vir;  // viral sensitivity
  real<upper=0> logit_fpr_vir;   // viral false-positive rate (1 - specificity)
  real<lower=0> logit_sens_bac;  // bacterial sensitivity (novip+Ib)
  real          delta_bac;       // logit-scale boost to bac sensitivity in VIP+Ib states
  real<upper=0> logit_fpr_bac;   // bacterial false-positive rate

}


// =============================================================================
transformed parameters {

  // Rates on probability scale
  real<lower=0, upper=1> gamma_v = inv_logit(logit_gamma_v);
  real<lower=0, upper=1> gamma_b = inv_logit(logit_gamma_b);
  real<lower=0, upper=1> rho     = inv_logit(logit_rho);
  real<lower=0, upper=1> lambda  = inv_logit(logit_lambda);

  // Emission probabilities
  real<lower=0.5, upper=1> sens_vir        = inv_logit(logit_sens_vir);
  real<lower=0,   upper=0.5> fpr_vir       = inv_logit(logit_fpr_vir);
  real<lower=0.5, upper=1> sens_bac        = inv_logit(logit_sens_bac);
  real<lower=0,   upper=1>   sens_bac_coinf = inv_logit(logit_sens_bac + delta_bac);
  real<lower=0,   upper=0.5> fpr_bac       = inv_logit(logit_fpr_bac);

  // Pre-computed person×time infection probability matrices.
  // cross_ih_trans / cross_ih_susc now apply to all VIP states (not just Iv).
  matrix[sum(hh_size), T_global] ih_prob_vir;
  matrix[sum(hh_size), T_global] ih_prob_bac;
  matrix[sum(hh_size), T_global] ih_prob_bac_trans;  // VIP transmitter
  matrix[sum(hh_size), T_global] ih_prob_bac_susc;   // VIP recipient
  matrix[sum(hh_size), T_global] ih_prob_bac_both;   // VIP transmitter + recipient
  matrix[sum(hh_size), T_global] eh_prob_vir;
  matrix[sum(hh_size), T_global] eh_prob_bac;
  matrix[sum(hh_size), T_global] eh_prob_bac_susc;   // VIP recipient EH

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

  // Emission probability matrices: obs_prob[k][row, state]
  //   row 1 = P(negative test | state)
  //   row 2 = P(positive test | state)
  //
  // Viral PCR (test 1):
  //   Positive iff in Iv: states 3 (Iv_Sb_vip) and 4 (Iv_Ib_vip)
  //
  // Bacterial PCR (test 2):
  //   VIP+Ib  (states 4, 6) -> sens_bac_coinf   [enhanced by viral inflammation]
  //   novip+Ib(states 2, 8) -> sens_bac
  //   Sb      (states 1,3,5,7) -> fpr_bac
  array[n_obs_type] matrix[n_unique_obs, 8] obs_prob;

  {
    // State order: 1=Sv_Sb, 2=Sv_Ib, 3=Iv_Sb_vip, 4=Iv_Ib_vip,
    //              5=Rv_Sb_vip, 6=Rv_Ib_vip, 7=Rv_Sb_novip, 8=Rv_Ib_novip

    row_vector[8] p_vir_pos = [fpr_vir,  fpr_vir,  sens_vir, sens_vir,
                                fpr_vir,  fpr_vir,  fpr_vir,  fpr_vir];
    obs_prob[1][2, ] = p_vir_pos;
    obs_prob[1][1, ] = 1 - p_vir_pos;

    row_vector[8] p_bac_pos = [fpr_bac,        sens_bac,       fpr_bac, sens_bac_coinf,
                                fpr_bac,        sens_bac_coinf, fpr_bac, sens_bac];
    obs_prob[2][2, ] = p_bac_pos;
    obs_prob[2][1, ] = 1 - p_bac_pos;
  }

}


// =============================================================================
model {

  // ---- Priors ----

  // Viral and bacterial recovery rates
  logit_gamma_v ~ normal(-2, 1);   // prior mean ~0.12 per day (~8 day infection)
  logit_gamma_b ~ normal(-2, 1);   // prior mean ~0.12 per day
  logit_rho     ~ normal(-3, 1);   // prior mean ~0.05 per day (~20 day immunity)

  // VIP waning: prior centred on ~3-week inflammatory window (informed by KM data)
  logit_lambda  ~ normal(-3, 1);   // prior mean ~0.05 per day (~20 day VIP window)

  // Baseline IH/EH log-odds
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

  // Emission model
  logit_sens_vir ~ normal(3, 1);   // prior sens_vir ~= 0.95
  logit_fpr_vir  ~ normal(-4, 1);  // prior fpr_vir  ~= 0.02
  logit_sens_bac ~ normal(3, 1);   // prior sens_bac ~= 0.95
  delta_bac      ~ normal(0, 1);   // logit-scale VIP boost to bac detectability
  logit_fpr_bac  ~ normal(-4, 1);  // prior fpr_bac  ~= 0.02

  // ---- Parallelised forward algorithm ----
  array[n_hh] int hh_indices;
  for (h in 1:n_hh) hh_indices[h] = h;

  target += reduce_sum(
    partial_log_lik,
    hh_indices,
    1,
    hh_size,
    n_obs_type, n_unique_obs,
    y, part_id, t_day,
    obs_per_hh, hh_start_ind, hh_end_ind, hh_tmin, hh_tmax,
    T_global,
    ih_prob_vir, ih_prob_bac, ih_prob_bac_trans, ih_prob_bac_susc, ih_prob_bac_both,
    eh_prob_vir, eh_prob_bac, eh_prob_bac_susc,
    gamma_v, gamma_b, rho, lambda,
    obs_prob, init_probs, epsilon
  );

}


// =============================================================================
generated quantities {

  // Per-household log-likelihoods for LOO-CV.
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
        gamma_v, gamma_b, rho, lambda,
        obs_prob, init_probs, epsilon
      );
    }
  }

}
