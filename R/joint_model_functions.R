# =============================================================================
# Joint SIRS-SIS model helpers
#
# Two model variants are supported:
#
# --- 6-state model (hmm_tv_cov_reduce_sum_joint.stan) ---
#   Joint state space:
#     1: (S_v, S_b)   2: (S_v, I_b)   3: (I_v, S_b)
#     4: (I_v, I_b)   5: (R_v, S_b)   6: (R_v, I_b)
#   Use: make_joint_obs_model() → make_joint_stan_data() → run_joint_model()
#
# --- 8-state VIP model (hmm_tv_cov_reduce_sum_joint_vip.stan) ---
#   Adds a Viral Inflammatory Phase (VIP) that persists beyond PCR positivity:
#     1: (S_v, S_b, novip)   2: (S_v, I_b, novip)
#     3: (I_v, S_b, vip)     4: (I_v, I_b, vip)
#     5: (R_v, S_b, vip)     6: (R_v, I_b, vip)
#     7: (R_v, S_b, novip)   8: (R_v, I_b, novip)
#   Emission model simplified to 5 parameters (priors fixed in Stan):
#     sens_vir, fpr_vir, sens_bac, delta_bac (VIP boost), fpr_bac
#   Use: make_joint_stan_data_vip() → run_joint_model_vip()
# =============================================================================


# -----------------------------------------------------------------------------
#' @title Make Observation Model for Joint SIRS-SIS Model
#'
#' @description
#' Specifies Beta distribution priors for the observation probabilities
#' of the joint SIRS-SIS model. There are 6 joint states and each test
#' type gets one probability per state.
#'
#' @param ... Named arguments, one per observation type (e.g. \code{viral_pcr},
#'   \code{bac_pcr}). Each argument must be a named list with exactly 6
#'   elements — one per joint state — in the order:
#'   \code{Sv_Sb}, \code{Sv_Ib}, \code{Iv_Sb}, \code{Iv_Ib}, \code{Rv_Sb},
#'   \code{Rv_Ib}. Each element is \code{c(alpha, beta)} giving the Beta
#'   prior hyperparameters; the prior mean is \code{alpha / (alpha + beta)}.
#'
#' @return A list with elements \code{$alpha}, \code{$beta}, \code{$lb} and
#'   \code{$ub} (each a matrix of dimension \code{[n_obs_type, 6]}) and a
#'   character vector \code{$test_names}.
#'
#'   \code{lb} and \code{ub} are auto-derived from the prior mean to break
#'   the HMM label-switching symmetry: states whose prior mean is \eqn{< 0.5}
#'   (a positive observation is unlikely) are bounded to \code{[0, 0.5]}, and
#'   states whose prior mean is \eqn{\geq 0.5} are bounded to \code{[0.5, 1]}.
#'   This forces each chain to land in the same labelling.
#'
#' @examples
#' \dontrun{
#' obs_joint <- make_joint_obs_model(
#'   viral_pcr = list(
#'     Sv_Sb = c(1, 99),   # FPR ~ 0.01
#'     Sv_Ib = c(1, 99),
#'     Iv_Sb = c(95, 5),   # TPR ~ 0.95
#'     Iv_Ib = c(95, 5),
#'     Rv_Sb = c(1, 99),
#'     Rv_Ib = c(1, 99)
#'   ),
#'   bac_pcr = list(
#'     Sv_Sb = c(1, 99),   # FPR ~ 0.01
#'     Sv_Ib = c(95, 5),   # TPR ~ 0.95
#'     Iv_Sb = c(1, 99),
#'     Iv_Ib = c(95, 5),   
#'     Rv_Sb = c(1, 99),
#'     Rv_Ib = c(95, 5)    
#'   )
#' )
#' }
#' @export
make_joint_obs_model <- function(...) {

  state_names <- c("Sv_Sb", "Sv_Ib", "Iv_Sb", "Iv_Ib", "Rv_Sb", "Rv_Ib")
  .dots <- list(...)
  n_obs_type <- length(.dots)

  alpha_mat <- matrix(NA_real_, nrow = n_obs_type, ncol = 6)
  beta_mat  <- matrix(NA_real_, nrow = n_obs_type, ncol = 6)

  for (i in seq_along(.dots)) {
    test_priors <- .dots[[i]]

    if (!is.list(test_priors))
      stop(sprintf(
        "Argument '%s' must be a named list of c(alpha, beta) pairs, one per joint state.",
        names(.dots)[i]
      ))

    if (length(test_priors) != 6)
      stop(sprintf(
        "Argument '%s' must have exactly 6 entries (one per joint state: %s).",
        names(.dots)[i], paste(state_names, collapse = ", ")
      ))

    # Allow unnamed if in the right order, but warn.
    if (is.null(names(test_priors))) {
      warning(sprintf(
        "Argument '%s' is unnamed; assuming order: %s.",
        names(.dots)[i], paste(state_names, collapse = ", ")
      ))
      names(test_priors) <- state_names
    }

    for (sn in state_names) {
      if (!(sn %in% names(test_priors)))
        stop(sprintf("State '%s' missing from argument '%s'.", sn, names(.dots)[i]))
      if (length(test_priors[[sn]]) != 2 || any(test_priors[[sn]] <= 0))
        stop(sprintf(
          "Entry '%s' of '%s' must be c(alpha, beta) with both values > 0.",
          sn, names(.dots)[i]
        ))
    }

    alpha_mat[i, ] <- sapply(test_priors[state_names], `[[`, 1)
    beta_mat[i, ]  <- sapply(test_priors[state_names], `[[`, 2)
  }

  rownames(alpha_mat) <- names(.dots)
  rownames(beta_mat)  <- names(.dots)
  colnames(alpha_mat) <- state_names
  colnames(beta_mat)  <- state_names

  # Derive per-element bounds from the prior mean.  Mean >= 0.5 → "positive
  # state" → bound [0.5, 1]; otherwise → "negative state" → bound [0, 0.5].
  mean_mat <- alpha_mat / (alpha_mat + beta_mat)
  lb_mat   <- ifelse(mean_mat >= 0.5, 0.5, 0)
  ub_mat   <- ifelse(mean_mat >= 0.5, 1,   0.5)
  dimnames(lb_mat) <- dimnames(alpha_mat)
  dimnames(ub_mat) <- dimnames(alpha_mat)

  list(
    alpha      = alpha_mat,
    beta       = beta_mat,
    lb         = lb_mat,
    ub         = ub_mat,
    test_names = names(.dots)
  )
}


# -----------------------------------------------------------------------------
#' @title Prepare Stan Data for Joint SIRS-SIS Model
#'
#' @description
#' Builds the list of inputs required by
#' \code{hmm_tv_cov_reduce_sum_joint.stan} from a combined viral+bacterial
#' observation data frame and time-varying covariate arrays.
#'
#' @param obs_model Output of \link{make_joint_obs_model}.
#' @param data Data frame with one row per (person, time-point). Must contain
#'   columns: \code{hh_id}, \code{part_id}, \code{t}, and one column per
#'   observation type in the same order as \code{obs_model$test_names}.
#'   Observations must be binary, coded as 0 (negative) or 1 (positive); use
#'   \code{NA} for missing. The Stan model assumes binary tests
#'   (\code{n_unique_obs = 2}).
#'   \code{t} is the global day index used to look up time-varying
#'   covariates; it is automatically converted to a per-household relative
#'   index for the HMM forward step.
#' @param obs_cols Character vector of column names in \code{data} for each
#'   observation type, in the same order as \code{obs_model$test_names}.
#' @param init_probs Length-6 numeric vector of initial state probabilities
#'   (will be normalised to sum to 1). Order: Sv_Sb, Sv_Ib, Iv_Sb, Iv_Ib,
#'   Rv_Sb, Rv_Ib.
#' @param ih_cov 3-D array \code{[T_global, N_people, k_ih]} of
#'   intra-household covariates. Used for both viral and bacterial IH
#'   transmission (separate coefficient vectors are estimated).
#' @param eh_cov 3-D array \code{[T_global, N_people, k_eh]} of
#'   extra-household covariates.
#' @param epsilon Small positive constant replacing exact zeros in the
#'   transition matrix. Default \code{1e-10}.
#'
#' @return A named list ready to pass to cmdstanr's \code{$sample()}.
#' @export
make_joint_stan_data <- function(obs_model,
                                 data,
                                 obs_cols,
                                 init_probs,
                                 ih_cov,
                                 eh_cov,
                                 epsilon = 1e-10) {

  # ---- Validate inputs -------------------------------------------------------
  if (!all(c("hh_id", "part_id", "t") %in% names(data)))
    stop("'data' must contain columns: hh_id, part_id, t.")

  if (length(obs_cols) != nrow(obs_model$alpha))
    stop("Length of obs_cols must match number of observation types in obs_model.")

  if (!all(obs_cols %in% names(data)))
    stop("Some obs_cols not found in data.")

  if (length(init_probs) != 6)
    stop("init_probs must be a length-6 vector.")

  if (length(dim(ih_cov)) != 3)
    stop("ih_cov must be a 3-D array [T_global, N_people, k_ih].")
  if (length(dim(eh_cov)) != 3)
    stop("eh_cov must be a 3-D array [T_global, N_people, k_eh].")

  # ---- Sort and index data ---------------------------------------------------
  dat <- data %>%
    dplyr::arrange(hh_id, t, part_id)

  dat$row_id <- seq_len(nrow(dat))

  # Convert global `t` to a per-household relative index (1, 2, ..., T_hh)
  # so the Stan forward-recursion check `t_day_hh[index] == tt` matches
  # observation rows correctly even when households have different start days.
  dat <- dat %>%
    dplyr::group_by(hh_id) %>%
    dplyr::mutate(t_rel = as.integer(t - min(t) + 1L)) %>%
    dplyr::ungroup()

  hh_sum <- dat %>%
    dplyr::group_by(hh_id) %>%
    dplyr::summarise(
      hh_size      = dplyr::n_distinct(part_id),
      hh_start_ind = min(row_id),
      hh_end_ind   = max(row_id),
      hh_tmin      = min(t),       # global day of first obs (covariate offset)
      hh_tmax      = max(t),       # global day of last obs
      obs_per_hh   = dplyr::n(),
      .groups      = "drop"
    )

  # ---- Recode observations (0/1 → 1/2; NA → -1) ----------------------------
  y_mat <- dat %>%
    dplyr::select(dplyr::all_of(obs_cols)) %>%
    as.matrix()

  y_mat <- apply(y_mat, 2, function(col) {
    col_int <- as.integer(col)
    col_int[is.na(col_int)] <- -1L
    col_int + 1L   # 0→1 (neg), 1→2 (pos); -1 stays 0 but checked as missing
  })

  # Missing (-1 before +1 becomes 0 after +1; reset to -1)
  y_mat[y_mat == 0L] <- -1L

  # ---- Dimension checks against covariates ----------------------------------
  T_global  <- dim(ih_cov)[1]
  N_people  <- dim(ih_cov)[2]
  k_ih      <- dim(ih_cov)[3]
  k_eh      <- dim(eh_cov)[3]

  if (dim(eh_cov)[1] != T_global)
    stop("ih_cov and eh_cov must have the same T_global (first dimension).")
  if (T_global < max(hh_sum$hh_tmax))
    stop("T_global must be >= max(hh_tmax).")
  if (N_people != sum(hh_sum$hh_size))
    stop("dim(ih_cov)[2] must equal sum(hh_size).")
  if (dim(eh_cov)[2] != N_people)
    stop("dim(eh_cov)[2] must equal sum(hh_size).")

  # ---- Normalise init_probs -------------------------------------------------
  init_probs <- init_probs / sum(init_probs)

  # ---- Assemble Stan data list -----------------------------------------------
  list(
    n_hh          = max(dat$hh_id),
    hh_size       = hh_sum$hh_size,

    n_obs         = nrow(dat),
    n_obs_type    = length(obs_cols),
    n_unique_obs  = 2L,
    y             = y_mat,
    part_id       = dat$part_id,
    t_day         = dat$t_rel,
    obs_per_hh    = hh_sum$obs_per_hh,
    hh_start_ind  = hh_sum$hh_start_ind,
    hh_end_ind    = hh_sum$hh_end_ind,
    hh_tmin       = hh_sum$hh_tmin,
    hh_tmax       = hh_sum$hh_tmax,

    k_ih          = k_ih,
    k_eh          = k_eh,
    T_global      = T_global,
    x_ih          = ih_cov,
    x_eh          = eh_cov,

    obs_prob_alpha = obs_model$alpha,
    obs_prob_beta  = obs_model$beta,
    obs_lb         = obs_model$lb,
    obs_ub         = obs_model$ub,

    init_probs    = init_probs,
    epsilon       = epsilon
  )
}


# -----------------------------------------------------------------------------
#' @title Run Joint SIRS-SIS Household Transmission Model
#'
#' @description
#' Compiles (if needed) and samples from the joint SIRS-SIS Stan model
#' \code{hmm_tv_cov_reduce_sum_joint.stan}.
#'
#' @param obs_model Output of \link{make_joint_obs_model}.
#' @param data See \link{make_joint_stan_data}.
#' @param obs_cols Character vector of observation column names in \code{data},
#'   in the same order as \code{obs_model$test_names}.
#' @param init_probs Length-6 vector of initial state probabilities.
#' @param ih_cov 3-D array \code{[T_global, N, k_ih]} of IH covariates.
#' @param eh_cov 3-D array \code{[T_global, N, k_eh]} of EH covariates.
#' @param epsilon Small positive constant for transition matrix stability.
#' @param file Path to the joint Stan model file. Defaults to the bundled
#'   \code{hmm_tv_cov_reduce_sum_joint.stan}.
#' @param iter Total number of MCMC iterations per chain (warmup = iter/2).
#' @param chains Number of chains.
#' @param parallel_chains Chains to run in parallel.
#' @param threads_per_chain Threads per chain for reduce_sum parallelism.
#' @param adapt_delta Target acceptance rate.
#' @param max_treedepth Maximum tree depth for NUTS.
#' @param init Optional list of initial values (length = chains). If
#'   \code{NULL} a sensible default is used.
#'
#' @return A \code{CmdStanMCMC} object.
#' @export
run_joint_model <- function(obs_model,
                             data,
                             obs_cols,
                             init_probs,
                             ih_cov,
                             eh_cov,
                             epsilon           = 1e-10,
                             file              = system.file(
                               "stan", "hmm_tv_cov_reduce_sum_joint.stan",
                               package = "hestia"
                             ),
                             iter              = 2000,
                             chains            = 4,
                             parallel_chains   = 4,
                             threads_per_chain = 4,
                             adapt_delta       = 0.9,
                             max_treedepth     = 12,
                             init              = NULL) {

  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("cmdstanr is required. See https://mc-stan.org/cmdstanr/")

  # Build data list
  dat_stan <- make_joint_stan_data(
    obs_model  = obs_model,
    data       = data,
    obs_cols   = obs_cols,
    init_probs = init_probs,
    ih_cov     = ih_cov,
    eh_cov     = eh_cov,
    epsilon    = epsilon
  )

  # Coerce any data frames to matrices (cmdstanr is strict it turns out....)
  dat_stan <- lapply(dat_stan, function(x) if (is.data.frame(x)) as.matrix(x) else x)

  # Default initialisations centred on plausible values
  if (is.null(init)) {
    k_ih <- dat_stan$k_ih
    k_eh <- dat_stan$k_eh
    n_obs_type <- dat_stan$n_obs_type

    init <- rep(list(list(
      logit_gamma_v   = qlogis(1/6),
      logit_gamma_b   = qlogis(1/12),
      logit_rho       = qlogis(1/30),
      beta0_ih_vir    = qlogis(0.03),
      beta0_eh_vir    = qlogis(0.02),
      beta_ih_vir     = rep(0, k_ih),
      beta_eh_vir     = rep(0, k_eh),
      beta0_ih_bac    = qlogis(0.01),
      beta0_eh_bac    = qlogis(0.04),
      beta_ih_bac     = rep(0, k_ih),
      beta_eh_bac     = rep(0, k_eh),
      cross_ih_trans  = 0,
      cross_ih_susc   = 0,
      # Init for obs_raw: back-transform the prior mean from
      # [obs_lb, obs_ub] back to [0, 1], then clamp interior.
      obs_raw         = {
        mean_mat <- obs_model$alpha / (obs_model$alpha + obs_model$beta)
        raw      <- (mean_mat - obs_model$lb) /
                    (obs_model$ub - obs_model$lb)
        matrix(pmin(0.98, pmax(0.02, raw)),
               nrow = n_obs_type, ncol = 6)
      }
    )), chains)
  }

  # Compile model
  mod <- cmdstanr::cmdstan_model(
    file,
    cpp_options = list(stan_threads = TRUE)
  )

  # Sample
  mod$sample(
    data              = dat_stan,
    init              = init,
    iter_warmup       = iter %/% 2,
    iter_sampling     = iter %/% 2,
    chains            = chains,
    parallel_chains   = parallel_chains,
    threads_per_chain = threads_per_chain,
    adapt_delta       = adapt_delta,
    max_treedepth     = max_treedepth,
    refresh           = 100
  )
}


# =============================================================================
# 8-state VIP model helpers
# =============================================================================

# -----------------------------------------------------------------------------
#' @title Prepare Stan Data for 8-State VIP Joint Model
#'
#' @description
#' Builds the list of inputs required by
#' \code{hmm_tv_cov_reduce_sum_joint_vip.stan}. The VIP model has a simplified
#' emission structure (5 parameters with priors fixed in Stan), so no
#' \code{obs_model} argument is needed.
#'
#' @param data Data frame with columns \code{hh_id}, \code{part_id}, \code{t},
#'   and one column per observation type. Observations coded 0/1/NA.
#' @param obs_cols Character vector of observation column names, in order
#'   (viral test first, bacterial test second).
#' @param init_probs Length-8 numeric vector of initial state probabilities
#'   (will be normalised). Order: Sv_Sb_novip, Sv_Ib_novip, Iv_Sb_vip,
#'   Iv_Ib_vip, Rv_Sb_vip, Rv_Ib_vip, Rv_Sb_novip, Rv_Ib_novip.
#' @param ih_cov 3-D array \code{[T_global, N_people, k_ih]} of IH covariates.
#' @param eh_cov 3-D array \code{[T_global, N_people, k_eh]} of EH covariates.
#' @param epsilon Small positive constant for transition matrix stability.
#'
#' @return A named list ready to pass to cmdstanr's \code{$sample()}.
#' @export
make_joint_stan_data_vip <- function(data,
                                     obs_cols,
                                     init_probs,
                                     ih_cov,
                                     eh_cov,
                                     epsilon = 1e-10) {

  # ---- Validate inputs -------------------------------------------------------
  if (!all(c("hh_id", "part_id", "t") %in% names(data)))
    stop("'data' must contain columns: hh_id, part_id, t.")

  if (!all(obs_cols %in% names(data)))
    stop("Some obs_cols not found in data.")

  if (length(init_probs) != 8)
    stop("init_probs must be a length-8 vector (one per VIP state).")

  if (length(dim(ih_cov)) != 3)
    stop("ih_cov must be a 3-D array [T_global, N_people, k_ih].")
  if (length(dim(eh_cov)) != 3)
    stop("eh_cov must be a 3-D array [T_global, N_people, k_eh].")

  # ---- Sort and index data ---------------------------------------------------
  dat <- data %>%
    dplyr::arrange(hh_id, t, part_id)

  dat$row_id <- seq_len(nrow(dat))

  dat <- dat %>%
    dplyr::group_by(hh_id) %>%
    dplyr::mutate(t_rel = as.integer(t - min(t) + 1L)) %>%
    dplyr::ungroup()

  hh_sum <- dat %>%
    dplyr::group_by(hh_id) %>%
    dplyr::summarise(
      hh_size      = dplyr::n_distinct(part_id),
      hh_start_ind = min(row_id),
      hh_end_ind   = max(row_id),
      hh_tmin      = min(t),
      hh_tmax      = max(t),
      obs_per_hh   = dplyr::n(),
      .groups      = "drop"
    )

  # ---- Recode observations (0/1 -> 1/2; NA -> -1) ---------------------------
  y_mat <- dat %>%
    dplyr::select(dplyr::all_of(obs_cols)) %>%
    as.matrix()

  y_mat <- apply(y_mat, 2, function(col) {
    col_int <- as.integer(col)
    col_int[is.na(col_int)] <- -1L
    col_int + 1L
  })
  y_mat[y_mat == 0L] <- -1L

  # ---- Covariate dimension checks -------------------------------------------
  T_global <- dim(ih_cov)[1]
  N_people <- dim(ih_cov)[2]
  k_ih     <- dim(ih_cov)[3]
  k_eh     <- dim(eh_cov)[3]

  if (dim(eh_cov)[1] != T_global)
    stop("ih_cov and eh_cov must have the same T_global (first dimension).")
  if (T_global < max(hh_sum$hh_tmax))
    stop("T_global must be >= max(hh_tmax).")
  if (N_people != sum(hh_sum$hh_size))
    stop("dim(ih_cov)[2] must equal sum(hh_size).")
  if (dim(eh_cov)[2] != N_people)
    stop("dim(eh_cov)[2] must equal sum(hh_size).")

  # ---- Normalise init_probs -------------------------------------------------
  init_probs <- init_probs / sum(init_probs)

  # ---- Assemble Stan data list ----------------------------------------------
  # Note: no obs_prob_alpha / beta / lb / ub — emission priors are in Stan.
  list(
    n_hh          = max(dat$hh_id),
    hh_size       = hh_sum$hh_size,

    n_obs         = nrow(dat),
    n_obs_type    = length(obs_cols),
    n_unique_obs  = 2L,
    y             = y_mat,
    part_id       = dat$part_id,
    t_day         = dat$t_rel,
    obs_per_hh    = hh_sum$obs_per_hh,
    hh_start_ind  = hh_sum$hh_start_ind,
    hh_end_ind    = hh_sum$hh_end_ind,
    hh_tmin       = hh_sum$hh_tmin,
    hh_tmax       = hh_sum$hh_tmax,

    k_ih          = k_ih,
    k_eh          = k_eh,
    T_global      = T_global,
    x_ih          = ih_cov,
    x_eh          = eh_cov,

    init_probs    = init_probs,
    epsilon       = epsilon
  )
}


# -----------------------------------------------------------------------------
#' @title Run 8-State VIP Joint Household Transmission Model
#'
#' @description
#' Compiles (if needed) and samples from the VIP joint Stan model
#' \code{hmm_tv_cov_reduce_sum_joint_vip.stan}. The emission model is
#' simplified to 5 parameters (sens_vir, fpr_vir, sens_bac, delta_bac,
#' fpr_bac) with priors fixed in Stan; no \code{obs_model} is required.
#'
#' @param data See \link{make_joint_stan_data_vip}.
#' @param obs_cols Character vector of observation column names (viral first,
#'   bacterial second).
#' @param init_probs Length-8 vector of initial state probabilities. Order:
#'   Sv_Sb_novip, Sv_Ib_novip, Iv_Sb_vip, Iv_Ib_vip, Rv_Sb_vip, Rv_Ib_vip,
#'   Rv_Sb_novip, Rv_Ib_novip.
#' @param ih_cov 3-D array \code{[T_global, N, k_ih]} of IH covariates.
#' @param eh_cov 3-D array \code{[T_global, N, k_eh]} of EH covariates.
#' @param epsilon Small positive constant for transition matrix stability.
#' @param file Path to the VIP Stan model file.
#' @param iter Total MCMC iterations per chain (warmup = iter/2).
#' @param chains Number of chains.
#' @param parallel_chains Chains to run in parallel.
#' @param threads_per_chain Threads per chain for reduce_sum parallelism.
#' @param adapt_delta Target acceptance rate.
#' @param max_treedepth Maximum NUTS tree depth.
#' @param init Optional list of initial values (length = chains). If
#'   \code{NULL} a sensible default is constructed.
#'
#' @return A \code{CmdStanMCMC} object.
#' @export
run_joint_model_vip <- function(data,
                                obs_cols,
                                init_probs,
                                ih_cov,
                                eh_cov,
                                epsilon           = 1e-10,
                                file              = system.file(
                                  "stan", "hmm_tv_cov_reduce_sum_joint_vip.stan",
                                  package = "hestia"
                                ),
                                iter              = 2000,
                                chains            = 4,
                                parallel_chains   = 4,
                                threads_per_chain = 4,
                                adapt_delta       = 0.9,
                                max_treedepth     = 12,
                                init              = NULL) {

  if (!requireNamespace("cmdstanr", quietly = TRUE))
    stop("cmdstanr is required. See https://mc-stan.org/cmdstanr/")

  # Build data list
  dat_stan <- make_joint_stan_data_vip(
    data       = data,
    obs_cols   = obs_cols,
    init_probs = init_probs,
    ih_cov     = ih_cov,
    eh_cov     = eh_cov,
    epsilon    = epsilon
  )

  dat_stan <- lapply(dat_stan, function(x) if (is.data.frame(x)) as.matrix(x) else x)

  # Default initialisations
  if (is.null(init)) {
    k_ih <- dat_stan$k_ih
    k_eh <- dat_stan$k_eh

    init <- rep(list(list(
      # Viral dynamics
      logit_gamma_v   = qlogis(1/6),     # ~6-day viral infection
      logit_rho       = qlogis(1/30),    # ~30-day immunity window
      beta0_ih_vir    = qlogis(0.03),
      beta0_eh_vir    = qlogis(0.02),
      beta_ih_vir     = rep(0, k_ih),
      beta_eh_vir     = rep(0, k_eh),

      # Bacterial dynamics
      logit_gamma_b   = qlogis(1/12),    # ~12-day bacterial carriage
      beta0_ih_bac    = qlogis(0.01),
      beta0_eh_bac    = qlogis(0.04),
      beta_ih_bac     = rep(0, k_ih),
      beta_eh_bac     = rep(0, k_eh),

      # VIP waning: ~21-day inflammatory window (informed by KM curve)
      logit_lambda    = qlogis(1/21),

      # Cross-immunity: start at no effect
      cross_ih_trans  = 0,
      cross_ih_susc   = 0,

      # Emission model (logit scale; constraints enforced in parameters block)
      # logit_sens_* >= 0 (sens > 0.5), logit_fpr_* <= 0 (fpr < 0.5)
      logit_sens_vir  = qlogis(0.90),    # ~0.90 viral sensitivity
      logit_fpr_vir   = qlogis(0.02),    # ~0.02 viral FPR
      logit_sens_bac  = qlogis(0.90),    # ~0.90 bacterial sensitivity
      delta_bac       = 0,               # no initial VIP detectability boost
      logit_fpr_bac   = qlogis(0.02)     # ~0.02 bacterial FPR
    )), chains)
  }

  # Compile model
  mod <- cmdstanr::cmdstan_model(
    file,
    cpp_options = list(stan_threads = TRUE)
  )

  # Sample
  mod$sample(
    data              = dat_stan,
    init              = init,
    iter_warmup       = iter %/% 2,
    iter_sampling     = iter %/% 2,
    chains            = chains,
    parallel_chains   = parallel_chains,
    threads_per_chain = threads_per_chain,
    adapt_delta       = adapt_delta,
    max_treedepth     = max_treedepth,
    refresh           = 100
  )
}
