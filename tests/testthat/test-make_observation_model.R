test_that("can't multiple have observations with same name", {
  expect_error(
    obs_process <- make_observation_model(
      pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      pcr = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    )
  )
})

test_that("length of output list equals number of observations", {
  expect_length(
    make_observation_model(
      a = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
      b = c("S" = 0.01, "I" = 0.01, "R" = 0.8),
      c = c("S" = 0.01, "I" = 0.01, "R" = 0.8)
    ),
    3
  )
})

test_that("entry validation rejects bad observation specifications", {
  # No observations supplied
  expect_error(
    make_observation_model(),
    "observation"
  )

  # Unnamed observation specification
  expect_error(
    make_observation_model(c("S" = 0.05, "I" = 0.95, "R" = 0.05)),
    "named"
  )

  # Probability outside [0, 1]
  expect_error(
    make_observation_model(pcr = c("S" = 0.05, "I" = 1.5, "R" = 0.05)),
    "pcr"
  )

  # Detection probabilities must themselves be named by compartment
  expect_error(
    make_observation_model(pcr = c(0.05, 0.95, 0.05)),
    "pcr"
  )
})

test_that("named entries are estimated and need valid priors", {
  pcr <- c("S" = "1-sp_pcr", "I" = "se_pcr", "R" = "1-sp_pcr")
  priors <- list(se_pcr = c(9, 1), sp_pcr = c(19, 1))

  obs <- make_observation_model(pcr = pcr, priors = priors)
  expect_equal(attr(obs, "priors"), priors)

  # Missing, unused and invalid priors
  expect_error(make_observation_model(pcr = pcr), "Missing priors")
  expect_error(
    make_observation_model(pcr = pcr, priors = c(priors, list(x = c(1, 1)))),
    "unknown"
  )
  expect_error(
    make_observation_model(
      pcr = pcr,
      priors = list(se_pcr = c(9, 0), sp_pcr = c(19, 1))
    ),
    "se_pcr"
  )

  # A name used both directly and as a complement
  expect_error(
    make_observation_model(
      pcr = c("S" = "1-se_pcr", "I" = "se_pcr", "R" = 0.05),
      priors = list(se_pcr = c(9, 1))
    ),
    "not both"
  )

  # A name shared between tests
  expect_error(
    make_observation_model(
      pcr = c("S" = 0.05, "I" = "se", "R" = 0.05),
      igg = c("S" = 0.01, "I" = 0.01, "R" = "se"),
      priors = list(se = c(9, 1))
    ),
    "single test"
  )

  # More than one specificity in a test
  expect_error(
    make_observation_model(
      pcr = c("S" = "1-sp1", "I" = "se", "R" = "1-sp2"),
      priors = list(se = c(9, 1), sp1 = c(19, 1), sp2 = c(19, 1))
    ),
    "at most one specificity"
  )

  # Fixed values given as strings must still be probabilities
  expect_error(
    make_observation_model(
      pcr = c("S" = "1-sp", "I" = "se", "R" = "1.5"),
      priors = list(se = c(9, 1), sp = c(19, 1))
    ),
    "pcr"
  )
})

test_that("stan data matches compartments by name", {
  inf_model <- sir_infection_model()
  data <- sir_subset(2)
  init_probs <- c(1 - 2 * 1e-10, 1e-10, 1e-10)

  ordered <- make_stan_data(
    inf_model,
    sir_observation_model(),
    data,
    init_probs
  )
  reordered <- make_stan_data(
    inf_model,
    make_observation_model(
      pcr = c("R" = 0.05, "S" = 0.05, "I" = 0.95),
      igg = c("I" = 0.01, "R" = 0.8, "S" = 0.01)
    ),
    data,
    init_probs
  )
  expect_equal(reordered$obs_fixed, ordered$obs_fixed)
  expect_equal(ordered$obs_fixed[1, ], c(0.05, 0.95, 0.05))

  # Compartment names must match the infection model
  expect_error(
    make_stan_data(
      inf_model,
      make_observation_model(pcr = c("S" = 0.05, "I" = 0.95, "X" = 0.05)),
      data,
      init_probs
    ),
    "each compartment"
  )
})

test_that("stan data indexes estimated observation parameters", {
  obs_model <- make_observation_model(
    pcr = c("S" = "1-sp_pcr", "I" = "se_pcr", "R" = "res_pcr"),
    igg = c("S" = 0.01, "I" = 0.01, "R" = "se_igg"),
    priors = list(
      se_pcr = c(9, 1),
      sp_pcr = c(19, 1),
      res_pcr = c(1, 1),
      se_igg = c(8, 2)
    )
  )
  dat_stan <- make_stan_data(
    sir_infection_model(),
    obs_model,
    sir_subset(2),
    c(1 - 2 * 1e-10, 1e-10, 1e-10)
  )

  expect_equal(dat_stan$n_spec, 1)
  expect_equal(dat_stan$n_sens, 3)
  expect_equal(dat_stan$obs_index, rbind(c(-1, 1, 2), c(0, 0, 3)))
  expect_equal(dat_stan$obs_fixed[2, ], c(0.01, 0.01, 0))
  # PCR parameters sit above PCR's false-positive rate, IgG has no specificity
  expect_equal(as.vector(dat_stan$sens_spec), c(1, 1, 0))
  expect_equal(rownames(dat_stan$sens_prior), c("se_pcr", "res_pcr", "se_igg"))
  expect_equal(unname(dat_stan$sens_prior[3, ]), c(8, 2))

  # Inits start inside the sens > 1 - spec region
  init <- init_obs_params(dat_stan)
  expect_equal(as.vector(init$spec), 0.95)
  expect_equal(as.vector(init$sens), c(0.9, 0.5, 0.8))
})
