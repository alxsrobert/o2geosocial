context("Test function custom_likelihoods")

## test timing_infections ##
test_that("test timing_infections: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  times <- 0:4
  alpha <- c(NA,rep(1,4))
  w <- c(.1, .2, .5, .2, .1)
  data <- outbreaker_data(dates = times, w_dens = w)
  config <- create_config(data = data, init_tree = alpha)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_timing_infections(data, param)
  out_few_cases <- cpp_ll_timing_infections(data, param, few_cases)
  out_rnd_cases <- cpp_ll_timing_infections(data, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -6.59584881763949)
  expect_equal(out_few_cases, -2.4932054526027)
})

## test timing_sampling ##
test_that("test timing_sampling: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  times <- 0:4
  alpha <- c(NA,rep(1,4))
  f <- c(.1, .2, .5, .2, .1)
  data <- outbreaker_data(dates = times +  c(1, 1, 2, 3, 4), f_dens = f)
  config <- create_config(data = data, init_t_inf = times, init_tree = alpha)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_timing_sampling(data, param)
  out_few_cases <- cpp_ll_timing_sampling(data, param, few_cases)
  out_rnd_cases <- cpp_ll_timing_sampling(data, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -8.300597)
  expect_equal(out_few_cases, -4.1979535)
})


## test age ##
test_that("test age: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  age_dens <- toy_outbreak_short$age_contact
  times <- 0:4
  alpha <- c(NA,rep(1,4))
  age <- c(1, 3, 3, 5, 1)
  data <- outbreaker_data(dates = times, age_group = age,
                          a_dens = age_dens)
  config <- create_config(data = data, init_tree = alpha)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_age(data, param)
  out_few_cases <- cpp_ll_age(data, param, few_cases)
  out_rnd_cases <- cpp_ll_age(data, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -11.9266839)
  expect_equal(out_few_cases, -6.8121909)
})


## test reporting ##
test_that("test reporting: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  times <- 0:4
  alpha <- c(NA,rep(1,4))
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .3, .3, .2, .1)
  kappa <- c(NA, 1, 1, 2, 2)
  data <- outbreaker_data(dates = times +  c(1, 1, 2, 3, 4), f_dens = f,
                          w_dens = w)
  config <- create_config(data = data, init_tree = alpha, init_kappa = kappa)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_reporting(data, param)
  out_few_cases <- cpp_ll_reporting(data, param, few_cases)
  out_rnd_cases <- cpp_ll_reporting(data, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -5.0266122)
  expect_equal(out_few_cases, -2.5133061)
})


## test space ##
test_that("test space: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  times <- 0:4
  alpha <- c(NA,rep(1,4))
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 60, 10, 60, 0, 15, 10, 15, 0), ncol = 3)

  a <- .7
  b <- .1
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  s_dens <- population ** b * exp(-b*distance)
  data <- outbreaker_data(dates = times, region = regions,
                          population = population,distance = distance,  
                          s_dens = s_dens)
  config <- create_config(data = data, init_tree = alpha, init_a = a, init_b = b,
                          move_a = FALSE, move_b = FALSE)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_space(data, config,param)
  out_few_cases <- cpp_ll_space(data, config, param, few_cases)
  out_rnd_cases <- cpp_ll_space(data, config, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -14.3956756)
  expect_equal(out_few_cases, -12.6518125)
})


