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
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)

  times <- 0:4
  alpha <- c(NA,rep(1,4))
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

## test all ##
test_that("test all: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  alpha <- c(NA,rep(1,4))
  
  times <- 0:4
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)
  
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 60, 10, 60, 0, 15, 10, 15, 0), ncol = 3)
  a <- .7
  b <- .1
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  s_dens <- population ** b * exp(-b*distance)
  
  data <- outbreaker_data(dates = times, region = regions,s_dens = s_dens,
                          population = population,distance = distance,
                          age_group = age, a_dens = age_dens, 
                          w_dens = w, f_dens = f)
  config <- create_config(data = data, init_tree = alpha, init_a = a, init_b = b,
                          move_a = FALSE, move_b = FALSE)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 3, replace = FALSE))
  
  ## tests
  out <- cpp_ll_all(data, config,param)
  out_few_cases <- cpp_ll_all(data, config, param, few_cases)
  out_rnd_cases <- cpp_ll_all(data, config, param, rnd_cases)
  
  
  expect_is(out, "numeric")
  expect_equal(out, -38.198228)
  expect_equal(out_few_cases, -25.4495926)
})

## test sum individual likelihoods ##
test_that("test indivs: ", {
  ## skip on CRAN
  skip_on_cran()
  
  ## generate data
  alpha <- c(NA,rep(1,4))
  
  times <- 0:4
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)
  
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 60, 10, 60, 0, 15, 10, 15, 0), ncol = 3)
  a <- .7
  b <- .1
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  s_dens <- population ** b * exp(-b*distance)
  
  data <- outbreaker_data(dates = times, region = regions,s_dens = s_dens,
                          population = population,distance = distance,
                          age_group = age, a_dens = age_dens, 
                          w_dens = w, f_dens = f)
  config <- create_config(data = data, init_tree = alpha, init_a = a, init_b = b,
                          move_a = FALSE, move_b = FALSE)
  param <- create_param(data = data, config = config)$current

  ## tests
  out_indiv_all <- sapply(1:data$N, function(X) cpp_ll_all(data, config, param, X))
  out_indiv_age <- sapply(1:data$N, function(X) cpp_ll_age(data, param, X))
  out_indiv_timing <- sapply(1:data$N, function(X) cpp_ll_timing(data, param, X))
  out_indiv_timing_inf <- sapply(1:data$N, function(X) 
    cpp_ll_timing_infections(data, param, X))
  out_indiv_timing_sam <- sapply(1:data$N, function(X) 
    cpp_ll_timing_sampling(data, param, X))
  out_indiv_space <- sapply(1:data$N, function(X) 
    cpp_ll_space(data, config, param, X))
  out_indiv_rep <- sapply(1:data$N, function(X) 
    cpp_ll_reporting(data, param, X))
  
  out_all <- cpp_ll_all(data, config,param)
  out_age <- cpp_ll_age(data, param)
  out_timing <- cpp_ll_timing(data, param)
  out_timing_inf <- cpp_ll_timing_infections(data, param)
  out_timing_sample <- cpp_ll_timing_sampling(data, param)
  out_space <- cpp_ll_space(data, config,param)
  out_rep <- cpp_ll_reporting(data, param)
  
  
  
  expect_is(out_all, "numeric")
  expect_is(out_age, "numeric")
  expect_is(out_timing, "numeric")
  expect_is(out_timing_inf, "numeric")
  expect_is(out_timing_sample, "numeric")
  expect_is(out_space, "numeric")
  expect_is(out_rep, "numeric")
  expect_is(out_indiv_all, "numeric")
  expect_is(out_indiv_age, "numeric")
  expect_is(out_indiv_timing, "numeric")
  expect_is(out_indiv_timing_inf, "numeric")
  expect_is(out_indiv_timing_sam, "numeric")
  expect_is(out_indiv_space, "numeric")
  expect_is(out_indiv_rep, "numeric")
  
  expect_equal(out_all, out_age + out_timing + out_space + out_rep)
  expect_equal(out_timing, out_timing_sample + out_timing_inf)
  expect_equal(out_all, sum(out_indiv_all))
  expect_equal(out_age, sum(out_indiv_age))
  expect_equal(out_timing, sum(out_indiv_timing))
  expect_equal(out_timing_inf, sum(out_indiv_timing_inf))
  expect_equal(out_timing_sample, sum(out_indiv_timing_sam))
  expect_equal(out_space, sum(out_indiv_space))

})

#Custom identical functions
test_that("Customisation with identical functions", {
  ## skip on CRAN
  skip_on_cran()
  
  ## check custom_likelihoods
  expect_identical(custom_likelihoods(),
                   custom_likelihoods(custom_likelihoods()))
  
  ## generate data
  alpha <- c(NA,rep(1,4))
  
  times <- 0:4
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)
  
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 60, 10, 60, 0, 15, 10, 15, 0), ncol = 3)
  a <- .7
  b <- .1
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  s_dens <- population ** b * exp(-b*distance)
  
  data <- outbreaker_data(dates = times, region = regions,s_dens = s_dens,
                          population = population,distance = distance,
                          age_group = age, a_dens = age_dens, 
                          w_dens = w, f_dens = f)
  config <- create_config(data = data, init_tree = alpha, init_a = a, init_b = b,
                          move_a = FALSE, move_b = FALSE)
  
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))
  
  
  ## generate custom functions with 2 arguments
  f_timing_infections  <-  function(data, param) cpp_ll_timing_infections(data, param)
  f_timing_sampling  <-  function(data, param) cpp_ll_timing_sampling(data, param)
  f_reporting  <-  function(data, param) cpp_ll_reporting(data, param)
  f_age  <-  function(data, param) cpp_ll_age(data, param)
  f_space  <-  function(data, param) cpp_ll_space(data, config, param)
  
  list_functions <- custom_likelihoods(age = f_age,
                                       space = f_space,
                                       timing_infections = f_timing_infections,
                                       timing_sampling = f_timing_sampling,
                                       reporting = f_reporting)
  
  ## tests
  expect_equal(cpp_ll_age(data, param, , list_functions$age),
               cpp_ll_age(data, param))

  expect_equal(cpp_ll_timing_infections(data, param, , list_functions$timing_infections),
               cpp_ll_timing_infections(data, param))

  expect_equal(cpp_ll_timing_sampling(data, param, , list_functions$timing_sampling),
               cpp_ll_timing_sampling(data, param))

  expect_equal(cpp_ll_space(data, config, param, , list_functions$space),
               cpp_ll_space(data, config, param))

  expect_equal(cpp_ll_reporting(data, param, , list_functions$reporting),
               cpp_ll_reporting(data, param))

  expect_equal(cpp_ll_timing(data, param, , list_functions),
               cpp_ll_timing(data, param))
  
  expect_equal(cpp_ll_all(data, config, param, , list_functions),
               cpp_ll_all(data, config, param))
  
})

#Test -inf
test_that("Function return -inf if incorrect parameters", {
  ## skip on CRAN
  skip_on_cran()
  
  ## check custom_likelihoods
  expect_identical(custom_likelihoods(),
                   custom_likelihoods(custom_likelihoods()))
  
  ## generate data
  alpha <- c(rep(5,4), NA)
  
  times <- 0:4
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  f_null <- function(data, param) return(0.0)
  
  data <- outbreaker_data(dates = times, 
                          w_dens = w, f_dens = f)
  config <- create_config(data = data, init_tree = alpha)
  likeli <- custom_likelihoods(reporting = f_null, space = f_null, age = f_null)
  likeli_all0 <- custom_likelihoods(reporting = f_null, space = f_null, age = f_null,
                                    timing_infections = f_null, 
                                    timing_sampling = f_null)
  param <- create_param(data = data, config = config)$current
  few_cases <- as.integer(c(1,3,4))
  rnd_cases <- sample(sample(seq_len(data$N), 5, replace = FALSE))
  
  out <- cpp_ll_timing(data, param)
  expect_equal(out, -Inf)
  out_all <- cpp_ll_all(data, config, param, ,likeli)
  expect_equal(out_all, -Inf)
  out_all0 <- cpp_ll_all(data, config, param, ,likeli_all0)
  expect_equal(out_all0, 0)
})

