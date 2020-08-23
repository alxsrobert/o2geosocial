context("Test function outbreaker")

## Test output format ##
test_that("Output have expected format", {
  
  ## get data
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
                          move_a = FALSE, move_b = FALSE, 
                          n_iter = 1000, n_iter_import = 500, burnin = 200)
  out <- outbreaker(data, config)
  
  print(out)
  print(out, n_row  = 2)
  print(out, n_col = 4)
  print(out, type = "cluster")
  print(out, type = "cluster", n_col = 2)
  print(out, type = "cluster", group_cluster = c(1,2,3))
  
  plot(out, type = "trace")
  plot(out, type = "trace", burnin = 200)
  plot(out, type = "hist")
  plot(out, type = "density")
  plot(out, type = "alpha")
  plot(out, type = "t_inf")
  plot(out, type = "kappa")
  plot(out, type = "network")
  plot(out, type = "cluster")
  plot(out, type = "cluster", group_cluster = c(1,2,3))
  
  
  summary(out, group_cluster = c(1,2,3))
  
  expect_error(summary(out, burnin = 2000),
               "burnin exceeds the number of steps in object")
  expect_error(plot(out, burnin = 2000),
               "burnin exceeds the number of steps in x")
  expect_error(plot(out, y = "error"), "error is not a column of x")
  
  out_df <- as.data.frame(out)
  
  ## check output
  expect_is(out, "outbreaker_chains")
  expect_is(out_df, "data.frame")
  expect_equal(nrow(out), 21)
  expect_true(!any(is.na(out_df$post)))
  expect_true(all(out_df[-1,"post"]> -1e30))
  
})

## Test convergence results ##
test_that("Results work, all component", {
  ## get data
  alpha <- c(NA,rep(1,4))
  
  times <- c(0, 3, 6, 7, 15)
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)
  
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 15, 10, 15, 0, 60, 10, 60, 0), ncol = 3)
  
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  
  genotype <- c("Not attributed", "B4", "Not attributed", "Not attributed", "B4")
  
  data <- outbreaker_data(dates = times, region = regions,
                          population = population,distance = distance,
                          age_group = age, a_dens = age_dens, 
                          w_dens = w, f_dens = f, genotype = genotype)
  config <- create_config(data = data, init_tree = alpha,
                          n_iter = 1000, n_iter_import = 500, burnin = 200)
  out <- outbreaker(data, config)
  
  out_summary <- summary(out, burnin = config$burnin)
  expect_true(all(out_summary$post > -50))
  expect_true(out_summary$tree[out_summary$tree$to == 2, "support"] > .9)
  expect_true(out_summary$tree[out_summary$tree$to == 2, "from"] == 1)
  expect_true(out_summary$tree[out_summary$tree$to == 5, "generations"] == 2)
  
})

## Test convergence results ##
test_that("Results work, 1 component at the time", {

  
  ## get data
  alpha <- c(NA,rep(1,4))
  
  times <- c(0, 3, 6, 7, 15)
  f <- c(.1, .2, .5, .2, .1)
  w <- c(.1, .2, .5, .2, .1)
  
  
  data(toy_outbreak_short)
  age_dens <- toy_outbreak_short$age_contact
  age <- c(1, 3, 3, 5, 1)
  
  regions <- c(1,1,2,2,3)
  population <- c(1e4, 5e4, 5e3)
  distance <- matrix(c(0, 15, 10, 15, 0, 60, 10, 60, 0), ncol = 3)
  names(population) <- colnames(distance) <- rownames(distance) <- 1:3
  
  genotype <- c("Not attributed", "B4", "Not attributed", "Not attributed", "B4")
  
  f_null <- function(data, param) {
    return(0.0)
  }
  
  data_time <- outbreaker_data(dates = times, w_dens = w, f_dens = f)
  config_time <- create_config(data = data_time, move_a = FALSE, move_b = FALSE, 
                               n_iter = 500, n_iter_import = 250, burnin = 100)
  like_time <- custom_likelihoods(space = f_null, age = f_null)
  out_time <- outbreaker(data = data_time, config = config_time, 
                         likelihoods = like_time)
  
  data_space <- outbreaker_data(dates = times, region = regions,
                                population = population,distance = distance)
  config_space <- create_config(data = data_space, 
                                n_iter = 500, n_iter_import = 250, burnin = 100)
  like_space <- custom_likelihoods(timing_sampling = f_null,
                                   timing_infections = f_null, age = f_null)
  out_space <- outbreaker(data = data_space, config = config_space, 
                          likelihoods = like_space)
  
  data_age <- outbreaker_data(dates = times, age_group = age, a_dens = age_dens)
  config_age <- create_config(data = data_age, move_a = FALSE, move_b = FALSE, 
                              n_iter = 500, n_iter_import = 250, burnin = 100)
  like_age <- custom_likelihoods(timing_sampling = f_null,
                                 timing_infections = f_null,
                                 space = f_null)
  out_age <- outbreaker(data = data_age, config = config_age, likelihoods = like_age)
  
  data_genotype <- outbreaker_data(dates = times, genotype = genotype,
                                   w_dens = w, f_dens = f)
  config_genotype <- create_config(data = data_genotype, move_a = FALSE, move_b = FALSE, 
                                   n_iter = 500, n_iter_import = 250, burnin = 100)
  like_genotype <- custom_likelihoods(space = f_null, age = f_null)
  out_genotype <- outbreaker(data = data_genotype, config = config_genotype, 
                             likelihoods = like_genotype)
  
  
  expect_is(out_time, "outbreaker_chains")
  expect_is(out_space, "outbreaker_chains")
  expect_is(out_age, "outbreaker_chains")
  expect_is(out_genotype, "outbreaker_chains")
  expect_equal(nrow(out_time), 11)
  expect_equal(nrow(out_space), 11)
  expect_equal(nrow(out_age), 11)
  expect_equal(nrow(out_genotype), 11)

})