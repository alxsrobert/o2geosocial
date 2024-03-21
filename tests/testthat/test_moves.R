context("Test function custom_moves")

## Test different moves ##
test_that("test: try different default movements", {
  ## get data
  data(toy_outbreak_short)
  x <- toy_outbreak_short
  dt_cases <- x$cases
  dt_cases <- dt_cases[order(dt_cases$Date), ]
  dt_regions <- x$dt_regions
  all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                                          rep(dt_regions$lat, nrow(dt_regions))), 
                                        ncol = 2), 
                                 matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                                          rep(dt_regions$lat, each = nrow(dt_regions))),
                                        ncol = 2))
  dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
  pop_vect <- dt_regions$population
  names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region
  
  w <- dnorm(x = 1:100, mean = 11.7, sd = 2.0)

  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect, 
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = dgamma(x = 1:100, scale = 0.43, shape = 27),
                          w_dens = w)
  config <- create_config(data = data)
  config_no_move <- create_config(data = data, move_alpha = FALSE, move_a = FALSE, 
                                  move_b = FALSE, move_t_inf = FALSE, 
                                  move_pi = FALSE, move_kappa = FALSE, 
                                  move_swap_cases = FALSE)
  pre_clust <- pre_clustering(data = data, config = config)
  data <- outbreaker_data(data = pre_clust[["data"]])
  config <- create_config(config = pre_clust[["config"]], 
                          data = data)
  data <- add_convolutions(data = data, config = config)
  if(!is.null(data$log_w_dens[config$max_kappa, ])){
    config$delta <- max(which(data$log_w_dens[config$max_kappa, ] > -20))
  }
  
  param <- create_param(data = data, config = config)$current
  likeli <- custom_likelihoods()
  priors <- custom_priors()
  
  moves <- bind_moves(config = config, data = data, likelihoods = likeli, 
                      priors = priors)
  moves_no <- bind_moves(config = config_no_move, data = data, likelihoods = likeli, 
                         priors = priors)
  expect_true(all(vapply(moves, is.function, logical(1))))
  expect_equal(length(moves), 8)
  expect_equal(length(moves_no), 0)
  
  ## Test format params after each movement
  for(i in seq_along(moves)){
    ## make moves
    set.seed(1)
    res <- moves[[i]](param = param)
    ## check that content in param after movements has identical shape
    expect_equal(length(param), length(res))
    expect_equal(length(unlist(param)), length(unlist(res)))
    expect_equal(names(param), names(res))

  }
})

## Test binding moves
test_that("test: try different default movements", {
  ## get data
  data(toy_outbreak_short)
  x <- toy_outbreak_short
  dt_cases <- x$cases
  dt_cases <- dt_cases[order(dt_cases$Date), ]
  dt_regions <- x$dt_regions
  all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                                          rep(dt_regions$lat, nrow(dt_regions))), 
                                        ncol = 2), 
                                 matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                                          rep(dt_regions$lat, each = nrow(dt_regions))),
                                        ncol = 2))
  dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
  pop_vect <- dt_regions$population
  names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region
  
  w <- dnorm(x = 1:100, mean = 11.7, sd = 2.0)
  
  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect, 
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = dgamma(x = 1:100, scale = 0.43, shape = 27),
                          w_dens = w)
  config <- create_config(data = data)
  moves <- custom_moves()
  likeli <- custom_likelihoods()
  priors <- custom_priors()
  expect_identical(moves, custom_moves(custom_moves()))
  expect_equal(length(moves), 8)
  expect_true(all(is.element(names(moves), c("a", "pi", "b", "alpha", "swap_cases", 
                                             "ancestors", "t_inf", "kappa"))))
  moves <- bind_moves(moves = moves, config = config, data = data, 
                      likelihoods = likeli, priors = priors)
  exp_names <- c("custom_prior", "custom_ll", "config", "data")
  expect_true(all(is.element(names(environment(moves$pi)), exp_names)))
  expect_true(all(is.element(names(environment(moves$a)), exp_names)))
  expect_true(all(is.element(names(environment(moves$b)), exp_names)))
  
  exp_names <- c("list_custom_ll", "config", "data")
  expect_true(all(is.element(names(environment(moves$alpha)), exp_names)))
  expect_true(all(is.element(names(environment(moves$swap_cases)), exp_names)))
  expect_true(all(is.element(names(environment(moves$ancestors)), exp_names)))
  expect_true(all(is.element(names(environment(moves$t_inf)), exp_names)))
  expect_true(all(is.element(names(environment(moves$kappa)), exp_names)))
})


## Test customization
test_that("test: try customization movements", {
  ## get data
  data(toy_outbreak_short)
  x <- toy_outbreak_short
  dt_cases <- x$cases
  dt_cases <- dt_cases[order(dt_cases$Date), ]
  dt_regions <- x$dt_regions
  all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                                          rep(dt_regions$lat, nrow(dt_regions))), 
                                        ncol = 2), 
                                 matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                                          rep(dt_regions$lat, each = nrow(dt_regions))),
                                        ncol = 2))
  dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
  pop_vect <- dt_regions$population
  names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region
  
  w <- dnorm(x = 1:100, mean = 11.7, sd = 2.0)
  
  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect, 
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = dgamma(x = 1:100, scale = 0.43, shape = 27),
                          w_dens = w, genotype = dt_cases$genotype)
  config <- create_config(data = data, n_iter = 200, burnin = 50,
                          n_iter_import = 100, gamma = 50, delta = 30)
  pre_clust <- pre_clustering(data = data, config = config)
  data <- outbreaker_data(data = pre_clust[["data"]])
  config <- create_config(config = pre_clust[["config"]])
  data <- add_convolutions(data = data, config = config)
  if(!is.null(data$log_w_dens[config$max_kappa, ])){
    config$delta <- max(which(data$log_w_dens[config$max_kappa, ] > -20))
  }
  
  param <- create_param(data = data, config = config)$current
  likeli <- custom_likelihoods()
  priors <- custom_priors()
  
  move_null <- function(data, param, config = NULL){
    return(param)
  }
  moves <- bind_moves(list(pi = move_null), data = data, config = config, 
                      likelihoods = likeli, priors = priors)
  expect_identical(moves$pi(param), param)
  expect_equal(length(moves), 8)
  
  out <- outbreaker(data, config, moves = list(pi = move_null))
  expect_true(all(out$pi == config$init_pi))
  
})

