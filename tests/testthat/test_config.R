context("Test outbreaker config")


## test settings ##
test_that("test: settings are processed fine", {
  ## skip on CRAN
  skip_on_cran()
  
  ## get data
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
  f <- dgamma(x = 1:100, scale = 0.43, shape = 27)
  
  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect, 
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = f, w_dens = w)
  expect_is(create_config(), "list")
  expect_is(create_config(data = data), "list")
  expect_equal(create_config(data = data, init_tree = c(NA, rep(1, data$N - 1)))$init_alpha,
               create_config(data = data, init_tree = c(NA, rep(1, data$N - 1)))$init_tree)
  expect_error(create_config(data = data, init_tree = rep(1, data$N)),
               "There should be an ancestor in the initial tree")
  expect_equal(create_config(data = data, init_kappa = c(NA, rep(1, data$N - 1)))$init_kappa,
               c(NA, rep(1, data$N - 1)))
  expect_error(create_config(fakearg = 2), "Additional invalid options: fakearg")
  expect_error(create_config(spatial_method = "invalid"), 
               "invalid value for spatial_method, spatial_method is either exponential, or power-law.")
  expect_error(create_config(gamma = -1), "gamma is below 0")
  expect_error(create_config(delta = -1), "delta is below 0")
  expect_error(create_config(init_kappa = 0), "init_kappa has values smaller than 1")
  expect_error(create_config(init_pi = -1), "init_pi is negative")
  expect_error(create_config(data = data, init_tree = c(-NA, 1)), 
               "length of init_alpha or init_tree incorrect")
  expect_warning(create_config(move_a = TRUE, max_kappa = 10), 
                 "If spatial kernel parameters are estimated, max_kappa is set to 2")
  expect_error(create_config(n_iter = 0),
               "n_iter is smaller than 2")
  expect_error(create_config(sample_every = 0),
               "sample_every is smaller than 1")
  
})


## Test init tree ##
test_that("test: initial tree does not mix genotypes", {
  ## skip on CRAN
  skip_on_cran()
  
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
  f <- dgamma(x = 1:100, scale = 0.43, shape = 27)

  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect,
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = f, w_dens = w, genotype = dt_cases$Genotype)
  config <- create_config(data = data)
  tree_ances <- config$init_alpha

  while(any(!is.na(tree_ances[tree_ances]))) 
    tree_ances[!is.na(tree_ances[tree_ances])] <- 
      tree_ances[tree_ances][!is.na(tree_ances[tree_ances])]
  
  tree_ances[is.na(tree_ances)] <- which(is.na(tree_ances))
  genotype_tree <- numeric(length(unique(tree_ances)))
  nb_gen_rep_per_tree <- sapply(unique(tree_ances), function(X) {
    gens <- unique(data$genotype[which(tree_ances == X)])
    return(length(gens[gens != "Not attributed"]))
  })

  expect_true(all(nb_gen_rep_per_tree < 2))

  expect_error(create_config(data = data, init_tree = c(NA, rep(1, data$N - 1))),
               "There should be one reported genotype per tree at most.")
})

