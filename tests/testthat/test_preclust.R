context("Test function pre_clustering")

## test pre_clustering ##
test_that("test: pre clustering gives expected output", {
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
  w <- w[which(w > 1e-15)]
  
  data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                          region = dt_cases$Cens_tract, population = pop_vect, 
                          distance = dist_mat, a_dens = x$age_contact,
                          f_dens = dgamma(x = 1:100, scale = 0.43, shape = 27),
                          w_dens = w)
  config <- create_config(data = data)
  
  test_clust <- pre_clustering(data = data, config = config)
  expect_equal(test_clust$data$is_cluster, data$is_cluster)
  expect_equal(test_clust$config$init_alpha, config$init_alpha)
  expect_equal(test_clust$config$init_kappa, config$init_kappa)
  expect_equal(test_clust$config$move_alpha, config$move_alpha)
  expect_equal(test_clust$config$move_kappa, config$move_kappa)
  
  config <- create_config(data = data, gamma = 1, delta = 1)
  test_clust <- pre_clustering(data = data, config = config)
  expect_true(all(is.na(test_clust$config$init_alpha)))
  
  })
