context("Test function bind_to_function")

## test data ##
test_that("test: ", {
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
  
  expect_error(bind_to_function(o2geosocial:::cpp_move_alpha),
               "'...' is empty but 'f' has more than one argument.")
  
  expect_error(bind_to_function(cpp_move_alpha, data), 
               "All arguments provided through '...' need to be named.")
  
  expect_error(bind_to_function(cpp_move_pi, data = data),
               paste("Arguments of cpp_move_pi missing from '...'",
                     "with no default: config"))
  
})