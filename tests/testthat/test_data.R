context("Test function outbreaker_data")

## test data ##
test_that("test: data are processed fine", {
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
  
  out <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                         region = dt_cases$Cens_tract, population = pop_vect, 
                         distance = dist_mat, a_dens = x$age_contact,
                         f_dens = dgamma(x = 1:100, scale = 0.43, shape = 27),
                         w_dens = w)
  
  max_range <- as.numeric(max(dt_cases$Date) - min(dt_cases$Date))
  ## check output
  expect_is(out, "list")
  expect_equal(out$max_range, max_range)
  expect_equal(out$N, nrow(dt_cases))
  expect_equal(out$import, rep(FALSE, nrow(dt_cases)))
  expect_equal(out$is_cluster, rep(1, nrow(dt_cases)))
  
  expect_is(out$log_w_dens, "matrix")
  expect_is(out$log_a_dens, "list")
  expect_is(out$log_f_dens, "numeric")
  
  expect_error(outbreaker_data(dates = 1, w_dens = c(0,-1)),
               "w_dens has negative entries")
  
  expect_error(outbreaker_data(dates = 1, w_dens = c(0,1), f_dens = c(0,-1)),
               "f_dens has negative entries")
  expect_error(outbreaker_data(dates = 1, region = 1, population = pop_vect[1:5], 
                               distance = dist_mat[11:15, 11:15]),
               "The vector population should have the same names as the distance matrix")
  expect_error(outbreaker_data(dates = 1, region = names(pop_vect[6]), 
                               population = pop_vect[1:5], 
                               distance = dist_mat[1:5, 1:5]),
               "Some regions are not in the population vector")
  expect_error(outbreaker_data(dates = 1, region = names(pop_vect[6]), 
                               population = pop_vect[1:5], 
                               distance = dist_mat[2:5, 2:5]),
               "The vector population should have the same names as the distance matrix")
  expect_error(outbreaker_data(dates = 1, region = 6, 
                               population = pop_vect[1:5], 
                               distance = dist_mat[1:5, 1:5]),
               "The length of the population vector is lower than the maximum value of the region vector")
  
})




