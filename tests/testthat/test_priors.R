context("Test prior functions")

test_that("Priors in cpp and R give identical results.", {
  ## set up param with parameters and basic config
  config <- create_config()
  param <- list(pi = .85, a = .8, b = .1)
  
  expect_equal(dbeta(x = param$pi, shape1 = config$prior_pi[1], 
                     shape2 = config$prior_pi[2], log = TRUE),
               cpp_prior_pi(param = param, config = config))
  expect_equal(dunif(x = param$a, min = config$prior_a[1], max = config$prior_a[2],
                     log = TRUE),
               cpp_prior_a(param = param, config = config))
  expect_equal(dunif(x = param$b, min = config$prior_b[1], max = config$prior_b[2],
                     log = TRUE),
               cpp_prior_b(param = param, config = config))
  
  expect_equal(cpp_prior_pi(param = param, config = config) + 
                 cpp_prior_a(param = param, config = config) + 
                 cpp_prior_b(param = param, config = config),
               cpp_prior_all(param = param, config = config))
  
})

test_that("Customization of prior distributions works well", {
  ## set up param with parameters and basic config
  config <- create_config(prior_pi = c(.8, .02),
                          prior_a = c(.9, .1), 
                          prior_b = c(.2, .02))
  param <- list(pi = .85, a = .7, b = .1)
  
  prior <- custom_priors()
  print(prior)
  expect_identical(prior, custom_priors(prior))
  
  normal_prior_pi <- function(param){
    return(dnorm(x = param$pi, mean = config$prior_pi[1], sd = config$prior_pi[2], log = T))
  }
  normal_prior_a <- function(param){
    return(dnorm(x = param$a, mean = config$prior_a[1], sd = config$prior_a[2], log = T))
  }
  normal_prior_b <- function(param){
    return(dnorm(x = param$b, mean = config$prior_b[1], sd = config$prior_b[2], log = T))
  }
  expect_error(custom_priors(pi = "error_pi"), "The following priors are not functions: pi")
  expect_error(custom_priors(pi = function(param1, param2) return(c(2))), 
               "The following priors don't have a single argument: pi")
  new_prior <- custom_priors(pi = normal_prior_pi, a = normal_prior_a, b = normal_prior_b)
  expect_identical(new_prior$pi, normal_prior_pi)
  expect_identical(new_prior$a, normal_prior_a)
  expect_identical(new_prior$b, normal_prior_b)
  expect_equal(cpp_prior_pi(param, config, new_prior$pi), -0.13191552)
  expect_equal(cpp_prior_a(param, config, new_prior$a), -0.61635344)
  expect_equal(cpp_prior_b(param, config, new_prior$b), -9.5069156)
  print(new_prior)
  
  expect_equal(cpp_prior_all(param, config, new_prior), 
               -0.13191552 + -0.61635344 + -9.5069156)
  
})