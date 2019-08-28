## #' Stores MCMC samples for outbreaker
## #'
## #' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
## #'
## #' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
## #'
## #' @inheritParams outbreaker
## #'
## #' @param param a list of output items as returned by \code{create_param}
## #'
## #' @param step an integer indicating the MCMC iteration being stored
## #'
## #' @export
## #'
outbreaker_mcmc_store <- function(param_current, param_store, data, config,
                                  likelihoods, priors, step) {
  ## UPDATE COUNTER
  counter <- param_store$counter <- param_store$counter + 1
  
  ## STORE STEP
  param_store$step[counter] <- step
  
  ## STORE LIKELIHOOD, PRIOR, POSTERIOR
  param_store$like[counter] <- cpp_ll_all(data = data, 
                                          config = config, 
                                          param = param_current, 
                                          i = NULL, 
                                          custom_functions = likelihoods)
  param_store$prior[counter] <- cpp_prior_all(param_current, config, priors)
  param_store$post[counter] <- param_store$like[counter] + param_store$prior[counter]
  
  ## PARAMETERS 
  param_store$pi[counter] <- param_current$pi
  param_store$alpha[[counter]] <- param_current$alpha
  param_store$t_inf[[counter]] <- param_current$t_inf
  param_store$kappa[[counter]] <- param_current$kappa
  param_store$a[[counter]] <- param_current$a
  param_store$b[[counter]] <- param_current$b
  
  return(param_store)
}

