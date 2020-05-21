## #' Initianilizes outputs for outbreaker
## #'
## #' This function creates initial outputs and parameter states for outbreaker.
## #'
## #' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
## #'
## #' @param param a list of data items as returned by \code{create_param}
## #'
## #' @param loglike a list of loglikelihood functions with enclosed data as returned by \code{custom_likelihood}
## #'
## #' @param priors a list of prior functions with enclosed parameters as returned by \code{custom_priors}
## #'
## #' @export
## #'
outbreaker_init_mcmc <- function(data, param_current, param_store, loglike, priors, config) {
  
  if(config$move_a == FALSE && config$move_b == FALSE && is.null(data$log_s_dens) &&
     is.null(loglike$space)){
    stop("Spatial likelihood will be -Inf, either define s_dens in outbreaker_data, or redefine the spatial component of the likelihood")
  }
  ## COMPUTE INITIAL LIKE/PRIOR/POST ##
  param_store$like[1] <- cpp_ll_all(data, config, param_current, NULL, loglike)
  param_store$prior[1] <- cpp_prior_all(param_current, config, priors)
  param_store$post[1] <- param_store$like[1] + param_store$prior[1]
  
  return(param_store)
}
