## #' Movements of augmented data and parameters for o2geosocial
## #'
## #' This function moves all parameters for o2geosocial
## #'
## #' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
## #'
## #'
## #' @param moves a list of movement functions as returned by \code{bind_moves}
## #'     (internal function)
## #'
## #' @param param a list of parameters as returned by
## #'     \code{create_param}
## #'
## #' @param densities a list containing lists of functions computing densities,
## #'     named: 'loglike' (log-likelihoods), 'priors' and 'posteriors'
## #'
## #' @return a potentially modified list of parameters as returned by
## #'     \code{create_param}
## #'
outbreaker_move <- function(moves, data, param_current,
                            param_store, config,
                            likelihoods, priors) {
  ## get number of moves ##
  J <- length(moves)
  if(config$find_import == TRUE){
    ## Define individual likelihood matrix
    n_measures <- sum(seq(0, config$n_iter, config$sample_every) > 0)
    if(config$sample_every != 1) n_measures = n_measures + 1
    influences <- matrix(0, ncol = data$N, nrow = n_measures)
    colnames(influences) <- 1:data$N
    influences[1,] <- - vapply(seq_len(data$N), function(case) 
      (cpp_ll_all(data = data,config = config,
                  param = param_current,
                  i = case, custom_functions = likelihoods)
      ),
      numeric(1))
  }
  
  ## RUN MCMC ##
  for (i in seq.int(2, config$n_iter, 1)) {
    ## move parameters / augmented data
    for (j in seq_len(J)) {
      ## move parameters
      param_current <- moves[[j]](param_current)
    }
    ## store outputs and influence if needed
    if ((i %% config$sample_every) == 0) {
      param_store <- outbreaker_mcmc_store(param_current, param_store, data,
                                           config, likelihoods, priors, i)
      influence_i <- - vapply(seq_len(data$N), function(case) 
        (cpp_ll_all(data = data,config = config, param = param_current, i = case,
                    custom_functions = likelihoods)
        ), numeric(1))
      if(config$verbatim == TRUE){
        message(paste0("Iteration number: ", i, "/", config$n_iter,
                       "|| likelihood = ", round(sum(influence_i), 2)))
      }
      if(config$find_import == TRUE){
        counter <- i / config$sample_every + 1
        if(config$sample_every == 1) counter <- i
        # Correct the influence for cases linked to an importation:
        ## Likelihoods of cases already classified as importations are not taken
        ## into account when adding new imports. Therefore, implausible latent
        ## periods may not result in the addition of an importation. To correct
        ## for this, we add the f_dens of importations to the influence score of 
        ## their infectee: 
        alpha <- param_current$alpha
        linked_to_imports <- which(is.na(alpha[alpha]) & !is.na(alpha))
        f_like_linked <- - vapply(
          linked_to_imports, function(i) (
            cpp_ll_timing_sampling(data = data, param = param_current, i = i,
                                   custom_function = likelihoods$timing_sampling)
          ), numeric(1))
        f_like_ances <- - vapply(
          alpha[linked_to_imports], function(i) (
            cpp_ll_timing_sampling(data = data, param = param_current, i = i,
                                   custom_function = likelihoods$timing_sampling)
          ), numeric(1))
        # Compute average f_dens for "linked_to_imports"
        average_f_link <- pmax(f_like_ances, f_like_linked)
        
        influence_i[linked_to_imports] <- influence_i[linked_to_imports] - 
          f_like_linked + average_f_link
        influences[counter,] <- influence_i
        
      }
    }
  } # end of the chain

  if(config$find_import == TRUE){
    ## Remove unlikely transmission links
    # Define threshold
    threshold <- -log(config$outlier_threshold)*5
    bad_ancestor <- rep(threshold, data$N)
    #Compare threshold to influence
    bad_ancestor_matrix <- (influences)<(bad_ancestor) | 
      matrix(!config$move_alpha, ncol = ncol(influences), nrow = nrow(influences), 
             byrow = T)
    # In cases where the likelihood in worst than the threshold, the transmission link is removed
    bad_ancestor_matrix[bad_ancestor_matrix == FALSE] <- NA
    bad_ancestor_list <- split(x = t(bad_ancestor_matrix), 
                               rep(1:(n_measures), each = data$N))
    # Update alpha and kappa
    param_store$alpha <- Map("*", param_store$alpha, bad_ancestor_list)
    param_store$kappa <- Map("*", param_store$kappa, bad_ancestor_list)
  }
  ## output is a list of saved chain states
  return(param_store)
}
