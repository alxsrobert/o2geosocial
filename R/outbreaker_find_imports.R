## This function is not exported and is meant for internal use.
## It is similar to the MCMC implemented by 'outbreaker.move'. 
## It has its own number of iterations ('config$n.iter.import') and sampling 
## frequency ('config$sample.every.import'), and stores individual likelihoods 
## for each case and saved iteration. These chains are then used to compute 
## the 'global influence' of each case, flag the number of outliers per 
## iteration based on these values and the threshold 
## ('config$outlier.threshold' and 'config$outlier.relative'). 
## We then compute the minimum number of outliers per iteration in the MCMC run
## and mark these cases as imported, i.e. for which the ancestor will be 'NA'.
#' @importFrom stats quantile
outbreaker_find_imports <- function(moves, data, param_current,
                                    param_store, config, likelihoods) {
  if (!config$find_import) {
    return(list(param_current = param_current,
                param_store = param_store))
  }
  ## get number of moves ##
  J <- length(moves)
  
  ## create matrix of individual influences ##
  n_measures <- sum(seq(0, config$n_iter_import, 
                        config$sample_every_import) > config$burnin)
  influences <- matrix(0, ncol = data$N, nrow = n_measures)
  colnames(influences) <- 1:data$N

  ## Function to run the short MCMC ##
  MCMC_imports <- function(config_imports, data_imports, moves_imports,
                           J_imports, likelihoods_imports,
                           influences_imports, 
                           param_current_imports){
    counter <- 1L
    for (i in seq.int(2, config_imports$n_iter_import, 1)) {
      ## move parameters / augmented data_imports
      for (j in seq_len(J_imports)) 
        param_current_imports <- moves_imports[[j]](param_current_imports)
      ## store outputs if needed
      if ((i %% config_imports$sample_every_import) == 0 && i>config$burnin) {
        influences_imports[counter,] <- 
          - vapply(seq_len(data_imports$N), 
                   function(i) 
                     (cpp_ll_all(data = data_imports, config = config_imports,
                                 param = param_current_imports, i = i,
                                 custom_functions = likelihoods_imports)
                      ), numeric(1))
        # Do not take the influence of importation into account
        # (we are only interested in connected cases)
        influences_imports[counter, is.na(param_current_imports$alpha)] <- -1
        if(config_imports$verbatim == TRUE) 
          message(paste0("Finding import, Iteration number: ", i, "/",
                         config_imports$n_iter_import, "|| likelihood = ", 
                         round(sum(influences_imports[counter,]), 2)))
        counter <- counter + 1L
      }
    } # end of the chain
    return(list(param_current = param_current_imports,
                influences = influences_imports))
  }
  outcome_MCMC <- MCMC_imports(config_imports = config,
                               data_imports = data,
                               moves_imports = moves,
                               likelihoods_imports = likelihoods,
                               J_imports = J,
                               influences_imports = influences,
                               param_current_imports = param_current)
  influences <- outcome_MCMC[["influences"]]
  param_current <- outcome_MCMC[["param_current"]]
  
  ## Influence matrix per cluster
  if(max(data$is_cluster) == 1) list_influences <- list(influences) else
    list_influences <- sapply(1:max(data$is_cluster), function(X)
      return(influences[, which(data$is_cluster == X)]))
  
  ## Compute the likelihood threshold. The likelihoods of connection lower than
  ## the threshold will be considered implausible
  threshold <- -log(config$outlier_threshold)*5
  
  # Do not take the influence of importation into account
  # (we are only interested in connected cases)
  influences_vect <- c(influences[influences >= 0])
  
  # Compute the relative threshold 
  if(config$outlier_relative == TRUE){
    threshold <- quantile(influences_vect, probs = config$outlier_threshold)
    message(
      "The import threshold computed from the likelihoods of connection is ",
      round(threshold, 1), ". It corresponds to an absolute threshold of ", 
      round(exp(-threshold/5), 3), 
      ". Therefore, a case with a probability of connection per component above ",
      round(exp(-threshold/5), 3), " will not be classified as an importation.") 
  }
  if(config$outlier_plot == TRUE){
    plot_importations(influences_vect = influences_vect, threshold = threshold, 
                      config = config)
  } else 
    message("Set config$outlier_plot = TRUE to plot the likelihoods of connection")
  
  # Set cases with no likely infector as import
  new_imports <- unlist(lapply(list_influences, function(X){
    if(!is.null(dim(X))){
      # Compute n_imports_iteration, the  number of import per iteration in each cluster
      bad_ancestor <- rep(threshold, dim(X)[2])
      bad_ancestor_matrix <- X>bad_ancestor
      n_imports_iteration <- apply(bad_ancestor_matrix, 1, sum)
      # Add the min(n_imports_iteration) new import in the cluster
      imports <- names(which(
        bad_ancestor_matrix[which.min(n_imports_iteration),] == TRUE))
      imports <- as.numeric(imports)
      return(imports)
    }
  }))
  
  message(paste0("Adding ", length(new_imports), " new importations"))
  
  ## store initial param values ##
  ini_param <- list(current = param_current, store = param_store)
  
  ## All outliers are considered as introductions, so that ancestries (alpha) are set to 'NA' and
  ## the number of generations between cases and their ancestor (kappa) is set to NA; the
  ## movements of alpha and kappa for these cases is also disabled; because the config has been
  ## altered in these cases, we systematically return the config as well as the initial
  ## parameters.
  ini_param$store$alpha[[1]][new_imports] <- ini_param$current$alpha[new_imports] <- NA
  ini_param$store$kappa[[1]][new_imports] <- ini_param$current$kappa[new_imports] <- NA
  ini_param$store$a[[1]] <- ini_param$current$a <- param_current$a
  ini_param$store$b[[1]] <- ini_param$current$b <- param_current$b
  ini_param$store$pi[[1]] <- ini_param$current$pi <- param_current$pi
  
  ini_param$store$like[1] <- cpp_ll_all(data, config, param_current, NULL, likelihoods)
  ini_param$store$post[1] <- ini_param$store$like[1] + ini_param$store$prior[1]
  
  return(list(param_current = ini_param$current,
              param_store = ini_param$store, threshold = as.numeric(threshold)))
}
