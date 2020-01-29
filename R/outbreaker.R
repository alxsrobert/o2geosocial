#' outbreaker: main function for reconstructing disease outbreaks
#'
#' The function \code{outbreaker} is the main function of the package. It runs
#' processes various inputs (data, configuration settings, custom priors,
#' likelihoods and movement functions) and explores the space of plausible
#' transmission trees of a densely sampled outbreaks.\cr
#'
#' @export
#'
#' @aliases outbreaker
#'
#' @rdname outbreaker
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @seealso \code{\link{outbreaker_data}} to process input data, and
#'     \code{\link{create_config}} to process/set up parameters
#'
#' @param data a list of named items containing input data as returned by
#'     \code{\link{outbreaker_data}}
#'
#' @param config a set of settings as returned by \code{\link{create_config}}
#'
#' @param likelihoods a set of log-likelihood functions as returned by
#'     \code{\link{custom_likelihoods}}
#'
#' @param priors a set of log-prior functions as returned by
#'     \code{\link{custom_priors}}
#'
#' @param moves a set of movement functions as returned by
#'     \code{\link{custom_moves}}

#' @seealso
#'
#' \itemize{
#'
#' \item \code{\link{outbreaker_data}}: function to process input data
#'
#' \item \code{\link{create_config}}: function to create default and customise
#' configuration settings
#'
#' \item \code{\link{custom_priors}}: function to specify customised prior
#' functions
#'
#' \item \code{\link{custom_likelihoods}}: function to specify customised likelihoods
#' functions
#'
#' \item \code{\link{custom_moves}}: function to create default and customise movement
#' functions
#'
#' }
#'
#' @examples
#'
#' ## get data
#' data(toy_outbreak)
#' dat <- toy_outbreak
#'
#' \dontrun{
#' ## run outbreaker
#' data("toy_outbreak")
#' dt_cases <- toy_outbreak$cases
#' dt_cases <- dt_cases[order(dt_cases$Date), ]
#' dt_regions <- toy_outbreak$dt_regions
#' all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
#'                                         rep(dt_regions$lat, nrow(dt_regions))), 
#'                                       ncol = 2), 
#'                                matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
#'                                         rep(dt_regions$lat, each = nrow(dt_regions))),
#'                                       ncol = 2))
#' 
#' dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
#' pop_vect <- dt_regions$population
#' names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region
#' 
#' data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
#'                         region = dt_cases$county, population = pop_vect, 
#'                         distance = dist_mat, a_dens = toy_outbreak$age_contact,
#'                         f_dens = dgamma(x = 1:300, scale = 0.43, shape = 27),
#'                         w_dens = dnorm(x = 1:300, mean = 11.7, sd = 2.0))
#' out <- outbreaker(data = data, config = list(n_iter = 200, sample_every = 5,
#'                                              n_iter_import = 100, sample_every_import = 5,
#'                                              gamma = 100, delta = 30, burnin = 20))
#' plot(out)
#'
#'
#' }
outbreaker <- function(data = outbreaker_data(),
                       config = create_config(),
                       priors = custom_priors(),
                       likelihoods = custom_likelihoods(),
                       moves = custom_moves()
) {
  
  ## CHECKS / PROCESS DATA ##
  data <- outbreaker_data(data = data)
  
  ## CHECK / PROCESS CONFIG ##
  config <- create_config(config, data = data)
  
  ## PRE CLUSTERING OF THE CASES
  pre_clust <- pre_clustering(data = data, config = config)
  data <- outbreaker_data(data = pre_clust[["data"]])
  config <- create_config(config = pre_clust[["config"]], 
                                    data = data)
  
  ## ADD CONVOLUTIONS TO DATA ##
  data <- add_convolutions(data = data, config = config)
  
  ## PROCESS CUSTOM FUNCTIONS FOR PRIORS AND LIKELIHOOD ##
  priors <- custom_priors(priors)
  loglike <- custom_likelihoods(likelihoods)
  
  
  ## CREATE AND INITIALIZE MCMC CHAIN ##
  temp <- create_param(data = data, config = config)
  param_store <- temp$store
  param_current <- temp$current
  param_store <- outbreaker_init_mcmc(data, param_current, param_store,
                                      loglike, priors, config)
  
  
  ## here we create a list of function for moving parameters
  moves <- bind_moves(moves = moves,
                      config = config,
                      data = data,
                      likelihoods = loglike,
                      priors = priors)
  
  
  ## IMPORTS
  
  ## preliminary run to detect imported cases this relies on a shorter run of
  ## the MCMC, then computing the average 'global influence' (-loglike) of
  ## each data point, identifying outliers (based on fixed threshold) and
  ## marking outliers down as 'imported cases'.
  if(config$find_import == TRUE){
    temp <- outbreaker_find_imports(moves = moves,
                                    data = data,
                                    param_current = param_current,
                                    param_store = param_store,
                                    config = config,
                                    likelihoods = loglike)
    param_current <- temp$param_current
    param_store <- temp$param_store
  }
  
  
  ## PERFORM MCMC
  
  ## procedure is the same as before, with some cases fixed as 'imported'
  param_store <- outbreaker_move(moves = moves,
                                 data = data,
                                 param_current = param_current,
                                 param_store = param_store,
                                 config = config,
                                 likelihoods = loglike,
                                 priors = priors)
  
  
  ## SHAPE RESULTS
  
  ## this takes the chains generated by 'outbreaker_move', stored as a list,
  ## and puts everything back together as a single data.frame; augmented data
  ## stored as vectors (e_g. 'alpha') become numbered columns of the
  ## data.frame (e_g. 'alpha_1', 'alpha_2' etc.)
  
  out <- outbreaker_mcmc_shape(param_store, data)
  
  return(out)
}

