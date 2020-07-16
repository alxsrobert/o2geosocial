#' Initializes outputs for outbreaker
#'
#' This function creates initial outputs and parameter states for outbreaker.
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @param data A list of data items as returned by \code{outbreaker_data}, or
#' arguments passed to this function.
#'
#' @param config A list of settings as returned by \code{create_config}, or
#' arguments passed to this function.
#'
#' @export
#'
#' @aliases outbreaker_store
#'
#' @aliases outbreaker_param
#'
#' @return
#'
#' A named list containing two components \code{$store} and
#' \code{$current}. \code{store} is a list with the class
#' \code{outbreaker_store}, used for storing 'saved' states of the
#' MCMC. \code{current} is a list with the class \code{outbreaker_param}, used
#' for storing 'current' states of the MCMC. \cr \cr
#'
#' \code{outbreaker_store} class content:
#' \itemize{
#'
#'  \item \code{size}{The length of the list, corresponding to the number of
#' samples saved from the MCMC.}
#'
#'  \item \code{step}{A vector of integers of length \code{size}, storing the
#' steps of the MCMC corresponding to the saved samples.}
#'
#'  \item \code{post}{A numeric vector of length \code{size}, storing
#' log-posterior values.}
#'
#'  \item \code{like}{A numeric vector of length \code{size}, storing
#' log-likelihood values.}
#'
#'  \item \code{prior}{A numeric vector of length \code{size},
#' storing log-prior values.}
#'
#'  \item \code{alpha}{A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing indices (from 1 to N) of
#' infectors for each case.}
#' 
#'  \item \code{t_inf}{A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing dates of infections for
#' each case.}
#'
#'  \item \code{kappa}{A list of length \code{size}. Each item of the list is
#' an integer vector of length \code{data$N}, storing the number of generations
#' before the last sampled ancestor for each case.}
#'
#'  \item \code{pi}{A numeric vector of length \code{size}, storing values of
#' the reporting probability.}
#' 
#'  \item \code{a}{A numeric vector of length \code{size}, storing values of
#' the first spatial parameter (population).}
#' 
#'  \item \code{b}{A numeric vector of length \code{size}, storing values of
#' the second spatial parameter (distance).}
#'
#'  \item \code{counter}{A counter used to keep track of the current iteration
#' of the MCMC (used internally).}
#'
#' }
#'
#'
#' \code{outbreaker_param} class content:
#' \itemize{
#'
#'  \item \code{alpha}{An integer vector of length \code{data$N}, storing
#' indices (from 1 to N) of infectors for each case.}
#'
#'  \item \code{t_inf}{An integer vector of length \code{data$N}, storing dates
#' of infections for each case.}
#'
#'  \item \code{kappa}{An integer vector of length \code{data$N}, storing the
#' number of generations before the last sampled ancestor for each case.}
#'
#'  \item \code{pi}{The value of the reporting probability.}
#'
#'  \item \code{a}{The value of the first spatial parameter (population).}
#'
#'  \item \code{b}{The value of the second spatial parameter (distance).}
#'  
#'  \item \code{log_s_dens}{The spatial likelihood matrix, calculated at each step
#'  from a and b if move_a == TRUE or move_b == TRUE.}
#'
#' }
#'
#' @importFrom geosphere distGeo
#' @examples
#'
#' ## load data
#' data("toy_outbreak_short")
#' dt_cases <- toy_outbreak_short$cases
#' dt_cases <- dt_cases[order(dt_cases$Date), ]
#' dt_regions <- toy_outbreak_short$dt_regions
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
#'                         region = dt_cases$Cens_tract, population = pop_vect, 
#'                         distance = dist_mat)
#'
#' ## modify config settings
#' config <- create_config(move_alpha = FALSE, n_iter = 2e5, sample_every = 1000)
#'
#' ## create param object
#' param <- create_param(data = data, config = config)
#'
create_param <- function(data = outbreaker_data(),
                         config = create_config()) {
  ## CREATE EMPTY OUTPUT VECTORS ##
  size <- round(config$n_iter/config$sample_every)
  step <- integer(size)
  post <- prior <- like <- pi <- a <- b <- double(size)
  alpha <- as.list(integer(size))
  t_inf <- as.list(integer(size))
  kappa <- as.list(integer(size))

  ## SET CURRENT VALUES AND COUNTER ##
  step[1] <- 1L
  current_alpha <- alpha[[1]] <- config$init_alpha
  current_kappa <- kappa[[1]] <- config$init_kappa
  current_pi <- pi[1] <- config$init_pi
  current_a <- a[1] <- config$init_a
  current_b <- b[1] <- config$init_b
  if (is.null(config$init_t_inf)) {
    current_t_inf <- t_inf[[1]] <- data$dates - which.max(data$f_dens) + 1L
  } else {
    current_t_inf <- t_inf[[1]] <- config$init_t_inf
  }
  counter <- 1L
  
  
  store <- list(
    size = size, step = step,
    post = post, like = like, prior = prior,
    alpha = alpha, t_inf = t_inf,
    kappa = kappa, pi = pi, a = a, b = b,
    counter = counter
  )
  class(store) <- c("outbreaker_store", "list")
  
  ## ADD SPATIAL LIKELIHOOD TO PARAM ##
  current_log_s_dens <- rep(list(matrix(-Inf, nrow = length(unique(data$region)), 
                                        ncol = length(unique(data$region)))),
                            config$max_kappa)
  if(config$move_a == FALSE && config$move_b == FALSE && !is.null(data$log_s_dens)){
    current_log_s_dens[[1]] <- data$log_s_dens
    if(config$max_kappa > 1){
      for(i in 2:config$max_kappa){
        s_dens <- data$s_dens %*% exp(current_log_s_dens[[i-1]])
        current_log_s_dens[[i]]  <- log(t(t(s_dens)/colSums(s_dens)))
      }
    } 
  }    
  current  <- list(
    alpha = current_alpha, 
    t_inf = current_t_inf, 
    kappa = current_kappa, pi = current_pi, 
    a = current_a, b = current_b,
    log_s_dens = current_log_s_dens
  )
  
  class(current) <- c("outbreaker_param", "list")

  ## SHAPE CHAIN ##
  out <- list(store = store,
              current = current)
  return(out)
}
