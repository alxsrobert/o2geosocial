#'
#' Customise likelihood functions for o2geosocial
#'
#' This function is used to specify customised likelihood functions for
#' o2geosocial Custom functions are specified as a named list or series of
#' comma-separated, named arguments, indicating which log-likelihood component
#' they compute. Values currently available are:
#'
#' \itemize{
#'
#' \item \code{timing_sampling}: the likelihood of sampling times; by default, the function
#' \code{cpp_ll_timing_sampling} is used.
#'
#' \item \code{timing_infections}: the likelihood of infection times; by default, the function
#' \code{cpp_ll_timing_infections} is used.
#'
#' \item \code{reporting}: the likelihood of the reporting process; by default,
#' the function \code{cpp_ll_reporting} is used.
#' 
#' \item \code{space}: the likelihood of spatial distances; by default,
#' the function \code{cpp_ll_space} is used.
#' 
#' \item \code{age}: the likelihood of the age contacts; by default,
#' the function \code{cpp_ll_age} is used.
#' 
#' }
#'
#' All log-likelihood functions should have the following arguments, in this
#' order:
#'
#' \itemize{
#'
#' \item \code{data}: a list of named items containing input data as returned by
#' \code{\link{outbreaker_data}}
#'
#' \item \code{param}: a list of parameters with the class
#' \code{\link{create_param}}
#'
#' }
#'
#'
#' @return
#' A named list of functions with the class \code{custom_likelihood}, each
#'     implementing a customised log-likelihood components of
#'     outbreaker. Functions which are not customised will result in a NULL
#'     component.
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @param ... a named list of functions, each computing a log-likelihood component.
#'
#' @export
#' 
#' @examples
#' 
#' ## specify a null model by disabling all likelihood components
#' f_null <- function(data, config = NULL, param, i) {
#'   return(0.0)
#' }
#' 
#' 
#' null_model <- custom_likelihoods(timing_sampling = f_null,
#'                                  timing_infections = f_null,
#'                                  reporting = f_null,
#'                                  space = f_null,
#'                                  age = f_null)
#'
#' null_config <- list(find_import = FALSE,
#'                     n_iter = 200, gamma = 100, delta = 30,
#'                     sample_every = 1)
#'
#' ## load data
#' data("toy_outbreak_short")
#' dt_cases <- toy_outbreak_short$cases
#' dt_cases <- dt_cases[order(dt_cases$Date), ][1:15,]
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
#' res_null <- outbreaker(data = data,
#'                        config = null_config,
#'                        likelihoods = null_model)
#'


## USING CUSTOM LIKELIHOOD FUNCTIONS

## Likelihood functions in o2geosocial are implemented using Rcpp. However,
## these functions can also be replaced by customized functions. These can be
## specified by the user, through the '...' argument of
## 'custom_likelihoods'. These functions must have 2 arguments:

## - data: a valid 'outbreaker_data' list

## - param: a list containing current parameter states, as returned by
## - create_param

custom_likelihoods <- function(...) {
  
  ll_functions <- list(...)
  
  if (length(ll_functions) == 1L && is.list(ll_functions[[1]])) {
    ll_functions <- ll_functions[[1]]
  }
  
  
  defaults <- list(reporting = NULL,
                   timing_infections = NULL,
                   timing_sampling = NULL,
                   space = NULL,
                   age = NULL
  )
  
  likelihoods <-  modify_defaults(defaults, ll_functions, FALSE)
  likelihoods_names <- names(likelihoods)
  
  
  
  ## check all likelihoods are functions
  
  function_or_null <- function(x) {
    is.null(x) || is.function(x)
  }
  
  is_ok <- vapply(likelihoods, function_or_null, logical(1))
  
  if (!all(is_ok)) {
    culprits <- likelihoods_names[!is_ok]
    msg <- paste0("The following likelihoods are not functions: ",
                  paste(culprits, collapse = ", "))
    stop(msg)
  }
  
  
  ## check they all have a single argument
  
  with_three_four_args <- function(x) {
    if(is.function(x)) {
      return (length(methods::formalArgs(x)) == 3L || 
                length(methods::formalArgs(x)) == 4L)
    }
    
    return(TRUE)
  }
  
  three_four_args <- vapply(likelihoods, with_three_four_args, logical(1))
  
  if (!all(three_four_args)) {
    culprits <- likelihoods_names[!three_four_args]
    msg <- paste0("The following likelihoods don't have three or four arguments: ",
                  paste(culprits, collapse = ", "))
    stop(msg)
  }
  
  
  class(likelihoods) <- c("custom_likelihoods", "list")
  return(likelihoods)
  
}







#' @rdname custom_likelihoods
#'
#' @export
#'
#' @aliases print.custom_likelihoods
#'
#' @param x an \code{outbreaker_config} object as returned by \code{create_config}.
#'

print.custom_likelihoods <- function(x, ...) {
  cat("\n\n ///// outbreaker custom likelihoods ///\n")
  cat("\nclass:", class(x))
  cat("\nnumber of items:", length(x), "\n\n")
  
  is_custom <- !vapply(x, is.null, FALSE)
  
  
  names_default <- names(x)[!is_custom]
  if (length(names_default) > 0) {
    cat("/// custom likelihoods set to NULL (default used) //\n")
    print(x[!is_custom])
  }
  
  
  names_custom <- names(x)[is_custom]
  if (length(names_custom) > 0) {
    cat("/// custom likelihoods //\n")
    print(x[is_custom])
  }
  
  return(invisible(NULL))
  
}

