
#' Customise samplers for outbreaker
#'
#' This function is used to specify customised movement functions
#' (a.k.a. samplers) for outbreaker. Custom functions are specified as a named
#' list or series of comma-separated, named arguments, indicating which type of
#' movement they implement. Values currently available are:
#'
#' \itemize{
#'
#' \item \code{pi}: movement of the reporting probability; by default, the function
#' \code{cpp_move_pi} is used.
#'
#' \item \code{a}: movement of the first spatial parameter; by default, the function
#' \code{cpp_move_a} is used.
#'
#' \item \code{b}: movement of the second  spatial parameter; by default, the function
#' \code{cpp_move_b} is used.
#'
#' \item \code{alpha}: movement of the transmission tree, by randomly proposing
#' infectors in the pool of cases infected before; by default, the function
#' \code{cpp_move_alpha} is used.
#'
#' \item \code{swap_cases}: movement of the transmission tree, by swapping
#' infectors and infected cases; by default, the function
#' \code{cpp_move_swap_cases} is used.
#'
#' \item \code{ancestors}: movement of the transmission tree, by changing
#' the ancestors of the different trees in a cluster; by default, the function
#' \code{cpp_move_ancestors} is used.
#'
#' \item \code{t_inf}: movement of the date of infection; by default, the
#' function \code{cpp_move_t_inf} is used.
#'
#' \item \code{kappa}: movement of the number generations between cases; by
#' default, the function \code{cpp_move_kappa} is used.
#'
#' }
#'
#'
#' Movement functions must have an argument \code{param}, which is a list of
#' parameters and augmented data of the class \code{\link{create_param}}.
#' Each movement function will be enclosed with its other arguments, so that the
#' resulting function will have a single argument 'param'. For non-standard
#' movements (i.e. none of the names specified above), the closure will contain:
#'
#' \itemize{
#'
#' \item \code{data}: a list of named items containing input data as returned by
#' \code{\link{outbreaker_data}}
#'
#' \item \code{config}:  a list of named items containing input data as returned by
#' \code{\link{create_config}}
#'
#' \item \code{likelihoods}: a list of named custom likelihood functions as returned by
#' \code{\link{custom_likelihoods}}
#'
#' \item \code{priors}: a list of named custom prior functions as returned by
#' \code{\link{custom_priors}}
#'
#' }
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @seealso See \href{http://www.repidemicsconsortium.org/outbreaker2/articles/customisation.html#customising-movements}{customization vignette} for detailed examples on how to customise movement functions.
#' 
#' @export
#'
#' @param ... A list or a series of named, comma-separated functions
#'     implementing movements of parameters or augmented data.
#'
#' @return A list of movement functions with a single argument 'param', with
#'     class \code{outbreaker_moves}.

custom_moves <- function(...) {
  
  move_functions <- list(...)
  
  if (length(move_functions) == 1L && is.list(move_functions[[1]])) {
    move_functions <- move_functions[[1]]
  }
  
  
  defaults <- list(pi = cpp_move_pi,
                   a = cpp_move_a,
                   b = cpp_move_b,
                   alpha = cpp_move_alpha,
                   swap_cases = cpp_move_swap_cases,
                   ancestors = cpp_move_ancestors,
                   t_inf = cpp_move_t_inf,
                   kappa = cpp_move_kappa
  )
  
  
  moves <-  modify_defaults(defaults, move_functions, FALSE)
  moves_names <- names(moves)
  
  
  
  ## check all moves are functions
  
  function_or_null <- function(x) {
    is.null(x) || is.function(x)
  }
  
  is_ok <- vapply(moves, function_or_null, logical(1))
  
  if (!all(is_ok)) {
    culprits <- moves_names[!is_ok]
    msg <- paste0("The following moves are not functions: ",
                  paste(culprits, collapse = ", "))
    stop(msg)
  }
  
  
  ## check they all have a 'param' argument
  
  param_is_arg <- function(x) {
    if(is.function(x)) {
      return ("param" %in% methods::formalArgs(x))
    }
    
    return(TRUE)
  }
  
  param_ok <- vapply(moves, param_is_arg, logical(1))
  
  if (!all(param_ok)) {
    culprits <- moves_names[!param_ok]
    msg <- paste0("The following moves dont' have a 'param' argument: ",
                  paste(culprits, collapse = ", "))
    stop(msg)
  }
  
  
  class(moves) <- c("outbreaker_moves", "list")
  return(moves)
}
