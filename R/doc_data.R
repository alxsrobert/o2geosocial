#' Simulated outbreak
#'
#' We generated a data containing 500 cases, from simulated outbreaks across the US,
#' used to illustrate \code{measlesoutbreaker}. This list contains the following:
#'
#' \itemize{
#'
#' \item \code{$cases}: A data table summarising the epidemiological features of the
#' 500 cases. It contains the State, onset date,  genotype,  county, age group,
#' import status, cluster and generation of the cases.
#'
#' \item \code{$distance}: A matrix indicating the distance between each county
#'
#' \item \code{$population}: A vector indicating the population per county
#'
#' \item \code{$age_contact}: A matrix indicating the number of contacts between 
#' age groups
#'
#' }
#'
#' @aliases toy_outbreak
#' @docType data
#' @author Alexis Robert \email{alexis.robert@lshtm.ac.uk}
#' @keywords datasets
#'
#' @examples
#' data("toy_outbreak")
#' names(toy_outbreak)
#' toy_outbreak
#' 
"toy_outbreak"
