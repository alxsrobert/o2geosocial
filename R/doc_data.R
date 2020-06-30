#' Simulated outbreaks
#'
#' We generated two datasets used to illustrate \code{o2geosocial}. 
#' The first one (toy_outbreak_long) contains 1,940 cases, from simulated outbreaks 
#' nationwide between 2010 and 2017. The list contains the following:
#'
#' \itemize{
#'
#' \item \code{$cases}: A data table summarising the epidemiological features of the
#' 1,940 cases. It contains the ID, State, onset date,  genotype,  county, age group,
#' import status, cluster, generation and infector of the cases.
#'
#' \item \code{$dt_regions}: A data table containing the ID, population, longitude 
#' and latitude of each region. Should be used to compute the distance matrix, using 
#' the package geosphere.
#' 
#' \item \code{$age_contact}: A matrix indicating the number of contacts between 
#' age groups
#'
#' }
#'
#' @aliases toy_outbreak_long
#' @docType data
#' @author Alexis Robert \email{alexis.robert@lshtm.ac.uk}
#' @keywords datasets
#'
#' @examples
#' data("toy_outbreak_long")
#' names(toy_outbreak_long)
#' toy_outbreak_long
#' 
"toy_outbreak_long"

#' Simulated outbreaks
#'
#' Second dataset used to illustrate \code{o2geosocial}. (toy_outbreak_short) is a 
#' smaller data set (75 cases), spread across different Census tracks in Ohio 
#' (population and location of each region taken from 
#' https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html).
#' The list contains the following:
#'
#' \itemize{
#'
#' \item \code{$cases}: A data table summarising the epidemiological features of the
#' 75 cases. It contains the ID, state, onset date,  genotype, Census tract, age group,
#' import status, cluster, generation and infector of the cases.
#'
#' \item \code{$dt_regions}: A data table containing the ID, population, longitude 
#' and latitude of each region. Should be used to compute the distance matrix, using 
#' the package geosphere.
#' 
#' \item \code{$age_contact}: A matrix indicating the number of contacts between 
#' age groups
#'
#' }
#'
#' @aliases toy_outbreak_short
#' @docType data
#' @author Alexis Robert \email{alexis.robert@lshtm.ac.uk}
#' @keywords datasets
#'
#' @examples
#' data("toy_outbreak_short")
#' names(toy_outbreak_short)
#' toy_outbreak_short
#' 
"toy_outbreak_short"
