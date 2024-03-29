% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_likelihoods.R
\name{custom_likelihoods}
\alias{custom_likelihoods}
\alias{print.custom_likelihoods}
\title{Customise likelihood functions for o2geosocial}
\usage{
custom_likelihoods(...)

\method{print}{custom_likelihoods}(x, ...)
}
\arguments{
\item{...}{a named list of functions, each computing a log-likelihood component.}

\item{x}{an \code{outbreaker_config} object as returned by \code{create_config}.}
}
\value{
A named list of functions with the class \code{custom_likelihood}, each
    implementing a customised log-likelihood components of
    outbreaker. Functions which are not customised will result in a NULL
    component.
}
\description{
This function is used to specify customised likelihood functions for
o2geosocial Custom functions are specified as a named list or series of
comma-separated, named arguments, indicating which log-likelihood component
they compute. Values currently available are:
}
\details{
\itemize{

\item \code{timing_sampling}: the likelihood of sampling times; by default, the function
\code{cpp_ll_timing_sampling} is used.

\item \code{timing_infections}: the likelihood of infection times; by default, the function
\code{cpp_ll_timing_infections} is used.

\item \code{reporting}: the likelihood of the reporting process; by default,
the function \code{cpp_ll_reporting} is used.

\item \code{space}: the likelihood of spatial distances; by default,
the function \code{cpp_ll_space} is used.

\item \code{age}: the likelihood of the age contacts; by default,
the function \code{cpp_ll_age} is used.

}

All log-likelihood functions should have the following arguments, in this
order:

\itemize{

\item \code{data}: a list of named items containing input data as returned by
\code{\link{outbreaker_data}}

\item \code{param}: a list of parameters with the class
\code{\link{create_param}}

}
}
\examples{

## specify a null model by disabling all likelihood components
f_null <- function(data, config = NULL, param, i) {
  return(0.0)
}


null_model <- custom_likelihoods(timing_sampling = f_null,
                                 timing_infections = f_null,
                                 reporting = f_null,
                                 space = f_null,
                                 age = f_null)

null_config <- list(find_import = FALSE,
                    n_iter = 200, gamma = 100, delta = 30,
                    sample_every = 1)

## load data
data("toy_outbreak_short")
dt_cases <- toy_outbreak_short$cases
dt_cases <- dt_cases[order(dt_cases$Date), ][1:15,]
dt_regions <- toy_outbreak_short$dt_regions
all_dist <- geosphere::distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                                        rep(dt_regions$lat, nrow(dt_regions))), 
                                      ncol = 2), 
                               matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                                        rep(dt_regions$lat, each = nrow(dt_regions))),
                                      ncol = 2))

dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
pop_vect <- dt_regions$population
names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region

data <- outbreaker_data(dates = dt_cases$Date, age_group = dt_cases$age_group,
                        region = dt_cases$Cens_tract, population = pop_vect, 
                        distance = dist_mat)

res_null <- outbreaker(data = data,
                       config = null_config,
                       likelihoods = null_model)

}
\author{
Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
}
