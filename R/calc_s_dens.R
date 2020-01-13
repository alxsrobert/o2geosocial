calc_s_dens <- function(data, config){
  log_s_dens <- cpp_log_like(data$population, data$distance, 
                             config$init_a, config$init_b,config$gamma, 
                             config$spatial_method, 
                             length(unique(data$postcode)))
  return(log_s_dens)
}