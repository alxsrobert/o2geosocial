calc_s_dens <- function(data, config){
  log_s_dens <- cpp_log_like(data$population, data$distance, data$can_be_ances_reg,
                             config$init_a, config$init_b, config$max_kappa, 
                             config$gamma, config$spatial_method, length(unique(data$region)))
  return(log_s_dens)
}