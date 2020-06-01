## add convolutions to data$log_w_dens
## rows = kapp avalue
## columns = time interval

log_sum <- function(u, v){
  return(max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v))))
}

log_sum_vec <- function(w){
  total <- w[1]
  if (length(w)<2) return(total)
  
  for (i in 2:length(w)){
    total <- log_sum(total, w[i]);
  }
  return(total)
}

convolve_log <- function(x, y) {
  n <- length(x)
  m <- length(y)
  
  r <- lapply(1:(n+m-1), function(k){
    i <- 1:max(m,n)
    i <- i[((i<=m) & ((k-m+i) <= n)) & ((k-m+i) > 0)]
    log_sum_vec(x[k-m+i]+y[i])
  })
  return(unlist(r))
}
## add convolutions to data$log_w_dens and data$log_a_dens
## rows = kappa value
## columns = time interval
add_convolutions <- function(data, config) {
  ## COMPUTE CONVOLUTIONS IF NEEDED ##
  if (config$max_kappa>1) {
    
    ## first compute convolutions on natural scale
    for (i in 2:config$max_kappa) {
      if (!is.null(data$log_a_dens)) {
        data$log_a_dens[[i]] <- log(exp(data$log_a_dens[[i-1]]) %*%
                                      exp(data$log_a_dens[[1]]))
      }
      if (!is.null(data$log_w_dens)) {
        data$log_w_dens <- rbind(data$log_w_dens,
                                 convolve_log(data$log_w_dens[i-1,],
                                              rev(data$log_w_dens)
                                 )[seq_len(ncol(data$log_w_dens))])
      }
    }
  }
  if(any(is.infinite(data$log_w_dens)))
    data$log_w_dens[is.na(data$log_w_dens)] <- 
      min(data$log_w_dens[is.finite(data$log_w_dens)])
  
  ## name rows/columns (useful if internal debugging needed)
  if (!is.null(data$log_w_dens)) {
    rownames(data$log_w_dens) <- paste("kappa", seq_len(nrow(data$log_w_dens)),
                                       sep="=")
    colnames(data$log_w_dens) <- seq_len(ncol(data$log_w_dens))
  }
  if (!is.null(data$log_a_dens)) {
    names(data$log_a_dens) <- 1:config$max_kappa
  }
  return(data)
}
