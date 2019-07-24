ralpha <- function (t_inf) {
  canBeAnces <- outer(t_inf, t_inf, FUN = "<")
  diag(canBeAnces) <- FALSE
  alpha <- apply(canBeAnces, 2, function(e) ifelse(length(which(e)) > 
                                                     0, sample(c(which(e),
                                                                 which(e)), 
                                                               1), NA))
  return(alpha)
}