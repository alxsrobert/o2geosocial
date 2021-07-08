## This function is not exported and is meant for internal use.
## It plots the sorted likelihood of connection for each individual, and
## shows the threshold (computed either as a quantile of the likelihoods, or 
## as an absolute value). 
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @importFrom graphics par
plot_importations <- function(influences_vect, threshold, config){
  par(mfrow = c(1,1), mar = c(5, 5, 1, 1), bty = "l")
  plot(sort(influences_vect), type = "l", xlab = "Connection", lwd = 2,
       ylab = "Likelihood of connection", 
       ylim = c(0, max(c(-log(0.05) * 5, influences_vect))))
  if(config$outlier_relative == T){
    abline(h = c(-log(0.05)*5, threshold), lty = c(2, 1), lwd = c(1, 1), 
           col = "red")
    legend("topleft", lwd = c(1,1,1), col = c("black", "red", "red"),  
           lty = c(1,1,2),bty = "n",
           legend = c("Likelihood of connection per case",
                      paste0("Threshold (", config$outlier_threshold * 100,
                             " centile)"),
                      "Value if absolute outlier_threshold = 0.05"))
  } else {
    abline(h = threshold, lty = 1, lwd = 1, col = "red")
    legend("topleft", lwd = c(1, 1), col = c("black", "red"),  
           lty = c(1, 1),bty = "n",
           legend = c("Likelihood of connection per case",
                      paste0("Threshold :", round(threshold, 1))))
  }
}
