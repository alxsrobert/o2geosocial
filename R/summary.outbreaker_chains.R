#' Basic methods for processing outbreaker results
#'
#' Several methods are defined for instances of the class
#' \code{outbreaker_chains}, returned by \code{\link{outbreaker}}, including:
#' \code{print}, \code{plot}
#'
#' @rdname outbreaker_chains
#'
#' @aliases outbreaker_chains print.outbreaker_chains plot.outbreaker_chains
#'   summary.outbreaker_chains
#'
#' @author Initial version by Thibaut Jombart, rewritten by Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @param x an \code{outbreaker_chains} object as returned by \code{outbreaker}.
#' @param n_row the number of rows to display in head and tail; defaults to 3.
#' @param n_col the number of columns to display; defaults to 8.
#' @param type a character string indicating the kind of data to be printed, if 
#' "chain", then the output of the outbreaker function will be printed ; 
#' if "cluster" then the cluster size distribution at every iteration will be printed.
#' @param ... further arguments to be passed to other methods
#'
#' @export
#' @importFrom utils head tail
#'
print.outbreaker_chains <- function(x, n_row = 3, n_col = 8, type = "chain", ...) {
  if(type == "chain"){
    cat("\n\n ///// outbreaker results ///\n")
    cat("\nclass: ", class(x))
    cat("\ndimensions", nrow(x), "rows, ", ncol(x), "columns")
    
    ## process names of variables not shown
    if (ncol(x) > n_col) {
      ori_names <- names(x)
      x <- x[, seq_len(min(n_col, ncol(x)))]
      
      not_shown <- setdiff(ori_names, names(x))
      
      alpha_txt <- paste(not_shown[range(grep("alpha", not_shown))], collapse=" - ")
      t_inf_txt <- paste(not_shown[range(grep("t_inf", not_shown))], collapse=" - ")
      kappa_txt <- paste(not_shown[range(grep("kappa", not_shown))], collapse=" - ")
      
      cat("\nancestries not shown:", alpha_txt)
      cat("\ninfection dates not shown:", t_inf_txt)
      cat("\nintermediate generations not shown:", kappa_txt)
    }
    
    ## heads and tails
    cat("\n\n/// head //\n")
    print(head(as.data.frame(x), n_row))
    cat("\n...")
    cat("\n/// tail //\n")
    print(tail(as.data.frame(x), n_row))
  } else if(type == "cluster"){
    ## Summary of cluster size distribution
    matrix_ances <- t(apply(x[, grep("alpha", colnames(x))], 
                            1, function(X){
                              while(any(!is.na(X[X]))) X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
                              X[!is.na(X)] <- names(X[X[!is.na(X)]])
                              X[is.na(X)] <- names(X[is.na(X)])
                              return(X)
                            }))
    max_clust_size <- max(unlist(apply(matrix_ances, 1, table)))
    # table_tot: Cluster size distribution
    table_tot <- t(apply(matrix_ances, 1, function(X){
      table_clust <- numeric(max_clust_size)
      names(table_clust) <- 1:max_clust_size
      table_clust[names(table(table(X)))] <- table(table(X))
      return(table_clust)
    }))
    if (ncol(table_tot) > n_col) {
      cat("\n Biggest cluster:", max_clust_size)
      cat("\n Cluster not shown:", n_col, "to", ncol(table_tot))
      table_tot <- table_tot[, seq_len(min(n_col, ncol(table_tot)))]
    }
    if((n_row * 2) < nrow(table_tot)){
      cat("\n\n/// head //\n")
      print(head(as.data.frame(table_tot), n_row))
      cat("\n...")
      cat("\n/// tail //\n")
      print(tail(as.data.frame(table_tot), n_row))
    } else
      print(as.data.frame(table_tot))
  } else stop("type should be chain or cluster")
}





#' @rdname outbreaker_chains
#'
#' @param y a character string indicating which element of an
#'   \code{outbreaker_chains} object to plot
#'
#' @param type a character string indicating the kind of plot to be used (see details)
#'
#' @param burnin the number of iterations to be discarded as burnin
#'
#' @param min_support a number between 0 and 1 indicating the minimum support of
#' ancestries to be plotted; only used if 'type' is 'network'
#'
#' @param labels a vector of length N indicating the case labels (must be
#'   provided in the same order used for dates of symptom onset)
#'   
#' @param group_cluster a numeric \code{vector} indicating the breaks to aggregate the
#' cluster size distribution.
#'
#' @export
#'
#' @details \code{type} indicates the type of graphic to plot:
#' 
#' \itemize{
#' 
#' \item \code{trace} to visualise MCMC traces for parameters or augmented data (plots the
#' log-likelihood by default)
#'
#' \item \code{hist} to plot histograms of quantitative values
#'
#' \item \code{density} to plot kernel density estimations of quantitative values
#'
#' \item \code{alpha} to visualise the posterior frequency of ancestries
#'
#' \item \code{network} to visualise the transmission tree; note that
#'  this opens up an interactive plot and requires a web browser with
#'  Javascript enabled; the argument `min_support` is useful to select only the
#'  most supported ancestries and avoid displaying too many links
#'
#' \item \code{kappa} to visualise the distributions generations between cases and their
#' ancestor/infector
#' 
#' \item \code{cluster} to visualise the cluster size distribution, grouped by
#' the value in group_cluster
#'
#' }
#'
#' @importFrom ggplot2 ggplot geom_line geom_point geom_histogram geom_density
#'   geom_violin aes aes_string coord_flip labs guides scale_size_area
#'   scale_x_discrete scale_y_discrete scale_color_manual scale_fill_manual 
#'   geom_bar geom_errorbar
#'
#' @importFrom grDevices xyTable
#' @importFrom graphics plot
#'
plot.outbreaker_chains <- function(x, y = "post",
                                   type = c("trace", "hist", "density", "cluster",
                                            "alpha", "t_inf", "kappa", "network"),
                                   burnin = 0, min_support = 0.1, labels = NULL, 
                                   group_cluster = NULL, ...) {
  
  ## CHECKS ##
  type <- match.arg(type)
  if (!y %in% names(x)) {
    stop(paste(y,"is not a column of x"))
  }
  
  frequency <- low <- med <- up <- categ <- NULL
  
  ## GET DATA TO PLOT ##
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in x")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  ## MAKE PLOT ##
  if (type == "trace") {
    out <- ggplot(x) + geom_line(aes_string(x="step", y = y)) +
      labs(x="Iteration", y = y, title = paste("trace:",y))
  }
  
  if (type == "hist") {
    out <- ggplot(x) + geom_histogram(aes_string(x = y)) +
      geom_point(aes_string(x = y, y = 0), shape="|", alpha = 0.5, size = 3) +
      labs(x = y, title = paste("histogram:",y))
  }
  
  if (type == "density") {
    out <- ggplot(x) + geom_density(aes_string(x = y)) +
      geom_point(aes_string(x = y, y = 0), shape="|", alpha = 0.5, size = 3) +
      labs(x = y, title = paste("density:",y))
  }
  
  if (type =="alpha") {
    alpha <- as.matrix(x[,grep("alpha", names(x))])
    colnames(alpha) <- seq_len(ncol(alpha))
    from <- as.vector(alpha)
    to <- as.vector(col(alpha))
    from[is.na(from)] <- 0
    out_dat <- data.frame(xyTable(from,to))
    names(out_dat) <- c("from", "to", "frequency")
    ## Calculate proportion among ancestries
    get_prop <- function(i) {
      ind <- which(out_dat$to == out_dat$to[i])
      out_dat[[3]][i]/sum(out_dat[[3]][ind])
    }
    ## Return labels, if provided
    get_alpha_lab <- function(axis, labels = NULL) {
      if(is.null(labels)) labels <- seq_len(ncol(alpha))
      if(axis == 'x') return(labels) else
        if(axis == 'y') return(c("Import", labels))
    }
    ## Return custom colors if provided
    get_alpha_color <- function(color = NULL) {
      if(is.null(color)) return(NULL)
      else return(scale_color_manual(values = color))
    }
    ## This joining function is needed so that the '...' argument can be passed
    ## to two functions with different arguments
    get_lab_color <- function(labels = NULL, color = NULL) {
      list(alpha_lab_x = get_alpha_lab('x', labels),
           alpha_lab_y = get_alpha_lab('y', labels),
           alpha_color = get_alpha_color(color))
    }
    
    tmp <- get_lab_color(labels, ...)
    out_dat[3] <- vapply(seq_along(out_dat[[3]]), get_prop, 1)
    out_dat$from <- factor(out_dat$from, levels = c(0, sort(unique(out_dat$to))))
    out_dat$to <- factor(out_dat$to, levels = sort(unique(out_dat$to)))
    out <- ggplot(out_dat) +
      geom_point(aes(x = to, y = from, size = frequency, color = to)) +
      scale_x_discrete(drop = FALSE, labels = tmp$alpha_lab_x) +
      scale_y_discrete(drop = FALSE, labels = tmp$alpha_lab_y) +
      labs(x = 'To', y = 'From', size = 'Posterior\nfrequency') +
      tmp$alpha_color +
      scale_size_area() +
      guides(colour = FALSE)
  }
  
  if (type =="t_inf") {
    get_t_inf_lab <- function(labels = NULL) {
      N <- ncol(t_inf)
      if(is.null(labels)) labels <- 1:N
      return(labels)
    }
    ## Return custom colors if provided
    get_t_inf_color <- function(color = NULL) {
      if(is.null(color)) return(NULL)
      else return(scale_fill_manual(values = color))
    }
    ## This joining function is needed so that the '...' argument can be passed
    ## to two functions with different arguments
    get_lab_color <- function(labels = NULL, color = NULL) {
      list(t_inf_lab_x = get_t_inf_lab(labels),
           t_inf_color = get_t_inf_color(color))
    }
    
    t_inf <- as.matrix(x[,grep("t_inf", names(x))])
    tmp <- get_lab_color(...)
    dates <- as.vector(t_inf)
    cases <- as.vector(col(t_inf))
    out_dat <- data.frame(cases = factor(cases), dates = dates)
    out <- ggplot(out_dat) +
      geom_violin(aes(x = cases, y = dates, fill = cases)) +
      coord_flip() + guides(fill = FALSE) +
      labs(y = 'Infection time', x = NULL) +
      tmp$t_inf_color +
      scale_x_discrete(labels = tmp$t_inf_lab)
  }
  
  if (type == "kappa") {
    get_kappa_lab <- function(labels = NULL) {
      N <- ncol(kappa)
      if(is.null(labels)) labels <- 1:N
      return(labels)
    }
    kappa <- as.matrix(x[,grep("kappa", names(x))])
    generations <- as.vector(kappa)
    cases <- as.vector(col(kappa))
    to_keep <- !is.na(generations)
    generations <- generations[to_keep]
    cases <- cases[to_keep]
    out_dat <- data.frame(xyTable(generations, cases))
    get_prop <- function(i) {
      ind <- which(out_dat$y == out_dat$y[i])
      out_dat[[3]][i]/sum(out_dat[[3]][ind])
    }
    out_dat[3] <- vapply(seq_along(out_dat[[3]]), get_prop, 1)
    names(out_dat) <- c("generations", "cases", "frequency")
    out <- ggplot(out_dat) +
      geom_point(aes(x = generations, y = as.factor(cases), size = frequency, color = factor(cases))) +
      scale_size_area() +
      scale_y_discrete(labels = get_kappa_lab(...)) +
      guides(colour = FALSE) +
      labs(title = "number of generations between cases",
           x = "number of generations to ancestor",
           y = NULL)
  }
  
  if (type == "network") {
    ## extract edge info: ancestries
    alpha <- x[, grep("alpha",names(x)), drop = FALSE]
    from <- unlist(alpha)
    to <- as.vector(col(alpha))
    N <- ncol(alpha)
    edges <- stats::na.omit(data.frame(xyTable(from, to)))
    edges[3] <- edges$number/nrow(alpha)
    names(edges) <- c("from", "to", "value")
    edges <- edges[edges$value > min_support,,drop = FALSE]
    edges$arrows <- "to"

    
    ## ## extract edge info: timing
    ## t_inf <- x[, grep("t_inf",names(x)), drop = FALSE]
    ## mean_time <- apply(t_inf, 2, mean)
    ## mean_delay <- mean_time[edges$to] - mean_time[edges$from]
    ## mean_delay[mean_delay<1] <- 1
    ## edges$label <- paste(round(mean_delay), "days")
    
    ## node info
    find_nodes_size <- function(i) {
      sum(from==i, na.rm = TRUE) / nrow(alpha)
    }
    get_node_lab <- function(labels = NULL) {
      if(is.null(labels)) labels <- 1:N
      return(labels)
    }
    nodes <- data.frame(id = seq_len(ncol(alpha)),
                        label = seq_len(ncol(alpha)))
    nodes$value <- vapply(nodes$id,
                          find_nodes_size,
                          numeric(1))
    nodes$shape <- rep("dot", N)
    nodes$label <- get_node_lab(...)
    
    smry <- summary(x, burnin = burnin)
    is_imported <- is.na(smry$tree$from)
    nodes$shaped[is_imported] <- "star"
    
    ## generate graph
    out <- visNetwork::visNetwork(nodes = nodes, edges = edges, ...)
    out <- visNetwork::visNodes(out, shadow = list(enabled = TRUE, size = 10),
                                color = list(highlight = "red"))
    out <- visNetwork::visEdges(out, arrows = list(
      to = list(enabled = TRUE, scaleFactor = 0.2)),
      color = list(highlight = "red"))
    
  }
  
  if (type == "cluster"){
    matrix_ances <- t(apply(x[x$step > burnin, grep("alpha", colnames(x))], 1, function(X){
      while(any(!is.na(X[X]))) X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
      X[!is.na(X)] <- names(X[X[!is.na(X)]])
      X[is.na(X)] <- names(X[is.na(X)])
      return(X)
    }))
    max_clust_size <- max(unlist(apply(matrix_ances, 1, table)))
    if(!is.null(group_cluster)) {
      max_clust_size <- max(max_clust_size, group_cluster)
      if(max_clust_size != max(group_cluster)) group_cluster = c(group_cluster, max_clust_size)
      if(min(group_cluster) != 0) group_cluster <- c(0, group_cluster)
    }
    # table_tot: Cluster size distribution
    table_tot <- t(apply(matrix_ances, 1, function(X){
      table_clust <- numeric(max_clust_size)
      names(table_clust) <- 1:max_clust_size
      table_clust[names(table(table(X)))] <- table(table(X))
      return(table_clust)
    }))
    if(!is.null(group_cluster)) {
      group_cluster[which.max(group_cluster)] <- group_cluster[which.max(group_cluster)] + 1
      aggreg_vector <- sapply(seq_len(max_clust_size), function(X) sum(X > group_cluster))
      aggreg_clust_size <- t(aggregate(t(table_tot), by = list(aggreg_vector), sum)[,-1])
      colnames(aggreg_clust_size) <- sapply(seq_len(length(group_cluster) - 1), function(X){
        if(group_cluster[X] == group_cluster[X+1] -1) return(as.character(group_cluster[X] + 1)) else
          return(paste0(group_cluster[X] + 1, " - ", group_cluster[X + 1]))
      })
      cluster_size_tot <- t(apply(aggreg_clust_size, 2, function(X) 
        quantile(X, probs = c(0.025, 0.5, 0.975))))
    } else 
      cluster_size_tot <- t(apply(table_tot, 2, function(X)
        quantile(X, probs = c(0.025, 0.5, 0.975))))
    out_dat <- cbind.data.frame(rownames(cluster_size_tot), cluster_size_tot, 
                                stringsAsFactors=FALSE)
    colnames(out_dat) <- c("categ", "low", "med", "up")
    out <- ggplot(out_dat) + geom_bar(aes(x = factor(categ, levels = categ), 
                                          y = med), stat = "identity") +
      geom_errorbar(aes(x = categ, ymin=low, ymax=up), width=.2) +
      guides(fill = FALSE) + labs(y = 'Number of clusters', x = "Cluster size")

  }
  
  return(out)
}


#' @rdname outbreaker_chains
#' @param object an \code{outbreaker_chains} object as returned by \code{outbreaker}.
#' @param group_cluster a numeric \code{vector} indicating the breaks to aggregate the
#' cluster size distribution.
#' @export
#' @importFrom stats median
#' @importFrom stats aggregate
summary.outbreaker_chains <- function(object, burnin = 0, group_cluster = NULL, ...) {
  ## check burnin ##
  x <- object
  if (burnin > max(x$step)) {
    stop("burnin exceeds the number of steps in object")
  }
  x <- x[x$step>burnin,,drop = FALSE]
  
  
  ## make output ##
  out <- list()
  
  
  ## summary for $step ##
  interv <- ifelse(nrow(x)>2, diff(tail(x$step, 2)), NA)
  out$step <- c(first = min(x$step),
                last = max(x$step),
                interval = interv,
                n_steps = length(x$step)
  )
  
  
  ## summary of post, like, prior ##
  out$post <- summary(x$post)
  out$like <- summary(x$like)
  out$prior <- summary(x$prior)
  
  
  ## summary for pi ##
  out$pi <- summary(x$pi)
  ## summary for pi ##
  out$a <- summary(x$a)
  ## summary for pi ##
  out$b <- summary(x$b)
  
  
  ## summary tree ##
  out$tree <- list()
  
  ## summary of alpha ##
  alpha <- as.matrix(x[,grep("alpha", names(x))])
  
  ## function to get most frequent item
  f1 <- function(x) {
    as.integer(names(sort(table(x, exclude = NULL), decreasing = TRUE)[1]))
  }
  out$tree$from <- apply(alpha, 2, f1)
  out$tree$to <- seq_len(ncol(alpha))
  
  ## summary of t_inf ##
  t_inf <- as.matrix(x[,grep("t_inf", names(x))])
  out$tree$time <- apply(t_inf, 2, median)
  
  ## function to get frequency of most frequent item
  f2 <- function(x) {
    (sort(table(x, exclude = NULL), decreasing = TRUE)/length(x))[1]
  }
  out$tree$support <- apply(alpha, 2, f2)
  
  ## summary of kappa ##
  kappa <- as.matrix(x[,grep("kappa", names(x))])
  out$tree$generations <- apply(kappa, 2, function(X) {
    return(median(X[!is.na(X)]))})
  
  ## shape tree as a data.frame
  out$tree <- as.data.frame(out$tree)
  rownames(out$tree) <- NULL
  
  ## Summary of cluster size distribution
  matrix_ances <- t(apply(x[x$step > burnin, grep("alpha", colnames(x))], 
                          1, function(X){
                            while(any(!is.na(X[X]))) X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
                            X[!is.na(X)] <- names(X[X[!is.na(X)]])
                            X[is.na(X)] <- names(X[is.na(X)])
                            return(X)
                          }))
  max_clust_size <- max(unlist(apply(matrix_ances, 1, table)))
  if(!is.null(group_cluster)) {
    max_clust_size <- max(max_clust_size, group_cluster)
    if(max_clust_size != max(group_cluster)) group_cluster = c(group_cluster, max_clust_size)
    if(min(group_cluster) != 0) group_cluster <- c(0, group_cluster)
  }
  # table_tot: Cluster size distribution
  table_tot <- t(apply(matrix_ances, 1, function(X){
    table_clust <- numeric(max_clust_size)
    names(table_clust) <- 1:max_clust_size
    table_clust[names(table(table(X)))] <- table(table(X))
    return(table_clust)
  }))
  if(!is.null(group_cluster)) {
    group_cluster[which.max(group_cluster)] <- group_cluster[which.max(group_cluster)] + 1
    aggreg_vector <- sapply(seq_len(max_clust_size), function(X) sum(X > group_cluster))
    aggreg_clust_size <- t(aggregate(t(table_tot), by = list(aggreg_vector), sum)[,-1])
    colnames(aggreg_clust_size) <- sapply(seq_len(length(group_cluster) - 1), function(X){
      if(group_cluster[X] == group_cluster[X+1] -1) return(as.character(group_cluster[X] + 1)) else
        return(paste0(group_cluster[X] + 1, " - ", group_cluster[X + 1]))
    })
    out$cluster <- apply(aggreg_clust_size, 2, summary)
  } else 
    out$cluster <- apply(table_tot, 2, summary)
  return(out)
}
