#' Pre cluster cases in groups according using the genotypes and the arbitrary
#' thresholds gamma (spatial) and delta (temporal).
#' 
#' This function updates the clusters and the  initial tree in the lists data 
#' and config
#' 
#' @author Alexis Robert (\email{alexis.robert@lshtm.ac.uk})
#'
#' @param data A list of data items as returned by \code{outbreaker_data}, or
#' arguments passed to this function.
#'
#' @param config A list of settings as returned by \code{create_config}, or
#' arguments passed to this function.
#'
#' @export
#' 
#' @return
#' 
#' A named list containing two components \code{$data} and
#' \code{$config}. \code{data} data items as returned by \code{outbreaker_data}.
#' \code{config} is a list of settings as returned by \code{create_config}. \cr \cr
#' 
pre_clustering <- function(data, config){
  ## List potential infectors for each case 
  find_infector <- function(i, data, config){
    pot_infect <- which(data$is_cluster == data$is_cluster[i] &
                          data$dates < data$dates[i])
    if(!is.null(config$delta)){
      pot_infect <- which(data$is_cluster == data$is_cluster[i] &
                            data$dates < data$dates[i] &
                            data$dates > data$dates[i]-config$delta)
    } else
      pot_infect <- which(data$is_cluster == data$is_cluster[i] &
                            data$dates < data$dates[i])
    if(data$genotype[i] != "Not attributed")
      pot_infect <- pot_infect[which(data$genotype[pot_infect] == "Not attributed" |
                                       data$genotype[pot_infect] == data$genotype[i])]
    if(!is.null(data$region) && !is.null(data$distance) && 
       !is.null(config$gamma)){
      dist_indiv <- data$distance[, data$region[i]]
      dist_indiv <- dist_indiv[data$region[pot_infect]]
      pot_infect <- pot_infect[which(dist_indiv < config$gamma)]
    }
    return(c(pot_infect))
  }
  infectors <- sapply(seq_along(data$dates), function(X) 
    return(find_infector(X, data = data,config = config)))

  new_clusters <- rep(NA, length(data$dates))
  count <- 1
  
  ## Cluster groups of potential infectors 
  for(i in seq_along(data$dates)){
    pot_inf <- infectors[[i]]
    if(any(!is.na(new_clusters[c(pot_inf, i)]))){
      clust_pot_inf <- new_clusters[c(pot_inf, i)]
      clust_pot_inf <- unique(clust_pot_inf[!is.na(clust_pot_inf)])
      clust_value <- min(clust_pot_inf)
      new_clusters[which(is.element(new_clusters, clust_pot_inf))] <- clust_value
      new_clusters[i] <- clust_value
    } else{
      new_clusters[i] <- count
      count <- count + 1
    }
    new_clusters[pot_inf] <- new_clusters[i]
  }
  
  clust_fact <- as.factor(new_clusters)
  new_clusters <- as.numeric(clust_fact)
  
  data$is_cluster <- new_clusters
  list_clust <- sapply(1:max(new_clusters), 
                       function(X) 
                         return(which(new_clusters == X)))
  data$cluster <- list_clust
  
  config$init_tree <- "star"
  config$init_kappa <- 1
  config <- create_config(config, data = data)
  
  # All cases in isolated chains are imports
  clust_isolated <- which(table(data$is_cluster) == 1)
  config$init_alpha[which(is.element(data$is_cluster, 
                                     clust_isolated))] <- NA
  config$init_kappa[which(is.element(data$is_cluster, 
                                     clust_isolated))] <- NA
  config$move_alpha[which(is.element(data$is_cluster, 
                                     clust_isolated))] <- FALSE
  config$move_kappa[which(is.element(data$is_cluster, 
                                     clust_isolated))] <- FALSE
  
  # Remove unused regions
  if((!is.null(data$region) && !is.null(data$distance) && 
      !is.null(config$gamma))){
    if(!is.infinite(config$gamma)){
      regions <- unique(data$region)
      exclude_regions <- which(colSums(data$distance[regions,] <= config$gamma) == 0)
      if(length(exclude_regions) > 0){
        data$distance <- data$distance[-exclude_regions, -exclude_regions]
        data$population <- data$population[-exclude_regions]
      }
    }
  }
  return(list(data = data, config = config))
}
