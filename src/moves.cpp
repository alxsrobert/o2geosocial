#include <Rcpp.h>
#include <Rmath.h>
#include "moves.h"
#include "internals.h"
#include "likelihoods.h"
#include "priors.h"



// IMPORTANT: ON INDEXING VECTORS AND ANCESTRIES

// Most of the functions implemented here are susceptible to be called from R
// via Rcpp, and are therefore treated as interfaces. This causes a number of
// headaches when using indices of cases defined in R (1:N) to refer to elements
// in Rcpp / Cpp vectors (0:N-1). By convention, we store all data on the
// original scale (1:N), and modify indices whenever accessing elements of
// vectors. In other words, in an expression like 'alpha[j]', 'j' should always
// be on the internal scale (0:N-1).

// In all these functions, 'SEXP i' is an optional vector of case indices, on
// the 1:N scale.

// ---------------------------

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_a(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                      Rcpp::RObject custom_ll = R_NilValue,
                      Rcpp::RObject custom_prior = R_NilValue) {
  // Import
  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector b = param["b"]; // these are just pointers
  
  Rcpp::String spatial = config["spatial_method"];
  Rcpp::IntegerVector region = data["region"];
  Rcpp::NumericMatrix distance = data["distance"];
  Rcpp::NumericVector population = data["population"];
  Rcpp::NumericVector limits = config["prior_a"];
  Rcpp::List new_log_s_dens = new_param["log_s_dens"];
  
  Rcpp::NumericVector new_a = new_param["a"]; // these are just pointers
  Rcpp::NumericMatrix probs = new_log_s_dens[0];
  double sd_a = static_cast<double>(config["sd_a"]);
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;
  
  // Move new_a
  // proposal (normal distribution with SD: config$sd_a)
  new_a[0] += R::rnorm(0.0, sd_a); // new proposed value
  
  if (new_a[0] < limits[0] || new_a[0] > limits[1]) {
    return param;
  }
  new_param["log_s_dens"] = cpp_log_like_s(population, distance, new_a[0], b[0], spatial);
  
  // compute likelihoods
  old_logpost = cpp_ll_space(data, config, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_space(data, config, new_param, R_NilValue, custom_ll);
  
  // compute priors
  old_logpost += cpp_prior_a(param, config, custom_prior);
  new_logpost += cpp_prior_a(new_param, config, custom_prior);
  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);

  // acceptance: the new value is already in a, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }

  return new_param;
}


// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_b(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                      Rcpp::RObject custom_ll = R_NilValue,
                      Rcpp::RObject custom_prior = R_NilValue) {
  // Import
  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector a = param["a"]; // these are just pointers
  Rcpp::NumericVector population = data["population"];
  Rcpp::NumericMatrix distance = data["distance"];
  Rcpp::NumericVector limits = config["prior_b"];
  
  Rcpp::String spatial = config["spatial_method"];
  Rcpp::IntegerVector region = data["region"];
  
  Rcpp::List new_log_s_dens = new_param["log_s_dens"];
  Rcpp::NumericMatrix probs = new_log_s_dens[0];
  Rcpp::NumericVector new_b = new_param["b"]; // these are just pointers
  
  
  double sd_b = static_cast<double>(config["sd_b"]);
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;
  
  // Move new_a
  // proposal (normal distribution with SD: config$sd_a)
  
  new_b[0] += R::rnorm(0.0, sd_b); // new proposed value
  
  
  if (new_b[0] < limits[0] || new_b[0] > limits[1]) {
    return param;
  }
  
  new_param["log_s_dens"] = cpp_log_like_s(population, distance, a[0], new_b[0], spatial);
  
  //compute likelihoods
  old_logpost = cpp_ll_space(data, config, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_space(data, config, new_param, R_NilValue, custom_ll);

  // compute priors
  
  old_logpost += cpp_prior_b(param, config, custom_prior);
  new_logpost += cpp_prior_b(new_param, config, custom_prior);

  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);

  
  // acceptance: the new value is already in b, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  return new_param;
}




// ---------------------------

// movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                       Rcpp::RObject custom_ll = R_NilValue,
                       Rcpp::RObject custom_prior = R_NilValue) {
  
  // deep copy here for now, ultimately should be an arg.
  
  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector pi = param["pi"]; // these are just pointers
  Rcpp::NumericVector new_pi = new_param["pi"]; // these are just pointers
  
  double sd_pi = static_cast<double>(config["sd_pi"]);
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;
  
  
  // proposal (normal distribution with SD: config$sd_pi)
  
  new_pi[0] += R::rnorm(0.0, sd_pi); // new proposed value
  
  
  // automatic rejection of pi outside [0;1]
  
  if (new_pi[0] < 0.0 || new_pi[0] > 1.0) {
    return param;
  }
  
  
  // compute likelihoods
  old_logpost = cpp_ll_reporting(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_reporting(data, new_param, R_NilValue, custom_ll);
  
  
  // compute priors
  
  old_logpost += cpp_prior_pi(param, config, custom_prior);
  new_logpost += cpp_prior_pi(new_param, config, custom_prior);
  
  
  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);
  
  
  // acceptance: the new value is already in pi, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  
  return new_param;
}



// ---------------------------

// Movement of infection dates are +/- 1 from current states. These movements
// are currently vectorised, i.e. a bunch of dates are proposed all together;
// this may not be sustainable for larger datasets. The non-vectorised option
// will be slower and speed-up with C/C++ will be more substantial then.

// This version differs from the initial R implementation in several points:

// 1. all cases are moved
// 2. cases are moved one by one
// 3. movement for each case is +/- 1 time unit

// Notes

// - when computing the timing log-likelihood, the descendents of each
// case are also affected.

// - we generate a new vector 'new_t_inf', which will replace the
// previous pointer defining param["t_inf"].

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data,
                          Rcpp::RObject list_custom_ll = R_NilValue) {
  
  // deep copy here for now, ultimately should be an arg.
  
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector t_inf = param["t_inf"];
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"]; // pointer to t_inf
  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::List cluster_list = data["cluster"];
  Rcpp::IntegerVector cluster_vec = data["is_cluster"];
  Rcpp::IntegerVector local_cases;
  int N = static_cast<int>(data["N"]);
  
  double old_loc_loglike = 0.0, new_loc_loglike = 0.0, p_loc_accept = 0.0;
  
  for (int i = 0; i < N; i++) {
    Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
    local_cases = cpp_find_descendents(param["alpha"], cluster_i, i+1);
    int n_loc = local_cases.size();
    // loglike with current value
    old_loc_loglike = cpp_ll_timing(data, param, i+1, list_custom_ll); // term for case 'i' with offset
    
    // term descendents of 'i'
    if (n_loc > 0) {
      old_loc_loglike += cpp_ll_timing(data, param, local_cases, list_custom_ll);
    }
    
    // proposal (+/- 1)
    new_t_inf[i] += unif_rand() > 0.5 ? 1 : -1; // new proposed value
    // Check the new proposed value is correct
    if(alpha[i] != NA_INTEGER){ 
      if(new_t_inf[i] < new_t_inf[alpha[i] - 1]) new_t_inf[i] = t_inf[i];
    }
    
    if (n_loc> 0){
      for(int k = 0; k < n_loc; k++)
        if(new_t_inf[local_cases[k] - 1] < new_t_inf[i]) new_t_inf[i] = t_inf[i];
    }
    // loglike with new value
    new_loc_loglike = cpp_ll_timing(data, new_param, i+1, list_custom_ll); // term for case 'i' with offset
    
    // term descendents of 'i'
    if (n_loc> 0) {
      new_loc_loglike += cpp_ll_timing(data, new_param, local_cases, list_custom_ll);
    }
    
    
    // acceptance term
    p_loc_accept = exp(new_loc_loglike - old_loc_loglike);
    
    // acceptance: the new value is already in t_inf, so we only act if the move
    // is rejected, in which case we restore the previous ('old') value
    if (p_loc_accept < unif_rand()) { // reject new values
      new_t_inf[i] = t_inf[i];
    } else t_inf[i] = new_t_inf[i];
  }
  
  return new_param;
  
}



// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data, Rcpp::List config, 
                          Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::StringVector genotype = data["genotype"]; // pointer to data$genotype
  Rcpp::StringVector all_gen = clone(genotype); // copy of data$genotype
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to kappa
  Rcpp::IntegerVector new_alpha = new_param["alpha"]; // pointer to new_param$alpha
  Rcpp::IntegerVector new_kappa = new_param["kappa"]; // pointer to new_param$alpha
  Rcpp::List cluster_list = data["cluster"];
  Rcpp::IntegerVector cluster_vec = data["is_cluster"];
  Rcpp::IntegerVector move_alpha = config["move_alpha"]; 
  int delta = config["delta"]; 
  int K = config["max_kappa"]; 
  
  int N = static_cast<int>(data["N"]);
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  Rcpp::IntegerVector possible_ancestors;
  Rcpp::IntegerVector t_inf_i;
  
  for (int i = 0; i < N; i++) {
    if (alpha[i] == NA_INTEGER ){
      Rcpp::IntegerVector tree;
      Rcpp::String gen_tree;
      int n_tree;
      Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
      
      tree = cpp_find_all_tree(alpha, t_inf, cluster_i, i+1);
      gen_tree = cpp_gen_tree(tree, cluster_i, genotype, i +1);
      n_tree = tree.size();
      for(int j = 0; j < n_tree; j++){
        all_gen[tree[j] - 1] = gen_tree;
      }
    }
  }
  
  for (int i = 0; i < N; i++) {
    Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
    t_inf_i = t_inf[cluster_i-1];
    // only non-NA ancestries are moved, if there is at least 1 option
    if (alpha[i] != NA_INTEGER && sum(t_inf_i < t_inf[i]) > 1 && 
        move_alpha[i] == TRUE) { 
      possible_ancestors = cpp_are_possible_ancestors(t_inf, alpha, genotype, 
                                                      all_gen, cluster_i, delta, i+1);
      
      if (possible_ancestors.size()>1){
        // loglike with current value
        old_loglike = cpp_ll_all(data, config, param, i+1, list_custom_ll); // offset
        // proposal (+/- 1)
        new_alpha[i] = possible_ancestors[unif_rand() * possible_ancestors.size()];
        new_kappa[i] = unif_rand() * K + 1;
        // loglike with current value
        new_loglike = cpp_ll_all(data, config, new_param, i+1, list_custom_ll);
        // acceptance term
        p_accept = exp(new_loglike - old_loglike);
        // which case we restore the previous ('old') value
        if (p_accept < unif_rand()) { // reject new values
          new_alpha[i] = alpha[i];
          new_kappa[i] = kappa[i];
        } else {
          int old_ances = alpha[i];
          alpha[i] = new_alpha[i];
          kappa[i] = new_kappa[i];
          
          // Update the genotype in the  tree alpha used to belong to
          Rcpp::IntegerVector tree_a = cpp_find_all_tree(alpha, t_inf, 
                                                         cluster_i, old_ances);
          Rcpp::String gen_tree_a = cpp_gen_tree(tree_a, cluster_i, genotype, 
                                                 old_ances);
          int n_tree_a = tree_a.size();
          int same_tree = 0;
          for(int j = 0; j < n_tree_a; j++){
            if(tree_a[j] == i+1) same_tree = 1;
            all_gen[tree_a[j] - 1] = gen_tree_a;
          }
          
          // Update the genotype in the  tree alpha now belongs to
          if(all_gen[i] != all_gen[alpha[i] - 1] && same_tree == 0){
            Rcpp::IntegerVector tree = cpp_find_all_tree(alpha, t_inf, cluster_i, i + 1);
            Rcpp::String gen_tree = cpp_gen_tree(tree, cluster_i, genotype, i + 1);
            int n_tree = tree.size();
            for(int j = 0; j < n_tree; j++){
              all_gen[tree[j] - 1] = gen_tree;
            }
          }
        }
      }
    }
  }
  
  return new_param;
}






// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_ancestors(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                              Rcpp::RObject list_custom_ll = R_NilValue) {
  // Rcpp::IntegerVector 
  Rcpp::List new_param = clone(param);
  
  Rcpp::StringVector genotype = data["genotype"]; // pointer to data$genotype
  Rcpp::StringVector all_gen = clone(genotype); // copy of data$genotype
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector move_alpha = config["move_alpha"];
  int delta = config["delta"]; 
  
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"];
  
  Rcpp::List cluster_list = data["cluster"];
  Rcpp::IntegerVector cluster_vec = data["is_cluster"];
  
  Rcpp::IntegerVector local_cases;
  Rcpp::IntegerVector desc_index;
  Rcpp::IntegerVector changes(2);
  int N = static_cast<int>(data["N"]);
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0, runif = 0.0;
  Rcpp::IntegerVector all_desc;
  int j_clust;
  
  for (int i = 0; i < N; i++) {
    if (alpha[i] == NA_INTEGER ){
      Rcpp::IntegerVector tree;
      Rcpp::String gen_tree;
      int n_tree;
      Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
      
      tree = cpp_find_all_tree(alpha, t_inf, cluster_i, i+1);
      gen_tree = cpp_gen_tree(tree, cluster_i, genotype, i +1);
      n_tree = tree.size();
      for(int j = 0; j < n_tree; j++){
        all_gen[tree[j] - 1] = gen_tree;
      }
    }
  }
  
  for (int i = 0; i < N; i++) {
    // Only NA ancestries are considered
    if (alpha[i] == NA_INTEGER && move_alpha[i] == TRUE){
      Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
      
      // Draw the new ancestor, can not already be an ancestor
      std::vector<int> possible_ances;
      for (int j = 0; j < cluster_i.size(); j++){
        j_clust = cluster_i[j]-1;
        if (alpha[j_clust] != NA_INTEGER){
          possible_ances.push_back(j_clust+1);
        }
      }
      size_t new_ances;
      Rcpp::IntegerVector possible_index;
      if (possible_ances.size()>0){
        new_ances = possible_ances[unif_rand()*possible_ances.size()];
        
        new_alpha[new_ances-1] = NA_INTEGER;
        
        // Draw i's new index 
        possible_index = cpp_are_possible_ancestors(t_inf, alpha, genotype, 
                                                    all_gen, cluster_i, delta, i+1);
        if (possible_index.size()>0){
          // Likelihood changed for 2 cases: i and new_ances
          changes[0] = i+1;// offset
          changes[1] = new_ances;
          old_loglike = cpp_ll_all(data, config, param, changes, list_custom_ll); 
          
          size_t new_index = possible_index[unif_rand() * possible_index.size()];
          
          new_alpha[i] = new_index;
          new_kappa[i] = new_kappa[new_ances-1];
          new_kappa[new_ances-1] = NA_INTEGER;
          // loglike with current value
          new_loglike = cpp_ll_all(data, config, new_param, changes, list_custom_ll);
          // acceptance term
          p_accept = exp(new_loglike - old_loglike);
          runif = unif_rand();
          // which case we restore the previous ('old') value
          if (p_accept < runif) { // reject new values
            new_alpha[i] = alpha[i];
            new_kappa[i] = kappa[i];
            new_alpha[new_ances-1] = alpha[new_ances-1];
            new_kappa[new_ances-1] = kappa[new_ances-1];
            
          } else {
            int old_ances = alpha[new_ances - 1];
            
            alpha[i] = new_alpha[i];
            kappa[i] = new_kappa[i];
            alpha[new_ances-1] = new_alpha[new_ances-1];
            kappa[new_ances-1] = new_kappa[new_ances-1];
            
            int same_tree = 0;
            
            // Update the genotype in the tree new_ances used to belong to
            Rcpp::IntegerVector tree_a = cpp_find_all_tree(alpha, t_inf, 
                                                           cluster_i, old_ances);
            Rcpp::String gen_tree_a = cpp_gen_tree(tree_a, cluster_i, genotype, 
                                                   old_ances);
            int n_tree_a = tree_a.size();
            for(int j = 0; j < n_tree_a; j++){
              if(tree_a[j] == i+1) same_tree = 1;
              all_gen[tree_a[j] - 1] = gen_tree_a;
            }
            
            // Update the genotype in the tree new_ances now belongs to
            Rcpp::IntegerVector tree_b = cpp_find_all_tree(alpha, t_inf, 
                                                           cluster_i, new_ances);
            Rcpp::String gen_tree_b = cpp_gen_tree(tree_b, cluster_i, genotype, 
                                                   new_ances);
            int n_tree_b = tree_b.size();
            
            for(int j = 0; j < n_tree_b; j++){
              if(tree_b[j] == i+1) same_tree = 1;
              all_gen[tree_b[j] - 1] = gen_tree_b;
            }
            
            // Update the genotype in the tree i now belongs to
            if(all_gen[i] != all_gen[alpha[i] - 1] && same_tree == 0){
              
              Rcpp::IntegerVector tree = cpp_find_all_tree(alpha, t_inf, cluster_i, 
                                                           i + 1);
              Rcpp::String gen_tree = cpp_gen_tree(tree, cluster_i, genotype, 
                                                   i + 1);
              int n_tree = tree.size();
              for(int j = 0; j < n_tree; j++){
                all_gen[tree[j] - 1] = gen_tree;
              }
            }            
          }
        } else{
          new_alpha[new_ances-1] = alpha[new_ances-1];
        }
      }
    }
  }
  return new_param;
}

// ---------------------------

// The basic movement of ancestries (picking an ancestor at random amongst in
// previous cases) makes swaps of ancestries (A->B) to (B->A) very
// difficult. This function addresses the issue. It is computer-intensive, but
// likely a determining factor for faster mixing. Unlike previous versions in
// the original 'outbreaker' package, all cases are 'moved' here. A swap is
// defined as:

// x -> a -> b  becomes a -> x -> b

// Obviously cases are moved one at a time. We need to use local likelihood
// changes for this move to scale well with outbreak size. The complicated bit
// is that the move impacts all descendents from 'a' as well as 'x'.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data, 
                               Rcpp::List config, 
                               Rcpp::RObject list_custom_ll = R_NilValue) {
  
  Rcpp::IntegerVector move_alpha = config["move_alpha"]; // pointer to config$move_alpha
  
  Rcpp::List cluster_list = data["cluster"];
  Rcpp::IntegerVector cluster_vec = data["is_cluster"];
  int N = static_cast<size_t>(data["N"]);
  
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to kappa
  Rcpp::StringVector genotype = data["genotype"]; // pointer to data$genotype
  
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector new_alpha = new_param["alpha"]; 
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"]; 
  Rcpp::IntegerVector new_kappa = new_param["kappa"]; 
  
  Rcpp::IntegerVector swap_alpha, swap_t_inf, swap_kappa;
  Rcpp::List swapinfo; // contains alpha, kappa and t_inf
  Rcpp::IntegerVector local_cases, descendents, desc_i_swap;
  Rcpp::IntegerVector tree_i, tree_j;
  Rcpp::String gen_tree, gen_tree_j;
  Rcpp::IntegerVector vec_swap(N);
  int vect_k, i_swap, size_descendents, n_loc, i, k, j_clust;
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  Rcpp::StringVector all_gen = clone(genotype); // copy of data$genotype
  for (int i = 0; i < N; i++) {
    if (alpha[i] == NA_INTEGER && all_gen[i] == "Not attributed"){
      Rcpp::IntegerVector tree;
      Rcpp::String gen_tree;
      Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
      
      tree = cpp_find_all_tree(alpha, t_inf, cluster_i, i+1);
      gen_tree = cpp_gen_tree(tree, cluster_i, genotype, i +1);
      all_gen[i] = gen_tree;
    }
  }
  
  for (i = 0; i < N; i++) {
    if(vec_swap[i] == 0){
      Rcpp::IntegerVector cluster_i = cluster_list[cluster_vec[i]-1];
      
      descendents = cpp_find_descendents(alpha, cluster_i, i+1);
      size_descendents = descendents.size();
      i_swap = descendents[unif_rand() * size_descendents];
      // The local likelihood is defined as the likelihood computed for the
      // cases affected by the swap; these include:
      // - 'i_swap' 
      // - the descendents of 'i_swap'
      // - 'i'
      // - the descendents of 'i' (other than 'i_swap')

      if(size_descendents > 0 && move_alpha[i_swap - 1] == TRUE && move_alpha[i] == TRUE){
        local_cases = cpp_find_local_cases(alpha, cluster_i, i_swap);
        n_loc = local_cases.size();
        
        // loglike with current parameters
        old_loglike = cpp_ll_all(data, config, param, local_cases, list_custom_ll); // offset
        
        
        // proposal: swap case 'i' and its ancestor
        swapinfo = cpp_swap_cases(param, cluster_i, move_alpha, i_swap);
        swap_alpha = swapinfo["alpha"];
        swap_t_inf = swapinfo["t_inf"];
        swap_kappa = swapinfo["kappa"];
        
        for (k = 0; k < n_loc; k++) {
          vect_k = local_cases[k];
          new_alpha[vect_k - 1] = swap_alpha[vect_k - 1];
          new_t_inf[vect_k - 1] = swap_t_inf[vect_k - 1];
          new_kappa[vect_k - 1] = swap_kappa[vect_k - 1];
        }
        // loglike with new parameters
        
        new_loglike = cpp_ll_all(data, config, new_param, local_cases, list_custom_ll);
        
        // acceptance term
        p_accept = exp(new_loglike - old_loglike);
        // acceptance: change param only if new values is accepted
        if (p_accept >= unif_rand()) { // accept new parameters
          for (int k = 0; k < n_loc; k++) {
            vect_k = local_cases[k];
            alpha[vect_k - 1] = new_alpha[vect_k - 1];
            t_inf[vect_k - 1] = new_t_inf[vect_k - 1];
            kappa[vect_k - 1] = new_kappa[vect_k - 1];
          }
          vec_swap[i_swap - 1] = 1;
        } else{
          for (int k = 0; k < n_loc; k++) {
            vect_k = local_cases[k];
            new_alpha[vect_k - 1] = alpha[vect_k - 1];
            new_t_inf[vect_k - 1] = t_inf[vect_k - 1];
            new_kappa[vect_k - 1] = kappa[vect_k - 1];
          }
        }
      } else if(alpha[i] == NA_INTEGER && move_alpha[i] == FALSE){
        // Special case to swap cases defined as importations (improve the mixing
        // of the trees)
        Rcpp::IntegerVector swap_ances;
        
        for (int j = 0; j < cluster_i.size(); j++){
          j_clust = cluster_i[j]-1;
          if (alpha[j_clust] == NA_INTEGER && j_clust != i &&
              (all_gen[i] == all_gen[j_clust] || all_gen[i] == "Not attributed" ||
              all_gen[j_clust] == "Not attributed")){
            swap_ances.push_back(j_clust+1);
          }
        }
        
        if(swap_ances.size() > 0){
          i_swap = swap_ances[unif_rand() * swap_ances.size()];
          desc_i_swap = cpp_find_descendents(alpha, cluster_i, i_swap);

          int t_inf_i = new_t_inf[i];
          int t_inf_i_swap = new_t_inf[i_swap - 1];
          new_t_inf[i_swap - 1] = t_inf_i;
          new_t_inf[i] = t_inf_i_swap;

          Rcpp::IntegerVector all_cases;
          
          for (k = 0; k < descendents.size(); k++) {
            vect_k = descendents[k];
            new_alpha[vect_k - 1] = i_swap;
            all_cases.push_back(vect_k);
          }
          for (k = 0; k < desc_i_swap.size(); k++) {
            vect_k = desc_i_swap[k];
            new_alpha[vect_k - 1] = i + 1;
            all_cases.push_back(vect_k);
          }
          all_cases.push_back(i_swap);
          all_cases.push_back(i + 1);
          
          old_loglike = cpp_ll_all(data, config, param, all_cases, list_custom_ll);
          new_loglike = cpp_ll_all(data, config, new_param, all_cases, list_custom_ll); // offset

          // acceptance term
          p_accept = exp(new_loglike - old_loglike);
          // acceptance: change param only if new values is accepted
          if (p_accept >= unif_rand()) { // accept new parameters
            for (k = 0; k < descendents.size(); k++) {
              vect_k = descendents[k];
              alpha[vect_k - 1] = new_alpha[vect_k - 1];
            }
            for (k = 0; k < desc_i_swap.size(); k++) {
              vect_k = desc_i_swap[k];
              alpha[vect_k - 1] = new_alpha[vect_k - 1];
            }
            t_inf[i_swap - 1] = new_t_inf[i_swap - 1];
            t_inf[i] = new_t_inf[i];
            vec_swap[i_swap - 1] = 1;
          } else{
            for (k = 0; k < descendents.size(); k++) {
              vect_k = descendents[k];
              new_alpha[vect_k - 1] = alpha[vect_k - 1];
            }
            for (k = 0; k < desc_i_swap.size(); k++) {
              vect_k = desc_i_swap[k];
              new_alpha[vect_k - 1] = alpha[vect_k - 1];
            }
            new_t_inf[i_swap - 1] = t_inf[i_swap - 1];
            new_t_inf[i] = t_inf[i];
          }
        }
      }
    }
  }
  return param;
}

// ---------------------------


// Movement of the number of generations on transmission chains ('kappa') is
// done for one ancestry at a time. As for infection times ('t_inf') we use a
// dumb, symmetric +/- 1 proposal. But because values are typically in a short
// range (e.g. [1-3]) we probably propose more dumb values here. We may
// eventually want to bounce back or use and correct for assymetric proposals.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                          Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  int N = static_cast<size_t>(data["N"]);
  int K = config["max_kappa"];
  size_t jump;
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  
  for (int i = 0; i < N; i++) {
    
    // only non-NA ancestries are moved
    if (alpha[i] != NA_INTEGER) {
      
      // propose new kappa
      jump = (unif_rand() > 0.5) ? 1 : -1;
      new_kappa[i] = kappa[i] + jump;
      
      
      // only look into this move if new kappa is positive and smaller than the
      // maximum value; if not, remember to reset the value of new_kappa to that
      // of kappa, otherwise we implicitely accept stupid moves automatically
      
      if (new_kappa[i] < 1 || new_kappa[i] > K) {
        new_kappa[i] = kappa[i];
      } else {
        
        // loglike with current parameters
        old_loglike = cpp_ll_all(data, config, param, i+1, list_custom_ll);
        
        
        // loglike with new parameters
        new_loglike = cpp_ll_all(data, config, new_param, i+1, list_custom_ll);
        
        // acceptance term
        p_accept = exp(new_loglike - old_loglike);

        // acceptance: change param only if new values is accepted
        if (p_accept >= unif_rand()) { // accept new parameters
          kappa[i] = new_kappa[i];
        }
      }
    }
  }
  
  return param;
}


