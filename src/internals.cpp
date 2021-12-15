#include "internals.h"



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

//   This function returns a vector of indices of cases which could be infector
//   of 'i' (i.e., their infection dates preceed that of 'i', and the genotype of the
//   tree they belong to is similar to "i"'s and its descendents' or unreported). 
//   Only tricky bit here is keep in mind that 't_inf' is indexed from 0 to N-1, 
//   while 'i' and 'alpha' (ancestors) are values from 1 to N.


// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, 
                                            Rcpp::IntegerVector alpha,
                                            Rcpp::StringVector genotype,
                                            Rcpp::StringVector gen_tree,
                                            Rcpp::IntegerVector cluster,
                                            int delta, size_t i) {
  size_t n = cluster.size();
  std::vector<int> out;
  out.reserve(n);
  int ref_t_inf = t_inf[i-1];
  
  // gen_ref: Genotype of i
  Rcpp::String gen_ref = genotype[i-1];
  int n_desc;
  Rcpp::IntegerVector all_descendents;
  // Find all descents from i
  all_descendents = cpp_find_all_descendents(alpha, t_inf, cluster, i);
  n_desc = all_descendents.size();
  int j = 0;
  int j_clust;
  while(gen_ref == "Not attributed" && j<n_desc){
    gen_ref = genotype[all_descendents[j]-1];
    ++j;
  }
  // If there was no genotype reported, then every case infected before i
  // and from the same group of cases ("cluster") can be considered an infector
  if(gen_ref == "Not attributed"){
    for (size_t j = 0; j < n; j++) {
      j_clust = cluster[j]-1;
      if (t_inf[j_clust] < ref_t_inf && t_inf[j_clust] > ref_t_inf - delta) { // offset
        out.push_back(j_clust+1);
      }
    }
  } else {
    // Otherwise, only cases whose tree had unreported genotype or the same 
    // genotype as i and i's descendents can be a cluster
    for (size_t j = 0; j < n; j++) {
      j_clust = cluster[j]-1;
      if (t_inf[j_clust] < ref_t_inf && t_inf[j_clust] > ref_t_inf - delta && 
          (gen_tree[j_clust] == gen_ref 
             || gen_tree[j_clust] == "Not attributed")) { // offset
        out.push_back(j_clust+1);
      }
    }
  }
  return out;
}

// ---------------------------

//  This function compute the spatial log likelihood distribution from parameters a 
//  and b

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::List cpp_log_like_s(Rcpp::NumericVector population, Rcpp::NumericMatrix distance,
                          double a, double b, Rcpp::String spatial) {
  int size_pop = population.size();
  Rcpp::NumericVector population_a(size_pop);
  Rcpp::NumericVector sum_pop(size_pop);
  Rcpp::NumericMatrix nb_move(size_pop, size_pop);
  
  Rcpp::NumericMatrix probs(size_pop, size_pop);
  int j, k;
  
  for(k = 0; k<size_pop; k++)
    population_a[k] = pow(population[k], a);
  if(spatial == "exponential"){
    for(k = 0; k<size_pop; k++){
      for(j = 0; j<size_pop; j++){
        if(j>k){
          nb_move(k, j) = population_a[k]*exp(-b*distance(k,j));
          nb_move(j, k) = nb_move(k,j) * population_a[j] / population_a[k];
          sum_pop[j] += nb_move(k, j);
          sum_pop[k] += nb_move(j, k);
        } else if(j == k){
          nb_move(k, j) = population_a[k];
          sum_pop[k] += nb_move(k, j);
        }
      }
      for(j = 0; j<size_pop; j++){
        probs(j, k) = nb_move(j, k) / sum_pop[k];
        probs(j, k) = (probs(j, k));
      }
    }
  } else if(spatial == "power-law"){
    for(k = 0; k<size_pop; k++){
      for(j = 0; j<size_pop; j++){
        if(j>k){
          nb_move(k, j) = population_a[k]*pow(1+distance(j,k), b);
          nb_move(j, k) = nb_move(k,j) * population_a[j] / population_a[k];
          sum_pop[j] += nb_move(k, j);
          sum_pop[k] += nb_move(j, k);
        } else if(j == k){
          nb_move(k, j) = population_a[k];
          sum_pop[k] += nb_move(k, j);
        }
      }
      for(j = 0; j<size_pop; j++){
        probs(j, k) = nb_move(j, k) / sum_pop[k];
        probs(j, k) = (probs(j, k));
      }
    }
  }
  Rcpp::List new_log_s_dens = Rcpp::List::create(probs);
  return(new_log_s_dens);
}


// [[Rcpp::export()]]
Rcpp::List cpp_log_like(Rcpp::NumericVector population, Rcpp::NumericMatrix distance,
                        Rcpp::NumericMatrix ances, double a, double b, int max_kappa,
                        double gamma, Rcpp::String spatial, int nb_cases) {
  Rcpp::List log_like = cpp_log_like_s(population, distance, a, b, spatial);
  return log_like;
}

// ---------------------------

// This function returns the descendents of a given case 'i' in the current
// ancestries; 'i' is on the scale 1:N. The output is also on the scale 1:N.

// Original R version:

// find.descendents <- function(param, i) {
//   ## find descendents
//     which(param.current$alpha==i)
//  }

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, 
                                         Rcpp::IntegerVector cluster, int i) {
  size_t counter = 0, n = 0;
  size_t length_cluster = cluster.size();
  
  // determine size of output vector and create it
  for (size_t j = 0; j < length_cluster; j++) {
    if (alpha[cluster[j]-1] == i) n++;
  }
  
  Rcpp::IntegerVector out(n);
  
  // fill in output vector
  for (size_t j = 0; j < length_cluster; j++) {
    if (alpha[cluster[j]-1] == i) {
      out[counter++] = cluster[j]; // offset
    }
  }
  return out;
}


// ---------------------------

// This function returns all the descendents of a given case 'i' in the current
// ancestries; 'i' is on the scale 1:N. The output is also on the scale 1:N.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
std::vector<int> cpp_find_all_descendents(Rcpp::IntegerVector alpha, 
                                          Rcpp::IntegerVector t_inf, 
                                          Rcpp::IntegerVector cluster,
                                          int i) {
  // Find how many descendents in total
  int n = alpha.size();
  std::vector<int> all_descents;
  all_descents.reserve(n);
  int k;
  int t_ref = t_inf[i-1];
  size_t length_cluster = cluster.size();
  int j_clust;
  // For every case in cluster, we go up the ancestry, until we find out if 
  // i is the ancestor. If so, j_clust is added to the list of descendents
  for (size_t j=0; j<length_cluster; j++){
    j_clust = cluster[j]-1;
    if(t_inf[j_clust] >= t_ref){
      k = j_clust+1;
      while(k != NA_INTEGER && k!=i){
        k = alpha[k-1];
      }
      if(k==i){
        all_descents.push_back(j_clust+1);
      }
    }
  }
  return all_descents;
}


// ---------------------------

// This function returns all the individuals in the same tree as a given case
// 'i' in the current  ancestries; 'i' is on the scale 1:N. The output is also
// on the scale 1:N.

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_all_tree(Rcpp::IntegerVector alpha,
                                      Rcpp::IntegerVector t_inf,
                                      Rcpp::IntegerVector cluster,
                                      size_t i) {
  // Find how many descendents in total
  Rcpp::IntegerVector all_tree;
  int k = i;
  // Go back to the ancestor of the tree i belongs to,
  while(alpha[k-1] != NA_INTEGER){
    k = alpha[k-1];
  }
  // Find all descendents from the ancestor k
  all_tree = cpp_find_all_descendents(alpha, t_inf, cluster, k);
  
  return all_tree;
}


// ---------------------------

// 

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::String cpp_gen_tree(Rcpp::IntegerVector tree,
                          Rcpp::IntegerVector cluster,
                          Rcpp::StringVector genotype,
                          size_t i){
  
  Rcpp::String gen = genotype[i-1];
  int j = 0;
  int n_tree = tree.size();
  
  while(j < n_tree && gen == "Not attributed"){
    gen = genotype[tree[j] - 1];
    ++j;
  }
  
  return(gen);
}

// ---------------------------

// This function returns a vector of indices of cases which are 'local' to a
// case 'i'. Locality is defined as the following set of cases:

// - 'i'
// - the descendents of 'i'
// - 'alpha[i-1]'
// - the descendents of 'alpha[i]' (excluding 'i')

// where 'alpha' is a IntegerVector storing ancestries. Note that 'i' and
// 'alpha' are on the scale 1:N. 

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha,
                                         Rcpp::IntegerVector cluster,
                                         int i) {
  // determine descendents of 'i':
  Rcpp::IntegerVector desc_i = cpp_find_descendents(alpha, cluster, i);
  size_t n = desc_i.size() + 1; // +1 is to count 'i' itself
  
  // determine descendents of 'alpha[i]':
  Rcpp::IntegerVector desc_alpha_i = cpp_find_descendents(alpha,cluster, alpha[i-1]);
  if (alpha[i-1] != NA_INTEGER) {
    n += desc_alpha_i.size();
  }
  
  // create output
  Rcpp::IntegerVector out(n);
  size_t counter = 0;
  
  // 'i'
  out[counter++] = i;
  
  // 'descendents of 'i'
  for (int j = 0; j < desc_i.size(); j++) {
    out[counter++] = desc_i[j];
  }
  
  if (alpha[i-1] != NA_INTEGER) {
    // alpha[i-1] ...
    out[counter++] = alpha[i-1];
    
    // ... and its descendents
    for (int j = 0; j < desc_alpha_i.size(); j++) {
      if ( desc_alpha_i[j] != i) {
        out[counter++] = desc_alpha_i[j];
      }
    }
  }
  
  return out;
}




// ---------------------------

// This function swaps cases in a transmission tree. The focus case is 'i', and
// is swapped with its ancestor 'x=alpha[i-1]'. In other words the change is
// from: x -> i to i -> x
// Involved changes are:

// - descendents of 'i' become descendents of 'x'
// - descendents of 'x' become descendents of 'i'
// - the infector if 'i' becomes the infector of 'x' (i.e. alpha[x-1])
// - the infector if 'x' becomes 'i'
// - infection time of 'i' becomes that of 'x'
// - infection time of 'x' becomes that of 'i'

// Note on indexing: 'i', 'x', and values of alpha are on the scale 1:N. The
// function's output is a list with new alpha and t_inf.

// Note on forbidden swaps: two types of swaps are excluded:
// - 'i' is imported, so that 'alpha[i-1]' is NA_INTEGER
// - 'x' is imported, so that 'alpha[x-1]' is NA_INTEGER

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export()]]
Rcpp::List cpp_swap_cases(Rcpp::List param, Rcpp::IntegerVector cluster, 
                          Rcpp::IntegerVector move_alpha, int i) {
  Rcpp::IntegerVector alpha_in = param["alpha"];
  Rcpp::IntegerVector t_inf_in = param["t_inf"];
  Rcpp::IntegerVector kappa_in = param["kappa"];
  Rcpp::IntegerVector alpha_out = clone(alpha_in);
  Rcpp::IntegerVector t_inf_out = clone(t_inf_in);
  Rcpp::IntegerVector kappa_out = clone(kappa_in);
  Rcpp::List out;
  
  size_t length_cluster = cluster.size();
  out["alpha"] = alpha_out;
  out["t_inf"] = t_inf_out;
  out["kappa"] = kappa_out;
  
  int j_clust;
  // escape if the case is imported, i.e. alpha[i-1] is NA
  
  if (alpha_in[i-1] == NA_INTEGER) {
    return out;
  }
  
  int x = alpha_in[i-1];

  
  // replace ancestries:
  // - descendents of 'i' become descendents of 'x'
  // - descendents of 'x' become descendents of 'i'
  
  for (size_t j = 0; j < length_cluster; j++) {
    j_clust = cluster[j]-1;
    if (alpha_in[j_clust] == i && move_alpha[j_clust] == TRUE) {
      alpha_out[j_clust] = x;
    } else if (alpha_in[j_clust] == x && move_alpha[j_clust] == TRUE) {
      alpha_out[j_clust] = i;
    }
  }
  
  
  // the ancestor of 'i' becomes an ancestor of 'x'
  
  alpha_out[i-1] = alpha_in[x-1];
  
  
  // 'i' is now the ancestor of 'x'
  alpha_out[x-1] = i;
  
  
  // swap infections times of 'i' and 'x'
  t_inf_out[i-1] =   t_inf_in[x-1];
  t_inf_out[x-1] =   t_inf_in[i-1];
  
  kappa_out[i-1] =   kappa_in[x-1];
  kappa_out[x-1] =   kappa_in[i-1];
  
  
  return out;
}


