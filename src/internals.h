#ifndef O2GEOSOCIAL_INTERNALS_H
#define O2GEOSOCIAL_INTERNALS_H

#include <Rcpp.h>

std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, Rcpp::IntegerVector alpha, Rcpp::StringVector genotype, Rcpp::StringVector gen_tree, Rcpp::IntegerVector cluster, size_t i);

Rcpp::List cpp_log_like(Rcpp::NumericVector population, Rcpp::NumericMatrix distance, Rcpp::NumericMatrix ances, double a, double b, int max_kappa, double gamma, Rcpp::String spatial, int nb_cases);

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, Rcpp::IntegerVector cluster, int i);

std::vector<int> cpp_find_all_descendents(Rcpp::IntegerVector alpha, Rcpp::IntegerVector t_inf, Rcpp::IntegerVector cluster, int i);

Rcpp::IntegerVector cpp_find_all_tree(Rcpp::IntegerVector alpha, Rcpp::IntegerVector t_inf, Rcpp::IntegerVector cluster, size_t i);

Rcpp::String cpp_gen_tree(Rcpp::IntegerVector tree, Rcpp::IntegerVector cluster, Rcpp::StringVector genotype, size_t i);

Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, Rcpp::IntegerVector cluster, int i);

Rcpp::List cpp_swap_cases(Rcpp::List param, Rcpp::IntegerVector cluster, int i);

#endif
