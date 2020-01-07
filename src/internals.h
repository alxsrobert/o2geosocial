#ifndef MEASLESOUTBREAKER_INTERNALS_H
#define MEASLESOUTBREAKER_INTERNALS_H

#include <Rcpp.h>

std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, Rcpp::IntegerVector alpha, Rcpp::StringVector genotype,Rcpp::IntegerVector cluster, size_t i);

size_t cpp_sample1(Rcpp::IntegerVector x);

Rcpp::List cpp_log_like(Rcpp::NumericVector population, Rcpp::NumericMatrix distance, double a, double b, double gamma, Rcpp::String spatial, int nb_cases);

size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, Rcpp::IntegerVector alpha, Rcpp::StringVector genotype, Rcpp::IntegerVector cluster, size_t i);

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, Rcpp::IntegerVector cluster, size_t i);

std::vector<int> cpp_find_all_descendents(Rcpp::IntegerVector alpha, Rcpp::IntegerVector t_inf, Rcpp::IntegerVector cluster, size_t i);

Rcpp::IntegerVector cpp_find_all_tree(Rcpp::IntegerVector alpha, Rcpp::IntegerVector cluster, size_t i);

Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, Rcpp::IntegerVector cluster, size_t i);

Rcpp::List cpp_swap_cases(Rcpp::List param, Rcpp::IntegerVector cluster, size_t i);

#endif
