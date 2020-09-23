
Rcpp::List cpp_move_a(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                      Rcpp::RObject custom_ll,
                      Rcpp::RObject custom_prior);

Rcpp::List cpp_move_b(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                      Rcpp::RObject custom_ll,
                      Rcpp::RObject custom_prior);

Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                       Rcpp::RObject custom_ll,
                       Rcpp::RObject custom_prior);

Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data,
                          Rcpp::RObject list_custom_ll);

Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                          Rcpp::RObject list_custom_ll);

Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                               Rcpp::RObject list_custom_ll);

Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                          Rcpp::RObject list_custom_ll);

Rcpp::List cpp_move_ancestors(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                              Rcpp::RObject list_custom_ll);