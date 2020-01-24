#include <Rmath.h>
#include <Rcpp.h>
#include "internals.h"
#include "likelihoods.h"


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

// This likelihood corresponds to the probability of observing infection dates
// of cases given the infection dates of their ancestors.
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i,
                                Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;
  
  if (custom_function == R_NilValue) {
    
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::NumericMatrix w_dens = data["log_w_dens"];
    int K = w_dens.nrow();
    size_t L = w_dens.ncol();
    size_t delay;
    size_t length_i;
    
    double out = 0.0;
    
    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
        if (alpha[j] != NA_INTEGER) {
          delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
          if (delay < 1 || delay > L) {
            return  R_NegInf;
          }
          if (kappa[j] < 1 || kappa[j] > K) {
            return  R_NegInf;
          }
          
          out += w_dens(kappa[j] - 1, delay - 1);
        }
      }
    } else {
      // only the cases listed in 'i' are retained
      length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
        size_t j = vec_i[k] - 1; // offset
        if (alpha[j] != NA_INTEGER) {
          delay = t_inf[j] - t_inf[alpha[j] - 1]; // offset
          if (delay < 1 || delay > L) {
            return  R_NegInf;
          }
          if (kappa[j] < 1 || kappa[j] > K) {
            return  R_NegInf;
          }
          
          out += w_dens(kappa[j] - 1, delay - 1);
        }
        
      }
    }
    
    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
    
    return Rcpp::as<double>(f(data, param));
  }
}
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, size_t i,
                                Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing_infections(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}





// ---------------------------

// This likelihood corresponds to the probability of reporting dates of cases
// given their infection dates.
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i,
                              Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;
  
  if (custom_function == R_NilValue) {
    
    Rcpp::IntegerVector dates = data["dates"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::NumericVector f_dens = data["log_f_dens"];
    size_t K = f_dens.size();
    double out = 0.0;
    
    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
        size_t delay = dates[j] - t_inf[j];
        if (delay < 1 || delay > K) {
          return  R_NegInf;
        }
        out += f_dens[delay - 1];
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
        size_t j = vec_i[k] - 1; // offset
        size_t delay = dates[j] - t_inf[j];
        if (delay < 1 || delay > K) {
          return  R_NegInf;
        }
        out += f_dens[delay - 1];
      }
    }
    
    return out;
  }  else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
    
    return Rcpp::as<double>(f(data, param));
  }
}
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, size_t i,
                              Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing_sampling(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}



// ---------------------------

// This likelihood corresponds to the probability of observing infection region
// of cases given the infection region of their ancestors.
double cpp_ll_space(Rcpp::List data, Rcpp::List config, 
                    Rcpp::List param, SEXP i,
                    Rcpp::RObject custom_function) {
  int N = static_cast<int>(data["N"]);
  if(N < 2) return 0.0;
  if (custom_function == R_NilValue) {
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector region = data["region"];
    Rcpp::String spatial = config["spatial_method"];
    Rcpp::NumericVector population = data["population"];
    Rcpp::NumericMatrix distance = data["distance"];
    int size_pop = population.size();
    double a = param["a"];
    double b = param["b"];
    double gamma = config["gamma"];
    Rcpp::List movement = cpp_log_like(population, distance, a, b, gamma, spatial);
    Rcpp::NumericMatrix nb_move = movement[0];
    Rcpp::NumericVector sum_pop = movement[1];
    double out = 0.0;
    int region_j;
    int region_index;
    double probs_gen_2;
    double probs_gen_1;
    if (i == R_NilValue) {
      for (int j = 0; j < N; j++) {
        if (alpha[j] != NA_INTEGER) {
          if (kappa[j] < 1 || kappa[j] > size_pop) {
            return  R_NegInf;
          }
          region_j = region[j];
          region_index = region[alpha[j]-1];
          if(kappa[j] == 1){
            probs_gen_1 = nb_move(region_index-1, region_j-1) / sum_pop[region_j-1];
            if(probs_gen_1 <= 0) out += -1000;
            else out += log(probs_gen_1);
          } else if(kappa[j] == 2){
            for(int k = 0; k < size_pop; k++)
              probs_gen_2 += nb_move(region_index-1, k) / sum_pop[k] * 
                nb_move(k, region_j - 1) / sum_pop[region_j-1];
            if(probs_gen_2 == 0) out += -1000;
            else out += log(probs_gen_2);
            probs_gen_2 = 0;
          }
        }
      }
    }
    else {
      int length_i = static_cast<int>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (int k = 0; k < length_i; k++) {
        int j = vec_i[k] - 1; // offset
        if (alpha[j] != NA_INTEGER) {
          if (kappa[j] < 1 || kappa[j] > size_pop) {
            return  R_NegInf;
          }
          region_j = region[j];
          region_index = region[alpha[j]-1];
          if(kappa[j] == 1){
            probs_gen_1 = nb_move(region_index-1, region_j-1) / sum_pop[region_j-1];
            if(probs_gen_1 == 0) out += -1000;
            else out += log(probs_gen_1);
          } else if(kappa[j] == 2){
            for(int k = 0; k < size_pop; k++)
              probs_gen_2 += nb_move(region_index-1, k) / sum_pop[k] * 
                nb_move(k, region_j - 1) / sum_pop[region_j-1];
            if(probs_gen_2 == 0) out += -1000;
            else out += log(probs_gen_2);
            probs_gen_2 = 0;
          }
        }
      } 
    }
    return out;
  }
  else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
    return Rcpp::as<double>(f(data, param));
  }
}
double cpp_ll_space(Rcpp::List data, Rcpp::List config, 
                    Rcpp::List param, int i,
                    Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_space(data, config, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}

// ---------------------------

// This likelihood corresponds to the probability of observing infection age
// of cases given the infection ages of their ancestors.

double cpp_ll_age(Rcpp::List data, Rcpp::List param, SEXP i,
                  Rcpp::RObject custom_function) {
  int N = static_cast<int>(data["N"]);
  if(N < 2) return 0.0;
  
  if (custom_function == R_NilValue) {
    
    Rcpp::IntegerVector alpha = param["alpha"];
    Rcpp::IntegerVector t_inf = param["t_inf"];
    Rcpp::IntegerVector kappa = param["kappa"];
    Rcpp::IntegerVector age_group = data["age_group"];
    Rcpp::List age_dens = data["log_a_dens"];
    Rcpp::NumericMatrix ageref = age_dens[0];
    int K = ageref.nrow();
    int L = ageref.ncol();
    double out = 0.0;
    
    if (i == R_NilValue) {
      for (int j = 0; j < N; j++) {
        if (alpha[j] != NA_INTEGER) {
          if (age_group[j] < 1 || age_group[j] > L) {
            return  R_NegInf;
          }
          if (kappa[j] < 1 || kappa[j] > K) {
            return  R_NegInf;
          }
          Rcpp::NumericMatrix contact = age_dens[kappa[j]-1];
          out += contact(age_group[alpha[j]-1]-1, age_group[j]-1);
        }
      }
    } else {
      // only the cases listed in 'i' are retained
      int length_i = static_cast<int>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (int k = 0; k < length_i; k++) {
        int j = vec_i[k] - 1; // offset
        if (alpha[j] != NA_INTEGER) {
          if (age_group[j] < 1 || age_group[j] > L) {
            return  R_NegInf;
          }
          if (kappa[j] < 1 || kappa[j] > K) {
            return  R_NegInf;
          }
          if(kappa[j]>1)
            Rcpp::NumericMatrix ageref = age_dens[kappa[j]-1];
          out += ageref(age_group[alpha[j]-1]-1, age_group[j]-1);
        }
      }
    }
    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
    
    return Rcpp::as<double>(f(data, param));
  }
}
double cpp_ll_age(Rcpp::List data, Rcpp::List param, int i,
                  Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_age(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}

// ---------------------------

// This likelihood corresponds to the probability of a given number of
// unreported cases on an ancestry.

// The likelihood is given by a geometric distribution with probability 'pi'
// to report a case

// - 'kappa' is the number of generation between two successive cases
// - 'kappa-1' is the number of unreported cases

double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i,
                        Rcpp::RObject custom_function) {
  size_t N = static_cast<size_t>(data["N"]);
  if(N < 2) return 0.0;
  
  double pi = static_cast<double>(param["pi"]);
  Rcpp::IntegerVector kappa = param["kappa"];
  
  // p(pi < 0) = p(pi > 1) = 0
  if (pi < 0.0 || pi > 1.0) {
    return R_NegInf;
  }
  
  if (custom_function == R_NilValue) {
    double out = 0.0;
    
    // all cases are retained
    if (i == R_NilValue) {
      for (size_t j = 0; j < N; j++) {
        if (kappa[j] != NA_INTEGER) {
          if (kappa[j] < 1) {
            return  R_NegInf;
          }
          out += R::dgeom(kappa[j] - 1.0, pi, 1); // first arg must be cast to double
        }
      }
    } else {
      // only the cases listed in 'i' are retained
      size_t length_i = static_cast<size_t>(LENGTH(i));
      Rcpp::IntegerVector vec_i(i);
      for (size_t k = 0; k < length_i; k++) {
        size_t j = vec_i[k] - 1; // offset
        if (kappa[j] != NA_INTEGER) {
          if (kappa[j] < 1) {
            return  R_NegInf;
          }
          out += R::dgeom(kappa[j] - 1.0, pi, 1); // first arg must be cast to double
        }
      }
    }
    
    return out;
  } else { // use of a customized likelihood function
    Rcpp::Function f = Rcpp::as<Rcpp::Function>(custom_function);
    
    return Rcpp::as<double>(f(data, param));
  }
}
double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, size_t i,
                        Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_reporting(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}


// ---------------------------

// This likelihood corresponds to the sums of the separate timing likelihoods,
// which include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling

double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i,
                     Rcpp::RObject custom_functions) {
  
  if (custom_functions == R_NilValue) {
    return cpp_ll_timing_infections(data, param, i) +
      cpp_ll_timing_sampling(data, param, i);
  } else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);
    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]);
    
  }
}
double cpp_ll_timing(Rcpp::List data, Rcpp::List param, size_t i,
                     Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_timing(data, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}




// ---------------------------

// This likelihood corresponds to the sums of the separate likelihoods, which
// include:

// - p(infection dates): see function cpp_ll_timing_infections
// - p(collection dates): see function cpp_ll_timing_sampling
// - p(missing cases): see function cpp_ll_reporting
// - p(space): see function cpp_ll_space 
// - p(age): see function cpp_ll_age 

double cpp_ll_all(Rcpp::List data, Rcpp::List config, 
                  Rcpp::List param, SEXP i,
                  Rcpp::RObject custom_functions) {
  if (custom_functions == R_NilValue) {
    
    return cpp_ll_timing_infections(data, param, i) +
      cpp_ll_timing_sampling(data, param, i) +
      cpp_ll_age(data, param, i) +
      cpp_ll_space(data, config, param, i) +
      cpp_ll_reporting(data, param, i);
    
  }  else { // use of a customized likelihood functions
    Rcpp::List list_functions = Rcpp::as<Rcpp::List>(custom_functions);
    
    return cpp_ll_timing_infections(data, param, i, list_functions["timing_infections"]) +
      cpp_ll_timing_sampling(data, param, i, list_functions["timing_sampling"]) +
      cpp_ll_age(data, param, i, list_functions["age"]) +
      cpp_ll_space(data, config, param, i, list_functions["space"]) +
      cpp_ll_reporting(data, param, i, list_functions["reporting"]);
    
  }
}
double cpp_ll_all(Rcpp::List data, Rcpp::List config, 
                  Rcpp::List param, size_t i,
                  Rcpp::RObject custom_function) {
  SEXP si = PROTECT(Rcpp::wrap(i));
  double ret = cpp_ll_all(data, config, param, si, custom_function);
  UNPROTECT(1);
  return ret;
}