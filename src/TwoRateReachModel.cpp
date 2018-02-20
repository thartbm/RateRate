#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @title Evaluate the Two-Rate Model.
//' @param par A list of named parameters.
//' @param schedule A sequence of feedback manipulations
//' @return A dataframe with three columns, the states of the \code{fast} and \code{slow} process and the \code{total} model output. Each row has the states and output for every trial in the rotation schedule.
//' @description Evaluate the Two-Rate Model Given Parameters and a Perturbation Schedule.
//' @seealso \code{\link{fitTwoRateReachModel}} and \code{\link{twoRateReachModelErrors}}
//' @export
// [[Rcpp::export]]
DataFrame twoRateReachModel(NumericVector par, NumericVector schedule) {

  // number of values to pre allocate for output
  // and number of loop iterations
  int n = schedule.size();

  // initialize the states:
  float Xs_t0 = 0;
  float Xf_t0 = 0;

  // initialize other variables:
  float e_t0 = 0;
  float Xs_t1 = 0;
  float Xf_t1 = 0;
  float X_t1 = 0;

  // check if the fast and slow state are to be initialized at some value:
  if(par.containsElementNamed("Is")) {
    Xs_t0 = par["Is"];
  }
  if(par.containsElementNamed("If")) {
    Xf_t0 = par["If"];
  }

  // only evaluate a process if both its retention and learning rate exist:
  bool doslow = TRUE;
  bool dofast = TRUE;
  if(par.containsElementNamed("Rs")==FALSE) {
    doslow = FALSE;
  }
  if(par.containsElementNamed("Ls")==FALSE) {
    doslow = FALSE;
  }
  if(par.containsElementNamed("Rf")==FALSE) {
    dofast = FALSE;
  }
  if(par.containsElementNamed("Lf")==FALSE) {
    dofast = FALSE;
  }

  // initialize total state, so that errors can be calculated:
  float X_t0 = Xs_t0 + Xf_t0;

  // vectors will be combined into a dataframe at the end:
  // IntegerVector a = df["a"];
  NumericVector slow(n);
  NumericVector fast(n);
  NumericVector total(n);


  for(int t = 0; t < n; ++t) {

    // calculate previous error:
    if( NumericVector::is_na(schedule[t]) ) {
      // if there is an NA, it was an error-calmped trial
      // so the error should be zero
      e_t0 = 0;
    } else {
      e_t0 = schedule[t] - X_t0;
    }

    // calculate current response:
    if (doslow) {
      Xs_t1 = (par["Rs"] * Xs_t0) + (par["Ls"] * e_t0);
    } else {
      Xs_t1 = 0;
    }
    if (dofast) {
      Xf_t1 = (par["Rf"] * Xf_t0) + (par["Lf"] * e_t0);
    } else {
      Xf_t1 = 0;
    }

    // total output is the sum of the two processes:
    X_t1 = Xs_t1 + Xf_t1;

    // store values in the vectors:
    slow[t] = Xs_t1;
    fast[t] = Xf_t1;
    total[t] = X_t1;

    // prepare for next iteration:
    Xs_t0 = Xs_t1;
    Xf_t0 = Xf_t1;
    X_t0 = X_t1;

  }

  // return a new data frame
  return DataFrame::create(_["slow"]= slow, _["fast"]= fast, _["total"]= total);

}

//' @title Return the Mean Squared Error Between a Two-Rate Model and a Dataset.
//' @param par A list of named parameters.
//' @param reaches A list of reach deviations.
//' @param schedule A sequence of feedback manipulations
//' @return The average squared error of the model's prediction of the reach deviations given the rotation schedule and model parameters.
//' @description Calculates the average squared error of the Two-Rate Model's prediction, given a set of parameters and a rotation schedule.
//' @seealso \code{\link{fitTwoRateReachModel}} and \code{\link{twoRateReachModel}}
//' @export
// [[Rcpp::export]]
double twoRateReachModelErrors(NumericVector par, NumericVector reaches, NumericVector schedule) {

  // first we check if the input parameters make sense
  // if not, we return infinity:
  double inf = std::numeric_limits<double>::infinity();


  // only evaluate rates if both its slow and fast version exist:
  // also check if each parameter is within bounds
  bool checkR = TRUE;
  bool checkL = TRUE;
  if(par.containsElementNamed("Rs")==FALSE) {
    checkR = FALSE;
  } else {
    if (par["Rs"] > 1.0) {
      return(inf);
    }
    if (par["Rs"] < 0.0) {
      return(inf);
    }
  }
  if(par.containsElementNamed("Rf")==FALSE) {
    checkR = FALSE;
  } else {
    if (par["Rf"] > 1.0) {
      return(inf);
    }
    if (par["Rf"] < 0.0) {
      return(inf);
    }
  }
  if(par.containsElementNamed("Ls")==FALSE) {
    checkL = FALSE;
  } else {
    if (par["Ls"] > 1.0) {
      return(inf);
    }
    if (par["Ls"] < 0.0) {
      return(inf);
    }
  }
  if(par.containsElementNamed("Lf")==FALSE) {
    checkL = FALSE;
  } else {
    if (par["Lf"] > 1.0) {
      return(inf);
    }
    if (par["Lf"] < 0.0) {
      return(inf);
    }
  }

  if (checkR) {
    // fast retention should not be larger than slow retention
    double Rf = par["Rf"];
    double Rs = par["Rs"];
    if (Rf > Rs) {
      return(inf);
    }
  }
  if (checkL) {
    // slow learning should not be larger than fast learning
    double Ls = par["Ls"];
    double Lf = par["Lf"];
    if (Ls > Lf) {
      return(inf);
    }
  }

  // parameters checked, we can now evaluate the model with the parameters
  DataFrame model = twoRateReachModel(par, schedule);
  // get only the total model output
  NumericVector total = model["total"];
  // these are the errors of the model in predicting the behavior for each trial:
  NumericVector errors = reaches - total;
  // now we square and then sum those errors:
  double sumerrors2 = std::inner_product(errors.begin(), errors.end(), errors.begin(), 0.0);
  // and return that, divided by the number of trials, the MSE:
  return(sumerrors2/schedule.size());

}
