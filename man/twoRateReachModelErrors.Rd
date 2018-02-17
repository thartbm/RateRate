twoRateReachModelErrors(NumericVector par, NumericVector reaches, NumericVector rotations)
% made by hand, fro Rcpp function
\name{twoRateReachModelErrors}
\alias{twoRateReachModelErrors}
\title{Average Square Errors of the Two-Rate Model.}
\usage{
twoRateReachModelErrors(par, reaches, rotations)
}
\arguments{
\item{par}{A list of named parameters.}

\item{reaches}{A list of reach deviations.}

\item{rotations}{A sequence of feedback manipulations.}

}
\value{
The average squared error of the model's prediction of the reach deviations given the rotation schedule and model parameters.
}
\description{
Calculates the average squared error of the Two-Rate Model's prediction, given a set of parameters and a rotation schedule.
}

