% made by hand, fro Rcpp function
\name{twoRateReachModel}
\alias{twoRateReachModel}
\title{Evaluate the Two-Rate Model.}
\usage{
twoRateReachModel(par, rot)
}
\arguments{
\item{par}{A list of named parameters.}

\item{rot}{A sequence of feedback manipulations.}

}
\value{
A dataframe with three columns, the states of the \code{fast} and \code{slow} process and the \code{total} model output. Each row has the states and output for every trial in the rotation schedule.
}
\description{
Evaluate the Two-Rate Model Given Parameters and a Rotation Schedule.
}

