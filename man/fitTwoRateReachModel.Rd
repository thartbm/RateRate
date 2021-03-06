% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoRateModel.R
\name{fitTwoRateReachModel}
\alias{fitTwoRateReachModel}
\title{Fit the Two-Rate Model To a Dataset.}
\usage{
fitTwoRateReachModel(reaches, schedule, oneTwoRates = 2, verbose = FALSE,
  grid = "uniform", gridsteps = 7, checkStability = TRUE, method = "NM",
  fnscale = 1)
}
\arguments{
\item{reaches}{A sequence of reach deviations.}

\item{schedule}{A sequence of feedback manipulations.}

\item{oneTwoRates}{How many processes to fit? (1 or 2)}

\item{verbose}{Should detailed information be outputted during the fitting?}

\item{grid}{How are parameters values for grid search distributed in [0,1]?
One of 'uniform' (default), 'restricted' or 'skewed'.}

\item{gridsteps}{How many values of each parameter are used in grid search?}

\item{checkStability}{Use additional stability constraints? (default=TRUE)}

\item{method}{Fitting method, currently one of "Nelder-Mead" (default and
very robust, also 'NM', a linear optimization method) or "BFGS" (a
quasi-Newton, non-linear method, also 'QN' or 'Quasi-Newton'). See
\code{\link{optim}} for details.}

\item{fnscale}{Fitting is done on function/fnscale, where fnscale is already
multiplied by -1 to make it an optimization where appropriate. By default it
is set to 1. To make different schedules more comparable, use the largest
deviation from zero in the schedule; fnscale=max(abs(schedule), na.rm=T).
See \code{\link{optim}} for details.}
}
\value{
The set of parameters that minimizes the difference between model
output and the \code{reaches} given the perturbation \code{schedule}.
}
\description{
Fit Smith et al's (2006) Two-Rate Model to a set of reach
deviations and a perturbation schedule.
}
\details{
This function runs a grid search first, and picks the best 5 are
fit with least square optimization after which the best fit is returned.
Mean squared error, as given by
\code{\link{twoRateReachModelErrors}} is used to determine quality of fit.

The sequences of \code{reaches} and the \code{schedule} should have the same
length. NAs in the \code{reaches} will be ignored, but in the
\code{schedule} they indicate error-clamp trials.

The model prediction base on the parameters can be retrieved by evaluating
them, based on the perturbation schedule, with
\code{\link{twoRateReachModel}}.

In the Two-Rate Model of motor learning, the motor output X on a trial t,
is the sum of the output of a slow and fast process:

X(t) = Xs(t) + Xf(t)

And each of these two processes retain part of their previous learning and
learn from previous errors:

Xs(t) = Rs . Xs(t-1) + Ls . E(t-1)

Xf(t) = Rf . Xf(t-1) + Lf . E(t-1)

The four parameters Rs, Ls, Rf and Lf are returned, except when a
one-process fit is requested, in which case only Rs and Ls are fit.
}
\examples{
data("RotAdapt")
param <- fitTwoRateReachModel(RotAdapt$reaches,RotAdapt$schedule, method='BFGS')
param

tworatemodel <- twoRateReachModel(param, RotAdapt$schedule)
str(tworatemodel)

plot(RotAdapt$reaches, ylim=c(-35,35), col='gray')
lines(tworatemodel$total)
lines(tworatemodel$slow, col='blue')
lines(tworatemodel$fast, col='red')
}
\seealso{
\code{\link{twoRateReachModelErrors}} and
\code{\link{twoRateReachModel}}

Smith MA, Ghazizadeh A, Shadmehr R (2006). Interacting Adaptive Processes
with Different Timescales Underlie Short-Term Motor Learning. PLoS Biol.
2006 Jun;4(6):e179. \url{https://doi.org/10.1371/journal.pbio.0040179}
}
