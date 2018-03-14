# this file only has a high-level function getting an optimal model fit
# the work-horse functions are written in Rcpp
#
# it uses a least squares methods that has actually been shown
# not to yield optimal fits in some respects, with a better
# alternative from Albert & Shadmehr (2017 or 2018?)

# the following lines added to organize the namespace:

#' @useDynLib RateRate, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# roxygen lines to create man pages:

#' @title Fit the Two-Rate Model To a Dataset.
#' @param reaches A sequence of reach deviations.
#' @param schedule A sequence of feedback manipulations.
#' @param oneTwoRates How many processes to fit? (1 or 2)
#' @param grid How are parameters values for grid search distributed in [0,1]?
#' One of 'uniform' (default), 'restricted' or 'skewed'.
#' @param gridsteps How many values of each parameter are used in grid search?
#' @param checkStability Use additional stability constraints? (default=TRUE)
#' @param method Fitting method, currently one of "Nelder-Mead" (default and
#' very robust, also 'NM', a linear optimization method) or "BFGS" (a
#' quasi-Newton, non-linear method, also 'QN' or 'Quasi-Newton'). See
#' \code{\link{optim}} for details.
#' @param fnscale Fitting is done on function/fnscale, where fnscale is already
#' multiplied by -1 to make it an optimization where appropriate. By default it
#' is set to 1. To make different schedules more comparable, use the largest
#' deviation from zero in the schedule; fnscale=max(abs(schedule), na.rm=T).
#' See \code{\link{optim}} for details.
#' @param verbose Should detailed information be outputted during the fitting?
#' @return The set of parameters that minimizes the difference between model
#' output and the \code{reaches} given the perturbation \code{schedule}.
#' @description Fit Smith et al's (2006) Two-Rate Model to a set of reach
#' deviations and a perturbation schedule.
#' @details This function runs a grid search first, and picks the best 5 are
#' fit with least square optimization after which the best fit is returned.
#' Mean squared error, as given by
#' \code{\link{twoRateReachModelErrors}} is used to determine quality of fit.
#'
#' The sequences of \code{reaches} and the \code{schedule} should have the same
#' length. NAs in the \code{reaches} will be ignored, but in the
#' \code{schedule} they indicate error-clamp trials.
#'
#' The model prediction base on the parameters can be retrieved by evaluating
#' them, based on the perturbation schedule, with
#' \code{\link{twoRateReachModel}}.
#'
#' In the Two-Rate Model of motor learning, the motor output X on a trial t,
#' is the sum of the output of a slow and fast process:
#'
#' X(t) = Xs(t) + Xf(t)
#'
#' And each of these two processes retain part of their previous learning and
#' learn from previous errors:
#'
#' Xs(t) = Rs . Xs(t-1) + Ls . E(t-1)
#'
#' Xf(t) = Rf . Xf(t-1) + Lf . E(t-1)
#'
#' The four parameters Rs, Ls, Rf and Lf are returned, except when a
#' one-process fit is requested, in which case only Rs and Ls are fit.
#'
#' @seealso \code{\link{twoRateReachModelErrors}} and
#' \code{\link{twoRateReachModel}}
#'
#' Smith MA, Ghazizadeh A, Shadmehr R (2006). Interacting Adaptive Processes
#' with Different Timescales Underlie Short-Term Motor Learning. PLoS Biol.
#' 2006 Jun;4(6):e179. \url{https://doi.org/10.1371/journal.pbio.0040179}
#'
#' @examples
#' data("RotAdapt")
#' param <- fitTwoRateReachModel(RotAdapt$reaches,RotAdapt$schedule)
#' param
#'
#' tworatemodel <- twoRateReachModel(param, RotAdapt$schedule)
#' str(tworatemodel)
#'
#' plot(RotAdapt$reaches, ylim=c(-35,35), col='gray')
#' lines(tworatemodel$total)
#' lines(tworatemodel$slow, col='blue')
#' lines(tworatemodel$fast, col='red')
#' @export
fitTwoRateReachModel <- function( reaches,
                                  schedule,
                                  oneTwoRates=2,
                                  verbose=FALSE,
                                  grid='uniform',
                                  gridsteps=7,
                                  checkStability=TRUE,
                                  method='NM',
                                  fnscale=1 ) {

  # find optimal starting parameters for this dataset, with grid search
  # first determine all combinations of parameters to test

  if (grid == 'uniform') {
    grid.steps <- seq(from=0.5/gridsteps, to=1.0-(0.5/gridsteps), by=1/gridsteps)
  }
  if (grid == 'restricted') {
    grid.steps <- seq(from=0.25/gridsteps, to=0.5-(0.25/gridsteps), by=0.5/gridsteps)
  }
  if (grid == 'skewed') {
    grid.steps <- (((exp(seq(from=0,to=2.5,by=2.5/(gridsteps+1))) - 1) / (exp(2.5)-1)))
  }

  stepmax <- min(diff(grid.steps))

  # print(grid.steps)
  # if (fitInitialState) {
  #   if (oneTwoRates == 2) {
  #     state.steps <- seq(from=-6,to=6, by=4)
  #     par.combos <- expand.grid(1-grid.steps,grid.steps,1-grid.steps,grid.steps,state.steps,state.steps)
  #     names(par.combos) <- c('Rs','Ls','Rf','Lf','Is','If')
  #   } else {
  #     state.steps <- seq(from=-6,to=6, by=4)
  #     par.combos <- expand.grid(grid.steps+0.5,grid.steps,state.steps)
  #     names(par.combos) <- c('Rs','Ls','Is')
  #   }
  # } else {

  if (oneTwoRates == 2) {
    par.combos <- expand.grid(1-grid.steps,grid.steps,1-grid.steps,grid.steps)
    names(par.combos) <- c('Rs','Ls','Rf','Lf')
  } else {
    par.combos <- expand.grid(1-grid.steps,grid.steps)
    names(par.combos) <- c('Rs','Ls')
  }
  # } // if you do want to fit initial states...

  # then test all combinations without fitting, collect MSEs
  MSEs <- c()
  N.combos <- dim(par.combos)[1]
  if (verbose) {
    cat(sprintf('evaluating %d parameter combinations, this may take some time\n',N.combos))
  }
  for (comboNo in 1:N.combos) {
    pc <- par.combos[comboNo,]
    par <- unlist(pc)
    # par <- c('Rs'=pc[[1]], 'Ls'=pc[[2]])
    # if (oneTwoRates == 2) {
    #   par <- c(par, c('Rf'=pc[[3]], 'Lf'=pc[[4]]))
    # }
    # if (fitInitialState) {
    #   par <- c(par, 'Is'=pc[[(oneTwoRates*2)+1]], 'If'=pc[[(oneTwoRates*2)+1]])
    # }
    # check if parameters are valid, if not: skip
    # this only makes sense for a two-rate model
    if (oneTwoRates == 2) {
      if (par['Rs'] <= par['Rf']) {
        MSEs <- c(MSEs, NA)
        next()
      }
      if (par['Ls'] >= par['Lf']) {
        MSEs <- c(MSEs, NA)
        next()
      }
      if (checkStability) {

        if ( (((par['Rf'] - par['Lf']) * (par['Rs'] - par['Ls'])) - (par['Lf'] * par['Ls'])) <= 0) {
          MSEs <- c(MSEs, NA)
          next()
        }
        p = par['Rf'] - par['Lf'] - par['Rs'] + par['Ls']
        q = (p^2) + (4 * par['Lf'] * par['Ls'])
        if (((par['Rf'] - par['Lf'] + par['Rs'] - par['Ls']) + q^0.5) >= 2) {
          MSEs <- c(MSEs, NA)
          next()
        }
      }
    }

    # get the error for this instance of the model (no fitting):
    MSEs <- c(MSEs, mean((RateRate::twoRateReachModel(par,schedule)$total - reaches)^2, na.rm=TRUE))
  }
  fitpars <- names(par.combos)
  par.combos$MSE <- MSEs

  # find the best 5 fits and fine-tune those:
  if (verbose) {
    cat('optimizing fit on the best 5\n')
  }
  runfull.idx <- order(par.combos$MSE, decreasing=FALSE)[1:5]
  # and optimize:
  MSEs <- c()
  models <- list()
  for (comboNo in 1:length(runfull.idx)) {
    pc <- par.combos[runfull.idx[comboNo],]
    par <- pc[fitpars]

    if (method %in% c('NM','Nelder-Mead')) {
      # this could be more efficient:
      control <- list('maxit'=10000, 'ndeps'=1e-16, 'fnscale'=1*fnscale )
      comboFit <- stats::optim(par=par, RateRate::twoRateReachModelErrors, gr=NULL, reaches, schedule, checkStability, control=control)
      models[[comboNo]] <- comboFit$par
    }
    if (method %in% c('BFGS','Quasi-Newton','QN')) {
      control <- list('maxit'=10000, 'ndeps'=rep(1e-16, length(par)), 'fnscale'=-1*fnscale )
      comboFit <- stats::optim(par=par, RateRate::twoRateReachModelErrors, gr=NULL, reaches, schedule, checkStability, method='BFGS', control=control)
      models[[comboNo]] <- comboFit$par
    }

    MSEs <- c(MSEs, mean((RateRate::twoRateReachModel(models[[comboNo]],schedule)$total - reaches)^2, na.rm=TRUE))
  }

  best.idx <- sort(MSEs, index.return=TRUE)$ix[1]
  fit <- models[[best.idx]]

  # Right here, we will just plot the output.
  if (verbose) {
    cat('parameter values yielding the best fit:\n')
    # Show the parameter fits:
    print(fit)
  }

  return(fit)

}

# The two functions below are replaced by Rcpp versions, speeding up the (fitting) process:

# twoRateReachModelErrors <- function(par,reaches,schedule,nonfitpar=c()) {
#
#   # This function is called by optim(), providing a set of parameter values, that are varied by optim.
#   # Optim() also passes the data, which we use here to calculate how far of the model with the given
#   # parameters is.
#   # print(par)
#   # First we check if all the parameters are within their constraints.
#   # If not, we don't even fit the model, but return an infinite error.
#
#   # A reminder of the what each parameter value is supposed to be:
#
#   # par[1] = slow retention
#   # par[2] = slow learning
#   # par[3] = fast retention
#   # par[4] = fast learning
#
#   # # if the fast process parameters are not provided, they are given a default value:
#   # if (all(c('Rf','Lf') %in% names(par))) {
#   #   # print('fast rates in par')
#   # } else {
#   #   # print('fast rates not in par')
#   #   if ('Rf' %in% names(nonfitpar)) {
#   #     par <- c(par, 'Rf'=nonfitpar['Rs'])
#   #   } else {
#   #     par <- c(par, 'Rf'=0)
#   #   }
#   #   if ('Lf' %in% names(nonfitpar)) {
#   #     par <- c(par, 'Lf'=nonfitpar['Ls'])
#   #   } else {
#   #     par <- c(par, 'Lf'=1)
#   #   }
#   #   print(par)
#   # }
#
#   # par[1] should be higher than par[3]
#   # par[2] should be lower than par[4]
#
#   # The slow process should retain more than the fast process:
#   if (all(c('Rs','Rf') %in% names(par))) {
#     if (par['Rs'] <= par['Rf']) {
#       # print('Slow process forgets faster than the fast process.')
#       return(Inf)
#     }
#   }
#   # The slow process should learn less than the fast process:
#   if (all(c('Ls','Lf') %in% names(par))) {
#     if (par['Ls'] >= par['Lf']) {
#       # print('Slow process learns faster than the fast process.')
#       return(Inf)
#     }
#   }
#   # All rate-parameters should be between 0 and 1:
#   if (any(par[c('Rs','Ls','Rf','Lf')[c('Rs','Ls','Rf','Lf') %in% names(par)]] < 0.0)) {
#     # print('some parameter values below 0')
#     # print(par)
#     return(Inf)
#   }
#   if (any(par[c('Rs','Ls','Rf','Lf')[c('Rs','Ls','Rf','Lf') %in% names(par)]] > 1.0)) {
#     # print('some parameter values above 1')
#     # print(par)
#     return(Inf)
#   }
#
#   # If we got this far, the parameters are valid.
#   # We can see what the model output is with those parameters.
#   errors <- (twoRateReachModel(c(par,nonfitpar),schdedule)$total - reaches)
#   # print(errors)
#   # print(mean(errors))
#
#   # In fact, we only need to know how far the model is off from the real data, or how large the errors are.
#   # The function twoRateReachModel() needs the parameters to generate model output, and returns on object
#   # with a property called 'output' which is the full model output.
#   # We immediately subtract the actual data from the model output. This should be zero for every trial
#   # if the model were perfect.
#
#   # We want to minimize the root mean squared errors, or RMSE:
#   # RMSE <- sqrt(mean(errors^2, na.rm=TRUE))
#   # No, the MSE, the fit will be better.
#   MSE <- mean(errors^2, na.rm=TRUE)
#
#   # Optim() needs to know what the MSE is:
#   return(MSE)
#
# }
#
#
# twoRateReachModel <- function(par,schedule) {
#
#   # This function generates model output, given a set of parameter values. It is used in fitting the model
#   # to data (finding the parameter values that best describe the data), but also to see what that fitted
#   # model actually does, so we can plot it, once it has been fitted.
#
#   # Set up vectors to store the two process outputs, and total output:
#   slow <- c()
#   fast <- c()
#   total <- c()
#
#   # initial states set to 0 at the start:
#   Xs_t0 <- 0
#   Xf_t0 <- 0
#
#   # # # # # # # # # # #
#   # EXCEPT WHEN THE INITIAL STATE IS TO BE FIT - OR NON-ZERO
#
#   if ('Is' %in% names(par)) {
#     Xs_t0 <- par['Is']
#   }
#   if ('If' %in% names(par)) {
#     Xf_t0 <- par['If']
#   }
#
#   X_t0 <- Xs_t0 + Xf_t0
#
#   doslow <- all(c('Rs','Ls') %in% names(par))
#   dofast <- all(c('Rf','Lf') %in% names(par))
#
#   # We loop through all trials:
#   for (trial in c(1:length(schedule))) {
#
#     # if (is.na()) {
#     #   slow <- c(slow, NA)
#     #   fast <- c(fast, NA)
#     #   output <- c(output, NA)
#     #   next()
#     # }
#
#     # determine error:
#     if (is.na(schedule[trial])) {
#       e_t0 <- 0 # error clamps: there is no learning, but there is forgetting...
#     } else {
#       e_t0 <- (schedule[trial] - X_t0)
#     }
#
#     # calculate what the output should be according to the slow and fast process:
#     if (doslow) {
#       Xs_t1 <- (par['Rs'] * Xs_t0) + (par['Ls'] * e_t0)
#     } else {
#       Xs_t1 <- 0
#     }
#     if (dofast) {
#       Xf_t1 <- (par['Rf'] * Xf_t0) + (par['Lf'] * e_t0)
#     } else {
#       Xf_t1 <- 0
#     }
#
#     # total output is their combined values:
#     X_t1 <- Xs_t1 + Xf_t1
#
#     # save this in the lists/vectors:
#     slow <- c(slow, Xs_t1)
#     fast <- c(fast, Xf_t1)
#     total <- c(total, X_t1)
#
#     # Now change the values for the "previous" model output, to match the current output.
#     # That way it is ready for the next iteration:
#     Xs_t0 <- Xs_t1
#     Xf_t0 <- Xf_t1
#     X_t0 <- X_t1
#
#   }
#
#   # Convert to a data frame:
#   model <- data.frame(slow,fast,total)
#   # With useful column names:
#   colnames(model) <- c('slow','fast','total')
#
#   # Return the data frame that contains the model output for the given parameter values:
#   return(model)
#
# }
