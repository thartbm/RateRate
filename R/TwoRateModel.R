# these are the workhorse functions that actually fit the model
#
# it uses a least squares methods that has actually been shown
# not to yield optimal fits in some respects, with a better
# alternative from Albert & Shadmehr (2017 or 2018?)
# we can still improve on the speed by coding in C++

#' Fit the Two-Rate Model To a Dataset.
#'
#' @param reaches A sequence of reach deviations.
#' @param rotations A sequence of feedback manipulations.
#' @param fitInitialState Boolean: is the initial state zero or does it have to be fit.
#' @param oneTwoRates How many processes to fit? (1 or 2)
#' @param verbose Should detailed information be outputted during the fitting?
#' @return The set of parameters that minimizes the difference between model output and the \code{reaches} given the \code{rotations}.
#' @seealso \code{\link{twoRateReachModelErrors}}, which generates the error that is minimized, and uses \code{\link{twoRateReachModel}} to evaluate the model
fitTwoRateReachModel <- function(reaches,rotations,fitInitialState=FALSE,oneTwoRates=2,verbose=FALSE) {

  # Parameters to fit:
  # 1) 'Rs' = slow retention rate
  # 2) 'Ls' = slow learning rate
  # 3) 'Rf' = fast retention rate
  # 4) 'Lf' = fast learning rate
  # 5) 'Is' = initial state of slow process
  # 6) 'If' = initial state of fast process
  # If the slow learning rate is higher than the fast learning rate, the fit should fail.
  # If any of the values learning or retention rates are outside the range 0-1 the fit should also fail.

  # use "optim" to fit the model to the data:
  control <- list('maxit'=10000, 'ndeps'=1e-9 )

  # find optimal starting parameters for this dataset, with grid search
  # first determine all combinations of parameters to test
  nsteps <- 7
  grid.steps <- seq(from=0.5/nsteps, to=1-(0.5/nsteps), by=1/nsteps)
  # print(grid.steps)
  if (fitInitialState) {
    if (oneTwoRates == 2) {
      state.steps <- seq(from=-6,to=6, by=4)
      par.combos <- expand.grid(grid.steps,grid.steps,grid.steps,grid.steps,state.steps,state.steps)
      names(par.combos) <- c('Rs','Ls','Rf','Lf','Is','If')
    } else {
      state.steps <- seq(from=-6,to=6, by=4)
      par.combos <- expand.grid(grid.steps,grid.steps,state.steps)
      names(par.combos) <- c('Rs','Ls','Is')
    }
  } else {
    if (oneTwoRates == 2) {
      par.combos <- expand.grid(grid.steps,grid.steps,grid.steps,grid.steps)
      names(par.combos) <- c('Rs','Ls','Rf','Lf')
    } else {
      par.combos <- expand.grid(grid.steps,grid.steps)
      names(par.combos) <- c('Rs','Ls')
    }
  }

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
    }
    # if (oneTwoRates == 1) {
    #   if (par['Rs'] <= 0) {
    #     MSEs <- c(MSEs, NA)
    #     next()
    #   }
    #   if (par['Ls'] <= 1) {
    #     MSEs <- c(MSEs, NA)
    #     next()
    #   }
    # }
    # get the error for this instance of the model (no fitting):
    MSEs <- c(MSEs, mean((twoRateReachModel(par,rotations)$total - reaches)^2, na.rm=TRUE))
  }
  fitpars <- names(par.combos)
  # print(length(MSEs))
  # print(length(which(is.na(MSEs))))
  # print(MSEs)
  par.combos$MSE <- MSEs

  # find the best 5 fits and finetune those:
  if (verbose) {
    cat('optimizing fit on the best 5\n')
  }
  runfull.idx <- order(par.combos$MSE, decreasing=FALSE)[1:5]
  # print(c('mean'=mean(par.combos$MSE, na.rm=TRUE),'sd'=sd(par.combos$MSE, na.rm=TRUE)))
  # print(par.combos$MSE[runfull.idx])
  # and optimize:
  MSEs <- c()
  models <- list()
  for (comboNo in 1:length(runfull.idx)) {
    pc <- par.combos[runfull.idx[comboNo],]
    # print(pc)
    par <- pc[fitpars]
    # print(par)
    models[[comboNo]] <- optim(par=par, twoRateReachModelErrors, gr=NULL, reaches, rotations, control=control)
    MSEs <- c(MSEs, mean((twoRateReachModel(models[[comboNo]]$par,rotations)$total - reaches)^2, na.rm=TRUE))
  }
  best.idx <- sort(MSEs, index.return=TRUE)$ix[1]
  fit <- models[[best.idx]]

  # # use "optim" to fit the model to the data:
  # control <- list('maxit'=10000, 'ndeps'=1e-9 )
  # fit <- optim(par=reachpar, twoRateReachModelErrors, gr=NULL, reaches, rotations, control=control)

  # The function optim fits the model for us (read it's help page).
  # Optim finds the parameter values, that minimize the error returned by the function we point it to.
  # We give optim some starting parameter values, a function to be minimized, a function for gradients,
  # that we don't use (it is set to NULL), stuff to pass to these functions (the reach data), and then we
  # provide some other arguments, including this 'control' thing. I set that to a value that allows optim
  # to search longer, so that the fit might be improved at the cost of longer computation.

  # Below, we will define the twoRateReachModelErrors function, which has the more important code.

  # Right here, we will just plot the output.
  if (verbose) {
    cat('parameter values yielding the best fit:\n')
    # Show the parameter fits:
    print(fit$par)
  }

  return(fit$par)

}

# The two functions below are replaced by Rcpp versions, speeding up the (fitting) process:

# twoRateReachModelErrors <- function(par,reaches,rotations,nonfitpar=c()) {
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
#   errors <- (twoRateReachModel(c(par,nonfitpar),rotations)$total - reaches)
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
# twoRateReachModel <- function(par,rot) {
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
#   for (trial in c(1:length(rot))) {
#
#     # if (is.na()) {
#     #   slow <- c(slow, NA)
#     #   fast <- c(fast, NA)
#     #   output <- c(output, NA)
#     #   next()
#     # }
#
#     # determine error:
#     if (is.na(rot[trial])) {
#       e_t0 <- 0 # error clamps: there is no learning, but there is forgetting...
#     } else {
#       e_t0 <- (rot[trial] - X_t0)
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
