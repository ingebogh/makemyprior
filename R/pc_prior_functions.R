

## function for computing rank
compute_rank <- function(mat, tol = 1e-12){
  eig <- base::eigen(mat)
  return(sum(Re(eig$values) > tol))
}


# calculates the KLD
calculate_KLD = function(logitW, w0, SigA, SigB, nA, nB, shift_l = -25, shift_u = 25, splitname = ""){

  w = 1/(1+exp(-logitW))

  # Precompute matrix describing shift
  Sig0 = w0*SigA + (1-w0)*SigB
  SigD = SigA-SigB

  # if (sum(SigA != SigB) == 0){
  #   warning("Not sure if this should be possible?")
  #   Sig0 <- SigD <- SigA
  # }

  # # removes last column to avoid singularities (used to be quickfix for breeding program only! not sure if this can be here always or not)
  # if (FALSE && rankMatrix(Sig0)[1] < nrow(Sig0) || rankMatrix(SigD)[1] < nrow(SigD)){
  #   Sig0_red <- Sig0[1:(nrow(Sig0)-1),1:(nrow(Sig0)-1)]
  #   SigD_red <- SigD[1:(nrow(Sig0)-1),1:(nrow(Sig0)-1)]
  #   A <- solve(Sig0_red, SigD_red)
  # } else {
  #   A = solve(Sig0, SigD)
  # }
  A = solve(Sig0, SigD)

  # Calculate eigenvalues
  epair = eigen(A)
  lam = Re(epair$values)

  # Precompute
  trA = sum(lam)
  if(w0 == 0){
    SigTmp = A%*%A
    f2 = sum(diag(SigTmp))
  }
  lamTmp = sort(lam)


  # Go through all values
  val = w
  der = w
  for (i in 1:length(val)){
    sgn = 1
    if (w0 != 0)
      if (logitW[i] < log(w0/(1-w0))) sgn = -1

      if (logitW[i] < shift_l){
        if (w0 == 0){
          wTmp = logitW[i]
          val[i] = log(f2/2)+2*wTmp
          der[i] = log(f2)+wTmp
        } else {
          lamLow = rev(lamTmp)[1:nB]
          lamOK = rev(lamTmp)[-(1:nB)]
          val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lamOK)/(w[i]-w0))) - nB*(logitW[i]-log(w0)))
          der[i] = log(sgn*(trA - sum(lamOK/(1+(w[i]-w0)*lamOK)) - nB*exp(-logitW[i])))
          if (is.infinite(der[i])) der[i] = log(nB)-logitW[i]
        }
      } else if (logitW[i] >= shift_u){
        lamHigh = lamTmp[1:nA]
        lamOK = lamTmp[-(1:nA)]

        val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lamOK)/(w[i]-w0))) + nA*(logitW[i]+log(1-w0)))
        der[i] = log(trA - sum(lamOK/(1+(w[i]-w0)*lamOK)) + nA*exp(logitW[i]))
        if (is.infinite(der[i])) der[i] = log(nA)+logitW[i]
      } else {
        val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lam)/(w[i]-w0))))
        der[i] = log(sgn*(trA - sum(lam/(1+(w[i]-w0)*lam))))
      }
      if (is.na(val[i])){
        stop(paste0(
          "Something went wrong when computing the KLD. Choose other parameters or use a Dirichlet prior for the ",
          splitname, " split.", sep = "", collapse = ""
        ), call. = FALSE)
      }
  }

  # Square root
  val = 0.5*val
  der = -log(2)-val+der

  if (any(c(Im(val), Im(der)) != 0)) {
    warning("Imaginary KLD, using only real value.", call. = TRUE)
  }

  return(list(d = Re(val), diff = Re(der)))
}

# calculates PC prior for when base model is intrinsic
prior_intrinsic_base = function(lW, shape, median, w0){

  if (w0 == 1) median <- 1-median # reversing the prior if the basemodel is at 1
  w = 1/(1+exp(-lW))
  d = sqrt(w)
  # optFun is a function of log(lambda)
  # TODO: find analytically for better stability ?
  optFun = function(x){
    return(abs(qgamma(0.5*pgamma(1, shape = shape, rate = exp(x)), shape = shape, rate = exp(x))-sqrt(median)))
  }
  # lambda = optimize(optFun, interval = c(0, 100))$minimum
  lambda = exp(optimize(optFun, interval = c(-20, 10))$minimum)

  piD = dgamma(d, shape = shape, rate = lambda, log = TRUE)-pgamma(1, shape = shape, rate = lambda, log.p = TRUE)
  piW = piD - log(2) - 0.5*log(w)
  piLW = piW - lW-2*log(1+exp(-lW))

  #return(piLW)

  return(
    if (w0 == 0) piLW else if (w0 == 1) rev(piLW) else stop("Something went wrong when computing the PC prior. Change prior or reformulate your model.", call. = FALSE)
  )

}

# logitW = logit weight values to find kld for
## w0 = amount of variance going to the above node
## shape = how much shrinkage we want on distance scale (this is the shape
## parameter in the gamma distribution, default = 1, cannot be changed by user (yet))
## median = where to put the median, only when base model is 0 or 1, else median is at base model
## wLeft = how much mass to put left of the base model when the base model is not at 0 or 1
## SigA = matrix for the alternative model
## SigB = matrix for the base model
## conc_param = concentation parameter
pc_dual_logit_stable = function(logitW, w0, shape = 1, median, wLeft = 0.5, SigA = NULL, SigB = NULL, conc_param = 0.5, splitname = ""){

  SigA <- as.matrix(SigA)
  SigB <- as.matrix(SigB)

  # TODO: we allow basemodel w0 = 1, we must fix this code
  # until then, we just swap the matrices
  changed_order <- FALSE
  if (w0 == 1){
    w0 <- 0
    tmp <- SigA
    SigA <- SigB
    SigB <- tmp
    changed_order <- TRUE
    median <- 1-median
  }

  # fake point(s) for prior calibration
  if (w0 == 0) {
    logitW = c(log(median/(1-median)), logitW)
  } else {
    #w0 <- 1-w0 # TODO: due to opposite of definition of what is basemodel matrix and not
    logit0 = log(w0/(1-w0))
    shVal = log(0.75/0.25)
    logitW = c(logit0-shVal, logit0+shVal, logitW)
  }

  # add quick-fix for basemodel 0.5 (which does not always work for some reason)
  if (w0 == 0.5) w0 <- 0.49999

  # Get distance and derivative
  if (!is.null(SigA)){
    nA = dim(SigA)[1]-(compute_rank(t(SigA)))
    nB = dim(SigB)[1]-(compute_rank(t(SigB)))
    # nA = dim(SigA)[1]-(base::qr(t(SigA)))$rank
    # nB = dim(SigB)[1]-(base::qr(t(SigB)))$rank
    shift_l <- -25
    shift_u <- 25
    if (nA == 0) shift_u <- 1000
    if (nB == 0 && w0 != 0) shift_l <- -1000
    KLD = calculate_KLD(logitW, w0, SigA, SigB, nA, nB, shift_l, shift_u, splitname)
    # print(shift)
  } else stop("Base model seems to not exist. Change prior or reformulate the model.", call. = FALSE)

  dist = KLD$d
  der = KLD$diff

  # Calculate hyperparameters and remove fake point(s)
  if (w0 == 0){
    medStd = qgamma(0.5, shape = shape, rate = 1)
    lambda = medStd/exp(dist[1])
    logitW = logitW[-1]
    dist = dist[-1]
    der = der[-1]
    #par(mfrow = c(1,2)); plot(dist); plot(der); par(mfrow = c(1,1))
  } else if (0 < w0 && w0 < 1){

    stopifnot(conc_param >= 0.5 && conc_param < 1)

    dValues = dist[1:2]
    dist = dist[-(1:2)]
    der  = der[-(1:2)]
    logitW = logitW[-(1:2)]
    #par(mfrow = c(1,2)); plot(dist); plot(der); par(mfrow = c(1,1))

    # Find correct lambda value
    # optFun = function(x) (conc_param-sum(c(wLeft, 1-wLeft)*pgamma(exp(dValues), shape = shape, rate = exp(x))))^2
    optFun = function(x) return((conc_param-sum(c(wLeft, 1-wLeft)*pgamma(exp(dValues), shape = shape, rate = exp(x))/pgamma(exp(c(dist[1], tail(dist,1))), shape = shape, rate = exp(x))))^2)

    # Grid search first
    xPot = matrix(seq(-10, 10, length.out = 100), nrow = 1)
    yPot = apply(xPot, c(2), optFun)
    xStart = xPot[which.min(yPot)]

    # Optimization
    lambdaNew = nlm(optFun, xStart)
    lambda = exp(lambdaNew$estimate)

    # Output information
    if(optFun(lambdaNew$estimate) > 0.01){
      x = lambdaNew$estimate
      stop(sprintf("Unable to achieve desired concentration %f, got %f instead for split %s. Choose another concentration or another prior.",
                   conc_param,
                   splitname,
                   sum(c(wLeft, 1-wLeft)*pgamma(exp(dValues), shape = shape, rate = exp(x))/pgamma(exp(c(dist[1], tail(dist,1))), shape = shape, rate = exp(x)))),
           call. = FALSE)
    } else {
      x = lambdaNew$estimate
      # print(sprintf("Achieved concentration was %f", sum(c(wLeft, 1-wLeft)*pgamma(exp(dValues), shape = shape, rate = exp(x))/pgamma(exp(c(dist[1], tail(dist,1))), shape = shape, rate = exp(x)))))
    }
    #lambdaNew = optimize(optFun, c(log(1e-10), 10)) #c(log(1e-10), log(1e10))) # larger than exp(10) = 22026.47 does not make sense anyway
    # print(lambdaNew); flush.console()
    # lambda = exp(lambdaNew$minimum)*lambda
    #lambda = exp(lambdaNew$minimum)
  }

  # Calculate prior
  lDens = dgamma(exp(dist), shape = shape, rate = lambda, log = TRUE);
  # print(lambda)

  # Jacobian of distance transformation
  lDens = lDens + der;

  # Jacobian of expit transformation
  lJac = -(logitW + 2*log(1+exp(-logitW)))
  lJac[logitW < -25] = logitW[logitW < -25]
  lDens = lDens + lJac;

  # Force desired split of probability
  idxL = which(logitW < log(w0/(1-w0)))
  if (w0 != 0){
    lDens[idxL] = lDens[idxL] - log(pgamma(exp(dist[1]), shape = shape, rate = lambda)) + log(wLeft)
    lDens[-idxL] = lDens[-idxL] - log(pgamma(exp(tail(dist, 1)), shape = shape, rate = lambda)) + log(1-wLeft)
  }

  # if we had base model at 1, we switched the matrices, and must return the prior in the "correct" (which is the opposite) order
  if (changed_order){
    lDens <- rev(lDens)
  }

  yy <- getSplinePrior(logitW, lDens)
  names(yy)[2] <- "coeffs"
  yy$n_knots <- 122
  # print(integrate(function(x) eval_spline_prior(x, yy), -Inf, Inf)); flush.console()

  return(lDens);

}


getSplinePrior = function(x, y){
  # Calculate interpolating spline
  hmm = splines::interpSpline(x, y)

  # Extract coefficients
  C = hmm$coefficients
  return(list(knots = x, C = C))

}























