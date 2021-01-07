
# these functions are just the old ones, slow and may crash often
# functions for calculating the spline coefficients for the PC prior
# TODO: Update so the functions are more stable/less messy code

# sigma2 <- bdiag(matrix(1, 3, 3), matrix(1, 3, 3), matrix(1, 3, 3))
# sigma1 <- diag(9)
# mm <- reduce_matrices(list(sigma1, sigma2))
# SigB <- sigma1 <- mm[[1]]
# SigA <- sigma2 <- mm[[2]]
# w0 <- 0.00
# w <- seq(0, 1, 0.1)
#
# n2 = dim(sigma2)[1]-(base::qr(t(sigma2)))$rank
# n1 = dim(sigma1)[1]-(base::qr(t(sigma1)))$rank
#
# logit_w <- seq(-20, 20, 1)
#
# plot(logit_w, calculate_KLD2(logit_w, w0, sigma1, sigma2)$d)
# points(logit_w, calculate_KLD(logit_w, w0, sigma2, sigma1, n2, n1, shift = if (w0 == 0) 25 else 1000)$d, pch = 4, col = "red")
#
#
# plot(logit_w, calculate_KLD2(logit_w, w0, sigma1, sigma2)$diff)
# points(logit_w, calculate_KLD(logit_w, w0, sigma2, sigma1, n2, n1, shift = if (w0 == 0) 25 else 1000)$diff, pch = 4, col = "red")

## calculating the KLD
## TODO: do not know what is done in the old code here, but I am trying to just write a new version
## note that the matrices going into this function should already be on the correct form, so we do not check that
## TODO: return error if matrices are not on correct form
# logit_w: logit weight values where KLD should be evaluated
# w0: value of the basemodel, the amount of sigma1 (so we use w0*sigma1 + (1-w0)*sigma2)
# sigma1: one of the matrices
# sigma2: the other matrix
calculate_KLD2 <- function(logit_w, w0, sigma1, sigma2){

  w <- 1/(1+exp(-logit_w))

  n <- nrow(sigma1)

  # basemodel-matrix
  sigma0 <- (1-w0)*sigma1 + w0*sigma2
  sigmaD <- sigma2 - sigma1

  if (w0 == 0){ # in this case sigma0 = sigma1, so we just get a diagonal matrix
    sigma0_inv_sigma1 <- diag(n)
  } else {
    sigma0_inv_sigma1 <- solve(sigma0, sigma1)
  }
  sigma0_inv_sigmaD <- solve(sigma0, sigmaD)
  sigmaD_inv_sigma1 <- solve(sigmaD, sigma1)

  # calculate eigenvalues to avoid finding determinant many times
  eig_0D <- Re(eigen(sigma0_inv_sigmaD)$values)
  eig_D1 <- Re(eigen(sigmaD_inv_sigma1)$values)

  ## solve(A, B) = inv(A) %*% B

  # find eigenvalues:
  logdist <- 0
  logderiv <- 0
  for (i in 1:length(w)){

    trace_tmp <- sum(diag(sigma0_inv_sigma1)) + w[i]*sum(diag(sigma0_inv_sigmaD))

    #logdet_tmp <- log(prod(eig_0D)) + log(prod(w[i] + eig_D1))
    logdet_tmp <- sum(log(abs(eig_0D))) + sum(log(abs(w[i] + eig_D1))) # TODO: ok to use absolute value here?

    trace_deriv <- sum(diag(sigma0_inv_sigmaD))

    logdet_deriv <- sum((eig_D1 + w[i])^(-1))

    kjerne <- trace_deriv - logdet_deriv

    logdist[i] <- 0.5*log(trace_tmp - n - logdet_tmp)
    logderiv[i] <- -log(2) - logdist[i] + log(abs(kjerne)) # TODO: ok to use absolute value here?

  }

  return(list(d = logdist, diff = logderiv))

}


prior_intrinsic_base2 <- function(logit_w, w_median, shape = 1){

  w <- 1/(1+exp(-logit_w))

  # finding the rate parameter (function of lambda)
  lambda_opt <- function(x) {
    lambda <- exp(x)
    # TODO: why is this wrong???????????
    #return(abs( 2*exp(-sqrt(w_median)*lambda) - exp(-lambda) -1 ))
    return(abs(qgamma(0.5*pgamma(1, shape = shape, rate = lambda), shape = shape, rate = lambda)-sqrt(w_median)))
  }

  lambda <- exp(optimize(lambda_opt, interval = c(-20, 10))$minimum)

  # TODO: why
  log_dens_d <- dgamma(sqrt(w), shape = shape, rate = lambda, log = TRUE) - pgamma(1, shape = shape, rate = lambda, log.p = TRUE)
  log_dens_w <- log_dens_d - log(2) - 0.5*log(w) # jacobian from distance to weight (d = sqrt(w))
  log_dens_logit_w <- log_dens_w - logit_w - 2*log1p(-exp(logit_w))

  return(log_dens_logit_w)

}

## calculates density for the PC prior (on dual split, for multisplits we use Dirichlet)
## for non-intrinsic basemodel
# logit_w: logit weight values where KLD should be evaluated
# w0: value of the basemodel, the amount of sigma1 (so we use w0*sigma1 + (1-w0)*sigma2)
# w_median: the value of the median (if basemodel != 0 or 1, this is the same as the basemodel)
# sigma1: one of the matrices
# sigma2: the other matrix
# concentration: concentration parameter around the basemodel when not at 0 or 1 (default = 0.5)
# shape: shape for gamma distribution (NOT WORKING YET) (default = 1)
pc_dual_logit <- function(logit_w, w0, w_median, sigma1, sigma2, concentration = 0.5, shape = 1){

  w_left <- 0.5 # how much mass to the left of the basemodel when basemodel not 0 or 1

  # make points for calibration of the prior
  if (w0 %in% c(0, 1)){
    logit_w_new <- c(log(w_median/(1-w_median)), logit_w)
  } else {
    logit_w_new <- c(log(w0/(1-w0)) + log(0.25/0.75), log(w0/(1-w0) + log(0.75/0.25)))
  }

  kld <- calculate_KLD(logit_w_new, w0, sigma2, sigma1, dim(sigma2)[1]-(base::qr(t(sigma2)))$rank, dim(sigma1)[1]-(base::qr(t(sigma1)))$rank, 25)

  # calculate hyperparameters
  if (w0 %in% c(0, 1)){
    kld_dist <- kld$d[-1]
    kld_der <- kld$diff[-1]
    kld_median <- kld$d[1]
    lambda <- log(2)/exp(kld_median) # the median at KLD-scale is the first value in the kld, median for exp.dist is = log(2)/lambda
  } else {
    kld_dist <- kld$d[-c(1,2)]
    kld_der <- kld$diff[-c(1,2)]
    kld_lower_upper <- kld$d[1:2]

    # find lambda (function of log(lambda):
    lambda_opt <- function(x) (concentration - sum(c(w_left, 1-w_left) * pgamma(exp(kld_lower_upper), shape = shape, rate = exp(x))))
    lambda <- exp(optimize(lambda_opt, c(-20, 10)))

  }

  log_density <- dgamma(exp(kld_dist), shape = shape, rate = lambda, log = TRUE)

  log_density <- log_density + kld_der # adding jacobian for distance transformation

  # jacobian of expit transformation (distance is on weight, not logit-weight, scale)
  log_jac <- -(logit_w + 2*log1p(exp(-logit_w)))
  # log_jac[logit_w < -25] <- logit_w[logit_w < -25] # not necessary I think
  log_density <- log_density + log_jac

  # if basemodel not 0 or 1, we want to get the correct amount of variance left of the median (=basemodel)
  if (!(w0 %in% c(0, 1))){
    lower_ind <- logit_w < log(w0/(1-w0))
    log_density[lower_ind] <- log_density[lower_ind] - log(pgamma(exp(kld_dist[1]), shape = shape, rate = lambda)) + log(w_left)
    log_density[!lower_ind] <- log_density[!lower_ind] - log(pgamma(tail(exp(kld_dist), 1), shape = shape, rate = lambda)) + log(1-w_left)
  }

  # may have to reverse the density, not sure

  return(log_density);

}


# calculates the KLD
calculate_KLD = function(logitW, w0, SigA, SigB, nA, nB, shift = 25){

  w = 1/(1+exp(-logitW))

  # Precompute matrix describing shift
  Sig0 = w0*SigA + (1-w0)*SigB
  SigD = SigA-SigB

  # TODO: THIS IS WRONG!!! Shold not be in this situation
  if (sum(SigA != SigB) == 0){
    warning("Not sure if this should be possible?")
    Sig0 <- SigD <- SigA
  }

  # removes last column to avoid singularities (used to be quickfix for breeding program only! not sure if this can be here always or not)
  if (FALSE && rankMatrix(Sig0)[1] < nrow(Sig0) || rankMatrix(SigD)[1] < nrow(SigD)){
    Sig0_red <- Sig0[1:(nrow(Sig0)-1),1:(nrow(Sig0)-1)]
    SigD_red <- SigD[1:(nrow(Sig0)-1),1:(nrow(Sig0)-1)]
    A <- solve(Sig0_red, SigD_red)
  } else {
    A = solve(Sig0, SigD)
  }

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
  for(i in 1:length(val)){
    sgn = 1
    if(w0 != 0)
      if(logitW[i] < log(w0/(1-w0)))
        sgn = -1

      if(logitW[i] < -shift){
        if(w0 == 0){
          wTmp = logitW[i]
          val[i] = log(f2/2)+2*wTmp
          der[i] = log(f2)+wTmp
        } else{
          lamLow = rev(lamTmp)[1:nB]
          lamOK = rev(lamTmp)[-(1:nB)]
          val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lamOK)/(w[i]-w0))) - nB*(logitW[i]-log(w0)))
          der[i] = log(sgn*(trA - sum(lamOK/(1+(w[i]-w0)*lamOK)) - nB*exp(-logitW[i])))
          if(is.infinite(der[i]))
            der[i] = log(nB)-logitW[i]
        }
      } else if(logitW[i] >= shift){
        lamHigh = lamTmp[1:nA]
        lamOK = lamTmp[-(1:nA)]

        val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lamOK)/(w[i]-w0))) + nA*(logitW[i]+log(1-w0)))
        der[i] = log(trA - sum(lamOK/(1+(w[i]-w0)*lamOK)) + nA*exp(logitW[i]))
        if(is.infinite(der[i]))
          der[i] = log(nA)+logitW[i]
      } else{
        val[i] = log((w[i]-w0)*(trA - sum(log1p((w[i]-w0)*lam)/(w[i]-w0))))
        der[i] = log(sgn*(trA - sum(lam/(1+(w[i]-w0)*lam))))
      }
      #if (is.na(val[i])) browser()
  }

  # Square root
  val = 0.5*val
  der = -log(2)-val+der

  return(list(d = val, diff = der))
}

# calculates PC prior for when base model is intrinsic
prior_intrinsic_base = function(lW, shape, median, w0){

  if (w0 == 1) median <- 1-median # reversing the prior if the basemodel is at 1
  w = 1/(1+exp(-lW))
  d = sqrt(w)
  # optFun is a function of log(lambda) TODO: find analytically for better stability
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
    if (w0 == 0) piLW else if (w0 == 1) rev(piLW) else stop("Should not be possible!")
  )

}

# logitW = logit weight values to find kld for
## w0 = amount of variance going to the alternative model (basemodel at SigB means w0 = 0, basemodel at SigA means w0=0)
## shape = how much shrinkage we want on distance scale (this is the shape
## parameter in the gamma distribution, default = 1)
## median = where to put the median, only when base model is 0 or 1, else median is at base model
## wLeft = how much mass to put left of the base model when the base model is not at 0 or 1
## SigA = matrix for the alternative model
## SigB = matrix for the base model
pc_dual_logit_stable = function(logitW, w0, shape = 1, median, wLeft = 0.5, SigA = NULL, SigB = NULL){

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
    w0 <- 1-w0 # TODO: due to opposite of definition of what is basemodel matrix and not
    logit0 = log(w0/(1-w0))
    shVal = log(0.75/0.25)
    logitW = c(logit0-shVal, logit0+shVal, logitW)
  }

  # add quick-fix for basemodel 0.5 (which does not always work for some reason)
  if (w0 == 0.5) w0 <- 0.49999

  # Get distance and derivative
  if (!is.null(SigA)){
    nA = dim(SigA)[1]-(base::qr(t(SigA)))$rank
    nB = dim(SigB)[1]-(base::qr(t(SigB)))$rank
    KLD = calculate_KLD(logitW, w0, SigA, SigB, nA, nB)
  } else stop("Base model seems to not exist")

  dist = KLD$d
  der = KLD$diff

  # Calculate hyperparameters and remove fake point(s)
  if (w0 == 0){
    medStd = qgamma(0.5, shape = shape, rate = 1)
    lambda = medStd/exp(dist[1])
    logitW = logitW[-1]
    dist = dist[-1]
    der = der[-1]
  } else if (0 < w0 && w0 < 1){
    dValues = dist[1:2]
    dist = dist[-(1:2)]
    der  = der[-(1:2)]
    logitW = logitW[-(1:2)]

    # Find correct lambda value
    optFun = function(x) (0.5-sum(c(wLeft, 1-wLeft)*pgamma(exp(dValues), shape = shape, rate = exp(x))))^2
    lambdaNew = optimize(optFun, c(log(1e-10), log(1e10)))
    # lambda = exp(lambdaNew$minimum)*lambda
    lambda = exp(lambdaNew$minimum)
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
  hmm = interpSpline(x, y)

  # Extract coefficients
  C = hmm$coefficients
  return(list(knots = x, C = C))

}






