library(Matrix)
library(MASS)
require(snowfall)
sfLibrary(glmnet)
library(e1071)
library(lattice)
library(snowfall)
library(stringr)
library(ggplot2)
library(glmnet)
library(energy)
library(msda)
library(mclust)
library(sparcl)
library(combinat)


#############0. Basic function#########
matpower <- function(a,alpha){
  small <- .00000001
  if (length(c(a))==1) {return(a^(alpha))} else {
    p1<-nrow(a)
    eva<-eigen(a)$values
    eve<-eigen(a)$vectors
    eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
    index<-(1:p1)[abs(eva)>small]
    evai<-eva
    evai[index]<-(eva[index])^(alpha)
    ai<-eve%*%diag(evai,length(evai))%*%t(eve)
    return(ai)}
}

proj<-function(v)
{return(v%*%matpower(t(v)%*%v,-1)%*%t(v))}

stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/(n-1)
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

deal_list<-function(x)
{
  return(as.vector(x)[-1])
}

allzero<-function(X)  
{
  return(all(X==0))
}

sumzero<-function(X)
{
  return(length(which(X!=0)))
}

thresholding<-function(X,threshold)
{
  return(all(sum(X^2)<=threshold))
}

slicing<-function(y,H) 
{
  n<-length(y)
  ord=order(y)
  y=y[ord]
  c=as.integer(n/H)
  if (length(levels(as.factor(y)))>H)
  {
    ytilde<-rep(0,H+1)
    ytilde[1]<-min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-y[h*c+1]
    }  
  }
  if (length(levels(as.factor(y)))<=H)
  {
    H <- length(levels(as.factor(y)))
    ytilde<-rep(0,H+1)
    ytilde[1]=min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-min(y[y>ytilde[h]])
    }
  } 
  ytilde[H+1]=max(y)+1
  prop<-rep(1,H)
  for (i in 1:H)
  {
    prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  return(res)
}

mhat_dr<-function(x,y,H){
  n<-nrow(x)
  p<-ncol(x)
  z<-stand(x)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  ind<-matrix(0,n,H)
  zbar<-matrix(0,p,H)
  for (j in 1:(H-1))
  {
    ind[,j]<-((y >= ytilde[j])&(y < ytilde[j+1]))
    zbar[,j]<- (t(z)%*%(ind[,j]))/sum(ind[,j])
  }
  ind[,H]<-(y >= ytilde[H])
  zbar[,H]<- (t(z)%*%(ind[,H]))/sum(ind[,H])
  A<-matrix(0,p,p)
  B<-matrix(0,p,p)
  C<-0
  for (q in 1:H)
  {
    Z<-(t(z))[,ind[,q]==1]-zbar[,q]  
    A<-A + prop[q]*((Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q]))%*%  
                      (Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q])) - diag(1,p))
    B<-B + sqrt(prop[j])*(zbar[,q]%*%t(zbar[,q]))
    C<-C + sqrt(prop[j])*(t(zbar[,q])%*%zbar[,q])
  }
  C<-as.vector(C)
  M<-2*A + 2*(B%*%B) + 2*B*C
  return(M)
}
##################################################################

###1.ADMM_algorithm############################################

AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# Cut small values in a matrix to zero. 
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}

# Evaluation based on distance correlation. (Szekely et al., 2007)
eval_dc <- function(Beta, x, y){
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      dcor(x %*% mat, y)
    }
  })
  return(result)
}

# Compute M and U matrices from observation data.
#########
# Input:
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in 1,...,H.
# type: Specifying the specific SEAS method. "sir" means SEAS-SIR, "intra" means SEAS-Intra and "pfc" means SEAS-PFC.
# FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3).
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
MU <- function(x, y, yclass=NULL, type='sir', FUN = NULL, categorical = FALSE, H = 5){
  if(is.null(yclass)){ # Construct the discretized response
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  cls <- sort(unique(yclass))
  nclass <- length(cls)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  mu <- colMeans(x)
  x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  M <- crossprod(x_c/sqrt(nobs)) # sample covariance of X
  
  if(type == 'sir'){
    U <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      U[, i] <- colMeans(x_c[yclass == cls[i],, drop=FALSE])
    }
  }else if(type == 'intra'){
    y_c <- y - mean(y)
    U <- matrix(0, nvars, nclass)
    lb <- quantile(y_c, 0.1)[[1]]
    ub <- quantile(y_c, 0.9)[[1]]
    y_c <- sapply(y_c, cut_func, lb = lb, ub = ub)
    for (i in 1:nclass){
      y_copy <- y_c
      y_copy[yclass!=cls[i]] <- 0
      U[, i] <- (1/nobs) * t(x_c) %*% (y_copy - mean(y_copy))
    }
  }else if(type == 'pfc'){
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) # the default function
    else Fmat <- t(sapply(y, FUN))
    Fmat_mean <- colMeans(Fmat)
    Fmat_c <- Fmat - matrix(Fmat_mean, NROW(Fmat), NCOL(Fmat), byrow = TRUE) # centered function f
    lb <- apply(Fmat_c, 2, quantile, 0.1)
    ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)){
      Fmat_c[,i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
    }
    U <- (1/nobs)*(t(x_c) %*% Fmat_c)
  }
  list(M = M, U = U, nclass = nclass, prior=prior)
}

# Cut extreme values in the samples. This function is used in MU function.
cut_func <- function(x, lb, ub){
  if(x < lb){
    return(lb)
  } else if(x > ub){
    return(ub)
  } else{
    return(x)
  }
}

# Estimate the rank of a matrix.
rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

# Subspace distance, defined in (19)
subspace <- function(A,B){
  if(is.vector(A)) A <- as.matrix(A)
  if(is.vector(B)) A <- as.matrix(B)
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

## ------------------------------------------------ ##
## The utility functions imported from R package 'msda'.
## These functions are used in 'msda' and 'cv.msda' functions.
formatoutput <- function(fit, maxit, pmax, p, H) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0){
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i){
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H)
      theta_i[ja,] <- a
      theta_i
    })
  }
  else{
    theta <- lapply(seq(nalam), function(x){matrix(0, p, H)})
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, " iterations; solutions for larger lambdas returned.\n", sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam){
  llam <- log(lam)
  if(length(llam) >= 3){lam[1] <- exp(2 * llam[2] - llam[3])}
  lam
}

seas <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', M = NULL, U = NULL, nobs = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), FUN = NULL, eps = 1e-3, maxit = 1e+3, ...){
  
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
        nclass <- as.integer(length(unique(yclass)))
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    if(is.null(gamma)){
      gamma <- c(10,30,50)
    }
    if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
      fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=10, lambda.factor=0.5, nfolds = 5, FUN = FUN, maxit=1e3)
      M <- fit_1$M  # The M matrix based on the full data.
      U <- fit_1$U  # The U matrix based on the full data.
      id_max_msda <- fit_1$id
      lam1_max_msda <- fit_1$lam_max  # The optimal lambda from msda
      beta_msda <- as.matrix(fit_1$beta)  # The optimal matrix from msda
      if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
      if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
      if (all(lam2 == 0)){
        lam2 <- 0
        warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
      }
    }else{
      MU_out <- MU(x, y, yclass, type, FUN)
      M <- MU_out$M
      U <- MU_out$U
    }
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(lam1) || is.null(lam2) || is.null(gamma)) stop("Sequences lam1, lam2 or gamma is missing.")
    if(is.null(nobs)) stop("Missing nobs.")
    nvars <- NCOL(M)
  }
  
  ## Error code
  code <- 0
  
  if(is.vector(lam1) && (length(lam1) == 1) && (lam1 == 0) && is.vector(lam2) && (length(lam2) == 1) && (lam2 == 0)){ # For degenerate case where lambda1 = lambda2 = 0, return B = M^{-1} U directly.
    B <- solve(M) %*% U
    if(is.null(d)) beta <- svd(B)$u
    else if(d == 0) beta <- matrix(0, nrow(Bnew), ncol(Bnew))
    else beta <- svd(B)$u[,1:d,drop=FALSE]
    vec <- as.vector(beta)
    vec[abs(vec) < 1e-3] <- 0
    beta <- matrix(vec, nrow(beta), ncol(beta))
    rank <- NCOL(beta)
    output <- list(beta = beta, B = B, rank = rank, lam1 = lam1, lam2 = lam2, code = code)
  }
  else{
    # Fit with admm function
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, ...)
    B_l <- fit$B
    beta_l <- fit$beta
    if (all(sapply(beta_l, is.null))){
      code <- 1
      warning("No converged results returned.")
      return(list(beta = beta_l, code = code))
    }
    rank_l <- fit$rank
    s_l <- fit$s
    step_l <- fit$step
    time_l <- fit$time
    if(length(B_l) == 1){
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}

cv.seas <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', lambda.factor=0.5, nlambda=10, nfolds = 5, foldid = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), plot = FALSE, FUN = NULL, eps = 1e-3, maxit = 1e+3, trace.it = FALSE, ...){
  # The inputs and outputs are similar to the ones in 'seas' functions. Only the different ones are listed below.
  # Inputs:
  # =======
  # nfolds: The number of folds in the cross-validation.
  # plot: If TRUE, (1) plot the evaluation for each tuning parameter in 'msda' function; (2) in each cross-validation data fold, plot the evaluation for each tuning parameter in 'seas' function.
  # trace.it: If TRUE, print the process of cross-validation.
  # 
  # Outputs:
  # ========
  # Refer to the outputs in 'seas' function.
  
  start_time <- Sys.time()
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  
  M <- U <- M_fold <- U_fold <- NULL
  nobs <- dim(x)[1]
  
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    if (nfolds < 3) stop("nfolds must be larger than 3")
    if (nfolds > nobs) stop("nfolds is larger than the sample size")
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }else{nfolds <- length(unique(foldid))}
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
    fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=nlambda, lambda.factor=lambda.factor, foldid = foldid, FUN = FUN, maxit=1e3, plot = plot)
    M <- fit_1$M  # The M matrix based on the full data.
    U <- fit_1$U  # The U matrix based on the full data.
    M_fold <- fit_1$M_fold  # The M matrix based on each four out of five folds.
    U_fold <- fit_1$U_fold  # The U matrix based on each four out of five folds.
    id_max_msda <- fit_1$id
    lam1_max_msda <- fit_1$lam_max
    beta_msda <- as.matrix(fit_1$beta)
    if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)){
      lam2 <- 0
      warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # The number of errors
  nerr <- 0
  code <- 0
  
  end_time <- Sys.time()
  time1 <- difftime(end_time, start_time, units = "secs")
  ## Record time1: estimate M and U matrices
  
  out_all <- lapply(1:nfolds, function(k){ # Cross-validation
    if(trace.it) cat(sprintf("Fold: %d/%d\n", k, nfolds))
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]
    
    if(is.null(M_fold) || is.null(U_fold)){
      x_train <- x[foldid!=k,,drop=FALSE]
      y_train <- y[foldid!=k]
      yclass_train <- yclass[foldid!=k]
      # Fit with seas function
      fit_fold <- seas(x_train, y_train, yclass = yclass_train, type = type, FUN = FUN, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    }
    else{
      fit_fold <- seas(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(foldid!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    }
    
    err <- 0
    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    
    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
    
    if(plot){ # If true, plot the evaluation for each tuning parameter
      dat <- data.frame(x = 1:length(eval_fold), y = eval_fold, rank = as.factor(rank_l))
      g <- ggplot(dat, aes(x = x, y = y, col = rank))+
        geom_point(size = 1)+
        labs(title=paste0("Fold ", k), x="", y="Distance correlation")+
        theme_bw()
      #print(g)
    }
    out <- list(eval_fold, err)
    out
  })
  
  # Combine the evaluations from each fold.
  eval_all <- do.call(rbind, lapply(out_all, "[[", 1))
  errs <- do.call(c, lapply(out_all, "[[", 2))
  nerr <- sum(errs)
  
  if((nerr != 0) && (nerr != nfolds)){
    code <- 3
    warning(paste0("No converged results returned in", nerr, "folds."))
  }else if(nerr == nfolds){
    code <- 4
    warning("No converged results returned in any fold.")
    return(list(beta = NULL, code = code))
  }
  
  if(is.vector(eval_all)){
    eval_all <- as.matrix(eval_all)
  }
  
  # Compute the cross-validation mean and standard error.
  cvm <- colMeans(eval_all, na.rm=TRUE)
  cvsd <- sqrt(colMeans(scale(eval_all, cvm, FALSE)^2, na.rm = TRUE)/(nfolds-1))
  
  # Select the best tuning parameter.
  id_max <- which.max(cvm)
  id_lam1 <- ceiling(id_max/(n2*n3))
  id_gamma <- ceiling((id_max-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_max-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_max <- lam1[id_lam1]
  gamma_max <- gamma[id_gamma]
  lam2_max <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  start_time <- Sys.time()
  # Refit with the selected tuning parameters.
  if(is.null(M) || is.null(U)){
    fit <- seas(x, y, yclass = yclass, type = type, FUN = FUN, lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }else{
    fit <- seas(M = M, U = U, nobs = NROW(x), lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }
  
  if(fit$code != 0){
    code <- 5
    warning("The estimated beta is null.")
    return(list(beta = NULL, code = code))
  }
  
  B <- fit$B
  beta <- fit$beta
  rank <- fit$rank
  
  end_time <- Sys.time()
  time2 <- difftime(end_time, start_time, units = "secs") # We do not include the time for tuning parameter selection.
  # Record time: one run with the selected tuning parameter
  
  time <- time1 + time2
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval_all, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_max = lam1_max, lam2_max = lam2_max, gamma_max = gamma_max, code = code, time = time)
  output
}

# admm algorithm function
admm <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, ...){
  # Inputs:
  # =======
  # M: The M matrix in optimization problem. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
  # U: The U matrix in optimization problem.
  # nobs: The number of observations.
  # nvars: The number of predictors.
  # lam1, lam2, gamma: The user-specified sequences of tuning parameter lambda_1, lambda_2 and gamma.
  # eps: The tolerance of convergence in ADMM algorithm. The value is passed to 'admm' function.
  # maxit: The maximal iterations in ADMM algorithm. The value is passed to 'admm' function.
  # d: The true structural dimension. The default is NULL.
  # 
  # Outputs:
  # ========
  # beta: A list containing the estimated basis matrices of central subspace.
  # B: A list containing estimated B matrices.
  # rank: A vector containing the estimated ranks.
  # s: A vector containing the estimated sparsity levels.
  # step: A vector containing the number of iterations to converge for each tuning parameter.
  # nlam: The number of converged matrices.
  
  # since the user is required to provide lam1, then set flmin=1
  if(is.null(dim(U)))
  {
    U<-cbind(U,rep(0,length(U)))
  }
  opts <- list(...)
  if(is.null(opts$nlam)) opts$nlam <- as.integer(1)
  if(is.null(opts$H)) opts$H <- as.integer(dim(U)[2])
  if(is.null(opts$nvars)) opts$nvars <- as.integer(nvars)
  if(is.null(opts$pf)) opts$pf <- as.double(rep(1, nvars))
  if(is.null(opts$dfmax)) opts$dfmax <- as.integer(nobs)
  if(is.null(opts$pmax)) opts$pmax <- as.integer(min(nobs*2+20, nvars))
  if(is.null(opts$flmin)) opts$flmin <- as.double(1)
  if(is.null(opts$eps_inner)) opts$eps_inner <- as.double(1e-04)
  if(is.null(opts$maxit_inner)) opts$maxit_inner <- as.integer(1e+6)
  if(is.null(opts$sml)) opts$sml <- as.double(1e-6)
  if(is.null(opts$verbose)) opts$verbose <- as.integer(FALSE)
  if(is.null(opts$nalam)) opts$nalam <- integer(1)
  if(is.null(opts$theta)) opts$theta <- double(opts$pmax * opts$H * opts$nlam)
  if(is.null(opts$itheta)) opts$itheta <- integer(opts$pmax)
  if(is.null(opts$ntheta)) opts$ntheta <- integer(opts$nlam)
  if(is.null(opts$alam)) opts$alam <- double(opts$nlam)
  if(is.null(opts$npass)) opts$npass <- integer(1)
  if(is.null(opts$jerr)) opts$jerr <- integer(1)
  
  M0 <- M
  U0 <- U
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2))
  n3 <- length(gam)
  nparams <- n1*n2*n3
  
  # The following lists save the corresponding objects for each tuning parameter, 
  B_l <- vector("list", nparams) # estimated B matrix
  beta_l <- vector("list", nparams) # estimated basis matrix
  step_l <- rep(NA_integer_, nparams) # iterations
  time_l <- rep(NA_real_, nparams) # running time
  rank_l <- rep(NA_integer_, nparams) # estimated rank
  s_l <- rep(NA_integer_, nparams) # estimated sparsity level
  
  # Count the number of converged matrices.
  nlam_cvg <- 0
  
  for(i in 1:n1){
    lambda1 <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){
        lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j,k])
        
        M <- M0 + gamma*diag(rep(1,ncol(M0)), ncol(M0),ncol(M0))
        
        # Initialize three matrices
        Bold <- matrix(0,dim(U0)[1], dim(U0)[2])
        Cold <- matrix(0,dim(U0)[1], dim(U0)[2])
        etaold <- matrix(0,dim(U0)[1], dim(U0)[2])
        
        # The MAIN loop of admm method
        step <- 0    
        start_time <- Sys.time()
        
        repeat{
          step <- step + 1
          
          # Update B
          U <- U0 - etaold + gamma * Cold
          out_B <- updateB(M, U, lambda1, opts)
          Bnew <- out_B$Bnew
          jerr <- out_B$jerr
          #if(jerr != 0) break
          
          # Update C
          Cnew <- updateC(Bnew, lambda2, gamma, etaold)
          
          # Update eta (omega in SEAS algorithm)
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          # Code 1: success
          if(max(abs(Bnew - Cnew)) < eps){
            jerr <- 1
            break
          }
          # Code 404: then maximal iteration is reached
          if(step > maxit){
            jerr <- 404
            warning('Maximal iteration is reached.')
            break
          }
          Bold <- Bnew
          Cold <- Cnew
          etaold <- etanew
        }# End of repeat 
        end_time <- Sys.time()  # The time for each repeat
        time <- difftime(end_time, start_time, units = "secs")
        # Code < -10000: non-sparse matrix
        #if(jerr < -10000){
        #  break
        #}
        jerr=1
        # Code 1: success, save the matrix and the related information.
        if(jerr==1){
          index <- (i-1)*n2*n3 + (j-1)*n2 + k
          nlam_cvg <- nlam_cvg + 1
          B_l[[index]] <- Bnew
          step_l[index] <- step
          time_l[index] <- time
          if(is.null(d)) rank <- rank_func(Cnew, thrd = eps)
          else rank <- d
          rank_l[index] <- rank
          # Cut and select the left singular vector of Bnew
          if(rank == 0){
            beta <- matrix(0, nrow(Bnew), ncol(Bnew))
          }else{
            tmp <- svd(Bnew)$u[,1:rank, drop = FALSE]
            vec <- as.vector(tmp)
            vec[abs(vec) < eps] <- 0
            beta <- matrix(vec, nrow(tmp), ncol(tmp))
          }
          beta_l[[index]] <- beta
          var_ind <- apply(beta, 1, function(x){any(x!=0)})
          s_l[index] <- sum(var_ind)
        }
      }# End of lambda2
      if(jerr < -10000) break
    }# End of gam
  }# End of lambda1
  return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg))
}

# Update B matrix in ADMM algorithm. Use the group-wise coordinate descent algorithm from 'msda' R package.
updateB <- function(M, U, lambda1, opts){
  U <- t(U)
  fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
  #if(fit$jerr != 0){return(list(Bnew = NULL, jerr = fit$jerr))} # Code: non-zero, abnormal result.
  outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
  Bnew <- as.matrix(outlist$theta[[1]])
  list(Bnew = Bnew, jerr = fit$jerr)
}

# Update C matrix in ADMM algorithm.
updateC <- function(Bnew, lambda2, gamma, etaold){
  Btemp <- Bnew + 1/gamma * etaold
  svd_B <- svd(Btemp)
  lamtemp <- pmax(0, svd_B$d-lambda2/gamma)
  Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
  Cnew
}

# ------------------ revised functions from 'msda' package ---------------------- #
# Revise 'msda' function to accommodate other forms of M and U matrices.
# Some inputs are similar to the ones in 'seas' function. Please refer to 'msda' package documentation for more details of the arguments in 'msda' function.
# 
# Outputs:
# ========
# lambda: The tuning parameter sequence.
# theta: The list of estimated matrix.
# M: The M matrix from samples. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix from samples, which depends on the argument type.
# rank: The list of estimated rank for each matrix.
msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', FUN = NULL, lambda.factor=NULL, nlambda=100, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, verbose=FALSE, perturb=NULL){
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, type, FUN)
    M <- MU_out$M
    U <- MU_out$U
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax*2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) 
    diag(M) <- diag(M) + perturb
  H <- as.integer(dim(U)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam),alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, H)
  rank <- rep(NA_integer_, length(outlist$theta))
  for (i in 1:length(outlist$theta)){
    if(!is.null(outlist$theta[[i]])){
      rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
    }
  }
  if(is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  outlist <- list(lambda = outlist$lambda, theta = outlist$theta, M = M, U = U, rank = rank)
  class(outlist) <- c("msda")
  outlist
}

# Revise 'cv.msda' function to accommodate other forms of M and U matrices. We also add the optional argument 'fold' to pass the user-specified folds index.
# Some inputs are similar to the ones in 'cv.seas' function. Please refer to 'msda' package documentation for more details of arguments used in 'cv.msda' function.
# 
# Outputs:
# ========
# beta: The optimal estimated matrix.
# id: The index of the optimal tuning parameter.
# lambda: The lambda sequence.
# lam_max: The optimal tuning parameter.
# rank: The rank of the optimal estimated matrix.
# M: The M matrix based on the full data. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix based on the full data, which depends on the argument type.
# M_fold: The M matrix list based on each cross-validation data fold.
# U_fold: The U matrix list based on each cross-validation data fold.
cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', lambda.factor=NULL, nlambda=100, nfolds=5, foldid = NULL, lambda = NULL, FUN = NULL, maxit = 1e3, plot = FALSE){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  nobs <- nrow(x)
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  # Fit the model on the full data, obtain the lambda sequence.
  fit <- msda(x, y, yclass = yclass, type = type, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  
  # Cross-validation
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(foldid))
  }
  
  cv_out <- lapply(1:nfolds, function(k){
    x_train <- x[foldid!=k,,drop=FALSE]
    x_val <- x[foldid==k,,drop=FALSE]
    y_train <- y[foldid!=k]
    y_val <- y[foldid==k]
    yclass_train <- yclass[foldid!=k]
    
    fit_fold <- msda(x_train, y_train, yclass_train, type = type, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc(beta_fold, x_val, y_val)
    if(length(eval_fold) != length(lambda)){
      eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    }
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2)
  U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)){
    eval_all <- t(as.matrix(eval_all))
  }
  
  ## No matrix is converged in any fold
  if(all(is.na(eval_all))) return(NULL)
  
  print(eval_all)
  print(typeof(eval_all))
  cvm <- colMeans(eval_all, na.rm = TRUE)
  # The optimal lambda1
  id_max <- which.max(cvm)
  lam_max <- lambda[id_max]
  beta <- as.matrix(beta_l[[id_max]])
  
  # Recalculate the rank
  rank <- rank_func(beta, thrd = 1e-3)
  
  if(plot){ # If TRUE, plot the cv evaluation for each tuning parameter.
    dat <- data.frame(x = 1:length(cvm), y = cvm)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 1)+
      xlab("")+
      ylab("Distance correlation")+
      theme_bw()
    #print(g)
  }
  
  list(beta = beta, id = id_max, lambda = lambda, lam_max = lam_max, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}

Sparse_vector<-function(Y,design_matrix,nfolds=5,lambda=NULL)
{
  if(is.null(lambda))
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,standardize=FALSE,type.measure="mse",nlambda=100,intercept=FALSE)
  }else
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,lambda=lambda,standardize=FALSE,type.measure="mse",intercept=FALSE)
  }
  return(as.numeric(coef(fit))[-1])
}

choose_weight<-function(x,rho=0.25)                
{
  x<-as.numeric(x)
  return(sum(x^2)^(-rho))
}

forward_column_selection <- function(X, threshold=1e-12, D=15) {
  n <- ncol(X)
  selected_cols <- c()
  col_norms <- apply(X, 2, function(col) sqrt(sum(col^2)))
  for (i in 1:D) {
    max_norm_index <- which.max(col_norms)
    max_norm_value <- col_norms[max_norm_index]
    if (max_norm_value < threshold) {
      break
    }
    selected_cols <- c(selected_cols, max_norm_index)
    selected_col <- X[, max_norm_index]
    for (j in 1:n) {
      if (!j %in% selected_cols) {
        projection <- sum(X[, j] * selected_col) / sum(selected_col^2) * selected_col
        X[, j] <- X[, j] - projection
        col_norms[j] <- sqrt(sum(X[, j]^2))
      }
    }
    col_norms[max_norm_index] <- 0
  }
  
  return(selected_cols)
}

U_TM<-function(x,index)
{
  X <- as.matrix(x)
  n <- nrow(X)
  p <- ncol(X)
  
  idx <- as.integer(index)
  
  Z2 <- X[, idx, drop = FALSE]^2
  crossprod(X, Z2) / n  
}

U_sin<-function(x, index) {
  X <- as.matrix(x)
  n <- nrow(X)
  p <- ncol(X)
  
  idx <- as.integer(index)

  S <- crossprod(scale(X, center = TRUE, scale = FALSE)) / n  # p x p
  Z <- X[, idx, drop = FALSE]
  SinZ <- sin(Z)
  CosZ <- cos(Z)

  term1 <- crossprod(X, SinZ) / n

  meanCos <- colMeans(CosZ) 
  term2 <- sweep(S[, idx, drop = FALSE], 2L, meanCos, `*`)
  
  term1 - term2
}


Model1<-function(p,n,cov_index)
{
  p=p
  n=n
  n1=as.integer(0.65*n)
  n2=as.integer(0.35*n)
  cov<-diag(p)
  mu1<-rep(0,p)
  mu1[1:4]<-c(1,1,-1,-1)
  mu1<-cov%*%mu1
  mu1<-mu1/norm(mu1,type="2")*1.5
  mu2<--mu1
  X1<-mvrnorm(n1,mu1,cov)
  X2<-mvrnorm(n2,mu2,cov)
  X<-rbind(X1,X2)
  X<-scale(X,scale=FALSE)
  return(X)
}


third_kernel<-function(x)
{
  x*x   
}


sin_kernel<-function(x)
{
  X <- as.matrix(x)
  p <- dim(X)[2]
  idx <- 1:p
  Z <- X[, idx, drop = FALSE]          # n x m
  meanCos <- colMeans(cos(Z))          # m
  V <- sin(Z) - sweep(Z, 2, meanCos, `*`)  # n x m
  V
}


WSS<-function(x)
{
  n<-nrow(as.matrix(x))
  x<-scale(x,scale=FALSE)
  return(log(sum(x^2)/n))
}

eval_lm <- function(Beta, M, U,n){
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      nnz_rows <- sum(rowSums(abs(mat)) > 0)
      reg<-0.5*sum((M %*% mat) * mat)-sum(mat * U)
      n*log(reg/n)+nnz_rows*log(n)
    }
  })
  return(result)
}

ADMM_estimation<-function(x, index, lam1, gam, method="TM",nfolds=5)
{
  p=ncol(as.matrix(x))
  n=nrow(as.matrix(x))
  lam2=0
  lam1=lam1
  gam=gam
  M=cov(x)
  if(method=="TM")
  {
    U=U_TM(x,index)
  }
  if(method=="SIN")
  {
    U=U_sin(x,index)
  }
  foldid <- c()
  count<-n
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))}
  out_all <- lapply(1:nfolds, function(k)
  { # Cross-validation
    x_val <- x[foldid==k,,drop=FALSE]
    x_train <- x[foldid!=k,,drop=FALSE]
    # Fit with seas function
    M_train=cov(x_train)
    M_val=cov(x_val)
    if(method=="TM")
    {
      U_train=U_TM(x_train,index)
      U_val=U_TM(x_val,index)
    }
    if(method=="sin")
    {
      U_train=U_sin(x_train,index)
      U_val=U_sin(x_val,index)
    }
    fit_fold <- admm(M_train, U_train, p,p, lam1=lam1, lam2=lam2, gam=gam, eps=1e-4, maxit=1e+3)
    B_l<-fit_fold$B
    eval_fold <- eval_lm(B_l, M_val, U_val,n)
    return(eval_fold)
  })
  eval_all <- do.call(rbind, lapply(out_all, deal_list))
  cvm <- colMeans(eval_all, na.rm=TRUE)
  lam1_max=lam1[which.min(cvm)]
  fit<-admm(M, U, p,p, lam1_max, 0, gam, eps=1e-3, maxit=1e+3)
  if(length(index)==1)
  {
    return((fit$B[[1]])[,1])
  }else
  {
    return(fit$B[[1]])
  }
}

ADMM_estimation_weight<-function(x,weight, index, lam1, gam, method="TM",nfolds=5)
{
  p=ncol(as.matrix(x))
  n=nrow(as.matrix(x))
  lam2=0
  lam1=lam1
  gam=gam
  M=cov(x)*outer(weight,weight)
  if(method=="TM")
  {
    U=weight*U_TM(x,index)
  }
  if(method=="SIN")
  {
    U=weight*U_sin(x,index)
  }
  foldid <- c()
  count<-n
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))}
  out_all <- lapply(1:nfolds, function(k)
  { # Cross-validation
    x_val <- x[foldid==k,,drop=FALSE]
    x_train <- x[foldid!=k,,drop=FALSE]
    # Fit with seas function
    M_train=cov(x_train)*outer(weight,weight)
    M_val=cov(x_val)*outer(weight,weight)
    if(method=="TM")
    {
      U_train=weight*U_TM(x_train,index)
      U_val=weight*U_TM(x_val,index)
    }
    if(method=="sin")
    {
      U_train=weight*U_sin(x_train,index)
      U_val=weight*U_sin(x_val,index)
    }
    fit_fold <- admm(M_train, U_train, p,p, lam1=lam1, lam2=lam2, gam=gam, eps=1e-4, maxit=1e+3)
    B_l<-fit_fold$B
    eval_fold <- eval_lm(B_l, M_val, U_val,n)
    return(eval_fold)
  })
  eval_all <- do.call(rbind, lapply(out_all, deal_list))
  cvm <- colMeans(eval_all, na.rm=TRUE)
  lam1_max=lam1[which.min(cvm)]
  fit<-admm(M, U, p,p, lam1_max, 0, gam, eps=1e-3, maxit=1e+3)
  if(length(index)==1)
  {
    return(weight*(fit$B[[1]])[,1])
  }else
  {
    return(weight*fit$B[[1]])
  }
}

Clustering<-function(X,nfolds=5, D=10, lambda=NULL,method="TM")
{
  p<-ncol(as.matrix(X))
  n<-nrow(as.matrix(X))
  if(method=="TM")
  {
  kernel_matrix<-third_kernel(X)
  }
  if(method=="SIN")
  {
  kernel_matrix<-sin_kernel(X)
  }
  sfInit(parallel = TRUE,cpus=15)
  sfLibrary(glmnet)
  sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=seq(0.025,0.2,0.025))
  #sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=NULL)
  sfStop()
  num=p
  index<-c()
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
      index=append(index,i)
    }
  }
  index<-c(1:num)[-index]
  for(i in index)
  {
    M=cov(X)
    if(method=="TM")
    {
      U=U_TM(X,i)
    }
    if(method=="SIN")
    {
      U=U_TM(X,i)
    }
    sparse_estimate[,i]=admm(M,U,p,p,0.3,0,5)$B[[1]]
  }
  index<-c()
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
      index=append(index,i)
    }
  }
  index<-c(1:num)[-index]
  if(length(index)==1)
  {
    sparse_index<-sample(1:p,1)
  }else
  {
    sparse_index<-forward_column_selection(sparse_estimate,D=10)
  }
  print(sparse_index)
  
  result=ADMM_estimation(X,sparse_index,c(0.04,0.08,0.04/49),5,method)
  weight=1/apply(result,MARGIN=1, choose_weight)
  result=ADMM_estimation_weight(X,weight,sparse_index,c(0.06,0.1,0.04/49),5,method)
  return(final_fit)
  #solution_equal<-Matrix_estimation(data=X,ytilde=kernel)
  #temp=do.call(cbind,(lapply(solution_equal,deal_list)))
  #weight<-apply(temp,MARGIN=1,choose_weight)
}


cal_acc <- function(true_labels, predicted_labels) {
  true_labels <- as.character(true_labels)
  predicted_labels <- as.character(predicted_labels)
  
  unique_true <- unique(true_labels)
  unique_pred <- unique(predicted_labels)
  if (length(unique_pred) > length(unique_true)) {
    unique_true_padded <- c(unique_true, rep("dummy_unmatched", length(unique_pred) - length(unique_true)))
    permutations <- permn(unique_true_padded)
  } else {
    permutations <- permn(unique_true)
  }
  
  accuracies <- sapply(permutations, function(perm) {
    perm_used <- perm[1:length(unique_pred)]
    mapping <- setNames(perm_used, unique_pred)
    mapped_pred <- mapping[predicted_labels]
    sum(mapped_pred == true_labels, na.rm = TRUE) / length(true_labels)
  })
  return(max(accuracies))
}


sample_tm <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  mean_X <- colMeans(X)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  result <- matrix(0, nrow = p, ncol = p^2)
  for (i in 1:n) {
    outer_prod <- X_centered[i,] %*% t(X_centered[i,])
    kronecker_prod <- as.vector(outer_prod)
    result <- result + X_centered[i,] %*% t(kronecker_prod)
  }
  result <- result / n
  return(result)
}

cov_decay<-function(p,rho)
{
  cov1<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      cov1[i,j]=rho^(abs(i-j))
    }
  }
  return(cov1)
}

cov_dense_biway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if((i<=a)||(j<=a)||(i>=(p-a+1))||(j>=(p-a+1)))
      {
        cov[i,j]=rho^abs(i-j)
      }else
        if(i==j)
        {
          cov[i,j]=1
        }
      else
      {
        cov[i,j]=rho
      }
    }
  }
  return(cov)
}

cov_dense_oneway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if((i<=a)||(j<=a))
      {
        cov[i,j]=rho^abs(i-j)
      }else
        if(i==j)
        {
          cov[i,j]=1
        }
      else
      {
        cov[i,j]=rho
      }
    }
  }
  return(cov)
}

Model1<-function(p,n,cov_index)
{
  p=p
  n=n
  n1=as.integer(0.65*n)
  n2=as.integer(0.35*n)
  cov<-diag(p)
  mu1<-rep(0,p)
  mu1[1:4]<-1.5*c(1,1,-1,-1)
  mu1<-cov%*%mu1
  mu1<-mu1/norm(mu1,type="2")*1.5
  mu2<--mu1
  X1<-mvrnorm(n1,mu1,cov)
  X2<-mvrnorm(n2,mu2,cov)
  X<-rbind(X1,X2)
  X<-scale(X,scale=FALSE)
  Y<-c(rep(1,as.integer(n*0.65)),rep(2,as.integer(n*0.35)))
  return(list(X=X,Y=Y))
}


Model2<-function(p,n,cov_index)
{
  r=p
  n=n
  n1=as.integer(0.55*n)
  n2=as.integer(0.25*n)
  n3=as.integer(0.2*n)
  sigma0=cov_decay(p,0.5)
  mu1[(as.integer(r/2-2)):(as.integer(r/2+3))] <- c(1,-1,1,-1,1,-1)
  mu1 <- sigma0 %*% mu1; mu1 <- mu1 / norm(mu1, type="2") * 2.5
  mu2 <- rep(0, r); mu3 <- -mu1
  gamma_v <- mu1 / 2.5
  cov1 <- sigma0 - exp(-1) * gamma_v %*% t(gamma_v)
  cov2 <- sigma0 - exp(-2) * gamma_v %*% t(gamma_v)
  cov3 <- sigma0 - exp(-3) * gamma_v %*% t(gamma_v)
  X1<-mvrnorm(n1,mu1,cov1)
  X2<-mvrnorm(n2,mu2,cov2)
  X3<-mvrnorm(n3,mu3,cov3)
  X<-rbind(X1,X2)
  X<-rbind(X,X3)
  X<-scale(X,scale=FALSE)
  Y<-c(rep(1,as.integer(n*0.55)),rep(2,as.integer(n*0.25)),rep(3,as.integer(n*0.2)))
  return(list(X=X,Y=Y))
}

Model3<-function(p,n,cov_index)
{
  r=p
  n=n
  n1=as.integer(0.5*n)
  n2=as.integer(0.3*n)
  n3=as.integer(0.2*n)
  sigma0 <- cov_decay(r, 0.5)
  mu1 <- rep(0, r)
  mu1[(as.integer(r/2-2)):(as.integer(r/2+3))] <- c(1,-1,1,-1,1,-1)
  mu1 <- sigma0 %*% mu1; mu1 <- mu1 / norm(mu1, type="2") * 3
  mu2 <- rep(0, r); mu3 <- -mu1
  gamma_v <- rep(0, r); gamma_v[1:4] <- c(1, 1, 1, 1)
  gamma_v <- sigma0 %*% gamma_v; gamma_v <- gamma_v / norm(gamma_v, type="2") * 3
  cov1 <- sigma0 + exp(-1) * gamma_v %*% t(gamma_v)
  cov2 <- sigma0 + exp(-2) * gamma_v %*% t(gamma_v)
  cov3 <- sigma0 + exp(-3) * gamma_v %*% t(gamma_v)
  X1<-mvrnorm(n1,mu1,cov1)
  X2<-mvrnorm(n2,mu2,cov2)
  X3<-mvrnorm(n3,mu3,cov3)
  X<-rbind(X1,X2)
  X<-rbind(X,X3)
  X<-scale(X,scale=FALSE)
  Y<-c(rep(1,as.integer(n*0.5)),rep(2,as.integer(n*0.3)),rep(3,as.integer(n*0.2)))
  return(list(X=X,Y=Y))
}



Model_X<-function(p,n,cov_index,model_index)
{
  if(model_index==1)
  {
    return(Model1(p,n,cov_index))
  }
  if(model_index==2)
  {
    return(Model2(p,n,cov_index))
  }
  if(model_index==3)
  {
    return(Model3(p,n,cov_index))
  }
}

beta_X<-function(model_index,p)
{
  if(model_index==1)
  {
    beta<-rep(0,p)
    beta[1:5]<-c(1,1,1,1,1)
    return(beta)
  }
  if(model_index==2)
  {
    beta1<-rep(0,p)
    beta1[1:5]<-c(1,1,1,1,1)/sqrt(5)
    beta2<-rep(0,p)
    beta2[1:5]<-c(1,1,-1,-1,-1)/sqrt(5)
    return(cbind(beta1,beta2))
  }
  if(model_index==3)
  {
    beta<-rep(0,p)
    beta[1:5]<-c(1,1,1,1,1)
    return(beta)
  }
}


loss_SDR<-rep(0,100) #the loss of estimating central subspace
mis_rate_SDR<-rep(0,100)
TIME_SDR<-0
mis_rate_sparsek<-rep(0,100)
TIME_sparsek<-0


assign_sparse_kmeans_vec <- function(X_test, centers, ws) {
  # 广播计算所有样本到所有中心的加权距离
  K <- nrow(centers)
  dist_mat <- sapply(1:K, function(k) {
    diff <- sweep(X_test, 2, centers[k, ])  # n_test x p
    rowSums(sweep(diff^2, 2, ws, `*`))      # n_test x 1
  })
  return(apply(dist_mat, 1, which.min))
}

assign_sparse_kmeans <- function(X_test, centers, ws) {
  K <- nrow(centers)
  dist_mat <- sapply(1:K, function(k) {
    diff <- sweep(X_test, 2, centers[k, ])
    rowSums(sweep(diff^2, 2, ws, `*`))
  })
  return(apply(dist_mat, 1, which.min))
}


#################k-means##################
p=500
n=500
model_index=1
cov_index=1
K=2
for(k in 1:200)
{
  data<-Model_X(p,n,cov_index,model_index)
  test_data<-Model_X(p,3000,cov_index,model_index)
  test_X<-test_data$X
  test_Y<-test_data$Y
  X=data$X
  real_label=data$Y
  X=scale(X,scale=FALSE)
  test_X=scale(test_X,scale=FALSE)
  time=as.numeric(Sys.time())
  km.perm <- KMeansSparseCluster.permute(X,K=K,wbounds=seq(1.5,7,len=50),nperms=5)
  km.out <- KMeansSparseCluster(X,K=2,wbounds=km.perm$bestw)
  time=as.numeric(Sys.time())-time
  ws <- km.out[[1]]$ws     
  Cs <- km.out[[1]]$Cs        
  centers <- t(sapply(1:K, function(k) {
    colMeans(X[Cs == k, , drop = FALSE])
  })) 
  test_pred <- assign_sparse_kmeans(test_X, centers, ws)
  mis_rate_sparsek[k]=1-cal_acc(test_pred,  test_Y)
  print(mis_rate_sparsek[k])
  TIME_sparsek=TIME_sparsek+time
}


TDRR <- function(eigenvalues, n, tau = 0.5, c_n = NULL, d_n = NULL) {
  p <- length(eigenvalues)
  if (is.null(c_n)) c_n <- 0.75 * log(n) / sqrt(n)
  if (is.null(d_n)) d_n <- 1.5 * log(n) / sqrt(n)
  lambda <- c(eigenvalues, 0)
  s_star <- (lambda[1:p] + c_n) / (lambda[2:(p+1)] + c_n)
  R <- (s_star[2:p] + d_n) / (s_star[1:(p-1)] + d_n)
  valid_indices <- which(R < tau)
  if (length(valid_indices) > 0) {
    q_hat <- max(valid_indices)
  } else {
    q_hat <- 0
  }
  return(list(q_hat = q_hat, R_ratios = R, s_star = s_star))
}


TM_SDR=0
##############DR-GMM##########################
method="TM"
method="SIN"
p=500
n=500
model_index=1
cov_index=1
for(k in 1:200)
{
  data<-Model_X(p,n,cov_index,model_index)
  test_data<-Model_X(p,3000,cov_index,model_index)
  X=data$X
  test_X<-test_data$X
  test_Y<-test_data$Y
  real_label=data$Y
  X=scale(X,scale=FALSE)
  test_X=scale(test_X,scale=FALSE)
  time=as.numeric(Sys.time())
  result=Clustering(X)
  index<-which(abs(as.matrix(result)[,1])>1e-6)
  result_new<-reulst[index,]
  d=TDRR(svd(result_new)$d,n,0.75*log(n/length(index))/sqrt(n/length(index)),1.5*log(n/length(index))/sqrt(n/length(index)))
  obj=svd(result)$u[,1:d]
  cls_result <- run_em(X%*%obj, G=K)
  time=as.numeric(Sys.time())-time
  pred <- predict(cls_result, newdata=test_X%*%obj)
  cls_label <- pred$classification
  mis_rate_SDR[k]=1-cal_acc(pred,  test_Y)
  TIME_sparsek=TIME_sparsek+time
}

