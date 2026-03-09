
library(Matrix)
library(MASS)
library(glmnet)
library(e1071)
library(lattice)
library(stringr)
library(ggplot2)
library(energy)
library(mclust)
library(combinat)
library(mvtnorm)
library(LiblineaR)
library(EMMIXmfa)     
source("clemm_em.R") 

######################## Utility Functions  ###################################

run_em <- function(data, G=3, replicate=5) {
  init_final <- list()
  lik <- -10000000
  for(k in 1:replicate) {
    init <- kmeans(data, centers=G)
    z <- matrix(0, nrow=nrow(data), ncol=G)
    z[cbind(1:nrow(data), init$cluster)] <- 1
    pi <- colMeans(z)
    mu <- init$centers
    sigma <- vector("list", G)
    for(kk in 1:G) {
      cluster_data <- data[init$cluster == kk, , drop=FALSE]
      sigma[[kk]] <- cov(as.matrix(cluster_data))
    }
    density <- matrix(0, nrow=nrow(data), ncol=G)
    for(kk in 1:G) {
      density[, kk] <- dmvnorm(data, mean=mu[kk, ], sigma=sigma[[kk]]) * pi[kk]
    }
    loglik <- sum(log(rowSums(density)))
    if(loglik > lik) {
      init_final <- init
      lik <- loglik
    }
  }
  result <- Mclust(data, G=G, initialization=list(classification=init_final$cluster), tol=1e-5)
  return(result)
}

run_init <- function(data, G=3, replicate=5) {
  init_final <- list()
  lik <- -10000000
  for(k in 1:replicate) {
    init <- kmeans(data, centers=G)
    z <- matrix(0, nrow=nrow(data), ncol=G)
    z[cbind(1:nrow(data), init$cluster)] <- 1
    pi <- colMeans(z)
    mu <- init$centers
    sigma <- vector("list", G)
    for(kk in 1:G) {
      cluster_data <- data[init$cluster == kk, , drop=FALSE]
      sigma[[kk]] <- cov(as.matrix(cluster_data))
    }
    density <- matrix(0, nrow=nrow(data), ncol=G)
    for(kk in 1:G) {
      density[, kk] <- dmvnorm(data, mean=mu[kk, ], sigma=sigma[[kk]]) * pi[kk]
    }
    loglik <- sum(log(rowSums(density)))
    if(loglik > lik) {
      init_final <- init
      lik <- loglik
    }
  }
  return(init_final)
}

cov_decay <- function(p, rho) {
  outer(1:p, 1:p, function(i, j) rho^abs(i - j))
}

matpower <- function(M, power) {
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(pmax(e$values, 1e-15)^power) %*% t(e$vectors)
}

proj <- function(B) {
  B <- as.matrix(B)
  if (ncol(B) == 0 || nrow(B) == 0) return(matrix(0, nrow(B), nrow(B)))
  B %*% solve(crossprod(B)) %*% t(B)
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

AR <- function(rho, p) {
  m <- matrix(0, p, p)
  for (i in 1:p) for (j in 1:p) m[i, j] <- rho^(abs(i - j))
  return(m)
}

cut_mat <- function(Beta, thrd, rank) {
  l <- length(Beta)
  for (i in 1:l) {
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat); nvars <- ncol(mat); r <- rank[i]
    if(r == 0) { Beta[[i]] <- matrix(0, nobs, nvars) }
    else { vec <- as.vector(mat); vec[abs(vec) < thrd] <- 0; Beta[[i]] <- matrix(vec, nobs, nvars) }
  }
  return(Beta)
}

eval_dc <- function(Beta, x, y) {
  if(!is.list(Beta)) Beta <- list(Beta)
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i) {
    if(is.null(Beta[[i]])) NA
    else { mat <- as.matrix(Beta[[i]]); dcor(x %*% mat, y) }
  })
  return(result)
}

MU <- function(x, y, yclass=NULL, type='sir', FUN = NULL, categorical = FALSE, H = 5) {
  if(is.null(yclass)) {
    if(categorical == FALSE) {
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    } else { yclass <- y }
  }
  cls <- sort(unique(yclass)); nclass <- length(cls)
  nobs <- as.integer(dim(x)[1]); nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i) mean(yclass == i))
  mu <- colMeans(x)
  x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE)
  M <- crossprod(x_c / sqrt(nobs))
  if(type == 'sir') {
    U <- matrix(0, nvars, nclass)
    for (i in 1:nclass) U[, i] <- colMeans(x_c[yclass == cls[i], , drop=FALSE])
  } else if(type == 'intra') {
    y_c <- y - mean(y); U <- matrix(0, nvars, nclass)
    lb <- quantile(y_c, 0.1)[[1]]; ub <- quantile(y_c, 0.9)[[1]]
    y_c <- sapply(y_c, cut_func, lb = lb, ub = ub)
    for (i in 1:nclass) {
      y_copy <- y_c; y_copy[yclass != cls[i]] <- 0
      U[, i] <- (1/nobs) * t(x_c) %*% (y_copy - mean(y_copy))
    }
  } else if(type == 'pfc') {
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) else Fmat <- t(sapply(y, FUN))
    Fmat_mean <- colMeans(Fmat)
    Fmat_c <- Fmat - matrix(Fmat_mean, NROW(Fmat), NCOL(Fmat), byrow = TRUE)
    lb <- apply(Fmat_c, 2, quantile, 0.1); ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)) Fmat_c[, i] <- sapply(Fmat_c[, i], cut_func, lb[i], ub[i])
    U <- (1/nobs) * (t(x_c) %*% Fmat_c)
  }
  list(M = M, U = U, nclass = nclass, prior=prior)
}

cut_func <- function(x, lb, ub) {
  if(x < lb) return(lb) else if(x > ub) return(ub) else return(x)
}

rank_func <- function(B, thrd) {
  d <- svd(B)$d; r <- sum(d >= thrd); return(r)
}

subspace <- function(A, B) {
  if(is.vector(A)) A <- as.matrix(A)
  if(is.vector(B)) B <- as.matrix(B)
  Pa <- qr.Q(qr(A)); Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B)); Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa - Pb, type="F") / sqrt(2 * d))
}

formatoutput <- function(fit, maxit, pmax, p, H) {
  nalam <- fit$nalam; ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta); lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0) {
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i) {
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H); theta_i[ja, ] <- a; theta_i
    })
  } else {
    theta <- lapply(seq(nalam), function(x) matrix(0, p, H))
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) msg <- ""
  if (n > 0) {
    if (n < 7777) msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) msg <- "All penalty factors are <= 0"
    n <- 1; msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    if (n > -10000) msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, " iterations; solutions for larger lambdas returned.\n", sep = "")
    if (n < -10000) msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    if (n < -20000) msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam) {
  llam <- log(lam)
  if(length(llam) >= 3) lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

seas <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', M = NULL, U = NULL, nobs = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), FUN = NULL, eps = 1e-3, maxit = 1e+3, ...) {
  if(is.null(M) || is.null(U)) {
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)) {
      if(categorical == FALSE) {
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
        nclass <- as.integer(length(unique(yclass)))
      } else { yclass <- y }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    if(is.null(gamma)) gamma <- c(10,30,50)
    if(is.null(lam1) || is.null(lam2)) {
      fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=10, lambda.factor=0.5, nfolds = 5, FUN = FUN, maxit=1e3)
      M <- fit_1$M; U <- fit_1$U
      id_max_msda <- fit_1$id; lam1_max_msda <- fit_1$lam_max
      beta_msda <- as.matrix(fit_1$beta)
      if(is.null(lam1)) lam1 <- (lam1_max_msda) * lam1_fac
      if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
      if (all(lam2 == 0)) { lam2 <- 0; warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.") }
    } else {
      MU_out <- MU(x, y, yclass, type, FUN); M <- MU_out$M; U <- MU_out$U
    }
    nobs <- as.integer(dim(x)[1]); nvars <- as.integer(dim(x)[2])
  } else {
    if(is.null(lam1) || is.null(lam2) || is.null(gamma)) stop("Sequences lam1, lam2 or gamma is missing.")
    if(is.null(nobs)) stop("Missing nobs.")
    nvars <- NCOL(M)
  }
  code <- 0
  if(is.vector(lam1) && (length(lam1) == 1) && (lam1 == 0) && is.vector(lam2) && (length(lam2) == 1) && (lam2 == 0)) {
    B <- solve(M) %*% U
    if(is.null(d)) beta <- svd(B)$u
    else if(d == 0) beta <- matrix(0, nrow(B), ncol(B))
    else beta <- svd(B)$u[, 1:d, drop=FALSE]
    vec <- as.vector(beta); vec[abs(vec) < 1e-3] <- 0
    beta <- matrix(vec, nrow(beta), ncol(beta))
    rank <- NCOL(beta)
    output <- list(beta = beta, B = B, rank = rank, lam1 = lam1, lam2 = lam2, code = code)
  } else {
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, ...)
    B_l <- fit$B; beta_l <- fit$beta
    if (all(sapply(beta_l, is.null))) {
      code <- 1; warning("No converged results returned.")
      return(list(beta = beta_l, code = code))
    }
    rank_l <- fit$rank; s_l <- fit$s; step_l <- fit$step; time_l <- fit$time
    if(length(B_l) == 1) {
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}

cv.seas <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', lambda.factor=0.5, nlambda=10, nfolds = 5, foldid = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), plot = FALSE, FUN = NULL, eps = 1e-3, maxit = 1e+3, trace.it = FALSE, ...) {
  start_time <- Sys.time()
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)) {
    if(categorical == FALSE) {
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    } else { yclass <- y }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  M <- U <- M_fold <- U_fold <- NULL; nobs <- dim(x)[1]
  if(is.null(foldid)) {
    ord <- order(y); y <- y[ord]; yclass <- yclass[ord]; x <- x[ord, ]
    if (nfolds < 3) stop("nfolds must be larger than 3")
    if (nfolds > nobs) stop("nfolds is larger than the sample size")
    count <- as.numeric(table(yclass)); foldid <- c()
    for(cnt in count) foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  } else { nfolds <- length(unique(foldid)) }
  if(is.null(gamma)) gamma <- c(10,30,50)
  if(is.null(lam1) || is.null(lam2)) {
    fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=nlambda, lambda.factor=lambda.factor, foldid = foldid, FUN = FUN, maxit=1e3, plot = plot)
    M <- fit_1$M; U <- fit_1$U; M_fold <- fit_1$M_fold; U_fold <- fit_1$U_fold
    id_max_msda <- fit_1$id; lam1_max_msda <- fit_1$lam_max; beta_msda <- as.matrix(fit_1$beta)
    if(is.null(lam1)) lam1 <- (lam1_max_msda) * lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)) { lam2 <- 0; warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.") }
  }
  n1 <- length(lam1); n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2]); n3 <- length(gamma)
  nerr <- 0; code <- 0
  end_time <- Sys.time(); time1 <- difftime(end_time, start_time, units = "secs")
  out_all <- lapply(1:nfolds, function(k) {
    if(trace.it) cat(sprintf("Fold: %d/%d\n", k, nfolds))
    x_val <- x[foldid==k, , drop=FALSE]; y_val <- y[foldid==k]
    if(is.null(M_fold) || is.null(U_fold)) {
      x_train <- x[foldid!=k, , drop=FALSE]; y_train <- y[foldid!=k]; yclass_train <- yclass[foldid!=k]
      fit_fold <- seas(x_train, y_train, yclass = yclass_train, type = type, FUN = FUN, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    } else {
      fit_fold <- seas(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(foldid!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    }
    err <- 0; beta_l <- fit_fold$beta; rank_l <- fit_fold$rank; step_l <- fit_fold$step; time_l <- fit_fold$time
    eval_fold <- eval_dc(beta_l, x_val, y_val)
    ind <- which(sapply(beta_l, is.null)); rank_l[ind] <- -1; eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
    list(eval_fold, err)
  })
  eval_all <- do.call(rbind, lapply(out_all, "[[", 1))
  errs <- do.call(c, lapply(out_all, "[[", 2)); nerr <- sum(errs)
  if((nerr != 0) && (nerr != nfolds)) { code <- 3; warning(paste0("No converged results returned in", nerr, "folds.")) }
  else if(nerr == nfolds) { code <- 4; warning("No converged results returned in any fold."); return(list(beta = NULL, code = code)) }
  if(is.vector(eval_all)) eval_all <- as.matrix(eval_all)
  cvm <- colMeans(eval_all, na.rm=TRUE); cvsd <- sqrt(colMeans(scale(eval_all, cvm, FALSE)^2, na.rm = TRUE)/(nfolds-1))
  id_max <- which.max(cvm); id_lam1 <- ceiling(id_max/(n2*n3))
  id_gamma <- ceiling((id_max-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_max-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_max <- lam1[id_lam1]; gamma_max <- gamma[id_gamma]
  lam2_max <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  start_time <- Sys.time()
  if(is.null(M) || is.null(U)) {
    fit <- seas(x, y, yclass = yclass, type = type, FUN = FUN, lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  } else {
    fit <- seas(M = M, U = U, nobs = NROW(x), lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }
  if(fit$code != 0) { code <- 5; warning("The estimated beta is null."); return(list(beta = NULL, code = code)) }
  B <- fit$B; beta <- fit$beta; rank <- fit$rank
  end_time <- Sys.time(); time2 <- difftime(end_time, start_time, units = "secs")
  time <- time1 + time2
  list(beta = beta, B = B, rank = rank, eval = eval_all, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_max = lam1_max, lam2_max = lam2_max, gamma_max = gamma_max, code = code, time = time)
}

admm <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, ...) {
  if(is.null(dim(U))) U <- cbind(U, rep(0, length(U)))
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
  M0 <- M; U0 <- U
  n1 <- length(lam1); n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2)); n3 <- length(gam)
  nparams <- n1 * n2 * n3
  B_l <- vector("list", nparams); beta_l <- vector("list", nparams)
  step_l <- rep(NA_integer_, nparams); time_l <- rep(NA_real_, nparams)
  rank_l <- rep(NA_integer_, nparams); s_l <- rep(NA_integer_, nparams)
  nlam_cvg <- 0
  for(i in 1:n1) {
    lambda1 <- as.double(lam1[i])
    for(j in 1:n3) {
      gamma <- gam[j]
      for(k in 1:n2) {
        lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j, k])
        M <- M0 + gamma * diag(rep(1, ncol(M0)), ncol(M0), ncol(M0))
        Bold <- matrix(0, dim(U0)[1], dim(U0)[2])
        Cold <- matrix(0, dim(U0)[1], dim(U0)[2])
        etaold <- matrix(0, dim(U0)[1], dim(U0)[2])
        step <- 0; start_time <- Sys.time()
        repeat {
          step <- step + 1
          U <- U0 - etaold + gamma * Cold
          out_B <- updateB(M, U, lambda1, opts); Bnew <- out_B$Bnew; jerr <- out_B$jerr
          if(jerr != 0) break
          Cnew <- updateC(Bnew, lambda2, gamma, etaold)
          etanew <- etaold + gamma * (Bnew - Cnew)
          if(max(abs(Bnew - Cnew)) < eps) { jerr <- 1; break }
          if(step > maxit) { jerr <- 404; warning('Maximal iteration is reached.'); break }
          Bold <- Bnew; Cold <- Cnew; etaold <- etanew
        }
        end_time <- Sys.time(); time <- difftime(end_time, start_time, units = "secs")
        if(jerr < -10000) break
        if(jerr == 1) {
          index <- (i-1)*n2*n3 + (j-1)*n2 + k; nlam_cvg <- nlam_cvg + 1
          B_l[[index]] <- Bnew; step_l[index] <- step; time_l[index] <- time
          if(is.null(d)) rank <- rank_func(Cnew, thrd = eps) else rank <- d
          rank_l[index] <- rank
          if(rank == 0) { beta <- matrix(0, nrow(Bnew), ncol(Bnew)) }
          else { tmp <- svd(Bnew)$u[, 1:rank, drop = FALSE]; vec <- as.vector(tmp); vec[abs(vec) < eps] <- 0; beta <- matrix(vec, nrow(tmp), ncol(tmp)) }
          beta_l[[index]] <- beta
          var_ind <- apply(beta, 1, function(x) any(x != 0)); s_l[index] <- sum(var_ind)
        }
      }
      if(jerr < -10000) break
    }
  }
  return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg))
}

updateB <- function(M, U, lambda1, opts) {
  U <- t(U)
  fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
  if(fit$jerr != 0) return(list(Bnew = NULL, jerr = fit$jerr))
  outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
  Bnew <- as.matrix(outlist$theta[[1]])
  list(Bnew = Bnew, jerr = fit$jerr)
}

updateC <- function(Bnew, lambda2, gamma, etaold) {
  Btemp <- Bnew + 1/gamma * etaold
  svd_B <- svd(Btemp); lamtemp <- pmax(0, svd_B$d - lambda2/gamma)
  Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
  Cnew
}

msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', FUN = NULL, lambda.factor=NULL, nlambda=100, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, verbose=FALSE, perturb=NULL) {
  if(is.null(M) || is.null(U)) {
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)) {
      if(categorical == FALSE) {
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      } else { yclass <- y }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, type, FUN); M <- MU_out$M; U <- MU_out$U
    nobs <- as.integer(dim(x)[1]); nvars <- as.integer(dim(x)[2])
  } else {
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass) <= nvars, 0.2, 1e-03)
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax * 2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) diag(M) <- diag(M) + perturb
  H <- as.integer(dim(U)[2])
  if (length(pf) != nvars) stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit); verbose <- as.integer(verbose); sml <- as.double(sml)
  pf <- as.double(pf); eps <- as.double(eps); dfmax <- as.integer(dfmax); pmax <- as.integer(pmax)
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor); ulam <- double(1)
  } else {
    flmin <- as.double(1)
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda))); nlam <- as.integer(length(lambda))
  }
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
  outlist <- formatoutput(fit, maxit, pmax, nvars, H)
  rank <- rep(NA_integer_, length(outlist$theta))
  for (i in 1:length(outlist$theta)) if(!is.null(outlist$theta[[i]])) rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
  if(is.null(lambda)) outlist$lambda <- lamfix(outlist$lambda)
  outlist <- list(lambda = outlist$lambda, theta = outlist$theta, M = M, U = U, rank = rank)
  class(outlist) <- c("msda"); outlist
}

cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', lambda.factor=NULL, nlambda=100, nfolds=5, foldid = NULL, lambda = NULL, FUN = NULL, maxit = 1e3, plot = FALSE) {
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)) {
    if(categorical == FALSE) {
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    } else { yclass <- y }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  nobs <- nrow(x); nvars <- ncol(x); nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass) <= nvars, 0.2, 1e-03)
  fit <- msda(x, y, yclass = yclass, type = type, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
  lambda <- fit$lambda; beta_l <- fit$theta; M <- fit$M; U <- fit$U
  rank_l <- fit$rank; beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  if(is.null(foldid)) {
    ord <- order(y); y <- y[ord]; yclass <- yclass[ord]; x <- x[ord, ]
    count <- as.numeric(table(yclass)); foldid <- c()
    for(cnt in count) foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
  } else { nfolds <- length(unique(foldid)) }
  cv_out <- lapply(1:nfolds, function(k) {
    x_train <- x[foldid!=k, , drop=FALSE]; x_val <- x[foldid==k, , drop=FALSE]
    y_train <- y[foldid!=k]; y_val <- y[foldid==k]; yclass_train <- yclass[foldid!=k]
    fit_fold <- msda(x_train, y_train, yclass_train, type = type, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
    M_fold <- fit_fold$M; U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta; rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    eval_fold <- eval_dc(beta_fold, x_val, y_val)
    if(length(eval_fold) != length(lambda)) eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2); U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)) eval_all <- t(as.matrix(eval_all))
  if(all(is.na(eval_all))) return(NULL)
  cvm <- colMeans(eval_all, na.rm = TRUE)
  id_max <- which.max(cvm); lam_max <- lambda[id_max]; beta <- as.matrix(beta_l[[id_max]])
  rank <- rank_func(beta, thrd = 1e-3)
  list(beta = beta, id = id_max, lambda = lambda, lam_max = lam_max, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}

################### TG-method & TM-method 核心函数 ###################################

tri_M <- function(X, B) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(B)) B <- as.matrix(B)
  n <- nrow(X); p <- ncol(X)
  if (n < 2) stop("Need at least two observations.")
  if (nrow(B) != p) stop("B must have p rows.")
  K <- ncol(B)
  if (K < 1) stop("B must have at least one column.")
  S <- cov(X)
  Tmat <- X %*% B; cosT <- cos(Tmat); sinT <- sin(Tmat)
  V_cos <- crossprod(X, cosT) / n; V_sin <- crossprod(X, sinT) / n
  m_sin <- colMeans(sinT); m_cos <- colMeans(cosT)
  G_cos <- -B %*% diag(m_sin, K); G_sin <- B %*% diag(m_cos, K)
  Est_sin <- admm(S, V_sin - S %*% G_sin, n, p, 0, 0, 1, maxit=100000)$B[[1]]
  Est_sin
}

U_beta <- function(x, B) {
  x <- as.matrix(x); B <- as.matrix(B)
  n <- nrow(x); p <- ncol(x); k <- ncol(B)
  mat <- matrix(0, p, k)
  for (j in 1:n) {
    tBx <- crossprod(B, x[j, ])
    w <- as.numeric(tBx)^2
    mat <- mat + (1/n) * (x[j, ] %o% w)
  }
  mat
}

###################################################
# Part III: 模型定义（系统化）
###################################################

define_models <- function() {
  models <- list()
  
  ## Model 1: K=2, d_true=1
  K <- 2; r <- 20
  sigma0 <- diag(r) * 0.5 + matrix(0.5, r, r)
  mu1 <- rep(0, r); mu1[1:5] <- c(1, 1, 0, -1, -1)
  mu1 <- sigma0 %*% mu1; mu1 <- mu1 / norm(mu1, type="2") * 1.5
  mu2 <- -mu1
  mu_mat <- cbind(mu1, mu2)
  Sig <- array(NA, c(r, r, K)); Sig[,,1] <- sigma0; Sig[,,2] <- sigma0
  models[["Model1"]] <- list(mu=mu_mat, Sigma=Sig, K=K, r=r, pi0=c(0.65, 0.35), d_true=1)
  
  ## Model 2: K=3, d_true=1
  K <- 3; r <- 20
  sigma0 <- cov_decay(r, 0.5)
  mu1 <- rep(0, r)
  mu1[(as.integer(r/2-2)):(as.integer(r/2+3))] <- c(1,-1,1,-1,1,-1)
  mu1 <- sigma0 %*% mu1; mu1 <- mu1 / norm(mu1, type="2") * 2.5
  mu2 <- rep(0, r); mu3 <- -mu1
  gamma_v <- mu1 / 2.5
  cov1 <- sigma0 - exp(-1) * gamma_v %*% t(gamma_v)
  cov2 <- sigma0 - exp(-2) * gamma_v %*% t(gamma_v)
  cov3 <- sigma0 - exp(-3) * gamma_v %*% t(gamma_v)
  mu_mat <- cbind(mu1, mu2, mu3)
  Sig <- array(NA, c(r, r, K)); Sig[,,1] <- cov1; Sig[,,2] <- cov2; Sig[,,3] <- cov3
  models[["Model2"]] <- list(mu=mu_mat, Sigma=Sig, K=K, r=r, pi0=c(0.5, 0.3, 0.2), d_true=1)
  
  ## Model 3: K=3, d_true=2
  K <- 3; r <- 20
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
  mu_mat <- cbind(mu1, mu2, mu3)
  Sig <- array(NA, c(r, r, K)); Sig[,,1] <- cov1; Sig[,,2] <- cov2; Sig[,,3] <- cov3
  models[["Model3"]] <- list(mu=mu_mat, Sigma=Sig, K=K, r=r, pi0=c(0.55, 0.25, 0.2), d_true=2)
  
  return(models)
}

###################################################
# Part IV: 数据生成
###################################################

generate_data <- function(model, N) {
  K <- model$K; r <- model$r; pi0 <- model$pi0
  dat <- matrix(NA, N, r)
  idx <- sample(1:K, N, replace=TRUE, prob=pi0)
  for (j in 1:K) {
    nj <- sum(idx == j)
    if (nj > 0) {
      x_tmp <- mvrnorm(nj, mu=model$mu[, j], Sigma=model$Sigma[,,j])
      if (nj == 1) x_tmp <- matrix(x_tmp, nrow=1)
      dat[idx == j, ] <- x_tmp
    }
  }
  list(dat=dat, idx=idx)
}

###################################################
# Part V: 各方法实现
###################################################

## ---------- TG 方法 ----------
run_TG <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  M1 <- tri_M(dat_c, diag(r))
  result <- svd(M1)$u[, 1:min(4, r), drop=FALSE]
  
  M2 <- tri_M(dat_c, result)
  result_matrix <- svd(M2)
  d_use <- max(TDRR(result_matrix$d, N)$q_hat, 1)
  result <- result_matrix$u[, 1:d_use, drop=FALSE]   # [修复] 原来写的 [,d_use]
  
  new_dat <- dat_c %*% result
  cls_result <- run_em(new_dat, G=K)
  time_elapsed <- proc.time()[3] - time_start
  
  test_proj <- test_c %*% result
  pred <- predict(cls_result, newdata=test_proj)
  cls_label <- pred$classification
  

  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed, d_est=d_use)
}

## ---------- Third Moment (TM) 方法 ----------

run_TM <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  M <- cov(dat_c)
  
  U <- U_beta(dat_c, diag(r))
  B1 <- admm(M, U, N, r, 0, 0, 1, maxit=20000)$B[[1]]
  result <- svd(B1)$u[, 1:min(4, r), drop=FALSE]
  
  U <- U_beta(dat_c, result)
  B2 <- admm(M, U, N, r, 0, 0, 1, maxit=20000)$B[[1]]
  result_matrix <- svd(B2)
  d_use <- max(TDRR(result_matrix$d, N)$q_hat, 1)
  result <- result_matrix$u[, 1:d_use, drop=FALSE]   # [修复] 原代码此处仍用旧 result
  
  new_dat <- dat_c %*% result
  cls_result <- run_em(new_dat, G=K)
  
  time_elapsed <- proc.time()[3] - time_start

  test_proj <- test_c %*% result
  pred <- predict(cls_result, newdata=test_proj)
  cls_label <- pred$classification
  

  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed, d_est=d_use)
}

## ---------- SVM 方法 ----------
run_SVM <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  
  time_start <- proc.time()[3]
  
  M_cov <- cov(dat)
  tildex <- mvrnorm(N, rep(0, r), M_cov)
  x_combined <- rbind(dat, tildex)
  y_binary <- c(rep(1, N), rep(-1, N))
  
  H <- 20
  step <- 1 / H
  pi.grid <- seq(step, 1 - step, by = step)
  w <- matrix(0, r, H - 1)
  
  z <- t(admm(matpower(cov(x_combined), 1/2), t(x_combined),
              2*N, r, 0, 0, 10, maxit=1000)$B[[1]])
  
  type_svm <- 3; cost <- 0.5; epsilon <- 1e-5
  for (h in 1:(H - 1)) {
    weights <- c("1" = 1 - pi.grid[h], "-1" = pi.grid[h])
    result_svm <- LiblineaR(data=z, target=as.factor(y_binary),
                            type=type_svm, cost=cost, epsilon=epsilon, wi=weights)
    coef <- result_svm$W
    w[, h] <- coef[1:r]
  }
  
  result_matrix <- svd(w)
  d_use <- max(TDRR(result_matrix$d, N)$q_hat, 1)  # [修复] 保护 d_use >= 1
  
  B_refined <- admm(matpower(cov(x_combined), 1/2),
                    result_matrix$u[, 1:d_use, drop=FALSE],
                    r, r, 0, 0, 10, maxit=10000)$B[[1]]
  result_matrix2 <- svd(B_refined)
  obj <- result_matrix2$u[, 1:d_use, drop=FALSE]
  
  new_dat <- dat %*% obj
  cls_result <- run_em(new_dat, G=K)
  
  test_proj <- test_dat %*% obj
  pred <- predict(cls_result, newdata=test_proj)
  cls_label <- pred$classification
  
  time_elapsed <- proc.time()[3] - time_start
  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed, d_est=d_use)
}

## ---------- Spectral (PCA) 方法 ----------
run_Spectral <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  M <- cov(dat_c)
  result <- eigen(M)$vectors[, 1:d_true, drop=FALSE]   # [修复] 用 d_true
  
  new_dat <- dat_c %*% result
  cls_result <- run_em(new_dat, G=K)
  
  # [修复] 在测试集上预测（原代码只看训练集）
  test_proj <- test_c %*% result
  pred <- predict(cls_result, newdata=test_proj)
  cls_label <- pred$classification
  
  time_elapsed <- proc.time()[3] - time_start
  mis_rate <- 1 - cal_acc(idx_test, cls_label)          # [修复] 用 idx_test
  list(mis_rate=mis_rate, time=time_elapsed)
}

## ---------- EM 方法（全维度） ----------
run_EM <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  cls_result <- run_em(dat_c, G=K)
  
  # 在测试集上预测
  pred <- predict(cls_result, newdata=test_c)
  cls_label <- pred$classification
  
  time_elapsed <- proc.time()[3] - time_start
  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed)
}

## ---------- MCFA 方法 ----------
predict_mcfa <- function(mcfa_fit, newdata) {
  # 根据 MCFA 模型参数计算测试集各分量后验概率
  g <- length(mcfa_fit$pivec)
  n <- nrow(newdata)
  log_post <- matrix(0, n, g)
  for (k in 1:g) {
    Sigma_k <- mcfa_fit$B %*% t(mcfa_fit$B) + diag(mcfa_fit$D[, k])
    log_post[, k] <- dmvnorm(newdata, mean=mcfa_fit$mu[, k],
                             sigma=Sigma_k, log=TRUE) + log(mcfa_fit$pivec[k])
  }
  apply(log_post, 1, which.max)
}

run_MCFA <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  # [修复] 用 K 和 d_true 代替硬编码的 2, 2
  mcfa_result <- mcfa(dat_c, g=K, q=max(d_true, 1), nkmeans=5)
  
  # [修复] 在测试集上预测
  cls_label <- predict_mcfa(mcfa_result, test_c)
  
  time_elapsed <- proc.time()[3] - time_start
  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed)
}

## ---------- CLEMM 方法 ----------
## 修复：初始化代码用 run_init() 代替有 bug 的拷贝粘贴
## 修复：用 d_true 代替硬编码；定义 opts
## 修复：在测试集上预测
## 注意：需要 source 加载 clemm_em() 及相关函数
predict_clemm <- function(clemm_fit, newdata, K) {
  # 根据 CLEMM 拟合结果的参数对测试数据进行预测
  # 假设 clemm_fit 包含 $mu (r x K), $cov (r x r x K), $wt (length K)
  # 请根据实际 clemm_em 的返回格式调整字段名
  n <- nrow(newdata)
  log_post <- matrix(0, n, K)
  for (k in 1:K) {
    log_post[, k] <- dmvnorm(newdata, mean=clemm_fit$mu[, k],
                             sigma=clemm_fit$cov[,,k], log=TRUE) + log(clemm_fit$wt[k])
  }
  apply(log_post, 1, which.max)
}

run_CLEMM <- function(dat, test_dat, idx_test, K, r, N, d_true) {
  train_mean <- colMeans(dat)
  dat_c  <- sweep(dat, 2, train_mean)
  test_c <- sweep(test_dat, 2, train_mean)
  
  time_start <- proc.time()[3]
  
  # [修复] 使用 run_init 代替有错误变量名的拷贝粘贴代码
  init_km <- run_init(dat_c, G=K, replicate=5)
  
  init <- list()
  init$centers <- t(init_km$centers)
  init$wt <- init_km$size / N
  init$cov <- array(NA, c(r, r, K))
  for (j in 1:K) {
    init$cov[,,j] <- cov(dat_c[init_km$cluster == j, , drop=FALSE])
  }
  
  # [修复] u 用 d_true；定义 opts
  opts <- list()    # 根据 clemm_em 的要求设置
  res_clemm <- clemm_em(dat_c, K, u=d_true, iter=500, opts=opts, init=init, stopping=1e-5)
  
  # [修复] 在测试集上预测（原代码用 clustering_err 只看训练集）
  cls_label <- predict_clemm(res_clemm, test_c, K)
  
  time_elapsed <- proc.time()[3] - time_start
  mis_rate <- 1 - cal_acc(idx_test, cls_label)
  list(mis_rate=mis_rate, time=time_elapsed)
}


###################################################
# Part VI: 系统化实验主循环
###################################################

run_experiment <- function(N_train = 1000, N_test = 3000, n_rep = 200,
                           method_names = c("TG", "TM", "Spectral", "EM")) {
  
  models <- define_models()
  
  # 存储所有结果
  all_results <- list()
  
  for (model_name in names(models)) {
    cat(sprintf("\n==================== %s ====================\n", model_name))
    model  <- models[[model_name]]
    K      <- model$K
    r      <- model$r
    pi0    <- model$pi0
    d_true <- model$d_true
    
    cat(sprintf("  K=%d, r=%d, d_true=%d, pi0=(%s)\n",
                K, r, d_true, paste(pi0, collapse=", ")))
    
    # 初始化结果矩阵
    mis_rates <- matrix(NA, n_rep, length(method_names),
                        dimnames=list(NULL, method_names))
    times     <- matrix(NA, n_rep, length(method_names),
                        dimnames=list(NULL, method_names))
    
    for (l in 1:n_rep) {
      if (l %% 20 == 0) cat(sprintf("  Replicate %d/%d\n", l, n_rep))
      
      # 每次实验所有方法共享同一份训练/测试数据
      train_data <- generate_data(model, N_train)
      dat <- train_data$dat; idx <- train_data$idx
      
      test_data <- generate_data(model, N_test)
      test_dat <- test_data$dat; idx_test <- test_data$idx
      
      for (method in method_names) {
        res <- tryCatch({
          switch(method,
                 "TG"       = run_TG(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "TM"       = run_TM(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "SVM"      = run_SVM(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "Spectral" = run_Spectral(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "EM"       = run_EM(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "MCFA"     = run_MCFA(dat, test_dat, idx_test, K, r, N_train, d_true),
                 "CLEMM"    = run_CLEMM(dat, test_dat, idx_test, K, r, N_train, d_true)
          )
        }, error = function(e) {
          cat(sprintf("    [%s] Error in rep %d: %s\n", method, l, e$message))
          return(NULL)
        })
        
        if (!is.null(res)) {
          idx_m <- which(method_names == method)
          mis_rates[l, idx_m] <- res$mis_rate
          times[l, idx_m]     <- res$time
        }
      }
    }
    
    # 汇总结果
    all_results[[model_name]] <- list(mis_rates=mis_rates, times=times)
    
    cat(sprintf("\n--- %s 结果汇总 (n_rep=%d) ---\n", model_name, n_rep))
    cat("平均误分类率 (Mean):\n")
    print(round(colMeans(mis_rates, na.rm=TRUE), 4))
    cat("误分类率标准差 (SD):\n")
    print(round(apply(mis_rates, 2, sd, na.rm=TRUE), 4))
    cat("平均计算时间/秒 (Mean Time):\n")
    print(round(colMeans(times, na.rm=TRUE), 4))
    cat("\n")
  }
  
  return(all_results)
}


###################################################
# Part VII: 运行实验
###################################################


methods_to_run <- c("TG", "TM", "Spectral", "EM","SVM","MCFA","CLEMM")

results <- run_experiment(
  N_train      = 1000,
  N_test       = 3000,
  n_rep        = 200,
  method_names = methods_to_run
)
