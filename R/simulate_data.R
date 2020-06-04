



.sim_10000_tcop_blockdiag <- function(n, cmat, df){
  
  chol_ <- t(chol(cmat))
  
  rn <- lapply(1:10, function(i) t(matrix(rnorm(n*1000), nrow =n)))
  R <- rchisq(n, df)
  w <- sqrt(df/R)
  
  V <- (cbind(t(chol_%*%rn[[1]])*w, t(chol_%*%rn[[2]])*w,
              t(chol_%*%rn[[3]])*w, t(chol_%*%rn[[4]])*w,
              t(chol_%*%rn[[5]])*w, t(chol_%*%rn[[6]])*w,
              t(chol_%*%rn[[7]])*w, t(chol_%*%rn[[8]])*w,
              t(chol_%*%rn[[9]])*w, t(chol_%*%rn[[10]])*w))
  return(apply(V, 2, pt, df=df))  
}

.sim_4000_tcop_blockdiag <- function(n, cmat, df){
  
  chol_ <- t(chol(cmat))
  
  rn <- lapply(1:4, function(i) t(matrix(rnorm(n*1000), nrow =n)))
  R <- rchisq(n, df)
  w <- sqrt(df/R)
  
  V <- (cbind(t(chol_%*%rn[[1]])*w, t(chol_%*%rn[[2]])*w,
              t(chol_%*%rn[[3]])*w, t(chol_%*%rn[[4]])*w))
  return(apply(V, 2, pt, df=df))  
}


simulate_X <- function(n, p, copula, margins, blockdiag_t, corr, df){
  
  inv_samps <- c(function(q) stats::qbinom(q, 1, 0.2),
                 function(q) stats::qbinom(q, 1, 0.4),
                 function(q) stats::qbinom(q, 1, 0.6),
                 function(q) stats::qbinom(q, 1, 0.8),
                 function(q) stats::qbeta(q, 1, 2),
                 function(q) stats::qbeta(q, 2, 1),
                 function(q) stats::qbeta(q, 2, 2),
                 function(q) stats::qgamma(q, 1, 3),
                 function(q) stats::qgamma(q, 3, 1),
                 function(q) stats::qgamma(q, 3, 3),
                 function(q) stats::qnorm(q),
                 function(q) stats::qt(q, 3),
                 function(q) stats::qt(q, 4),
                 function(q) stats::qt(q, 6),
                 function(q) stats::qpois(q, lambda = 1),
                 function(q) stats::qpois(q, lambda = 3),
                 function(q) stats::qpois(q, lambda = 5))
  
  if (blockdiag_t){
    if (p == 10000){
      X <- .sim_10000_tcop_blockdiag(n=n, corr, df)
    } else if (p == 4000){
      X <- .sim_4000_tcop_blockdiag(n=n, corr, df)
    } else {
      stop("If blockdiag_t==TRUE, p must be either 10000 or 4000.")
    }
    
  } else if (is.null(copula)) {
    X <- copula::rCopula(copula=copula::claytonCopula(param=c(2), dim=p),
                         n=n)
  } else {
    X <- copula::rCopula(copula=copula, n=n)
  }
  
  if (is.null(margins)) {
    col_dists <- sample(inv_samps, p, replace=TRUE)
    
    for (j in 1:p) {
      X[, j] <- col_dists[[j]](X[, j])
    }
  } else {
    for (j in 1:p) {
      X[, j] <- inv_samps[[margins[j]]](X[, j])
    }
  }
  
  return(X)
  
}


simulate_data <- function(n, p, marg_p, spars=NULL, noise=0, link="logit",
                          copula=NULL, betas=NULL, margins=NULL, blockdiag_t=F,
                          corr=NULL, df=2){


  if (is.null(spars) & is.null(betas)){
    stop("Must specify either spars or betas")
  }
  
  if (blockdiag_t & is.null(corr)){
    stop("Must specify 1000x1000 correlation matrix when blockdiag_t == TRUE.")
  }
    
  X <- simulate_X(n, p, copula, margins, blockdiag_t, corr, df)
  
  if (is.null(betas)){

    betas <- stats::rnorm(p, 0, 1)
    B <- stats::rbinom(p, 1, spars)

    while (sum(B) == 0 & spars > 0){
      B <- stats::rbinom(p, 1, spars)
    }

    betas <- betas*B

  }

  noise_term <- stats::rnorm(n, 0, sd(X%*%betas)*noise)

  if (length(which(c("logit", "invNormCDF", "cloglog") == link)) >= 1) {
    p_func <- switch(which(c("logit", "invNormCDF", "cloglog") == link),
                     stats::plogis, stats::pnorm, function(x) 1 - exp(-exp(x)))
    q_func <- switch(which(c("logit", "invNormCDF", "cloglog") == link),
                     stats::qlogis, stats::qnorm, function(x) log(-log(1 - x)))

  } else {
    stop(paste0(c("Link function ", link, " not implemented.\n")))
  }

  p_alpha <- function(alpha, X, betas, noise_term){
    mean(p_func(alpha + X%*%betas + noise_term))
  }

  alpha <- stats::optim(fn= function(alpha){(p_alpha(alpha, X, betas,
                                                        noise_term) - marg_p)**2},
                        par=-mean(X%*%betas) + q_func(marg_p), method="BFGS",
                        control=list(maxit=1000, fnscale=0.01))$par

  score <- alpha + X%*%betas + noise_term
  Y <- stats::rbinom(n, 1, prob = p_func(score))

  return(list(alpha=alpha, beta=betas, X=X, Y=Y, score=score))

}

make_xgb_rule <- function(n_make_rule, p, nonzero_cols, copula, margins,
                          max_depth, nrounds, blockdiag_t=F,
                          corr=NULL, df=2){

  if (blockdiag_t & is.null(corr)){
    stop("Must specify 1000x1000 correlation matrix when blockdiag_t == TRUE.")
  }
  
  X <- simulate_X(n_make_rule, p, copula, margins, blockdiag_t, corr, df)

  B_lat <- rbinom(n_make_rule, 1, 0.5)
  latent <- B_lat*rexp(n_make_rule, 0.2) - (1 - B_lat)*rexp(n_make_rule, 0.1)

  rulemd <- xgboost::xgb.train(params=list(eta=0.1, max_depth=max_depth),
                               nrounds=nrounds,
                               data=xgboost::xgb.DMatrix(X[1:n_make_rule,
                                                           nonzero_cols],
                                                         label=latent))

  return(rulemd)

}

simulate_data_from_xgb <- function(n, p, marg_p, copula=NULL, margins=NULL,
                                   rulemd=NULL, nonzero_cols=NULL, blockdiag_t=F,
                                   corr=NULL, df=2){

  if (is.null(nonzero_cols)){
    nonzero_cols <- 1:p
  }
  
  if (blockdiag_t & is.null(corr)){
    stop("Must specify 1000x1000 correlation matrix when blockdiag_t == TRUE.")
  }
  
  X <- simulate_X(n, p, copula, margins, blockdiag_t, corr, df)
  
  sc <- predict(rulemd,
                newdata=xgboost::xgb.DMatrix(
                  X[, nonzero_cols]))

  p_alpha <- function(alpha, sc){
    mean(stats::plogis(alpha + sc))
  }
  
  alpha <- stats::optim(fn= function(alpha){(p_alpha(alpha, sc) - marg_p)**2},
                        par=-mean(sc) + stats::qlogis(marg_p), method="BFGS",
                        control=list(maxit=1000, fnscale=0.01))$par
  
  score <- alpha + sc
  Y <- stats::rbinom(n, 1, prob=stats::plogis(score))
  
  return(list(alpha=alpha, X=X,
              Y=Y, score=score))
    
}
