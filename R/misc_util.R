
fit_xgb_md <- function(x, y, par, nrounds){
  
  md <- xgboost::xgb.train(params=par,
                           data=xgboost::xgb.DMatrix(x, label=y),
                           nrounds=nrounds)
  return(md)
}

score_xgb <- function(x, md){
  D <- xgboost::xgb.DMatrix(x)
  return(
    sapply(1:md$niter, function(ntree) predict(md, D, ntreelimit=ntree))
  )
}

c_err_from_score <- function(sc, y){
  
  c_err <- function(yhat, y){
    mean(abs(yhat - y))
  }
  
  return(
    sapply(1:99, function(i) apply(apply(sc, 2, .pred_fun, q=i/100),
                                   2, c_err, y=y))
  )
}

f_err_from_score <- function(sc, y){
  
  f_err <- function(yhat, y){
    sum(yhat*(1 - y))/sum(yhat)
  }
  
  return(
    sapply(1:99, function(i) apply(apply(sc, 2, .pred_fun, q=i/100),
                                   2, f_err, y=y))
  )
}

comp_auc <- function(y, x){

  #
  # Copied from glmnet after 'auc' no longer exported from 'namepace:glmnet'
  #

  rprob <- rank(x)
  n1 <- sum(y)
  n0 <- length(y)-n1
  u <- sum(rprob[y==1])-n1*(n1+1)/2

  return(
    exp(log(u) - log(n1) - log(n0))
  )
}