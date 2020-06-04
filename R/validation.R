
.pred_fun <- function(x, q){
  k <- max(1, round(q*length(x)))
  thresh <- x[order(-x)[k]]
  return(1*(x >= thresh))
}

.cmp_agg <- function(lambda, type, y, nfolds, folds, prop_select, error, fraud,
                     auc, lab) {

  # compute aggregates
  if (type == "boot") {

    I_ib <- matrix(0, nrow=length(y), ncol=nfolds)
    for (f in 1:nfolds){
      I_ib[(1:length(y))[-folds[[f]]], f] <- 1
    }

    multfact1 <- apply(I_ib, 2, sum)
    multfact2 <- lapply(1:length(prop_select),
                        function(v) sapply(1:nfolds, 
                                           function(f) apply(lab[[v]][[f]],
                                                             2, sum)))
    multfact3 <- sapply(1:nfolds,
                        function(f) sum(I_ib[, f]*y)*sum(I_ib[, f]*(1 - y)))

    class_agg <- lapply(1:length(prop_select),
                        function(v) (error[[v]]%*%multfact1)/sum(multfact1))
    fraud_agg <- lapply(1:length(prop_select),
                        function(v) apply(fraud[[v]]*multfact2[[v]], 1, sum) /
                          apply(multfact2[[v]], 1, sum))
    auc_agg <- (auc%*%multfact3)/sum(multfact3)

  } else {
    class_agg <- lapply(1:length(prop_select),
                        function(v) apply(error[[v]], 1, mean))
    fraud_agg <- lapply(1:length(prop_select),
                        function(v) apply(fraud[[v]], 1, mean))
    auc_agg <- apply(auc, 1, mean)

  }

  return(list(lambda=lambda, class=error, fraud=fraud, auc=auc,
              class_agg=class_agg, fraud_agg=fraud_agg, auc_agg=auc_agg))

}


validate_fraud_glmnet <- function(x, y, prop_select, lambda=NULL,
                                  lambda_min_ratio=10**(-5), nlambda=50,
                                  alpha=0, folds=NULL, nfolds=10, stratified=F,
                                  parallel_val=F, val_type="cv", cv_repeat=1) {

  if(any(is.null(prop_select) | (prop_select > 1 | prop_select < 0))){
    stop("You must specify prop_select as a number between 0 and 1, or
         a list of numbers between 0 and 1.")
  }

  # Check if folds are supplied, compute folds otherwise
  if (is.null(folds)){
    folds <- .draw_folds(y, nfolds, stratified, val_type, cv_repeat)
    nfolds <- folds$nfolds
    folds <- folds$folds
  } else {
    nfolds <- length(folds)
  }

  # Check if values of lambda are supplied
  if (is.null(lambda)) {
    lambda <- glmnet::glmnet(x=x, y=y, alpha=alpha, family="binomial",
                             nlambda=nlambda,
                             lambda.min.ratio=lambda_min_ratio)$lambda
  }

  # Set up validation datasets
  sets <- list(train=lapply(1:nfolds,
                            function(f) list(x=x[folds[[f]], ],
                                                   y=y[folds[[f]]])),
               test=lapply(1:nfolds,
                           function(f) list(x=x[-folds[[f]], ],
                                                  y=y[-folds[[f]]])))

  # Fit models to the split datasets
  if (parallel_val){
    no_cores <- parallel::detectCores()
    if(is.nan(no_cores)){
      no_cores <- 4
    }
    cl <- parallel::makeCluster(no_cores)
    models <- parallel::parLapply(cl, 1:nfolds,
                                  function(f)
                                    glmnet::glmnet(x=sets$train[[f]]$x,
                                                   y=sets$train[[f]]$y,
                                                   family="binomial",
                                                   lambda=lambda, alpha=alpha))

    parallel::stopCluster(cl)

  } else {
    models <- lapply(1:nfolds,
                     function(f) glmnet::glmnet(x=sets$train[[f]]$x,
                                                y=sets$train[[f]]$y,
                                                family="binomial",
                                                lambda=lambda, alpha=alpha))
  }

  # Compute the scores X%*%beta for each fold and each value of lambda
  score <- lapply(1:nfolds,
                  function(f) sets$test[[f]]$x %*% models[[f]]$beta)

  # Compute the predicted labels with the prediction rule of classifing any
  # observation with a score over the quantile prop_select as 1.
  lab <- lapply(1:length(prop_select),
                function(v) lapply(1:nfolds,
                                   function(f) apply(score[[f]], 2, .pred_fun,
                                                     q=prop_select[v])))

  # Compute the classification error for these predictions
  test_error <- lapply(1:length(prop_select),
                       function(v)
                         sapply(1:nfolds,
                                function(f)
                                  apply(lab[[v]][[f]], 2,
                                        function(yhat)
                                          mean(abs(yhat - sets$test[[f]]$y))
                                        )
                                )
                       )

  # Compute the "fraud loss", the error in terms of the number of false
  # positives for these predictions
  test_fraud <- lapply(1:length(prop_select),
                       function(v)
                         sapply(1:nfolds,
                                function(f)
                                  apply(lab[[v]][[f]], 2,
                                        function(yhat)
                                          sum((1-sets$test[[f]]$y)*yhat) /
                                          sum(yhat)
                                        )
                                )
                       )

  # Also compute the auc for the predicted scores
  test_auc <- sapply(1:nfolds, function(f) apply(score[[f]], 2, comp_auc,
                                                 y=sets$test[[f]]$y))

  return(.cmp_agg(lambda, val_type, y, nfolds, folds, prop_select,
                  test_error, test_fraud, test_auc, lab))

}


validate_fraud_xgboost <- function(x, y, param, nrounds, prop_select,
                                   folds=NULL, nfolds=10, stratified=F,
                                   parallel_val=F, val_type="cv",
                                   cv_repeat=1) {

  if(any(is.null(prop_select) | (prop_select > 1 | prop_select < 0))){
    stop("You must specify prop_select as a number between 0 and 1, or a list of
         numbers between 0 and 1.")
  }

  # Check if folds are supplied, compute folds otherwise
  if (is.null(folds)){
    folds <- .draw_folds(y, nfolds, stratified, val_type, cv_repeat)
    nfolds <- folds$nfolds
    folds <- folds$folds
  } else {
    nfolds <- length(folds)
  }

  # Set up validation datasets
  sets <- list(train=lapply(1:nfolds,
                            function(f)
                              xgboost::xgb.DMatrix(data=x[folds[[f]], ],
                                                   label=y[folds[[f]]])),
               test=lapply(1:nfolds,
                           function(f)
                             xgboost::xgb.DMatrix(data=x[-folds[[f]], ],
                                                  label=y[-folds[[f]]]))
               )

  # Fit models to the split datasets
  if (parallel_val){
    no_cores <- parallel::detectCores()  
    #cl <- parallel::makeCluster(no_cores)
    #parallel::clusterExport(cl, varlist=c("param", "nrounds", "sets"))
    #models <- parallel::parLapply(cl, 1:nfolds,
    #                              function(f)
    #                                xgboost::xgb.train(params=param,
    #                                                   data=sets$train[[f]],
    #                                                   nrounds=nrounds))
    models <- lapply(1:nfolds,
                     function(f) xgboost::xgb.train(params=param,
                                                    data=sets$train[[f]],
                                                    nrounds=nrounds, nthread=no_cores))
    #parallel::stopCluster(cl)

  } else {
    models <- lapply(1:nfolds,
                     function(f) xgboost::xgb.train(params=param,
                                                    data=sets$train[[f]],
                                                    nrounds=nrounds))
  }

  # Compute the scores X%*%beta for each fold and each value of lambda
  score <- lapply(1:nfolds,
                  function(f) sapply(1:nrounds,
                                     function(ntree)
                                       predict(models[[f]],
                                               newdata=sets$test[[f]],
                                               ntreelimit=ntree)
                  )
  )

  # Compute the predicted labels with the prediction rule of classifing any
  # observation with a score over the quantile prop_select as 1.
  lab <- lapply(1:length(prop_select),
                function(v) lapply(1:nfolds,
                                   function(f) apply(score[[f]], 2, .pred_fun,
                                                     q=prop_select[v]))
                )

  # Compute the classification error for these predictions
  error <- lapply(1:length(prop_select),
                  function(v)
                    sapply(1:nfolds,
                           function(f)
                             apply(lab[[v]][[f]], 2,
                                   function(yhat)
                                     mean(abs(yhat -
                                                xgboost::getinfo(
                                                  sets$test[[f]], "label"))))))

  # Compute the "fraud loss", the error in terms of false positives for these
  # predictions
  fraud <- lapply(1:length(prop_select),
                  function(v)
                    sapply(1:nfolds, 
                           function(f)
                             apply(lab[[v]][[f]], 2,
                                   function(yhat)
                                     sum((1-xgboost::getinfo(sets$test[[f]],
                                                             "label")) * 
                                           yhat) / sum(yhat))))

  # Also compute the auc for the predicted scores
  auc <- sapply(1:nfolds,
                function(f) apply(score[[f]],
                                  2, comp_auc,
                                  y=xgboost::getinfo(sets$test[[f]], "label")))

  return(.cmp_agg(1:nrounds, val_type, y, nfolds, folds, prop_select,
                  error, fraud, auc, lab))

}
