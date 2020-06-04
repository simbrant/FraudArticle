.Kfold <- function (n, K) {

  foldids <- c()
  if (n > K) {

    for (k in 1:(ceiling(n/K) - 1)) {
      foldids <- c(foldids, 1:K)
    }
    foldids <- c(foldids, sample(x=1:K, size=n - (ceiling(n/K)*K - K),
                                 replace=FALSE))
  } else if (n == K) {
    foldids <- 1:n
  }

  return(sample(x=foldids, size=n, replace=FALSE))

}

.draw_folds <- function(y, nfolds, stratified, val_type, cv_repeat){

  if (val_type=="boot") {

    if (stratified) {
      folds <- lapply(1:nfolds,
                      function(fold_no)
                        c(sample(which(y==0),
                                 size=sum(y==0),
                                 replace=T),
                          sample(which(y==1),
                                 size=sum(y==1),
                                 replace=T)))
    } else {
      folds <- lapply(1:nfolds,
                      function(fold_no)
                        sample(1:length(y),
                               size=length(y),
                               replace=T))
    }

  } else if (val_type=="cv"){

    if(stratified){
      cv_inds <- lapply(1:cv_repeat,
                        function(dummy) rep(NA, length(y)))

      for (rep in 1:cv_repeat){
        cv_inds[[rep]][y==1] <- .Kfold(sum(y==1), nfolds)
        cv_inds[[rep]][y==0] <- .Kfold(sum(y==0), nfolds)
      }

    } else{
      cv_inds <- lapply(1:cv_repeat,
                        function(dummy) .Kfold(length(y), nfolds))
    }

    folds <- lapply(1:(nfolds*cv_repeat),
                    function(r) which(cv_inds[[1 + (r-1)%/%nfolds]]
                                      != (1 + (r-1)%%nfolds)))
    nfolds <- nfolds*cv_repeat

  } else {

    stop(paste(c("val_type=", val_type,
                 " not implemented.",
                 " Specify val_type as 'cv' or 'boot'.",
                 " For repeated cv specify val_boot='cv', and",
                 " cv_repeat as an",
                 " integer that is larger than 1."),
               collapse = ""))
  }

  return(list(folds=folds, nfolds=nfolds))

}
