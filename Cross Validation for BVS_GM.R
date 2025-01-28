source("EM Algorithm for BVS_GM.R")
source("Helper Functions.R")
library(pbmcapply)

# Due to the nature of our hyperparameters, we only need to optimize one of them 
# at a time, assuming the others are fixed. For this purpose, we performed 
# cross-validation for v_1. 

# We include two cross validation methods here - the first uses prediction error 
# as the measure to choose the best v_1, the second chooses likelihood value. 

############### Cross Validation for v1 ########################################

# using prediction error as measure
cv_our_v2 = function(X,Y,v1_vals){
  cv_id = sample(1:5,N,replace= T,prob = c(0.2,0.2,0.2,0.2,0.2))
  
  #Defining the function to do cross validation for a particular fold k.
  cv_fold = function(k){
    X_CV_train = X[cv_id != k,]
    X_CV_test = as.matrix(X[cv_id == k,],nrow = sum(cv_id == k),ncol = ncol(X))
    Y_CV_train = Y[cv_id != k,]
    Y_CV_test = as.matrix(Y[cv_id == k,],nrow = sum(cv_id == k),ncol = ncol(X))

    Theta_init = NULL
    Lambda_init = NULL
    error = numeric(length(v1_vals))

    for(i in 1:length(v1_vals)){
      v1 = v1_vals[i]
      ans = BVS.GM(X_CV_train,Y_CV_train,v_1 = v1,
                   Theta_init = Theta_init,
                   Lambda_init = Lambda_init)
      Theta_est = ans$Theta
      Lambda_est  = ans$Lambda
      B_est = -chol2inv(chol(Lambda_est))%*%t(Theta_est)

      Y_pred = matrix(nrow = nrow(Y_CV_test),ncol = ncol(Y_CV_test))
      for(k in 1:nrow(Y_CV_test)){
        Y_pred[k,] = B_est%*%X_CV_test[k,]
      }

      error[i] = mean(apply(Y_CV_test-Y_pred,1,l2norm))

      Theta_init = Theta_est
      Lambda_init = Lambda_est

    }
    return(error)
  }
  cv_results = pbmclapply(1:5,cv_fold,mc.cores = 5)
  cv_results = Reduce("+",cv_results)/5
  best = which(cv_results == min(cv_results))[1]

  return(list("v1_best" = v1_vals[best],"cv_results" = cv_results))
}

# ---------------------------------------------------------------------------- #
# using likelihood value as measure
cv_our_v1 = function(X,Y,v1_vals){
  cv_id = sample(1:5,N,replace= T,prob = c(0.2,0.2,0.2,0.2,0.2))
  
  #Defining the function to do cross validation for a particular fold k.
  cv_fold = function(k){
    X_CV_train = X[cv_id != k,]
    X_CV_test = matrix(X[cv_id == k,],nrow = sum(cv_id == k),ncol = ncol(X))
    Y_CV_train = Y[cv_id != k,]
    Y_CV_test = matrix(Y[cv_id == k,],nrow = sum(cv_id == k),ncol = ncol(X))

    Theta_init = NULL
    Lambda_init = NULL
    likval = numeric(length(v1_vals))

    for(i in 1:length(v1_vals)){
      v1 = v1_vals[i]
      ans = BVS.GM(X_CV_train,Y_CV_train,
                   v_1 = v1,Theta_init = Theta_init,
                   Lambda_init = Lambda_init)
      Theta_est = ans$Theta
      Lambda_est  = ans$Lambda
      B_est = -chol2inv(chol(Lambda_est))%*%t(Theta_est)
      pred = B_est
      
      # Defining the covariance matrices required to compute likelihood
      Z_CV_test =cbind(X_CV_test,Y_CV_test)
      
      S_CV_test = cov(Z_CV_test)
      S_xx_CV_test = S_CV_test[1:q,1:q]
      S_xy_CV_test= S_CV_test[1:q,(q+1):(p+q)]
      S_yy_CV_test = S_CV_test[(q+1):(p+q),(q+1):(p+q)]

      likval[i] = lik(Lambda_est,Theta_est,S_yy_CV_test,S_xy_CV_test,
                      S_xx_CV_test,nrow(X_CV_test))

      Theta_init = Theta_est
      Lambda_init = Lambda_est

    }
    return(likval)
  }
  cv_results = pbmclapply(1:5,cv_fold,mc.cores = 5)
  cv_results = Reduce("+",cv_results)/5
  best = which(cv_results == max(cv_results))[1]

  return(list("v1_best" = v1_vals[best],"cv_results" = cv_results))
}
