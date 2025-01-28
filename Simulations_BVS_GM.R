## Loading the function
library(ggplot2)
library(data.table)
library(pbmcapply)
library(pROC)
library(MASS)
source("EM Algorithm for BVS_GM.R")
source("Helper Functions.R")
source("Cross Validation for BVS_GM.R")

set.seed(2609)
p = 50
q = 100

N_choices = c(100,500,1000,2000)

#### -------------------------- OUR --------------------------------------------
results_our = function(N_choice){
  N = N_choice
  
  # l2 error
  error_lambda = 0
  error_theta = 0
  error_B = 0
  
  # 0 elements in parameter correctly estimated as 0
  selcons_lambda = 0 
  selcons_theta = 0
  selcons_B = 0
  
  # non-0 elements in parameter correctly estimated as non-0
  selecActive_lambda = 0
  selecActive_theta = 0
  selecActive_B = 0
  
  # Selection consistency for parameters via MCC
  MCC_lambda = 0
  MCC_theta = 0
  MCC_B = 0
  
  # Column selection consistency via MCC
  colMCC = 0
  
  for(k in 1:10){
    #Generating random matrices for Theta and Lambda using binom
    Lambda = gen_Lambda(p,100)
    #Theta = gen_Theta(q,p,100)
    #Theta = gen_Theta_fixnorm(q,p,16)
    Theta = gen_Theta_fixnorm_sparse(q,p,4)
    B = -chol2inv(chol(Lambda))%*%t(Theta)
    
    Y = matrix(0, N, p)
    X = mvrnorm(N, rep(0,q), chol2inv(chol(toeplitz(c(1,0.3,rep(0,q-2))))))
    
    for(i in 1:N){
      Y[i,] = mvrnorm(1,B%*%X[i,],chol2inv(chol(Lambda)))  
    }
    
    #Doing a CV to choose the best hyperparameter 
    # if(k ==1){
    #   v1_vals = sqrt(1/(log(max(p,q))*N))*c(1,5,10,100,500,1000,10000)
    #   v1_best = cv_our_v2(X,Y,v1_vals)
    # }
    # ans = BVS.GM(X,Y,v_1 = v1_best$v1_best)
    
    ans = BVS.GM(X,Y)
    
    Theta_est = ans$Theta
    Lambda_est = ans$Lambda
    B_est = -chol2inv(chol(Lambda_est))%*%t(Theta_est)
    
    
    error_lambda = error_lambda + norm(Lambda_est - Lambda,"F")/10
    error_theta = error_theta+norm(Theta_est - Theta,"F")/10
    error_B = error_B+norm(B_est - B,"F")/10
    
    
    #Defining the elements of the confidence matrix
    FP_lambda = sum(!(which(Lambda == 0) %in% which(Lambda_est == 0)))
    TN_lambda = sum(Lambda == 0) - FP_lambda
    FN_lambda = sum(!(which(Lambda != 0) %in% which(Lambda_est!= 0)))
    TP_lambda = sum(Lambda != 0) - FN_lambda
    
    FP_theta = sum(!(which(Theta == 0) %in% which(Theta_est == 0)))
    TN_theta = sum(Theta == 0) - FP_theta
    FN_theta = sum(!(which(Theta != 0) %in% which(Theta_est!= 0)))
    TP_theta = sum(Theta != 0) - FN_theta
    
    FP_B = sum(!(which(B == 0) %in% which(B_est == 0)))
    TN_B = sum(B == 0) - FP_B
    FN_B = sum(!(which(B != 0) %in% which(B_est!= 0)))
    TP_B = sum(B != 0) - FN_B
    
    selcons_lambda = selcons_lambda + 
      (FP_lambda/sum(Lambda == 0))/10
    selcons_theta = selcons_theta +
      (FP_theta/sum(Theta == 0))/10
    selcons_B = selcons_B+
      (FP_B/sum(B == 0))/10
    
    selecActive_lambda = selecActive_lambda +
      (FN_lambda/sum(Lambda!= 0))/10
    selecActive_theta = selecActive_theta +
      (FN_theta/sum(Theta!= 0))/10
    selecActive_B = selecActive_B +
      (FN_B/sum(B!= 0))/10
    
    MCC_lambda = MCC_lambda + MCC(TP_lambda,TN_lambda,FP_lambda,FN_lambda)/10
    MCC_theta = MCC_theta + MCC(TP_theta,TN_theta,FP_theta,FN_theta)/10
    MCC_B = MCC_B + MCC(TP_B,TN_B,FP_B,FN_B)/10
    
    act_colind = (colSums(B)  != 0)
    pred_colind =(colSums(B_est)  != 0)
    confmat = confusionMatrix(data= factor(pred_colind),reference = factor(act_colind))
    colMCC = colMCC + 
      MCC(confmat$table[2,2],confmat$table[1,1],confmat$table[2,1],confmat$table[1,2])/10
  }
  return(list("error_lambda" = error_lambda,"error_theta"=error_theta,
              "error_B"=error_B,"selcons_lambda"=selcons_lambda,
              "selcons_theta"=selcons_theta,"selcons_B"=selcons_B,
              "selecActive_lambda" = selecActive_lambda,
              "selecActive_theta" = selecActive_theta,
              "selecActive_B" = selecActive_B,
              "MCC_lambda"=MCC_lambda,
              "MCC_theta"=MCC_theta,
              "MCC_B"=MCC_B,
              "colMCC" = colMCC))
}

results = pbmclapply(N_choices,results_our,mc.cores = 8)


error_lambda_N = unlist(sapply(results,function(t){t[1]}))
error_theta_N = unlist(sapply(results,function(t){t[2]}))
error_B_N = unlist(sapply(results,function(t){t[3]}))
selcons_lambda_N = unlist(sapply(results,function(t){t[4]}))
selcons_theta_N = unlist(sapply(results,function(t){t[5]}))
selcons_B_N = unlist(sapply(results,function(t){t[6]}))
selecActive_lambda_N = unlist(sapply(results,function(t){t[7]}))
selecActive_theta_N = unlist(sapply(results,function(t){t[8]}))
selecActive_B_N = unlist(sapply(results,function(t){t[9]}))
MCC_lambda_N = unlist(sapply(results,function(t){t[10]}))
MCC_theta_N = unlist(sapply(results,function(t){t[11]}))
MCC_B_N = unlist(sapply(results,function(t){t[12]}))
colMCC_N = unlist(sapply(results,function(t){t[13]}))

# Dataframe for error
dfE_our = data.frame("N" = N_choices,
                     "error_tht" = error_theta_N,"error_lam" = error_lambda_N,
                     "error_B" = error_B_N)
dfE_our  = melt(data = data.table(dfE_our) ,id.vars = "N")
dfE_our$method = rep("Our",nrow(dfE_our))

# Dataframe for selection consistency
dfS_our = data.frame("N" = N_choices,
                     "selcons_lambda" = selcons_lambda_N,
                     "selcons_theta" = selcons_theta_N,"selcons_B" = selcons_B_N)
dfS_our  = melt(data = data.table(dfS_our) ,id.vars = "N")

dfS_our$method = rep("Our",nrow(dfS_our))

# Dataframe for selection consistency of active set
dfA_our = data.frame("N" = N_choices,
                     "selecActive_lambda" = selecActive_lambda_N,
                     "selecActive_theta" = selecActive_theta_N,
                     "selecActive_B" = selecActive_B_N)
dfA_our  = melt(data = data.table(dfA_our) ,id.vars = "N")

dfA_our$method = rep("Our",nrow(dfA_our))

#Dataframe for MCC

dfMCC_our = data.frame("N" = N_choices,
                       "MCC_lambda" = MCC_lambda_N,
                       "MCC_theta" = MCC_theta_N,
                       "MCC_B" = MCC_B_N)
dfMCC_our  = melt(data = data.table(dfMCC_our) ,id.vars = "N")
dfMCC_our$method = rep("Our",nrow(dfMCC_our))

# Dataframe for colMCC

dfCS_our = data.frame("N" = N_choices,
                      "colMCC" = colMCC_N)
dfCS_our$method = rep("Our",nrow(dfCS_our))

#write.csv(dfE_our,"dfE_our  LN HDim.csv")
# write.csv(dfS_our,"dfS_our.csv")
# write.csv(dfA_our,"dfA_our.csv")
#write.csv(dfMCC_our,"dfMCC_our  LN HDim.csv")
#write.csv(dfCS_our,"dfCS_our.csv")




