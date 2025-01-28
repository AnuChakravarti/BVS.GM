library(pbmcapply)
source("Helper Functions.R")

########################### EM Algorithm #######################################

# MODEL:

# r_row follows Bern(rho)
# r_tht follow Bern(eta_tht)
# Theta_ij|r_theta,r_row = 1 follows Spike and slab v0 v1
# Theta_ij|r_theta,r_row = 0 follows Spike v_0
# r_lam follows Bern(eta_lam)
# Lambda_ij|r_lam (off diagonal) follows Spike and slab v0 v1
# Lambda (diagonal) follows Unif(a,b)

# Estimate of r_lam - plam (also posterior prob)
# Estimate of r_tht - ptht (also posterior prob)

# ------------ Writing the main function for the EM Algorithm ------------------

# PARAMETERS:
# X - The covariate matrix
# y - The response matrix
# Prior parameters for Theta - v_0,v_1,eta_tht,rho
# Prior parameters for Lambda - v_0,v_1,eta_lam
# Armijo's rule parameter - sig
# Iteration cutoffs - miniter, maxiter
# Lambda_init - Initialiser for Lambda
# Theta_init - Initialiser for Theta
# verbose - Whether iteration number should be printed

BVS.GM = function(X,y,v_0 = NULL,v_1= NULL,eta_tht = 0.5,
                  eta_lam = 0.5,rho = 0.8,sig = 0.1,maxiter=500,
                  miniter = 0, Lambda_init = NULL,Theta_init = NULL,
                  verbose = F){
  
  N = nrow(X)
  q = ncol(X)
  p = ncol(y)
  
  Z =cbind(X,y)
  
  S = cov(Z)
  S_xx = S[1:q,1:q]
  S_xy= S[1:q,(q+1):(p+q)]
  S_yy = S[(q+1):(p+q),(q+1):(p+q)]
  
  #Setting default values for v_0 and v_1 if nothing specified
  if(is.null(v_0)) {
    v_0 = 0.5*sqrt(1/(log(max(p,q))*N))
  } else {
    v_0 = v_0
  }
  
  if(is.null(v_1)) {
    v_1 = 10*sqrt(1/(log(max(p,q))*N))
  } else {
    v_1 = v_1
  }
  
  #Initialising the values: Theta = 0, Lambda = I if nothing specified
  if(is.null(Lambda_init)) {
    Lambda = diag(1,nrow=p)
  } else {
    Lambda = Lambda_init
  }
  
  if(is.null(Theta_init)) {
    Theta = matrix(0,nrow=q,ncol=p)
  } else {
    Theta = Theta_init
  }
  
  # Values at iteration 0 (before algoritm updates)
  Theta0=10
  Lambda0=1000
  
  # iteration counter 
  iter=0
  while((max(abs(Theta-Theta0))>0.0001|(max(abs(Lambda-Lambda0))>0.0001)|iter<miniter)&(iter<=maxiter)){
    iter = iter + 1
    
    # Defining some matrices for ease of computation 
    # We use chol2inv(chol(Lambda)) to get lambda_inv instead of solve
    Lambda_inv=chol2inv(chol(Lambda))
    A = Lambda_inv %*% t(Theta) %*% S_xx %*% Theta %*% Lambda_inv
    B = Lambda_inv + 2*A
    C = Lambda_inv %*% t(Theta) %*% S_xx
    
    # Defining the active set for the next iteration:
    if(iter == 1){
      act_lambda = matrix(TRUE,p,p)
      act_theta = matrix(TRUE,q,p)
    }
    if(iter != 1){
      del_lambda = abs(N/2*(S_yy - Lambda_inv - A))
      act_lambda = (del_lambda > abs(del_pen(plam,Lambda,v_0,v_1)))|(Lambda!=0)
      
      del_theta = abs(N/2*(2*S_xy +2*t(C)))
      act_theta = (del_theta > abs(del_pen(ptht,Theta,v_0,v_1)))|(Theta!=0)
    }
    
    #### E-Step 
    
    #For lambda
    plam = eta_1(Lambda,eta_lam,v_1,v_0)
    diag(plam) = 1
    
    #For theta
    eta2_ans = eta_2(Theta,rho,eta_tht,v_1,v_0)
    ptht = eta_1(Theta,eta_tht,v_1,v_0)*matrix(rep(eta2_ans,p),q,p)
    
    #Values at time 0 for convergence of M Step
    Lambda0=Lambda
    Theta0=Theta
    Theta2=1000
    Lambda2=1000
    
    # M-Step
    count = 0
    while((max(abs(Lambda-Lambda2))>0.0001|max(abs(Theta-Theta2))>0.0001)&count<=200){
      Lambda2=Lambda
      Theta2=Theta
      count=count+1
      
      # V = DelTheta*Lambda^-1, U = Lambda^-1*DelLambda
      U=matrix(0,p,p)
      V=matrix(0,q,p)
      D_y=matrix(0,p,p)
      D_xy=matrix(0,q,p)
      
      # Update off diagonal of Lambda
      for(i in 1:(p-1)){
        for(j in (i+1):p){
          if(act_lambda[i,j]){
            a = N/2*(Lambda_inv[i,j]^2 + (Lambda_inv[i,i]*Lambda_inv[j,j]) +
                       (Lambda_inv[i,i]*A[j,j])+ (2*Lambda_inv[i,j]*A[i,j])+
                       (Lambda_inv[j,j]*A[i,i]))
            b = -N*(-S_yy[i,j]+Lambda_inv[i,j]+A[i,j]+(C[i,]%*%V[,j])
                    +(C[j,]%*%V[,i])
                    - (U[i,]%*%Lambda_inv[,j])
                    -(U[i,]%*%A[,j])
                    - (U[j,]%*%A[,i]))
            c = Lambda[i,j] + D_y[i,j]
            tau = (plam[i,j]/v_1)+ ((1-plam[i,j])/v_0)
            
            delta_lam=soft_thresh(a,b,c,tau) 
            D_y[i,j]=D_y[i,j]+delta_lam
            D_y[j,i]=D_y[j,i]+delta_lam
            U[,j]=U[,j]+Lambda_inv[,i]*delta_lam
            U[,i]=U[,i]+Lambda_inv[,j]*delta_lam
          }
          else{
            U[,j]=U[,j]-Lambda_inv[,i]*D_y[i,j]
            U[,i]=U[,i]-Lambda_inv[,j]*D_y[j,i]
            D_y[i,j]=0
            D_y[j,i]=0
          }
        }
      }
      
      #Update diagonal entries of Lambda
      for(i in 1:p){
        a = N*(Lambda_inv[i,i]*B[i,i])
        b = -N*(-S_yy[i,i]+Lambda_inv[i,i]+A[i,i]+(2*C[i,]%*%V[,i])
                - (U[i,]%*%Lambda_inv[,i])-(2*U[i,]%*%A[,i]))
        c = Lambda[i,i] + D_y[i,i]
        tau = 0
        delta_diag=soft_thresh(a,b,c,tau)
        U[,i]=U[,i]+Lambda_inv[,i]*delta_diag
        D_y[i,i]=D_y[i,i]+delta_diag
      }
      
      #Update entries of Theta
      for(i in 1:q){
        for(j in 1:p){
          if(act_theta[i,j]){
            a = N*Lambda_inv[j,j]*S_xx[i,i]
            b = N*(S_xy[i,j]+C[j,i]
                   +(S_xx[i,]%*%V[,j])-(U[j,]%*%C[,i])) 
            c = Theta[i,j] + D_xy[i,j]
            tau = (ptht[i,j]/v_1)+ ((1-ptht[i,j])/v_0)
            delta_tht=soft_thresh(a,b,c,tau)
            D_xy[i,j]=D_xy[i,j]+delta_tht
            V[i,]=V[i,]+Lambda_inv[j,]*delta_tht
          }
          else{
            V[i,]=V[i,]-Lambda_inv[j,]*D_xy[i,j]
            D_xy[i,j]=0
          }
        }
      }
      
      # For finding alpha we check 3 conditions:
      # (1) Positive-Definiteness of Lambda Matrix
      # (2) Armijos rule
      # (3) Bound of resultant Lambda Matrix
      
      # Note: lik = -f
      
      alpha = 1
      likhood_0 = -lik(Lambda0,Theta0,S_yy,S_xy,S_xx,N)
      
      # Checking that Lambda update is PD
      while(min(eigen(Lambda0+alpha*D_y)$value)<0){
        alpha=alpha/2
      }
      
      # Checking Armijo's rule
      likhood=-lik(Lambda0+alpha*D_y,Theta0+alpha*D_xy,S_yy,S_xy,S_xx,N)
      inv_lambda0=chol2inv(chol(Lambda0))
      
      #delta_t = del_{wrt to Lambda}*D_y + del_{wrt to Theta}*D_xy 
      # i.e. its the Armijo's rule RHS check
      
      delta_t = N/2*(trace_mat(
        t(-inv_lambda0+S_yy
          -inv_lambda0%*%t(Theta0)%*%S_xx%*%Theta0%*%inv_lambda0)%*%D_y)
        +2*trace_mat(t(S_xy+S_xx%*%Theta0%*%inv_lambda0)%*%D_xy))
      
      #Note: Only for off - diagonal upper triangular part of Lambda
      # A - diag(diag(A)) removes diagonal elements of matrix A 
      
      
      diff=(pen(plam,Lambda-diag(diag(Lambda)),v_0,v_1)-
              pen(plam,Lambda0-diag(diag(Lambda0)),v_0,v_1))/2+
        pen(ptht,Theta,v_0,v_1)-
        pen(ptht,Theta0,v_0,v_1) 
      
      extra_term = sum(((plam/
                           v_1 + (1-plam)/v_0)*sign(Lambda0-diag(diag(Lambda0)))*D_y)/2)+
        sum((ptht/v_1+(1-ptht)/v_0)*sign(Theta0)*D_xy)
      
      iter1 = 0
      
      while(iter1<20){
        #Check Armijo's Rule
        if((likhood<likhood_0+alpha*sig*delta_t-diff+alpha*sig*extra_term)&
           min(eigen(Lambda0+alpha*D_y)$value)>0){
          break
        }
        alpha=alpha*0.5
        likhood=-lik(Lambda0+alpha*D_y,Theta0+alpha*D_xy,S_yy,S_xy,S_xx,N)
        iter1=iter1+1
      }
      if(iter1<20){
        Lambda = Lambda0 + alpha*D_y
        Theta=Theta0+alpha*D_xy
      }else{
        Lambda=Lambda0
        Theta=Theta0
      }
    }
    if(verbose == T){
      print(iter)
    }
  }
  return(list(Lambda=Lambda,plam=plam,Theta=Theta,ptht=ptht))
}
