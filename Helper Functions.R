library(caret)
library(MASS)
library(expint)

#---------------------- Creating helper functions -------------------------------
# Data generating function for Lambda
gen_Lambda = function(p,num_non0)
{
  # Making 5 non-zero observations in off diagonal Lambda
  # Randomly generate the elements in the lower triangular Lambda which will be 
  # non-sparse
  
  rows = sample(2:p,num_non0, replace = T)
  cols= rows
  for ( i in 1:length(rows)){
    cols[i] = sample((rows[i]-1),1) 
  }
  
  # Create a symmetric Lambda
  Lambda = matrix(0,p,p)
  Lambda[(cols-1)*p + rows] = runif(num_non0,2,3)*sample(c(-1,1),num_non0,replace = T)
  # Making the signal weaker
  #Lambda[(cols-1)*p + rows] = runif(num_non0,0.1,0.2)*sample(c(-1,1),num_non0,replace = T)
  Lambda =  (Lambda+t(Lambda))/2
  
  # Randomly generate values for off-diagonal of Lambda
  Lambda[0:(p-1)*p + 1:p] = runif(p,2,3)*sample(c(-1,1),p,replace = T)
  # Making the signal weaker
  #Lambda[0:(p-1)*p + 1:p] = runif(p,0.1,0.2)*sample(c(-1,1),p,replace = T)
  diag(Lambda) = apply(Lambda,1,function(x)sum(abs(x)))+0.2
  
  return(Lambda)
}

# Data generating function for Theta
gen_Theta = function(q,p,num_non0){
  ### Making 5 non-zero values in Theta 
  # only consider 1 -prow_act rows for putting random non sparse values in 
  prow_act = 0.7
  Theta = matrix(0,q,p)
  non0_rows = rbinom(q,1,1-prow_act)
  dummy = matrix(0,sum(non0_rows),p)
  
   # dummy[sample(sum(non0_rows)*p,num_non0)] = runif(num_non0,4,6)*
   #   sample(c(-1,1),num_non0,replace = T)
  # Making the signal weaker
  #dummy[sample(sum(non0_rows)*p,num_non0)] = runif(num_non0,0.1,0.2)*
  #  sample(c(-1,1),num_non0,replace = T)
  # if condition for the case when we want to make too many non-zero, 
  # max non-zero possible is p*(0.3 prob selection of row) - should be around 150 
  # but can be less
  if(sum(non0_rows)*p < num_non0){
    dummy[sample(sum(non0_rows)*p,num_non0,replace= T)] = runif(num_non0,4,6)*
      sample(c(-1,1),num_non0,replace = T)
  }
  else{
    dummy[sample(sum(non0_rows)*p,num_non0)] = runif(num_non0,4,6)*
      sample(c(-1,1),num_non0,replace = T)
  }
  Theta[as.logical(non0_rows),] = dummy
  return(Theta)
}

# Function finds trace of matrix
trace_mat = function(mat){
  sum(diag(mat))
}

# Function to do soft thresholding
soft_thresh = function(a,b,c,tau)
{
  d = -c + sign(c-b/a)*max(0,abs(c-b/a)-tau/a)
  return(d)
}

# Function to calculate likelihood 
lik = function(Lambda,Theta,S_yy,S_xy,S_xx,N)
{
  N/2*(log(det(Lambda)) - trace_mat(S_yy%*%Lambda+ 2*t(S_xy)%*%Theta + 
                                      solve(Lambda)%*%t(Theta)%*%S_xx%*%Theta))
}

# Function to calculate penalty for Lambda and Theta given the value of p
pen<-function(pval,Mat,v_0,v_1){
  sum((pval/v_1+(1-pval)/v_0)*abs(Mat))
}

# Function to calculate derivative of penalty for Lambda and Theta
del_pen<-function(pval,Mat,v_0,v_1){
  pval/v_1+(1-pval)/v_0
}

# Function for Laplace density
lap = function(X,scale){
  1/(2*scale)*exp(-abs(X)/scale)
}

# Function to calculate eta_1 - same as plam
eta_1 = function(Mat,eta,v_1,v_0){
  eta*lap(Mat,v_1)/(eta*lap(Mat,v_1)+(1-eta)*lap(Mat,v_0))
}

# Function to calculate eta_2
eta_2 = function(Mat,rho,eta,v_1,v_0){
  S_1 = apply(eta*lap(Mat,v_1)+(1-eta)*lap(Mat,v_0),1,prod)
  S_2 = apply(lap(Mat,v_0),1,prod)+1.476829e-160
  rho*S_1/(rho*S_1+ (1-rho)*S_2)
}

# Function to calculate the Matthews Correlation Coefficient
MCC = function(TP,TN,FP,FN){
  if((as.numeric(TP+FP)*as.numeric(TP+FN)*as.numeric(TN+FP)*as.numeric(TN+FN))!= 0){
    ans = (TP*TN - FP*FN)/sqrt((as.numeric(TP+FP)*as.numeric(TP+FN)*as.numeric(TN+FP)*as.numeric(TN+FN)))
  }
  if((TP == 0 & FP == 0 & FN == 0 & TN != 0)|(TP !=0 & FP == 0 &FN == 0 &TN == 0)){
    ans = 1
  }
  if((TP == 0 & FP == 0 & FN != 0 & TN != 0)|(TP ==0 & FP != 0 &FN == 0 &TN != 0)|
     (TP!= 0 & FP == 0 & FN != 0& TN == 0)|(TP != 0 & FP != 0 & FN == 0 & TN == 0)){
    ans = 0
  }
  if((TP == 0 & FP == 0 & FN != 0 & TN == 0)|(TP ==0 & FP != 0 &FN == 0 &TN == 0)){
    ans = -1
  }
  return(ans)
}

# Function to calculate l2 norm
l2norm = function(x) sqrt(sum(x^2))

# Function to generate points uniformly over a n-dim sphere
sample_sphere = function(num_pts,dim,r){
  x = mvrnorm(num_pts,mu = rep(0,dim),Sigma = diag(1,dim))
  if(num_pts == 1){
    ssq = sum(x^2)
  }
  else{
    ssq = rowSums(x^2)
  }
  fr = r*(pgamma(ssq/2,shape = dim/2,scale = 1))^(1/dim)/sqrt(ssq)
  frtiled = matrix(rep(fr, dim), nrow = num_pts, byrow = FALSE)
  p = x*frtiled
  return(p)
}

# Function to generate non-sparse Theta with fixed norm
gen_Theta_fixnorm = function(q,p,norm){
  # only consider 1 -prow_act rows for putting random non sparse values in 
  prow_act = 0.7
  Theta = matrix(0,q,p)
  non0_rows = rbinom(q,1,1-prow_act)
  dummy = sample_sphere(sum(non0_rows),dim = p,norm)
  Theta[as.logical(non0_rows),] = dummy
  return(Theta)
}

# Function to generate sparse Theta with fixed norm

gen_Theta_fixnorm_sparse = function(q,p,norm){
  # only consider 1 -prow_act rows for putting random non sparse values in 
  prow_act = 0.7
  Theta = matrix(0,q,p)
  non0_rows = rbinom(q,1,1-prow_act)
  for(i in 1:length(non0_rows)){
    if(as.logical(non0_rows[i])){
      dummy = sample_sphere(1,dim = floor(p*runif(1,0.1,0.5)),norm)
      Theta[i,] = sample(c(rep(0,(p-length(dummy))),dummy))
    }
  }
  return(Theta)
}
