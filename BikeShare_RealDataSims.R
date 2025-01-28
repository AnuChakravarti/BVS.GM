library(data.table)
library(dplyr)
library(tidyverse)
library(glasso)
library(capme) # install.packages("capme_1.3.tar.gz", repos = NULL, type = "source") 
source("Cross Validation for BVS_GM.R")
source("BayesCRF Efficient.R")
source("EM Algorithm for BVS_GM.R")
source("Helper Functions.R")
source("5fold-graph-cv.R")


#Reading 2023 January-July bikeshare data
bikes = NULL
for(i in 1:3){
  path = paste("Capital_Bikeshare_data/20230"
               ,i,"-capitalbikeshare-tripdata.csv",sep = "")
  bikes = rbind(bikes,fread(path))
}

#bikes = fread("Capital_Bikeshare_data/202301-capitalbikeshare-tripdata.csv")


#Creating our data
bikes_mod = bikes %>%
  # Considering 5 stations
  dplyr::filter(start_station_name %in% c("Metro Center / 12th & G St NW",
                                   "14th St & New York Ave NW",
                                   "14th & G St NW",
                                   "13th St & New York Ave NW",
                                   "11th & F St NW")) %>% 
  # Considering only classic bikes
  dplyr::filter(rideable_type == "classic_bike") %>% 
  
  # time- Morning: < 12 am, Afternoon: 12pm to 5 pm, Evening >= 5pm 
  dplyr::mutate(date_trip = substr(started_at,1,10),
         time_trip = ifelse(substr(started_at,12,13) <12, "morning",
                            ifelse(substr(started_at,12,13) >=17,"evening",
                                   "afternoon"))) %>% 
  # Removing unnecessary variables
  dplyr::select(-c("ride_id","rideable_type","started_at",
            "ended_at","start_station_id","end_station_name",
            "end_station_id","start_lat","start_lng","end_lat",
            "end_lng","member_casual")) %>% 
  # Finding count of bikes used
  group_by(start_station_name,date_trip,time_trip) %>%
  summarise(n = n()) %>% 
  # Making the table in wide format
  pivot_wider(names_from = c(start_station_name,time_trip), values_from = n, names_sep="_") %>% 
  # Replacing NAs with 0
  replace(is.na(.), 0)

  # Creating the Y and X dataframes
  Y = as.data.frame(bikes_mod[,-1])
  X = NULL
  for(t in 1:3){
    temp = apply(Y,2,lag,n= t)
    colnames(temp) = paste(colnames(Y),".",t,sep = "")
    X = cbind(temp,X) 
}
rm(temp)


X = X[-(1:3),]
Y = Y[-(1:3),]

#Train - 80% Test - 20%
num_train = floor(dim(X)[1]*0.8)

X_train = X[1:num_train,]
Y_train = Y[1:num_train,]
X_test = X[(num_train+1):dim(X)[1],]
Y_test = Y[(num_train+1):dim(Y)[1],]


# Using the train data to estimate Lambda and Theta
# ------------------ Our Method -------------------------------------

### With CV
# v1_vals = sqrt(1/(log(max(p,q))*N))*c(10000,1000,500,100,50,10,5,3)
# v1_best = cv_our_v2(X_train,Y_train,v1_vals)
# ans = BVS.GM(X_train,Y_train,v_1 = v1_best$v1_best)

ans = BVS.GM(X_train,Y_train,v_1 = sqrt(1/(log(max(p,q))*N))*10000)

Theta_est = ans$Theta
Lambda_est = ans$Lambda
B_est = -chol2inv(chol(Lambda_est))%*%t(Theta_est)

Y_pred = matrix(nrow = nrow(Y_test),ncol = ncol(Y_test))
for(i in 1:nrow(X_test)){
  Y_pred[i,] = B_est%*%X_test[i,]
}

our_method_error = mean(apply(Y_test-Y_pred,1,l2norm))

# ------------------ GLASSO  -------------------------------------
glasso_rho = graph.cv(X_train, as.matrix(Y_train),
                      c(0.001,0.005,0.01,0.05),method ="glasso" )
res_glasso = glasso(S,glasso_rho)$wi
glasso.Lambda = res_glasso[(q+1):(p+q),(q+1):(p+q)]
glasso.Theta = res_glasso[1:q,(q+1):(p+q)]
glasso.B = -solve(glasso.Lambda)%*%t(glasso.Theta)

Y_pred = matrix(nrow = nrow(Y_test),ncol = ncol(Y_test))
for(i in 1:nrow(X_test)){
  Y_pred[i,] = glasso.B%*%X_test[i,]
}

glasso_method_error = mean(apply(Y_test-Y_pred,1,l2norm))

# ------------------ CAPME -------------------------------------
capme_rho = graph.cv(X_train, as.matrix(Y_train),
                     c(0.05,0.1,0.3,0.5,0.7,1),method ="capme" )
res_capme = capme(X_train,Y_train,capme_rho,capme_rho)
capme.B = t(res_capme$Gammalist[[1]])

Y_pred = matrix(nrow = nrow(Y_test),ncol = ncol(Y_test))
for(i in 1:nrow(X_test)){
  Y_pred[i,] = capme.B%*%X_test[i,]
}

capme_method_error =mean(apply(Y_test-Y_pred,1,l2norm))

# ------------------ BCRF -------------------------------------
ans = BayesCRF(S_xx = S_xx, S_yy = S_yy,S_xy = S_xy,
               N = N,v_0 = sqrt(1/(log(max(p,q))*N))*0.5,v_1 = sqrt(1/(log(max(p,q))*N))*3)

BCRF.Theta = ans$Theta
BCRF.Lambda = ans$Lambda

BCRF.B = -chol2inv(chol(BCRF.Lambda))%*%t(BCRF.Theta)

Y_pred = matrix(nrow = nrow(Y_test),ncol = ncol(Y_test))
for(i in 1:nrow(X_test)){
  Y_pred[i,] = BCRF.B%*%X_test[i,]
}
bcrf_method_error =mean(apply(Y_test-Y_pred,1,l2norm))

# ------------------ Results -------------------------------------

data.frame("Method" = c("Our","Our_diag","GLASSO","CAPME","BCRF"),
           "Error" = c(our_method_error,ourdiag_method_error,glasso_method_error,
                       capme_method_error,bcrf_method_error))


################################################################################
############# Predicting half entries of Y in the test set given rest known ####
################################################################################


# Half of test data Y known as well
set.seed(2910)
Y_test = as.matrix(Y_test)
Y_test_known = matrix(NA,nrow = nrow(Y_test),ncol = ncol(Y_test))
Y_test_unknown = matrix(0,nrow = nrow(Y_test),ncol = ncol(Y_test))

known_index = sample(1:length(Y_test),length(Y_test)/2)
unknown_index = setdiff(1:length(Y_test),known_index)

Y_test_known[known_index] = Y_test[known_index]
Y_test_unknown[unknown_index] = Y_test[unknown_index]


# Using the train data to estimate Lambda and Theta
# ------------------ Our Method -------------------------------------
ans = BVS.GM(X_Train,Y_train,v_1 = sqrt(1/(log(max(p,q))*N))*10000)

Theta_est = ans$Theta
Lambda_est = ans$Lambda

prec_mat = matrix(0,nrow = q+p,ncol = q+p)
prec_mat[1:q,(q+1):(q+p)] = Theta_est
prec_mat[(q+1):(q+p),1:q] = t(Theta_est)
prec_mat[(q+1):(q+p),(q+1):(q+p)] = Lambda_est

Y_pred_unknown = matrix(0,nrow = nrow(Y_test),ncol = ncol(Y_test))

for(i in 1:nrow(X_test)){
  K_row = prec_mat[c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  H_row = prec_mat[!c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  
  Y_pred_unknown[i,is.na(Y_test_known[i,])] = (-solve(K_row)%*%t(H_row))%*%c(X_test[i,],na.omit(Y_test_known[i,]))
}

result.our = mean(apply(Y_test_unknown-Y_pred_unknown,1,l2norm))

# ------------------ BCRF Method -------------------------------------


ans = BayesCRF(S_xx = S_xx, S_yy = S_yy,S_xy = S_xy,
               N = N,v_0 = sqrt(1/(log(max(p,q))*N))*0.5,v_1 = sqrt(1/(log(max(p,q))*N))*3)

BCRF.Theta = ans$Theta
BCRF.Lambda = ans$Lambda


prec_mat = matrix(0,nrow = q+p,ncol = q+p)
prec_mat[1:q,(q+1):(q+p)] = BCRF.Theta
prec_mat[(q+1):(q+p),1:q] = t(BCRF.Theta)
prec_mat[(q+1):(q+p),(q+1):(q+p)] =BCRF.Lambda

Y_pred_unknown = matrix(0,nrow = nrow(Y_test),ncol = ncol(Y_test))

for(i in 1:nrow(X_test)){
  K_row = prec_mat[c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  H_row = prec_mat[!c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  
  Y_pred_unknown[i,is.na(Y_test_known[i,])] = (-solve(K_row)%*%t(H_row))%*%c(X_test[i,],na.omit(Y_test_known[i,]))
}

result.BCRF = mean(apply(Y_test_unknown-Y_pred_unknown,1,l2norm))

# ------------------ GLASSO -------------------------------------

glasso_rho = graph.cv(X_train,Y_train, c(0.0001,0.0005,0.001,0.005,0.01),method ="glasso")
res_glasso = glasso(S,glasso_rho)$wi
glasso.Lambda = res_glasso[(q+1):(p+q),(q+1):(p+q)]
glasso.Theta = res_glasso[1:q,(q+1):(p+q)]

prec_mat = matrix(0,nrow = q+p,ncol = q+p)
prec_mat[1:q,(q+1):(q+p)] = glasso.Theta
prec_mat[(q+1):(q+p),1:q] = t(glasso.Theta)
prec_mat[(q+1):(q+p),(q+1):(q+p)] =glasso.Lambda

Y_pred_unknown = matrix(0,nrow = nrow(Y_test),ncol = ncol(Y_test))

for(i in 1:nrow(X_test)){
  K_row = prec_mat[c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  H_row = prec_mat[!c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  
  Y_pred_unknown[i,is.na(Y_test_known[i,])] = (-solve(K_row)%*%t(H_row))%*%c(X_test[i,],na.omit(Y_test_known[i,]))
}


result.GLASSO = mean(apply(Y_test_unknown-Y_pred_unknown,1,l2norm))

# ------------------ CAPME -------------------------------------
capme_rho = graph.cv(X_train, as.matrix(Y_train), c(0.05,0.1,0.3,0.5,0.7,1),method ="capme" )
res_capme = capme(X_train,Y_train,capme_rho,capme_rho)
capme.B = t(res_capme$Gammalist[[1]])
capme.Lambda = capme.Lambda = res_capme$Omegalist[[1]]
capme.Theta = t(-capme.Lambda%*%capme.B)


prec_mat = matrix(0,nrow = q+p,ncol = q+p)
prec_mat[1:q,(q+1):(q+p)] = capme.Theta
prec_mat[(q+1):(q+p),1:q] = t(capme.Theta)
prec_mat[(q+1):(q+p),(q+1):(q+p)] =capme.Lambda

Y_pred_unknown = matrix(0,nrow = nrow(Y_test),ncol = ncol(Y_test))

for(i in 1:nrow(X_test)){
  K_row = prec_mat[c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  H_row = prec_mat[!c(rep(FALSE,q),is.na(Y_test_known[i,])),c(rep(FALSE,q),is.na(Y_test_known[i,]))]
  
  Y_pred_unknown[i,is.na(Y_test_known[i,])] = (-solve(K_row)%*%t(H_row))%*%c(X_test[i,],na.omit(Y_test_known[i,]))
}

result.CAPME = mean(apply(Y_test_unknown-Y_pred_unknown,1,l2norm))

data.frame("Method" = c("our","our_diag","BCRF","GLASSO","CAPME"),
           "Result" = c(result.our,result.ourdiag,result.BCRF,result.GLASSO,result.CAPME))
























