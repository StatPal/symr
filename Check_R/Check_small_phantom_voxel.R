###### Reading and cropping the data ##########

library(oro.nifti)
(ffd <- readNIfTI("~/R/Researchs/MRI/Read_Data/new_phantom.nii.gz"))
image(ffd)

# tmp <- ffd[124:126, 124:127, 1, ]
# tmp <- ffd[50:55, 57:60, 1, ]
tmp <- ffd[74:83, 124:133, 1, ]
# tmp <- ffd[74:83+1, 124:133+1, 1, ]
dim(tmp) <- c(dim(tmp)[1], dim(tmp)[2], 1, dim(tmp)[3])
(image1 <- nifti(tmp, datatype = datatype(ffd)))
dim(image1)
image(image1)



#### Preprocess_data ######
image2 <- array(dim=dim(image1))
for(i in 1:(dim(image1)[1])){
  for(j in 1:(dim(image1)[2])){
    image2[i,j,,] = ifelse(image1[i,j,,]==0, 0.5, image1[i,j,,])
  }
}


(n_x <- dim(image2)[1])
(n_y <- dim(image2)[2])
(n_z <- dim(image2)[3])
(M <- dim(image2)[4])
n <- n_x*n_y*n_z

TE_example <- c(0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10)
TR_example <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)
sigma_example <- c(20.6334, 18.3911, 20.5477, 17.6867, 17.3431, 14.9071, 23.2545, 
                   21.6556, 23.1835, 21.3345, 23.8697, 18.8937, 25.8299, 24.1174, 
                   26.5838, 23.8338, 25.7627, 20.9439)

(TE_scale = 2.01/min(TE_example))
(TR_scale = 2.01/min(TR_example))
(TE_example = TE_example * TE_scale)
(TR_example = TR_example * TR_scale)

lb <- c(0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale)))
ub <- c(450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale)))

(W_1_init = exp(-1/(2.0*TR_scale)))
(W_2_init = exp(-1/(0.1*TE_scale)))

(lb <- ifelse(lb < 1e-10, 1e-10, lb))
lb; ub



r <- array(dim = c(n_x*n_y*n_z, length(TE_example)))  # n x m matrix (Check the order)
for(i in 1:length(TE_example)){
  r[,i] <- as.vector(image2[,,,i])
}


## normalize:
(normalizer <- max(r))
normalizer <- 10
r <- r/normalizer
(sigma_example <- sigma_example/normalizer)
# ub[1] <- ub[1]/normalizer





###### Indices:

train_ind <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
test_ind <- c(17, 18)

train_ind <- c(1, 7, 14)
test_ind <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)

train <- r[,train_ind]
test <- r[,test_ind]
TE_train <- TE_example[train_ind]; TE_test <- TE_example[test_ind] 
TR_train <- TR_example[train_ind]; TR_test <- TR_example[test_ind]
sigma_train <- sigma_example[train_ind]; sigma_test <- sigma_example[test_ind]






### Make the initial W matrix:

mean_img <- rowMeans(train)
W_init <- array(dim=c(length(mean_img),3))
W_init[,1] <- ifelse(mean_img>450/normalizer, 425/normalizer, mean_img)
W_init[,2] = W_1_init
W_init[,3] = W_2_init
dim(W_init)




# source("scheme_new_numeric.R")
source("~/R/Researchs/MRI/Check_R_cleaned/scheme_new_numeric.R")








##### LS #####

W_LS <- least_sq_solve(W_init, TE_train, TR_train, 
                       train, TE_scale, TR_scale, if_full = T)
head(W_LS)


## C and R gives diff results: 
# R stops if diff in growth is very small (even though grad is quite high). 
# But C continues and stops when grad is ~0
# numDeriv::grad(temp_lsq_star, c(1.7411657301, 0.0000000001, 0.9280897239), side = c(NA,1,1), TE = TE_train, TR = TR_train, i = 1)
# numDeriv::grad(temp_lsq_star, c(404.043, 1e-08, 0.054274), side = c(NA,1,1), TE = TE_train, TR = TR_train, i = 1)





## W_LS_2 is the modified version of W_LS
W_LS_2 <- W_LS
W_LS_2[,1] <- (W_LS[,1]+0.001)/2
W_LS_2[,3] <- (W_LS[,3]+ub[3])/2







#### MRF check: ########
update_n_reverse(n_x, n_y, n_z)



####### OLS functions #################
Q_OSL_per_voxel <- function(W_i, TE, TR, sigma, r_i, W_old_i, c_i, i){

  v_i = Bloch_vec(W_i, TE, TR)
  v_old_i = Bloch_vec(W_old_i, TE, TR)
  m = length(TE); likeli_sum = 0
  
  ## Rice part
  ## Q_i = \sum_j ( -v_ij^2/2 + r_ij * v_ij besselI1_I0(r_ij*v_ij_old/sigma_j^2) )/sigma_j^2 
  for(j in 1:m) {
    tmp = besselI1_I0(r_i[j]*v_old_i[j]/SQ(sigma[j]))
    likeli_sum = likeli_sum + v_i[j]*(- 0.5*v_i[j] + r_i[j]*tmp)/SQ(sigma[j])
  }
  
  ## MRF part: (c_i \cdot W_i would be there)
  for(k in 1:3){
    likeli_sum = likeli_sum - c_i[k] * W_i[k]
  }
  return (-likeli_sum)
}



## Negative Gradient of Penalised Q function per voxel - grad of J evaluated at old parameter ##
Q_OSL_grad_per_voxel <- function(W_i, TE, TR, sigma, r_i, W_old_i, c_i, i){
  
  v_i = Bloch_vec(W_i, TE, TR)
  v_old_i = Bloch_vec(W_old_i, TE, TR)
  W_grad = W_i; m = length(TE)
  
  ## Likelihood part
  for(k in 1:3){
    temp_sum = 0
    for(j in 1:m){
      ## dQ_i/dW_ik = \Sum_j dQ_i/dv_ij dv_ij/dW_ik 
      ## dQ_i/dv_ij = -v_ij/sigma_j^2 + r_ij besselI1_I0(...)/sigma_j^2
      
      tmp = - v_i[j]/SQ(sigma[j]) + besselI1_I0(v_old_i[j]*r_i[j]/SQ(sigma[j])) * r_i[j]/SQ(sigma[j])
      temp_sum = temp_sum + tmp * simple_dee_v_ij_dee_W_ik(W_i, TE, TR, j, k)
    }
    W_grad[k] = temp_sum - c_i[k]
  }
  return (-W_grad)
}


## Negative penalised log likelihood
neg_log_l <- function(W, TE, TR, sigma, r, Psi_inv, beta1, just_MLE = F){
  v = v_mat(W, TE, TR)
  m = length(TE)
  n = nrow(W)
  likeli_sum <- 0
  
  ## Rice part
  for(i in 1:n){
    for(j in 1:m) {
      tmp = (r[i, j]/SQ(sigma[j]))
      # likeli_sum = likeli_sum + log(tmp) -
      #   (r_i[j]^2 + v_i[j]^2)/(2*(sigma[j]^2)) + logBesselI0(tmp * v_i[j])
      likeli_sum = likeli_sum + 
        log(r[i,j]/SQ(sigma[j])) - 
        ((r[i,j]^2+v[i,j]^2)/(2*sigma[j]^2)) + 
        logBesselI0(r[i,j]*v[i,j]/SQ(sigma[j]))
    }
  }
  
  ## MRF part: 
  if(!just_MLE){
    likeli_sum = likeli_sum + MRF_log_likeli(W, Psi_inv, beta1)
    ## Check sign again 
  }
  
  return (-likeli_sum)
}






## Q _ MRF :
param <- c(0.5, 0.1, 0.9, 0.1, 0.2, 0.3, 0.4)
# param <- c(0.52, 0.15, 0.85, 0.12, 0.17, 0.27, 0.36)
# lb_MRF <- c(1e-5, -Inf, -Inf, 1e-5, -Inf, 1e-5, 1e-5)
# ub_MRF <- rep(Inf, 7)
# lb_MRF <- c(1e-2, -Inf, -Inf, 1e-2, -Inf, 1e-2, 1e-2)*2
# ub_MRF <- rep(Inf, 7)
lb_MRF <- c(1e-5, -255, -255, 1e-5, -255, 1e-5, 1e-5)
ub_MRF <- rep(255, 7)






### Basic MLE/MPLE OSL function:  #####
OSL_optim <- function(W1, Psi_inv, beta1, TE1, TR1, sigma1, r1,  
                      TE_scale, TR_scale, maxiter = 10, 
                      just_MLE = F, abs_diff = 1e-6, rel_diff = 1e-3, 
                      verbose = 0, if_full = T){
  
  ### MRF based initial values ###
  L <- t(chol(Psi_inv))
  x_MRF <- array(dim=7)
  x_MRF[1:6] <- from_L_mat(L); x_MRF[7] <- beta1[1]
  
  old_val <- Q_star_other_parameters(x_MRF, W = W1)
  optim_obj_MRF <- optim(x_MRF, Q_star_other_parameters, Q_grad_vec_other_parameter, W = W1, 
                         method="L-BFGS-B", lower=lb_MRF, upper=ub_MRF,
                         control = list(maxit = 10000))
  
  temp_L <- optim_obj_MRF$par[1:6]
  L_mat <- to_L_mat(temp_L)
  Psi_inv <- L_mat %*% t(L_mat)

  if(n_z>1){
    beta1[1] <- optim_obj_MRF$par[7]; beta1[2] <- optim_obj_MRF$par[8]; beta1[3] <- 0.1
  } else {
    beta1[1] <- optim_obj_MRF$par[7]; beta1[2] <- 1; beta1[3] <- 0
  }  

  new_val <- Q_star_other_parameters(optim_obj_MRF$par, W = W1)
  cat("\nMRF log likeli new:", -new_val, "\n")
  
  
  
  
  
  
  
  ### Voxel based optimization ###
  iter <- 0
  W_old <- W1
  old_log_likeli <- neg_log_l(W1, TE1, TR1, sigma1, r1, Psi_inv, beta1, just_MLE)
  loop_vec <- 1:n
  
  
  while(iter < maxiter){
    print(paste0(rep("-", 40), collapse="")); print(paste("Iteration", iter))
    iter = iter + 1
    
    Gamma_inv = Lambda(beta1)
    MRF_grad = Gamma_inv %*% W_old %*% Psi_inv
    ## Check +- signs
    
    
    for(i1 in loop_vec){
      if(i1 %% 1000 == 1) cat("i: ", i1, "\n")
      
      c_i <- MRF_grad[i1,]
      if(just_MLE){
        c_i <- c(0,0,0)  ## Corresponds to MLE
      }
      x <- W1[i1,]
      
      
      old_val <- Q_OSL_per_voxel(W1[i1,], TE = TE1, TR = TR1, sigma = sigma1, 
                                 r_i = r1[i1,], W_old_i = W_old[i1,], c_i = c_i, i = i1)
      optim_obj <- optim(W1[i1,], Q_OSL_per_voxel, Q_OSL_grad_per_voxel, method="L-BFGS-B",
                         lower=lb, upper=ub, TE = TE1, TR = TR1, sigma = sigma1,
                         r_i = r1[i1,], W_old_i = W_old[i1,], c_i = c_i, i = i1,
                         control = list(maxit = 10000))
      # optim_obj <- optim(W_LS_2[i1,], Q_OSL_per_voxel, Q_OSL_grad_per_voxel, method="L-BFGS-B",
      #                    lower=lb, upper=ub, TE = TE1, TR = TR1, sigma = sigma1, 
      #                    r_i = r1[i1,], W_old_i = W_old[i1,], c_i = c_i, i = i1,
      #                    control = list(maxit = 10000)) 
      x <- optim_obj$par; fx <- optim_obj$value
      
      if(fx >= old_val){ 
        print(optim_obj$message)
        print(fx)
        print(old_val)
        print(i1)
      } else {
        W1[i1,] <- x
      }
    }
    
    cat("Abs diff(W):", mean(abs(W_old - W1)), "\n")
    new_log_likeli <- neg_log_l(W1, TE1, TR1, sigma1, r1, 
                                Psi_inv, beta1, just_MLE)
    cat("new_log_likeli: ", new_log_likeli, ", old_log_likeli: ", old_log_likeli, "\n")
    cat("rel_diff(lld): " , abs(new_log_likeli - old_log_likeli)/abs(new_log_likeli), "\n")
    if(abs(new_log_likeli - old_log_likeli)/abs(new_log_likeli) < rel_diff){
      break
    }
    W_old <- W1
    old_log_likeli <- new_log_likeli
  }
  
  print("OSL ends!")
  print(paste0(rep("-", 40), collapse=""))
  
  return(list(W = W1, Psi_inv = Psi_inv, beta = beta1))
}

  

Psi_inv <- Diagonal(3)
beta1 <- c(0.1, 1, 0)
MRF_log_likeli(W_LS, Psi_inv, beta1)


## Change W_LS to W_LS_2 for the changed W
W_MLE_obj <- OSL_optim(W_LS, Psi_inv, beta1, TE_train, TR_train, sigma_train, train,  
                         TE_scale, TR_scale, maxiter = 100, just_MLE = T,
                         abs_diff = 1e-2, rel_diff = 1e-4,
                         verbose = 0)
head(W_MLE_obj$W)
sum(is.nan(W_MLE_obj$W))


## MPLE, but MRF parameters are not updated.
# W_final_obj <- OSL_optim(W_LS, Psi_inv, beta1, TE_train, TR_train, sigma_train, train,  
#                          TE_scale, TR_scale, maxiter = 100, just_MLE = F,
#                          rel_diff = 1e-4,
#                          verbose = 0)
# 
# head(W_final_obj$W)






### i = 3 has some problem: 4.500000e+02 1.000906e-10 0.03423947
i_bad <- 3
train[i_bad,]

W1_i <- W_LS[i_bad,]
W_old_i <- W1_i
iter <- 1
maxiter <- 20
while(iter < maxiter){
  print(paste0(rep("-", 40), collapse="")); print(paste("Iteration", iter))
  iter = iter + 1
  c_i <- c(0,0,0)  ## Corresponds to MLE
  
  i1 <- i_bad
  x <- W1_i
  old_val <- Q_OSL_per_voxel(W1_i, TE = TE_train, TR = TR_train, sigma = sigma_train, 
                             r_i = train[i_bad,], W_old_i = W_old_i, c_i = c_i, i = i1)
  optim_obj <- optim(W1_i, Q_OSL_per_voxel, Q_OSL_grad_per_voxel, method="L-BFGS-B",
                     lower=lb, upper=ub, TE = TE_train, TR = TR_train, sigma = sigma_train,
                     r_i = train[i_bad,], W_old_i = W_old_i, c_i = c_i, i = i1,
                     hessian = T)
  x <- optim_obj$par; fx <- optim_obj$value
  
  
  cat("Grad: ")
  print(Q_OSL_grad_per_voxel(x, TE = TE_train, TR = TR_train, sigma = sigma_train, 
                    r_i = train[i_bad,], W_old_i = W_old_i, c_i = c_i, i = i1))
  print(optim_obj$hessian)
  print(optim_obj$message)
  print(fx)
  print(old_val)
  print(i1)
  W1_i <- x
  W_old_i <- W1_i
  cat(W1_i)
}

W1_i
Bloch_vec(W_LS[i_bad,], TE_train, TR_train)
Bloch_vec(W1_i, TE_train, TR_train)











####### All matrices and perfromance measures:###########


head(W_init)
head(W_LS)
head(W_MLE_obj$W)
head(W_final_obj_2$W)

head(round(W_init, 4))
head(round(W_LS, 4))
head(round(W_MLE_obj$W, 4))
head(round(W_final_obj_2$W, 4))


## test - Subrata: Performance matrix - for the mean image
Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 1)
Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 1)
Performance_test(W_init, test, TE_test, TR_test, sigma_test, 1, 2)
Performance_test(W_init, test, TE_test, TR_test, sigma_test, 3, 2)

Performance_test(W_LS, test, TE_test, TR_test, sigma_test, 1, 1)
Performance_test(W_LS, test, TE_test, TR_test, sigma_test, 3, 1)
Performance_test(W_LS, test, TE_test, TR_test, sigma_test, 1, 2)
Performance_test(W_LS, test, TE_test, TR_test, sigma_test, 3, 2)


Performance_test(W_MLE_obj$W, test, TE_test, TR_test, sigma_test, 1, 1)
Performance_test(W_MLE_obj$W, test, TE_test, TR_test, sigma_test, 3, 1)
Performance_test(W_MLE_obj$W, test, TE_test, TR_test, sigma_test, 1, 2)
Performance_test(W_MLE_obj$W, test, TE_test, TR_test, sigma_test, 3, 2)


Performance_test(W_final_obj_2$W, test, TE_test, TR_test, sigma_test, 1, 1)
Performance_test(W_final_obj_2$W, test, TE_test, TR_test, sigma_test, 3, 1)
Performance_test(W_final_obj_2$W, test, TE_test, TR_test, sigma_test, 1, 2)
Performance_test(W_final_obj_2$W, test, TE_test, TR_test, sigma_test, 3, 2)


# summary(W_LS - W_MLE_obj$W)
# summary(W_LS - W_final_obj_2$W)





####### CHECK THE IMAGES #########


### Check the W's
tmp_img <- W_LS[,1]
dim(tmp_img) <- c(10,10)
tmp_img_2 <- test[,1]
dim(tmp_img_2) <- c(10,10)
tmp_img_3 <- W_MLE_obj$W[,1]
dim(tmp_img_3) <- c(10,10)
tmp_img_4 <- W_final_obj_2$W[,1]
dim(tmp_img_4) <- c(10,10)

library(fields)
image.plot(tmp_img_2)
image.plot(tmp_img)
image.plot(tmp_img_3)
image.plot(tmp_img_4)

## Differences
image.plot(tmp_img - tmp_img_3, zlim = c(-20, 10))
image.plot(tmp_img - tmp_img_4, zlim = c(-20, 10))





### Check the actual test vs estimated v_ij's
v_LS <- v_mat(W_LS, TE_test, TR_test)
v_MLE <- v_mat(W_MLE_obj$W, TE_test, TR_test)
v_MPLE <- v_mat(W_final_obj_2$W, TE_test, TR_test)



tmp_img <- test[,1]
tmp_img_2 <- v_LS[,1]
tmp_img_3 <- v_MLE[,1]
tmp_img_4 <- v_MPLE[,1]

dim(tmp_img) <- c(10,10)
dim(tmp_img_2) <- c(10,10)
dim(tmp_img_3) <- c(10,10)
dim(tmp_img_4) <- c(10,10)


image.plot(tmp_img)
image.plot(tmp_img_2)
image.plot(tmp_img_3)
image.plot(tmp_img_4)


## difference
image.plot(tmp_img - tmp_img_2, zlim = c(-10, 10))
image.plot(tmp_img - tmp_img_3, zlim = c(-10, 10))
image.plot(tmp_img - tmp_img_4, zlim = c(-10, 10))

