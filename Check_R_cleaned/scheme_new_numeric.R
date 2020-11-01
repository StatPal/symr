###########################################
## Changes: n_x, n_y,n_z removed from 
## lambda, sp_log_det_specific, sp_log_inv_specific


library(Matrix)

SQ = function(x) x^2
tr <- function(A) sum(diag(A))

## G is needed for one of the derivatives for symmetric matrix
G <- matrix( 
  c(1, 0, 0,  0, 0, 0,  0, 0, 0,
    0, 1, 0,  1, 0, 0,  0, 0, 0,
    0, 0, 1,  0, 0, 0,  1, 0, 0,
    0, 0, 0,  0, 1, 0,  0, 0, 0,
    0, 0, 0,  0, 0, 1,  0, 1, 0,
    0, 0, 0,  0, 0, 0,  0, 0, 1),
  6, 9, byrow = T)


J_n <- function(n_x){
  temp <- array(0, dim = c(n_x, n_x))
  if(n_x == 1){
    temp[1, 1] = 0.0
  } else if(n_x == 2){
    temp[1, 1] = 1
    temp[1, 2] = -1
    temp[2, 1] = -1
    temp[2, 2] = 1
  } else {
    for(i in 1:(n_x-1)){
      temp[i, i] = 2.0
      temp[i,i+1] = -1
      temp[i+1,i] = -1
    }
    temp[1,1] = 1
    temp[n_x,n_x] = 1
  }
  return(as(temp, "dgCMatrix"))
}

# eigenvals_J_n = function(n){
#   A = as.matrix(J_n(n))
#   return (eigen(A)$values)
# }

eigenvals_J_n <- function(k, if_narural_order = T){
  d_vec <- array(dim = k)
  for(i in 1:k){
    d_vec[i] <- 2*(1-cos(pi*(i-1)/k))
  }
  rev(d_vec)
}




####### MRF parameters: ############

update_n_reverse <- function(n_x, n_y, n_z){
  
  n <<- n_x*n_y*n_z
  H_1 <<- Diagonal(n_z*n_y) %x% J_n(n_x)
  H_2 <<- (Diagonal(n_z) %x% J_n(n_y)) %x% Diagonal(n_x)
  H_3 <<- J_n(n_z) %x% Diagonal(n_y*n_x)
  
  eigenval_1_small <<- eigenvals_J_n(n_x)
  eigenval_2_small <<- eigenvals_J_n(n_y)
  eigenval_3_small <<- eigenvals_J_n(n_z)
  
  one_1 <<- rep(1, n_x)
  one_2 <<- rep(1, n_y)
  one_3 <<- rep(1, n_z)
  one_23 <<- rep(1, n_y*n_z)
  one_12 <<- rep(1, n_x*n_y)
  
  eigenval_1 <<- one_23 %x% eigenval_1_small
  eigenval_2 <<- one_3 %x% (eigenval_2_small %x% one_1)
  eigenval_3 <<- eigenval_3_small %x% one_12
  
  final_vec <<- rep(0, n)
  eigens    <<- rep(0, n)
  
  ### Checks: 
  # all.equal(sort(eigen(H_1)$values), sort(as.numeric(eigenval_1)))
  # all.equal(sort(eigen(H_2)$values), sort(as.numeric(eigenval_2)))
  # all.equal(sort(eigen(H_3)$values), sort(as.numeric(eigenval_3)))
  
  return(1)
}

## Add possible check here? 



## Lambda(beta_vec) = beta_x * H_1 + beta_y * H_2 + beta_3 * H_3  
Lambda <- function(beta_vec){
  return(beta_vec[1]*H_1 + beta_vec[2]*H_2 + beta_vec[3]*H_3)
}


## |Lambda(beta1)|
sp_log_det_specific <-function(beta1, thres = 0.000001){
  eigens = beta1[1] * eigenval_1 + beta1[2]* eigenval_2 + beta1[3] * eigenval_3
  sum(log(eigens[-length(eigens)]))
}


## Needed for Derivative of MRF wrt beta_x / beta_y 
sp_log_inv_specific <- function(beta1, k, thres = 0.000001){
  
  ## Eigenvalue of Lambda
  eigens <-  beta1[1] * eigenval_1 + beta1[2] * eigenval_2 + beta1[3] * eigenval_3
  ## Check
  # all.equal(sort(eigen(Lambda(beta1))$values), sort(as.numeric(eigens)))
  
  if(k == 1) {
    final_vec = eigenval_1
  } else if(k == 2) {
    final_vec = eigenval_2
  } else if(k == 3) {
    final_vec = eigenval_3
  }
  for(i in 1:(n-1)){
    final_vec[i] = final_vec[i]/eigens[i]
  }
  temp = sum(final_vec)
  return (temp)
}


##### log likelihood due to MRF part: (Lambda = Gamma^{-1})
# (-0.5*trace(Psi_inv*W.transpose()*Lambda(beta1)*W)) 
# - 1.5*n*log(2*M_PI) + 1.5*log|Gamma^{-1}| + (n/2)*log|Psi^{-1}|
MRF_log_likeli <- function(W, Psi_inv, beta1) {
  
  likeli_sum = -0.5 * tr(Psi_inv %*% t(W) %*% Lambda(beta1) %*% W)
  
  likeli_sum = likeli_sum + 
    ( 3 * sp_log_det_specific(beta1) + 
        n * determinant(Psi_inv, logarithm = T)$modulus - 3*n*log(2*pi) )/2
  
  return (likeli_sum)
}


## gradient of likelihood w.r.t. i-th row of W. ##
MRF_grad_fn <- function( W, Psi_inv, beta1, i){
  tmp_W_Psi_inv = W %*% Psi_inv
  tmp1 = H_1[i,] %*% tmp_W_Psi_inv
  tmp2 = H_2[i,] %*% tmp_W_Psi_inv
  tmp3 = H_3[i,] %*% tmp_W_Psi_inv
  MRF_grad = beta1[1]*tmp1 + beta1[2]*tmp2 + beta1[3]*tmp3
  return (MRF_grad)
}

## W' Lambda(beta1) W
Wt_L_W <- function(W, beta1){
  Lambda_init <<- Lambda(beta1)
  tmp_Lambda_W <- Lambda_init%*%W
  return (t(W)%*%tmp_Lambda_W)
}



## Gradient of MRF likelihood w.r.t. all MRF parameters
MRF_log_likeli_grad <- function(W, Psi_inv, beta1, if_ginv = T) {
  
  Lambda_init <- Lambda(beta1)
  
  if(if_ginv){
    tmp_mat <- as.matrix(n * MASS::ginv(Psi_inv) - t(W) %*% Lambda_init %*% W)
  } else{
    tmp_mat <- as.matrix(n * solve(Psi_inv) - t(W) %*% Lambda_init %*% W)
  }
  ## Psi_inv is computaionally singular at lb_MRF
  
  Psi_grad = 0.5 * G %*% c(tmp_mat)
  
  tmp_W_Psi_inv = W %*% Psi_inv
  beta_x_grad = 1.5*sp_log_inv_specific(beta1, 1) - 0.5 * tr(t(W) %*% H_1 %*% tmp_W_Psi_inv)
  beta_y_grad = 1.5*sp_log_inv_specific(beta1, 2) - 0.5 * tr(t(W) %*%  H_2 %*% tmp_W_Psi_inv)
  
  
  grad = array(dim=8)
  grad[1:6] = Psi_grad				# 6x1 matrix I guess
  grad[7] = beta_x_grad
  grad[8] = beta_y_grad
  
  return (grad)
}





## Bloch Vector: 
# v_i = W_i0 * (1 - W_i1^TR1) * W_i2^TE1
Bloch_vec <- function(W_row, TE1, TR1){
  W_row <- as.numeric(W_row)
  W_row[1] * (1 - exp(TR1*log(W_row[2]))) * exp(TE1*log(W_row[3]))
}


## Bloch vector for whole W
v_mat <- function(W, TE, TR){
  nCol = length(TE)
  nRow = NROW(W)
  tmp = array(0, dim=c(nRow, nCol))
  for(i in 1:nRow) {
    for(j in 1:nCol) {
      tmp[i,j] = W[i,1] * exp(TE[j]*log(W[i,3])) * (1-exp(TR[j]*log(W[i,2])))
    }
  }
  return(tmp)
}



## dv_ij/dW_ik
simple_dee_v_ij_dee_W_ik = function(W, TE, TR, j, k){
  if(k == 1){
    return( exp(TE[j]*log(W[3])) * (1-exp(TR[j]*log(W[2]))) )
  } else if(k == 2){
    return( -W[1] * TR[j] * exp(TE[j]*log(W[3])) * exp((TR[j]-1)*log(W[2])) )
  } else if(k == 3){
    return( W[1] * TE[j] * exp((TE[j]-1)*log(W[3])) * (1-exp(TR[j]*log(W[2]))) )
  } else {
    return (-10000)
  }
}


## Corresponding double derivative (for variance only)
simple_dee_2_v_ij_dee_W_ik_dee_W_ik1 <- function(W, TE, TR, j, k, k1){
  
  if( k != 1 && k != 2 && k != 3){
    warning("k is not 1/2/3:")
  }
  if( k1 != 1 && k1 != 2 && k1 != 3){
    warning("k1 is not 1/2/3:")
  }
  
  if(k == 1 & k1 == 1){
    return(0)
  } else if((k == 1 & k1 == 2) | (k == 2 & k1 == 1)){
    return( -TR[j] * exp(TE[j]*log(W[3])) * exp((TR[j]-1)*log(W[2])) )
  } else if((k == 1 & k1 == 3) | (k == 3 & k1 == 1)){
    return( TE[j] * exp((TE[j]-1)*log(W[3])) * (1-exp(TR[j]*log(W[2]))) )
  } else if(k == 2 & k1 == 2){
    return( -W[1]*TR[j]*(TR[j]-1) * exp(TE[j]*log(W[3])) * exp((TR[j]-2)*log(W[2])) )
  } else if((k == 2 & k1 == 3) | (k == 3 & k1 == 2)){
    return( -W[1]*TR[j]*TE[j] * exp((TE[j]-1)*log(W[3])) * exp((TR[j]-1)*log(W[2])) )
  } else if(k == 3 & k1 == 3){
    return( W[1]*TE[j]*(TE[j]-1) * exp((TE[j]-2)*log(W[3])) * (1-exp(TR[j]*log(W[2]))) )
  }
}



## Not used  now, would be useful if one shot EM was done
to_param <- function(W, Psi_inv, beta1){
  n <- nrow(W)
  temp <- array(dim=3*n+6+2)
  temp[1:n] <- W[,1]
  temp[n+1:n] <- W[,2]
  temp[2*n+1:n] <- W[,3]
  temp[3*n+1] <- Psi_inv[1,1]
  temp[3*n+2] <- Psi_inv[1,2]
  temp[3*n+3] <- Psi_inv[1,3]
  temp[3*n+4] <- Psi_inv[2,2]
  temp[3*n+5] <- Psi_inv[2,3]
  temp[3*n+6] <- Psi_inv[3,3]
  temp[3*n+7] <- beta1[1]
  temp[3*n+8] <- beta1[2]
  return(temp)
}



### Gradient corresponding too previous vector
to_param_vec_grad <- function(W_grad, Psi_grad, beta_x_grad, beta_y_grad){
  n <- nrow(W_grad)
  temp <- array(dim=3*n+6+2)
  temp[1:n] <- W_grad[,1]
  temp[n+1:n] <- W_grad[,2]
  temp[2*n+1:n] <- W_grad[,3]
  temp[3*n+(1:6)] <- Psi_grad
  temp[3*n+7] <- beta_x_grad
  temp[3*n+8] <- beta_y_grad
  return(temp)
}



############# Init val ###############
# log(I_0(x))
logBesselI0 <- function(x, process_type = 1){
  if(process_type){
    log(besselI(x, 0, TRUE)) + x
  } else{
    log(besselI(x, 0))
  }
}

# h <- Vectorize(function(x)logBesselI0(x, process_type =1))
# curve(h, xlim=c(0,15))
# h <- Vectorize(function(x)logBesselI0(x, process_type =0))
# curve(h, xlim=c(0,15), col=2)



# I_1(x)/I_0(x)
besselI1_I0 <- function(x){
  if(x < 100000){
    return(besselI(x, 1, TRUE)/besselI(x, 0, TRUE))    
  } else {
    return(1.0 - 0.5/x)
  }
  
}

## Rice mean with parameters nu and sigma: 
mean_rice <- function(nu, sigma){
  x = - SQ(nu)/(2*SQ(sigma))
  # return (sigma * sqrt(pi/2) * exp(x/2)*( (1-x)*besselI(-x/2, 0) - x * besselI(-x/2, 1)))
  return (sigma * sqrt(pi/2) * ( (1-x)*besselI(-x/2, 0, T) - x * besselI(-x/2, 1, T)))
}


check_nan_W <- function(W, W_old){
  W[is.nan(W)] <- W_old[is.nan(W)]
}



###### Least Square #########


## ||v_i - r_i||^2
## To be minimised
temp_lsq_star <- function(W_i, TE, TR, i = 1){
  if(any(TE<1)){
    warning("TE smaller than 1")
  }
  if(any(TR<1)){
    warning("TR smaller than 1")
  }
  v_new <- Bloch_vec(W_i, TE, TR)
  
  return(sum((v_new - t(train[i,]))^2))
}


## Corresponding grad
temp_lsq_grad_vec <- function(W_i, TE, TR, i = 1) {
  if(any(TE<1)){
    warning("TE smaller than 1")
  }
  if(any(TR<1)){
    warning("TR smaller than 1")
  }
  v_new <- Bloch_vec(W_i, TE, TR)
  
  grad = array(0, dim=3)
  for(j in 1:length(TE)){
    grad[1] = grad[1] -  2*(train[i,j] - v_new[j]) * 
      (1 - exp(TR[j]* log(W_i[2]))) * 
      exp(TE[j]*log(W_i[3]))
    
    grad[2] = grad[2] -  2*(v_new[j] - train[i,j]) * 
      W_i[1]*TR[j] * exp((TR[j]-1)*log(W_i[2])) *
      exp(TE[j]*log(W_i[3]))
    
    grad[3] = grad[3] - 2*(train[i,j] - v_new[j]) * 
      (1 - exp(TR[j]*log(W_i[2]))) * 
      W_i[1]*TE[j] * exp((TE[j]-1)*log(W_i[3]))
  }
  return(grad)
}


## Checks: 

Bloch_vec(c(0.5, 0.1, 0.9), TE_train, TR_train)   # 0.5*(1-(0.1^1.01))*(0.9^1.01)
temp_lsq_star(c(0.5, 0.1, 0.9), TE_train, TR_train)
temp_lsq_grad_vec(c(0.5, 0.1, 0.9), TE_train, TR_train)
numDeriv::grad(temp_lsq_star, c(0.5, 0.1, 0.9), TE = TE_train, TR = TR_train)
# numDeriv::grad(temp_lsq_star,c(0.5, 0, 1), side=rep(1,3), TE = TE_train, TR = TR_train)
temp_lsq_grad_vec(c(0.5, 1e-100, 1), TE = TE_train, TR = TR_train)
temp_lsq_grad_vec(c(0.5, 0, 1), TE = TE_train, TR = TR_train)


# optim(c(0.5, 0.1, 0.9), temp_lsq_star, method="L-BFGS-B", 
#       lower=c(0,0,0), upper=c(255, 1, 1),  TE = TE_train, TR = TR_train)
# optim(c(0.5, 0.1, 0.9), temp_lsq_star, temp_lsq_grad_vec, method="L-BFGS-B",
#       lower=c(0,0,0), upper=c(255, 1, 1),  TE = TE_train, TR = TR_train)
# # boundary in temp gradient
# optim(c(0.5, 0.1, 0.9), temp_lsq_star, temp_lsq_grad_vec, method="L-BFGS-B", 
#       lower=c(0,0,0), upper=c(255, 1, 1), control=list(trace=6),  TE = TE_train, TR = TR_train)
# # temp_lsq_star(c(1.40556, 0.0, 0.91088))
# ## Concern: problem percists in new situation also 
# ## - problem in temp_lsq_star
# ## there was an additional minus
# 
# optim(c(0.5, 0.1, 0.9), temp_lsq_star, temp_lsq_grad_vec, method="L-BFGS-B",
#       lower=lb, upper=ub,  TE = TE_train, TR = TR_train, i = 1)






####### Least Sq solve for whole part:


# W1 <- W_init; TE1 <- TE_train; TR1 <- TR_train; r1 <- train; i1 <- 4068
# optim(W1[i1,], temp_lsq_star, temp_lsq_grad_vec, method="L-BFGS-B",
#        lower=lb, upper=ub, TE = TE1, TR = TR1, i = i1)
# Bloch_vec(c(2.205954e+01,-5.551115e-17,8.418582e-01), TE1, TR1)

least_sq_solve <- function(W1, TE1, TR1, r1, TE_scale, TR_scale, if_full = T){
  
  old_val = 0.0
  bad_count_o = 0; bad_count_o_2 = 0; bad_bound_1 = 0; bad_bound_2 = 0;
  n = nrow(r1)
  
  # Loop over voxels:
  loop_vec <- 1:n
  for(i1 in loop_vec){
    old_val <- temp_lsq_star(W1[i1,], TE = TE1, TR = TR1, i = i1)
    optim_obj <- optim(W1[i1,], temp_lsq_star, temp_lsq_grad_vec, method="L-BFGS-B",
                         lower=lb, upper=ub, TE = TE1, TR = TR1, i = i1) 
    x <- optim_obj$par; fx <- optim_obj$value
    all(x <= ub); all(x >= lb)
      
    if(fx >= old_val){
      print(paste("Value not decreased!! old x:", toString(W1[i1,]), " val: ", old_val, 
                  ";\t x: ", toString(x), " val: ", fx, " i1:", i1))
      bad_count_o = bad_count_o + 1
      if(fx > old_val){
        bad_count_o_2 = bad_count_o_2 + 1
      }
    } else {
      W1[i1,] <- x
    }
  }
  return(W1)
}







########## Other parameters: ########

## 6 dim vector to lower triangular matrix for LLt
to_L_mat <- function(temp_L){
  L <- array(0, dim = c(3, 3))
  
  L[1,1] <- temp_L[1]
  L[2,1] <- temp_L[2]; L[2,2] <- temp_L[4]
  L[3,1] <- temp_L[3]; L[3,2] <- temp_L[5]; L[3,3] <- temp_L[6]
  
  return(L)
}

## The opposite to the previous case
from_L_mat <- function(L){
  temp_L <- array(0, dim = 6)
  temp_L[1] <- L[1,1]
  temp_L[2] <- L[2,1]
  temp_L[3] <- L[3,1]
  temp_L[4] <- L[2,2]
  temp_L[5] <- L[3,2]
  temp_L[6] <- L[3,3]
  
  temp_L
}

from_Cholesky <- function(L){
  return (L %*% t(L))
}

# to_L_mat(1:6) %*% t(to_L_mat(1:6))
# t(chol(to_L_mat(1:6) %*% t(to_L_mat(1:6))))
## In R, chol be default gives Upper triangular matrix


# /*
# * Chain rule for converting gradient with Cholesky parametrization
# * Input: 6x1 vector (l = nonzero-vec(L))
# * Output: 6x6 Lower Traingular(?) matrix: d LL'/d vec(l) ?
# * -- Recheck:
# * / 2L_0	L_1		L_2		0		0		0    \
# * | 0		L_0		0		2L_1	L_2		0    |
# * | 0		0		L_0		0		L_1		2L_2 |
# * | 0		0		0		2L_3	L_4		0    |
# * | 0		0		0		0		L_3		2L_4 |
# * \ 0		0		0		0		0		2L_5 /
# * 
# * 
# * We have calculated dl_star/d vec_symm(A)
# * To caculate: dl_star/d vec_chol(L) = dl_star/d vec_symm(A) * d vec_symm(A)/ d vec_chol(L)
# * A = LL'
# * vec_symm(A) = [a_00, a_10, a_20, a_11, a_12, a_22]
# * vec_chol(L) = [l_00, l_10, l_20, l_11, l_12, l_22]				// changed due to change in symmetry??
# * The parameter is: 
# */

to_grad_Cholesky <- function(L){
  
  D = array(0, c(6, 6))
  
  D[1,1] = 2*L[1];	D[1,2] = L[2];		D[1,3] = L[3];
  D[2,2] = L[1];		D[2,4] = 2*L[2];	D[2,5] = L[3];
  D[3,3] = L[1];		D[3,5] = L[2];		D[3,6] = 2*L[3];
  D[4,4] = 2*L[4];	D[4,5] = L[5];
  D[5,5] = L[4];		D[5,6] = 2*L[5];
  D[6,6] = 2*L[6];
  
  return (D)
}

## Test - Subrata
# to_grad_Cholesky(1:6)

## Check:
# tmp <- function(temp_L, i){
#   L_mat <- to_L_mat(temp_L)
#   A = L_mat %*% t(L_mat)
#   
#   switch(i, A[1,1], A[1,2], A[1,3], A[2,2], A[2,3], A[3,3])
# }
# 
# tmp(1:6, 1)
# 
# to_grad_Cholesky(1:6)
# 
# numDeriv::grad(tmp, 1:6, i = 1)   ## 
# numDeriv::grad(tmp, 1:6, i = 2)
# numDeriv::grad(tmp, 1:6, i = 3)
# numDeriv::grad(tmp, 1:6, i = 4)
# numDeriv::grad(tmp, 1:6, i = 5)
# numDeriv::grad(tmp, 1:6, i = 6)

### Missed a transpose while multiplying the matrices ?????

tmp2 <- function(temp_L, i) {
  L_mat <- to_L_mat(temp_L)
  A = L_mat %*% t(L_mat)
  
  c(A[1,1], A[1,2], A[1,3], A[2,2], A[2,3], A[3,3])
}

## Test - Subrata
# round(numDeriv::jacobian(tmp2, 1:6) - t(to_grad_Cholesky(1:6)), 6)


update_n_reverse(n_x, n_y, n_z)





### The MRF  part wrt the other parameters(7/8 dim vector):
### 
# param contains the Psi_inv-Cholesky and beta vector
Q_star_other_parameters <- function(param, W){
  
  temp_L <- param[1:6]
  L_mat <- to_L_mat(temp_L)
  Psi_inv <- L_mat %*% t(L_mat)
  
  beta1 <- array(0, 3)
  if(n_z > 1){
    beta1[1] <- param[7]; beta1[2] <- param[8]; beta1[3] <- 0.1
  } else {
    beta1[1] <- param[7]; beta1[2] <- 1; beta1[3] <- 0.0
  }
  
  # cat("Inside Q_star_other_parameters:\n")
  # print(Psi_inv); print(beta1)
  
  return(-MRF_log_likeli(W, Psi_inv, beta1))
}


## Corresponding grad: 
Q_grad_vec_other_parameter <- function(param, W){
  
  temp_L <- param[1:6]
  L_mat <- to_L_mat(temp_L)
  Psi_inv <- L_mat %*% t(L_mat)
  
  beta1 <- array(0, 3)
  if(n_z > 1){
    beta1[1] <- param[7]; beta1[2] <- param[8]; beta1[3] <- 0.1
  } else {
    beta1[1] <- param[7]; beta1[2] <- 1; beta1[3] <- 0.0
  }
  
  grad = MRF_log_likeli_grad(W, Psi_inv, beta1)
  if(n_z==1){
    grad = grad[1:7]
  }
  
  
  chain = grad[1:6]
  # grad[1:6] = to_grad_Cholesky(temp_L) * chain
  ## Is this element wise? Check
  grad[1:6] = (to_grad_Cholesky(temp_L)) %*% chain
  ## transpose added now - no, without transpose seems fine 
  ## Posibly it was a BUG - no BUG, but check Cpp code to be sure
  
  return(-grad)
  
}

## Test - Subrata (change for 2D)
# tmp_vec <- c(1, 0.1, 0.1, 1, 0.1, 1, 0.2, 3)
# 
# Q_star_other_parameters(tmp_vec, W = W_init)
# Q_grad_vec_other_parameter(tmp_vec, W = W_init)
# numDeriv::grad(Q_star_other_parameters, tmp_vec, W = W_init)
# 
# round(Q_grad_vec_other_parameter(tmp_vec, W = W_init) -
#         numDeriv::grad(Q_star_other_parameters, tmp_vec, W = W_init), 5)
# round(numDeriv::grad(Q_star_other_parameters, tmp_vec, W = W_init)/
#         Q_grad_vec_other_parameter(tmp_vec, W = W_init), 2)









Symmetrize <- function(temp_L){
  L <- array(0, dim = c(3, 3))
  
  L[1,1] <- temp_L[1]; L[1,2] <- temp_L[2]; L[1,3] <- temp_L[3]; 
  L[2,1] <- temp_L[2]; L[2,2] <- temp_L[4]; L[2,3] <- temp_L[5];
  L[3,1] <- temp_L[3]; L[3,2] <- temp_L[5]; L[3,3] <- temp_L[6];
  
  L
}


## New tests: 
# # param contains the Psi_inv-Cholesky and beta vector
# Q_test_MRF <- function(param, W){
#   
#   Psi_inv <- Symmetrize(param[1:6])
#   beta1 <- array(0, 3)
#   beta1[1] <- param[7]; beta1[2] <- param[8]; beta1[3] <- 0.1
#   
#   MRF_log_likeli(W, Psi_inv, beta1)
# }
# 
# Q_test_MRF_grad <- function(param, W){
#   
#   Psi_inv <- Symmetrize(param[1:6])
#   beta1 <- array(0, 3)
#   beta1[1] <- param[7]; beta1[2] <- param[8]; beta1[3] <- 0.1
#   
#   grad = MRF_log_likeli_grad(W, Psi_inv, beta1)
#   
#   return(grad)
# }



## Test - Subrata
# tmp_vec <- c(1, 0.1, 0.1, 1, 0.1, 1, 0.2, 3)
# 
# Q_test_MRF(tmp_vec, W = W_init)
# Q_test_MRF_grad(tmp_vec, W = W_init)
# numDeriv::grad(Q_test_MRF, tmp_vec, W = W_init)
# 
# round(Q_test_MRF_grad(tmp_vec, W = W_init) -
#         numDeriv::grad(Q_test_MRF, tmp_vec, W = W_init), 5)
# ## There are problems here also
# ## 1, Problem was, G had not byrow = T   -- 1 BUG resolved here
# ## Psi_inv is solved
# ## But, there is still numerical instabilities in betas - check whether this is a bug or not
# round(numDeriv::grad(Q_test_MRF, tmp_vec, W = W_init)/
#         Q_test_MRF_grad(tmp_vec, W = W_init), 2)





####### Performance test ########



#  Performance matrices:
#   W: 		parameter matrix:
#   test: 	test set image matrix
#   TE_test:			
#   TR_test:			
#   sigma_test:		
#   v_type 
#     1 -- compared with v
#     2 -- compared with mode of rice distn (not implemented yet - need another set of optimization)
#     3 -- compared with rice mean
# 
#   measure_type: 
#     1 -- abs deviation from mean(not median?)
#     2 -- squared deviation from mean
# 
#   Scale: Scaled measure vs Unscaled measure 
#     0 -- no
#     1 -- yes
#   Scaling factor is sd / MAD w.r.t. mean

Performance_test <- function(W, test, TE_test, TR_test, sigma_test, 
                             v_type = 1, measure_type = 1, scale = 1, verbose = F){
  
  n_test = length(TE_test)
  stopifnot(length(sigma_test) == n_test)
  Performance_test_vec = array(0, dim=n_test)
  # Perf_mat = array(0, dim = c(nrow(W), n_test))
  
  for(i in 1:nrow(W)){
    v_new = Bloch_vec(W[i,], TE_test, TR_test)			# v_{ij}
    v_star = array(dim = n_test)
    
    
    if(v_type == 1){
      v_star = v_new
    } else if (v_type == 2){
      # Need to be done - mode of rice distn to be calculated with 
      # a seperate optimizing(maybe NR) method
    } else if (v_type == 3){
      for(j in 1:n_test){
        v_star[j] = mean_rice(v_new[j], sigma_test[j])
      }
    }
    
    tmp1 = abs(v_star - test[i,])
    
    if(measure_type == 2){
      for(j in 1:n_test){
        tmp1[j] = SQ(tmp1[j])
      }
    }
    
    
    
    
    # Perf_mat[,i] = t(v_star) - test[i,]
    
    if(verbose){
      if(i <= 20){
        cat("i: ", i, "\nW[i,]: ", W[i,], "\nv_new: ", v_new, "\n")
        cat("v_star: ", v_star, "\ntmp1: ", tmp1, "\n", test[i,], "\n")
      }
    }
    Performance_test_vec = Performance_test_vec + array(tmp1)
  }
  
  Performance_test_vec = Performance_test_vec/nrow(W)
  
  
  if(measure_type == 2){
    for(j in 1:n_test){
      Performance_test_vec[j] = sqrt(Performance_test_vec[j])
    }
  }
  
  if(verbose){
    print(Performance_test_vec)
    cat("\n\n")
  }
  
  
  if(scale){
    scale_factor = 1.0
    for(j in 1:n_test){
      if(measure_type == 1){
        scale_factor = mean(abs(test[,j]-mean(test[,j])))    # abs_dev_mean(test[,j])) 
      } else if (measure_type == 2){
        scale_factor = sqrt(var(test[,j])*((nrow(W)-1)/nrow(W)))
      }
      Performance_test_vec[j] = Performance_test_vec[j]/scale_factor
    }
  }
  
  
  return (Performance_test_vec)
}








############### Checks #################
# log_likelihood_R = function(r, v, sigma1){
#   tmp2 = r/sigma1^2
#   tmp3 = (r^2+v^2)/(sigma1^2)
#   tmp1 = logBesselI0(tmp2*v, 0)
#   return(log(tmp2) + log(tmp1) - 0.5*tmp3)
# }
# 
# Q_fn_R <- function(r, v, sigma1, v_old){
#   tmp2 = r*v_old/(sigma1^2)
#   tmp3 = logBesselI0(tmp2, 0)
#   return(v*(- 0.5*v + r*tmp3)/(sigma1^2))
# }
# 
# # plot(function(x) log_likelihood_R(r = 5.2, v=x, sigma1 = 2), xlim=c(0,10))
# best_opt <- optim(2.5, fn = function(x) {-log_likelihood_R(r = 5.2, v=x, sigma1 = 2)},
#                   method="Brent", lower = 0.0,upper=20)
# best_opt$par
# 
# vec_test <- seq(0, 20, 0.08)
# plot(vec_test, 
#      sapply(vec_test,function(x) log_likelihood_R(r = 5.2, v=x, sigma1 = 2)), type="l", ylim=c(-30,20))
# abline(v=best_opt$par, lty=1) # best value
# 
# lines(vec_test,
#       sapply(vec_test,function(x) Q_fn_R(r = 5.2, v=x, sigma1 = 2, v_old = 2)),type="l",col=2, lwd=0.5)
# abline(v=2, lty=2, col = 2)
# 
# next_opt <- optim(2.5, fn = function(x) {-Q_fn_R(r = 5.2, v=x, sigma1 = 2, v_old = 2.5)},
#                   method="Brent", lower = 0.0,upper=20)
# abline(v=next_opt$par, lty=2, col=2)
# 
# ##########################################
# vec_test <- seq(0, 20, 0.05)
# plot(vec_test, 
#      sapply(vec_test,function(x) log_likelihood_R(r = 5.2, v=x, sigma1 = 2)), type="l", ylim=c(-50,200))
# lines(vec_test,
#       sapply(vec_test,function(x) Q_fn_R(r = 5.2, v=x, sigma1 = 2, v_old = 2.5)),type="l",col=2)
# lines(vec_test,
#       sapply(vec_test,function(x) Q_fn_R(r = 5.2, v=x, sigma1 = 2, v_old = 9.308)),type="l",col=2)
# 
# ## We expect 
# # log(p(X|theta))-log(p(X|theta*)) >= Q(theta|theta*) - Q(theta*|theta*)
# 
# 
# 
# # log_likelihood_R(0.5,0.5,1)
# # log_likelihood_R(image2[150,1,1,], mean_img[150,1,1],1)
# 
# 
# log_likelihood_sum <- function(r, v, sigma){
#   v_mat_1 <- v_mat(as.vector(W), TE_example, TR_example)
#   for(i in 1:(dim(image)[1])){
#     for(j in 1:(dim(image)[2])){
#       
#     }
#   }
# }




######## Calculations for Hessian ###########

## I_2(x)/I_0(x)
besselI2_I0 <- function(x){
  return(besselI(x, 2, TRUE)/besselI(x, 0, TRUE))
}

# d(I_1(x)/I_0(x))/dx
h <- function(x){
  return(0.5*(1 + besselI2_I0(x) -2*SQ(besselI1_I0(x))))
}



## test purpose: 
## Hesian of Q function - needs to be cleaned
# Q_hess_per_voxel <- function(W1_row, TE, TR, sigma1, r_i,  W_old_i, i, 
#                              verbose = 0){
#   
#   m <- length(sigma1)
#   obs_I <- matrix(0, 3, 3)
#   v_i = Bloch_vec(W1_row, TE, TR)
#   for(k in 1:3){
#     for(k1 in 1:3){  ## make 1:k
#       
#       for(j in 1:m){
#         if(verbose){
#           print("k: ");print(k)
#           print("k1: ");print(k1)
#           print("j: ");print(j)
#         }
#         
#         tmp <- r_i[j]/SQ(sigma1[j])
#         if(verbose){
#           print(tmp)
#         }
#         
#         tmp_2 <- (-v_i[j]/SQ(sigma1[j])) + tmp * besselI1_I0(tmp*v_i[j])
#         if(verbose){
#           print(tmp_2)
#         }
#         tmp_2 <- tmp_2 * 
#           simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W1_row, TE, TR, j, k, k1)
#         if(verbose){
#           print(simple_dee_2_v_ij_dee_W_ik_dee_W_ik1(W1_row, TE, TR, j, k, k1))
#           print(tmp_2)   ## 0 when k and k1
#           # This is little different
#         }
#         
#         
#         
#         
#         tmp_1 <- (-1)/SQ(sigma1[j]) + SQ(tmp) * h(tmp*v_i[j]/SQ(sigma1[j]))
#         if(verbose){
#           print(tmp)
#           print(v_i[j])
#           print(SQ(sigma1[j]))
#           print(tmp*v_i[j]/SQ(sigma1[j]))
#           print(h(tmp*v_i[j]/SQ(sigma1[j])))  # NaN
#           print(tmp_1)  ## NaN
#           # little diff - creating the whole difference
#         }
#         tmp_1 <- tmp_1 * 
#           simple_dee_v_ij_dee_W_ik(W1_row, TE, TR, j, k) *
#           simple_dee_v_ij_dee_W_ik(W1_row, TE, TR, j, k1)
#         if(verbose){
#           print(simple_dee_v_ij_dee_W_ik(W1_row, TE, TR, j, k))
#           print(simple_dee_v_ij_dee_W_ik(W1_row, TE, TR, j, k)*
#                   simple_dee_v_ij_dee_W_ik(W1_row, TE, TR, j, k1))
#           print(tmp_1)  ## Very Different for some - NaN now!!!
#         }
#         
#         
#         obs_I[k, k1] <- obs_I[3*(i-1)+k, 3*(i-1)+k1] + tmp_1 + tmp_2
#         ## Double counting????
#         if(verbose){
#           print("added:")
#           print(tmp_1)
#           print(tmp_2)
#           print(tmp_1 + tmp_2)
#           cat("\n")
#           # print(obs_I)
#           cat("\n\n")
#         }
#       }
#     }
#   }
#   obs_I
# }











