## phantom

## Concern that I just noticed: 
## If 1st 3-4 images are taken, all TR is 1!!

TE_example <- c(0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10)
TR_example <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)
sigma_example <- c(20.6334, 18.3911, 20.5477, 17.6867, 17.3431, 14.9071, 23.2545, 
                   21.6556, 23.1835, 21.3345, 23.8697, 18.8937, 25.8299, 24.1174, 
                   26.5838, 23.8338, 25.7627, 20.9439)

(TE_scale = 2.01/min(TE_example))
(TR_scale = 2.01/min(TR_example))
(TE_example = TE_example * TE_scale)
(TR_example = TR_example * TR_scale)
(lb <- c(0.0001, exp(-1/(0.01*TR_scale)), exp(-1/(0.001*TE_scale))))
(ub <- c(450.0, exp(-1/(4.0*TR_scale)), exp(-1/(0.2*TE_scale))))
lb[2] <- 1e-8

(W_1_init = exp(-1/(2.0*TR_scale)))
(W_2_init = exp(-1/(0.1*TE_scale)))


TE_train <- TE_example[1:3]
TR_train <- TR_example[1:3]
sigma_train <- sigma_example[1:3]




SQ <- function(x){x^2}

Bloch_vec <- function(W_row, TE1, TR1){
  W_row <- as.numeric(W_row)
  W_row[1] * (1 - exp(TR1*log(W_row[2]))) * exp(TE1*log(W_row[3]))
}

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

besselI1_I0 <- function(x){
  if(x < 100000){
    return(besselI(x, 1, TRUE)/besselI(x, 0, TRUE))    
  } else {
    return(1.0 - 0.5/x)
  }
}





## Q functioin
Q_OSL_per_voxel <- function(W_i, TE, TR, sigma, r_i, W_old_i, c_i, i){
  
  v_i = Bloch_vec(W_i, TE, TR)
  v_old_i = Bloch_vec(W_old_i, TE, TR)
  print(W_i)
  # print(v_i)
  # print(v_old_i)
  m <- length(TE)
  likeli_sum <- 0
  
  ## Rice part ##
  ## Q_i = \sum_j ( -v_ij^2/2 + r_ij * v_ij besselI1_I0(r_ij*v_ij_old/sigma_j^2) )/sigma_j^2 
  for(j in 1:m) {
    tmp = besselI1_I0(r_i[j]*v_old_i[j]/SQ(sigma[j]))
    # cat("r_i[j]: ", r_i[j], "v_old_i[j]: ", v_old_i[j], 
    #     "SQ(sigma[j]): ", SQ(sigma[j]),
    #     "tmp2: ", r_i[j]*v_old_i[j]/SQ(sigma[j]), 
    #     "tmp3:", tmp, "\n")
    likeli_sum = likeli_sum + v_i[j]*(- 0.5*v_i[j] + r_i[j]*tmp)/SQ(sigma[j])
    # cat( v_i[j]*(- 0.5*v_i[j] + r_i[j]*tmp)/SQ(sigma[j]), " ", likeli_sum, "\n")
    cat(likeli_sum, "\n")
  }
  cat("-----\n")
  return (-likeli_sum)
}


## Negative Gradient of Penalised Q function per voxel - grad of J evaluated at old parameter ##
Q_OSL_grad_per_voxel <- function(W_i, TE, TR, sigma, r_i, W_old_i, c_i, i){
  m = length(TE)
  v_i = Bloch_vec(W_i, TE, TR)
  v_old_i = Bloch_vec(W_old_i, TE, TR)
  W_grad = c(0,0,0)
  
  ## Likelihood part: ##
  for(k in 1:3){
    temp_sum = 0
    for(j in 1:m){
      ## dQ_i/dW_ik = \Sum_j dQ_i/dv_ij dv_ij/dW_ik 
      ## dQ_i/dv_ij = -v_ij/sigma_j^2 + r_ij besselI1_I0(...)/sigma_j^2
      tmp = - v_i[j]/SQ(sigma[j]) + 
        besselI1_I0(v_old_i[j]*r_i[j]/SQ(sigma[j])) * r_i[j]/SQ(sigma[j])
      temp_sum = temp_sum + tmp * simple_dee_v_ij_dee_W_ik(W_i, TE, TR, j, k)
    }
    W_grad[k] = temp_sum - c_i[k]
  }
  return (-W_grad)
}



#### Specific value
W_i <- c(26.3589, 0.422502, 0.92809)

r_i <- c(7,36,12)
train_i <- r_i


Q_OSL_per_voxel(W_i = W_i, TE = TE_train, TR = TR_train, sigma = sigma_train, 
                r_i = train_i, W_old_i = W_i, c_i = c(0,0,0), i = 1)
Q_OSL_grad_per_voxel(W_i = W_i, TE = TE_train, TR = TR_train, sigma = sigma_train, 
                     r_i = train_i, W_old_i = W_i, c_i = c(0,0,0), i = 1)
numDeriv::grad(Q_OSL_per_voxel, W_i, 
               TE = TE_train, TR = TR_train, sigma = sigma_train, 
               r_i = train_i, W_old_i = W_i, c_i = c(0,0,0), i = 1)

W_old_i <- c(0,0,0)

optim(W_i, Q_OSL_per_voxel, Q_OSL_grad_per_voxel,
      method="L-BFGS-B", lower=lb, upper=ub, 
      TE = TE_train, TR = TR_train, sigma = sigma_train, 
      r_i = train_i, W_old_i = W_i, c_i = c(0,0,0), i = 1)


### Iteration 
while(sum(abs(W_i - W_old_i)) > 1e-5){
  W_old_i <- W_i
  optim_new <- optim(W_i, Q_OSL_per_voxel, Q_OSL_grad_per_voxel,
                     method="L-BFGS-B", lower=lb, upper=ub, 
                     TE = TE_train, TR = TR_train, sigma = sigma_train, 
                     r_i = train_i, W_old_i = W_old_i, c_i = c(0,0,0), i = 1)
  W_i <- optim_new$par
  print(W_i)
  print(optim_new$value)
  cat("Another iter stops\n\n\n\n")
}

