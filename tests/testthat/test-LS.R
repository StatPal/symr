##### LS #####

source("test-common.R")




###### Least Square #########


## ||v_i - r_i||^2
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



#test_that("Performance measures are okay",{
#  
#  expect_equal(least_sq_solve(), Init_val_least_sq_R())
#}
#)


