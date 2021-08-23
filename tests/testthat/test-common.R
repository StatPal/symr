SQ <- function(x){x*x}

Bloch_vec <- function(W_row, TE1, TR1){
  W_row <- as.numeric(W_row)
  W_row[1] * (1 - exp(TR1*log(W_row[2]))) * exp(TE1*log(W_row[3]))
}

### dv_ij/dW_ik
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


eigenvals_J_n <- function(k, if_narural_order = T){
  d_vec <- array(dim = k)
  for(i in 1:k){
    d_vec[i] <- 2*(1-cos(pi*(i-1)/k))
  }
  d_vec
}


mean_rice <- function(nu, sigma){
  x = - SQ(nu)/(2*SQ(sigma))
  return (sigma * sqrt(pi/2) * ( (1-x)*besselI(-x/2, 0, T) - x * besselI(-x/2, 1, T)))
}




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


