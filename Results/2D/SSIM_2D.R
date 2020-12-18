### Calculating SSIM: 


###### Reading and cropping the data ##########
library(oro.nifti)
(ffd <- readNIfTI("~/R/Researchs/MRI/Read_Data/new_phantom.nii.gz"))
# (ffd <- readNIfTI("~/Researches/MRI/data/new_phantom.nii"))
image(ffd)

tmp <- ffd
(image1 <- nifti(tmp, datatype = datatype(ffd)))
dim(image1)
image(image1)
rm(tmp)



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


## phantom
TE_example <- c(0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10)
TR_example <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)

(TE_scale = 2.01/min(TE_example))
(TR_scale = 2.01/min(TR_example))
(TE_example = TE_example * TE_scale)
(TR_example = TR_example * TR_scale)




r <- array(dim = c(n_x*n_y*n_z, length(TE_example)))  # n x m matrix (Check the order)
for(i in 1:length(TE_example)){
  r[,i] <- as.vector(image2[,,,i])
}


## normalize:
(normalizer <- max(r))
normalizer <- 10
r <- r/normalizer






###### Indices:

## phantom
train_ind <- c(1, 10, 15)
test_ind <- setdiff(1:18, train_ind)
test <- r[,test_ind]




SSIM_val <- array(dim=c(length(test_ind),4))



v_mat_LS <- as.matrix(read.table("~/R/Researchs/MRI/Git/synthetic_mri_1/Results/2D/v_predicted_LS_29.txt"))
v_mat_MLE <- as.matrix(read.table("~/R/Researchs/MRI/Git/synthetic_mri_1/Results/2D/v_predicted_Likeli_29.txt"))
v_mat_OSL <- as.matrix(read.table("~/R/Researchs/MRI/Git/synthetic_mri_1/Results/2D/v_predicted_OSL_26.txt"))
v_mat_MPLE <- as.matrix(read.table("~/R/Researchs/MRI/Git/synthetic_mri_1/Results/2D/v_predicted_MPLE_29.txt"))

# v_mat_LS <- as.matrix(read.table("~/Researches/MRI/TRY_EIGEN_6_phantom/result/v_predicted_LS_29.txt"))
# v_mat_MLE <- as.matrix(read.table("~/Researches/MRI/TRY_EIGEN_6_phantom/result/v_predicted_Likeli_29.txt"))
# v_mat_OSL <- as.matrix(read.table("~/Researches/MRI/TRY_EIGEN_6_phantom/result/v_predicted_OSL_26.txt"))
# v_mat_MPLE <- as.matrix(read.table("~/Researches/MRI/TRY_EIGEN_6_phantom/result/v_predicted_MPLE_29.txt"))


for(i in 1:length(test_ind)){
  temp_LS <- v_mat_LS[,i];  temp_MLE <- v_mat_MLE[,i]; temp_OSL <- v_mat_OSL[,i]; temp_MPLE <- v_mat_MPLE[,i];
  temp_test <- test[,i]
  
  dim(temp_LS) <- c(n_x, n_y); dim(temp_MLE) <- c(n_x, n_y); dim(temp_OSL) <- c(n_x, n_y);dim(temp_MPLE) <- c(n_x, n_y);
  dim(temp_test) <- c(n_x, n_y)
  
  SSIM_val[i, 1] <- SpatialPack::SSIM(temp_LS, temp_test)$SSIM
  SSIM_val[i, 2] <- SpatialPack::SSIM(temp_MLE, temp_test)$SSIM
  SSIM_val[i, 3] <- SpatialPack::SSIM(temp_OSL, temp_test)$SSIM
  SSIM_val[i, 4] <- SpatialPack::SSIM(temp_MPLE, temp_test)$SSIM
}

# SSIM_val; plot(SSIM_val)
xtable::xtable(rbind(colMeans(SSIM_val), apply(SSIM_val, 2, sd)), digits=c(0, 6, 6, 6, 6))

