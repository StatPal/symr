DEEP_LEARNING_COMPARE <- TRUE

# install.packages("reticulate")
library(reticulate)
# py_install("nibabel")
# py_install("scikit-image")

dir.create("W_values")

nib <-  import("nibabel")
skim <- import("skimage")
ssim <-  skim$metrics$structural_similarity
mad_new <- function(x){mean(abs(x-median(x)))}

# Mask
mask_all <-  nib$load('./DeepSynMRI/data/mask/subject47_crisp_v.mnc.gz')
mask_all <-  mask_all$get_fdata()
mask <-  (mask_all == 0)
mask_reshaped <-  mask[(1:36)*10, (1:217)*2, (1:181)*2]
mask <- mask_reshaped
mask_vec <- 1 - c(mask)


img <-  nib$load(paste0('./DeepSynMRI/data/noise-1-INU-00/brainweb_0.mnc.gz'))
shapes <-  unlist(img$shape)


## Other noises etc
library(symR)
train_ind <- c(1, 9, 10)
n <- prod(shapes)
phantom <- array(dim = c(n, 12))
test_ind <- setdiff(1:ncol(phantom), train_ind)

for (i in 1:3) {
    img <-  nib$load(paste0('./DeepSynMRI/data/noise-1-INU-00/brainweb_', i-1, '.mnc.gz'))    ## BUG, this was 0, not i
    data <-  img$get_fdata()
    phantom[, train_ind[i]] <- c(data)
}

for (i in 1:9) {
    img <-  nib$load(paste0('./DeepSynMRI/data/test-noise-0-INU-00/brainweb_', i-1, '.mnc.gz'))
    data <-  img$get_fdata()
    phantom[, test_ind[i]] <- c(data)
}

train_scale_factor <- 400 / max(phantom, na.rm = T)
phantom <- phantom * train_scale_factor
apply(phantom, 2, max)

phantom[phantom == 0.0] <- 0.5 ## Pre-processing to remove the -Inf issue in likelihood.

## Other input parameters:
TE_values <- c(0.01, 0.015, 0.02, 0.01, 0.03, 0.04, 0.01, 0.04, 0.08, 0.01, 0.06, 0.1)
TR_values <- c(0.6,    0.6,  0.6,    1,    1,    1,    2,    2,    2,    3,    3,   3)
sigma_values <- c(1.99146, 1.81265, 1.82837, 2.30221, 1.63414, 1.71876, 3.13695, 1.77141, 1.55651, 2.72191, 1.63068, 1.4359)
if(!exists("sigma_values")) {
    sigma_values <- array(dim=12)
    for (i in 1:3) {
        print(i)
        sigma_values[train_ind[i]] <- symR::estimate.sigma.j(phantom[,train_ind[i]])
    }
}


TE.scale <- 2.01 / min(TE_values); TR.scale <- 2.01 / min(TR_values); r_scale <- 1.0
TE_values <- TE_values * TE.scale; TR_values <- TR_values * TR.scale
phantom <- phantom / r_scale; sigma_values <- sigma_values / r_scale

## Divide into train and test with 3 train images:
train <- phantom[, train_ind];sigma.train <- sigma_values[train_ind];TE.train <- TE_values[train_ind]; TR.train <- TR_values[train_ind]
test <- phantom[, test_ind]; sigma_test <- sigma_values[test_ind]; TE_test <- TE_values[test_ind]; TR_test <- TR_values[test_ind]

test[mask==1,] <- 0
dimen <- c(4, shapes, 1) ## First element correspond to dim+1
















print("Methods will start")



## LS - C
if(!file.exists("W_values/W_LS_1.rds")){
    t1<- Sys.time()
    W_LS <- symr(NULL,
                method = "LS", dimen, TE.train, TR.train, sigma.train, train,
                r_scale, TE.scale, TR.scale, mask,
                maxiter.LS = 100)
    saveRDS(W_LS, "W_values/W_LS_1.rds")
    t2 <-  Sys.time()
    print(t2-t1)
}
W_LS <- readRDS("W_values/W_LS_1.rds")


W_LS[mask==1,1] <- 0.597891
W_LS[mask==1,2] <- 1e-8        # controvertial
W_LS[mask==1,3] <- 0.975431

SSIM_vals <- array(dim=9)
pred_LS <- bloch.image(W_LS, TE_test, TR_test)
pred_LS[mask==1,] <- 0
for(i in 1:9){
    pred_3D = array(pred_LS[,i], shapes)   ## BUG shapes
    test_3D = array(test[,i], shapes)
    SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
}
print(SSIM_vals)
tmp_diff <- abs(pred_LS - test)
tmp_pred_val <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                        sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)



## OSL after usual LS
if(!file.exists("W_values/W_OSL_1.rds")){
    t1 <- Sys.time()
    W_OSL <- symr(W_LS,
    method = "OSL-EM", dimen, TE.train, TR.train, sigma.train, train,
    r_scale, TE.scale, TR.scale, as.numeric(mask),
                maxiter = 100
    )$W
    saveRDS(W_OSL, "W_values/W_OSL_1.rds")
    t2 <- Sys.time()
    print(t2-t1)
}
W_OSL <- readRDS("W_values/W_OSL_1.rds")


W_OSL[mask==1,1] <- 0.597891
W_OSL[mask==1,2] <- 1e-8        # controvertial
W_OSL[mask==1,3] <- 0.975431

SSIM_vals <- array(dim=9)
pred_OSL <- bloch.image(W_OSL, TE_test, TR_test)
pred_OSL[mask==1,] <- 0
for(i in 1:9){
    pred_3D = array(pred_OSL[,i], shapes)    ## Bug shapes
    test_3D = array(test[,i], shapes)
    SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
}
print(SSIM_vals)
tmp_diff <- abs(pred_OSL - test)
tmp_pred_val_OSL <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                        sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)



## AECM after usual LS
if(!file.exists("W_values/W_AECM_1.rds")){
    t1 <- Sys.time()
    W_AECM <- symr(W_LS,
    method = "AECM", dimen, TE.train, TR.train, sigma.train, train,
    r_scale, TE.scale, TR.scale, as.numeric(mask),
                maxiter = 100
    )$W
    saveRDS(W_AECM, "W_values/W_AECM_1.rds")
    t2 <- Sys.time()
    print(t2-t1)
}
W_AECM <- readRDS("W_values/W_AECM_1.rds")


W_AECM[mask==1,1] <- 0.597891
W_AECM[mask==1,2] <- 1e-8        # controvertial
W_AECM[mask==1,3] <- 0.975431

SSIM_vals <- array(dim=9)
pred_AECM <- bloch.image(W_AECM, TE_test, TR_test)
pred_AECM[mask==1,] <- 0
for(i in 1:9){
    pred_3D = array(pred_AECM[,i], shapes)   ## Bug shapes
    test_3D = array(test[,i], shapes)
    SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
}
print(SSIM_vals)
tmp_diff <- abs(pred_AECM - test)
tmp_pred_val_AECM <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                        sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)







if(DEEP_LEARNING_COMPARE){

    ## From python: 
    W_LS_py <- as.matrix( read.csv('./DeepSynMRI/LS-py/W_1.csv.gz', header=F) )
    # py orientation to R orientation:
    for(i in 1:3){
        tmp <- W_LS_py[,i]
        tmp <- array_reshape(tmp, shapes)
        W_LS_py[,i] <- c(tmp)
    }


    W_LS_py[,1] <- ifelse(W_LS_py[,1] < 0.59789100, 0.59789100, W_LS_py[,1])
    W_LS_py[,2] <- ifelse(W_LS_py[,2] < 0.00000001, 0.00000001, W_LS_py[,2])
    W_LS_py[,3] <- ifelse(W_LS_py[,3] < 0.09120034, 0.09120034, W_LS_py[,3])
    W_LS_py[mask==1,1] <- 0.597891
    W_LS_py[mask==1,2] <- 1e-8        # controvertial
    W_LS_py[mask==1,3] <- 0.975431

    SSIM_vals <- array(dim=9)
    pred_LS_py <- bloch.image(W_LS_py, TE_test, TR_test)
    pred_LS_py[mask==1,] <- 0
    for(i in 1:9){
        pred_3D = array(pred_LS_py[,i], shapes)
        test_3D = array(test[,i], shapes)
        SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
    }
    print(SSIM_vals)
    tmp_diff <- abs(pred_LS_py - test)
    tmp_pred_val_py <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                            sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)




    ## OSL after py LS
    t1 <- Sys.time()
    W_OSL_py <- symr(W_LS_py,
    method = "OSL-EM", dimen, TE.train, TR.train, sigma.train, train,
    r_scale, TE.scale, TR.scale, as.numeric(mask),
                maxiter = 100
    )$W
    saveRDS(W_OSL_py, "W_values/W_OSL_py_1.rds")
    t2 <- Sys.time()
    print(t2-t1)
    W_OSL_py <- readRDS("W_values/W_OSL_py_1.rds")

    W_OSL_py[mask==1,1] <- 0.597891
    W_OSL_py[mask==1,2] <- 1e-8        # controvertial
    W_OSL_py[mask==1,3] <- 0.975431

    print(SSIM_vals)
    pred_OSL_py <- bloch.image(W_OSL_py, TE_test, TR_test)
    pred_OSL_py[mask==1,] <- 0
    for(i in 1:9){
        pred_3D = array(pred_OSL_py[,i], shapes)     ## Bug shapes
        test_3D = array(test[,i], shapes)
        SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
    }
    print(SSIM_vals)
    tmp_diff <- abs(pred_OSL_py - test)
    tmp_pred_val_OSL_py <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                            sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)




    ## AECM after py LS
    t1 <- Sys.time()
    W_AECM_py <- symr(W_LS_py,
    method = "AECM", dimen, TE.train, TR.train, sigma.train, train,
    r_scale, TE.scale, TR.scale, as.numeric(mask),
    maxiter = 100
    )$W
    saveRDS(W_AECM_py, "W_values/W_AECM_py_1.rds")
    t2 <- Sys.time()
    print(t2-t1)
    W_AECM_py <- readRDS("W_values/W_AECM_py_1.rds")

    W_AECM_py[mask==1,1] <- 0.597891
    W_AECM_py[mask==1,2] <- 1e-8        # controvertial
    W_AECM_py[mask==1,3] <- 0.975431

    SSIM_vals <- array(dim=9)
    pred_AECM_py <- bloch.image(W_AECM_py, TE_test, TR_test)
    pred_AECM_py[mask==1,] <- 0
    for(i in 1:9){
        pred_3D = array(pred_AECM_py[,i], shapes)
        test_3D = array(test[,i], shapes)
        SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
    }
    print(SSIM_vals)
    tmp_diff <- abs(pred_AECM_py - test)
    tmp_pred_val_AECM_py <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                            sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)




    ## DL 
    W_LS_DL_py <- as.matrix( read.csv('./DeepSynMRI/LS-py/W_DL_1.csv.gz', header=F) )
    # py orientation to R orientation:
    for(i in 1:3){
        tmp <- W_LS_DL_py[,i]
        tmp <- array_reshape(tmp, shapes)
        W_LS_DL_py[,i] <- c(tmp)
    }

    W_LS_DL_py[mask==1,1] <- 0.597891
    W_LS_DL_py[mask==1,2] <- 1e-8        # controvertial
    W_LS_DL_py[mask==1,3] <- 0.975431

    SSIM_vals <- array(dim=9)
    pred_LS_DL_py <- bloch.image(W_LS_DL_py, TE_test, TR_test)
    pred_LS_DL_py[mask==1,] <- 0
    for(i in 1:9){
        pred_3D = array(pred_LS_DL_py[,i], shapes)
        test_3D = array(test[,i], shapes)
        SSIM_vals[i] = ssim(pred_3D, test_3D, data_range = max(c(pred_3D, test_3D)))
    }
    print(SSIM_vals)
    tmp_diff <- abs(pred_LS_DL_py - test)
    tmp_pred_val_DL_py <- cbind(colMeans(tmp_diff[mask_vec==1,])/apply(test[mask_vec==1,], 2, mad_new), 
                            sqrt(colMeans(tmp_diff[mask_vec==1,]^2))/apply(test[mask_vec==1,], 2, sd), SSIM_vals)


}



# print(tmp_pred_val)
# print(tmp_pred_val_OSL)
# print(tmp_pred_val_AECM)
# print(tmp_pred_val_py)
# print(tmp_pred_val_OSL_py)
# print(tmp_pred_val_AECM_py)
# print(tmp_pred_val_DL_py)


print(mean(tmp_pred_val))
print(mean(tmp_pred_val_OSL))
print(mean(tmp_pred_val_AECM))
if(DEEP_LEARNING_COMPARE){
    print(mean(tmp_pred_val_py))
    print(mean(tmp_pred_val_OSL_py))
    print(mean(tmp_pred_val_AECM_py))
    print(mean(tmp_pred_val_DL_py))
}



## Having the results in a transposed format
# tmp_all_1 <- rbind(tmp_pred_val[,1],tmp_pred_val_OSL[,1], tmp_pred_val_AECM[,1], tmp_pred_val_py[,1], tmp_pred_val_OSL_py[,1],tmp_pred_val_AECM_py[,1], tmp_pred_val_DL_py[,1])
# tmp_all_2 <- rbind(tmp_pred_val[,2],tmp_pred_val_OSL[,2], tmp_pred_val_AECM[,2], tmp_pred_val_py[,2], tmp_pred_val_OSL_py[,2],tmp_pred_val_AECM_py[,2], tmp_pred_val_DL_py[,2])
# tmp_all_3 <- 100 * rbind(tmp_pred_val[,3],tmp_pred_val_OSL[,3], tmp_pred_val_AECM[,3], tmp_pred_val_py[,3], tmp_pred_val_OSL_py[,3],tmp_pred_val_AECM_py[,3], tmp_pred_val_DL_py[,3])


# rowMeans(tmp_all_1)
# rowMeans(tmp_all_2)
# rowMeans(tmp_all_3)



# all_val <- rbind(tmp_all_1, tmp_all_2, tmp_all_3)
# colnames(all_val) <- paste("test", 1:9)
# rownames(all_val) <- rep(c("MAPE", "RMSPE", "SSIM"), each=7)
# saveRDS(all_val, "values/SE_1.rds")

# (ord_1 <- order(rowMeans(tmp_all_1)))
# (ord_2 <- order(rowMeans(tmp_all_2)))
# (ord_3 <- order(rowMeans(tmp_all_3)))

# method <- c("LS-C", "OSL-C","AECM-C","LS-py","OSL-py","AECM-py", "DL")

# method[ord_1]
# method[ord_2]
# method[ord_3]



## Plots, with DL
# pdf('images_1.pdf', height=4*2, width=16*2)
# par(mfrow=c(1,4))

# for(i in 1:9) {
#     test_3D <- array(test[,i], shapes)
#     pred_3D <- array(pred_LS_py[,i], shapes)
#     pred_3D_AECM <- array(pred_AECM_py[,i], shapes)
#     pred_3D_DL <- array(pred_LS_DL_py[,i], shapes)

#     range_val <- range(c(test_3D[18,,], pred_3D[18,,], pred_3D_AECM[18,,], pred_3D_DL[18,,]))

#     image(test_3D[18,,], zlim=range_val)
#     image(pred_3D[18,,], zlim=range_val)
#     image(pred_3D_AECM[18,,], zlim=range_val)
#     image(pred_3D_DL[18,,], zlim=range_val)
# }
# dev.off()


