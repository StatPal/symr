library(reticulate)
# py_install("nibabel")
# py_install("skimage")

nib <-  import("nibabel")
skim <- import("skimage")
ssim <-  skim$metrics$structural_similarity
mad_new <- function(x){mean(abs(x-median(x)))}

# Mask
mask_all <-  nib$load('../DeepSynMRI/data/mask/subject47_crisp_v.mnc.gz')
mask_all <-  mask_all$get_fdata()
mask <-  (mask_all == 0)
mask_reshaped <-  mask[(1:36)*10, (1:217)*2, (1:181)*2]
mask <- mask_reshaped
mask_vec <- 1 - c(mask)





img <-  nib$load(paste0('../DeepSynMRI/data/noise-1-INU-00/brainweb_0.mnc.gz'))
shapes <-  unlist(img$shape)


## Other noises etc
library(symR)
train_ind <- c(1, 9, 10)
n <- prod(shapes)
phantom <- array(dim = c(n, 12))
test_ind <- setdiff(1:ncol(phantom), train_ind)

for (i in 1:3) {
    img <-  nib$load(paste0('../DeepSynMRI/data/noise-1-INU-00/brainweb_', i-1, '.mnc.gz'))    ## BUG, this was 0, not i
    data <-  img$get_fdata()
    phantom[, train_ind[i]] <- c(data)
}

for (i in 1:9) {
    img <-  nib$load(paste0('../DeepSynMRI/data/test-noise-0-INU-00/brainweb_', i-1, '.mnc.gz'))
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
### sigma_values <- c(1.99146, 1.81265, 1.82837, 2.30221, 1.63414, 1.71876, 3.13695, 1.77141, 1.55651, 2.72191, 1.63068, 1.4359)
if(!exists("sigma_values")) {
    sigma_values <- array(dim=12)
    for (i in 1:3) {
        print(i)
        sigma_values[train_ind[i]] <- symR::estimate.sigma.j(phantom[,train_ind[i]])
    }
}
# sigma_values <- array(dim=12)
# sigma_values[1] <- 5.688486; sigma_values[9] <- 3.451482; sigma_values[10] <- 6.014322; # pre computed values

print(sigma_values)


TE.scale <- 2.01 / min(TE_values); TR.scale <- 2.01 / min(TR_values); r_scale <- 1.0
TE_values <- TE_values * TE.scale; TR_values <- TR_values * TR.scale
phantom <- phantom / r_scale; sigma_values <- sigma_values / r_scale

## Divide into train and test with 3 train images:
train <- phantom[, train_ind];sigma.train <- sigma_values[train_ind];TE.train <- TE_values[train_ind]; TR.train <- TR_values[train_ind]
test <- phantom[, test_ind]; sigma_test <- sigma_values[test_ind]; TE_test <- TE_values[test_ind]; TR_test <- TR_values[test_ind]

test[mask==1,] <- 0
dimen <- c(4, shapes, 1) ## First element correspond to dim+1
















print("Methods will start")

if (!dir.exists('W_values')) {dir.create('W_values')}



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


## AECM after usual LS
if(!file.exists("W_values/W_AECM_1.rds")){
    t1 <- Sys.time()
    W_AECM_all <- symr(W_LS,
    method = "AECM", dimen, TE.train, TR.train, sigma.train, train,
    r_scale, TE.scale, TR.scale, as.numeric(mask),
                maxiter = 100
    )
    saveRDS(W_AECM_all, "W_values/W_AECM_1.rds")
    t2 <- Sys.time()
    print(t2-t1)
}

W_AECM_all <- readRDS("W_values/W_AECM_1.rds")
W_AECM <- W_AECM_all$W


W_AECM[mask==1,1] <- 0.597891
W_AECM[mask==1,2] <- 1e-8      
W_AECM[mask==1,3] <- 0.975431




# Contrast_GM
Contrast_GM <- nib$load('phantom_1.0mm_msles1_gray_matter.mnc.gz')$get_fdata()
Contrast_GM <- Contrast_GM[(1:36)*5, (1:217)*1, (1:181)*1]
Contrast_GM <- c(Contrast_GM)
Contrast_GM <- Contrast_GM/sum(Contrast_GM)


# sigma_test is NA, putting some temporary values
sigma_test <- runif(9, 3, 6)

library(Matrix)

prediction.var(W_AECM, W_AECM_all$Psi_inv, W_AECM_all$beta, Contrast_GM,
      dimen, TE.train, TR.train, sigma.train, train,
      TE_test, TR_test, sigma_test, test, r_scale, TE.scale, TR.scale,
      as.numeric(mask), 75L, 1e-06, 1L)




