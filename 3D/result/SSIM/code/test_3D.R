dyn.load("~/R/Researchs/MRI/Git/synthetic_mri_1/3D/result/SSIM/code/src/SSIM.so")
source("~/R/Researchs/MRI/Git/synthetic_mri_1/3D/result/SSIM/code/SSIM.R")

library(SpatialPack)
data(texmos2)
detach("package:SpatialPack", unload=TRUE)

y <- SpatialPack::imnoise(texmos2, type = "gaussian")
SSIM(texmos2, y)
mean(y)


dim(texmos2) <- c(512, 32, 16)
dim(y) <- c(512, 32, 16)

SSIM(texmos2, y)
mean(y)

