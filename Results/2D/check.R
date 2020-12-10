abc <- read.csv("Performances_29_all_3_2D.csv", header=F)
bcd <- read.csv("~/R/Researchs/MRI/Git/synthetic_mri_1/TRY_EIGEN_6_phantom/result/Performances_29_all_3.csv", header=F)

summary(abc - bcd)
