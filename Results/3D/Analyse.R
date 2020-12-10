setwd("/home/subrata/R/Researchs/MRI/Git/synthetic_mri_1/TRY_EIGEN_CODE_Rowmajor/result")
# setwd("/home/subrata/R/Researchs/MRI/Git/synthetic_mri_1/TRY_EIGEN_6_phantom/result")
test_dim <- 4

abc <- read.csv(paste0("Performances_29_all_", test_dim,  ".csv"), header=F)
abc[,1:test_dim] <- abc[,1:test_dim]+1
head(abc)
our_order <- order(abc[,test_dim+11]+abc[,test_dim+12]) ## Avg RMPSE


abc[our_order[1],test_dim+9+0:3]
rowMeans(abc[our_order[1],test_dim+9+0:1])


head(abc[our_order,], 10)


round(head(abc[our_order,], 10), 4)

head(abc[our_order,c(1:test_dim, test_dim+9+0:3)], 10)

xtable::xtable(head(abc[our_order,c(1:test_dim, test_dim+9+0:3)], 10), 
               digits = c(0,0,0,0,4,4,4,4))

# our_order_MLE <- order(abc$V8+abc$V9)
# head(abc[our_order_MLE,], 10)
