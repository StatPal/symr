#abc <- read.csv("../examples/2D/Hessian_Matrix_26.csv", header=F)

#n <- length(abc$V1)

#library(Matrix)
#A <- new("dgTMatrix", i = abc$V1, j = abc$V2, x = as.numeric(abc$V3), Dim=c(n,n))

#val <- band(A,0,0)@x
#summary(val)
#min(which(val < 0))


