
abc <- read.csv("Performances_29_all_3.csv", header=F)
abc[,1:3] <- abc[,1:3]+1
head(abc)
our_order <- order(abc$V14+abc$V15)
head(abc[our_order,], 10)

head(abc[our_order,1:3], 10)


our_order_MLE <- order(abc$V8+abc$V9)
head(abc[our_order_MLE,], 10)
