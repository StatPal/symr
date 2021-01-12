library(RColorBrewer)

TE <- c(0.03, 0.08, 0.01)
TR <- c(3, 3, 0.035)

n_x <- 181; n_y <- 217; n_z <- 181

seq_x <- 1:181; seq_y <- 1:217; seq_z <- 1:181
s_x <- length(seq_x)/100
s_y <- length(seq_y)/100
s_z <- length(seq_z)/100



pdf("brainweb_0.pdf", width=s_y*3, height=s_x*1)
par(mai = c(0,0,0,0), mfrow = c(1,3), mar=rep(0.01, 4))

W_orig_1 <- read.table("W_orig_29_brainweb_1.0.txt")
for(i in 1:3){
	temp_test <- W_orig_1[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()



pdf("brainweb_1.pdf", width=s_y*3, height=s_x*3)
par(mai = c(0,0,0,0), mfrow = c(3,3), mar=rep(0.01, 4))

W_LS_1 <- read.table("W_LS_29_brainweb_1.0.txt")
W_OSL_1 <- read.table("W_OSL_29_brainweb_1.0.txt")
W_AECM_1 <- read.table("W_AECM_29_brainweb_1.0.txt")

for(i in 1:3){
	temp_test <- W_LS_1[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_OSL_1[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_AECM_1[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()







pdf("brainweb_3.pdf", width=s_y*3, height=s_x*3)
par(mai = c(0,0,0,0), mfrow = c(3,3), mar=rep(0.01, 4))

#W_orig_3 <- read.table("W_orig_29_brainweb_3.0.txt")
W_LS_3 <- read.table("W_LS_29_brainweb_3.0.txt")
W_OSL_3 <- read.table("W_OSL_29_brainweb_3.0.txt")
W_AECM_3 <- read.table("W_AECM_29_brainweb_3.0.txt")

for(i in 1:3){
	temp_test <- W_LS_3[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_OSL_3[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_AECM_3[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()






pdf("brainweb_5.pdf", width=s_y*3, height=s_x*3)
par(mai = c(0,0,0,0), mfrow = c(3,3), mar=rep(0.01, 4))

#W_orig_5 <- read.table("W_orig_29_brainweb_5.0.txt")
W_LS_5 <- read.table("W_LS_29_brainweb_5.0.txt")
W_OSL_5 <- read.table("W_OSL_29_brainweb_5.0.txt")
W_AECM_5 <- read.table("W_AECM_29_brainweb_5.0.txt")

for(i in 1:3){
	temp_test <- W_LS_5[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_OSL_5[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_AECM_5[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()






pdf("brainweb_10.pdf", width=s_y*3, height=s_x*3)
par(mai = c(0,0,0,0), mfrow = c(3,3), mar=rep(0.01, 4))

#W_orig_10 <- read.table("W_orig_29_brainweb_10.0.txt")
W_LS_10 <- read.table("W_LS_29_brainweb_10.0.txt")
W_OSL_10 <- read.table("W_OSL_29_brainweb_10.0.txt")
W_AECM_10 <- read.table("W_AECM_29_brainweb_10.0.txt")

for(i in 1:3){
	temp_test <- W_LS_10[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_OSL_10[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_AECM_10[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()





pdf("brainweb_20.pdf", width=s_y*3, height=s_x*3)
par(mai = c(0,0,0,0), mfrow = c(3,3), mar=rep(0.01, 4))

#W_orig_20 <- read.table("W_orig_29_brainweb_20.0.txt")
W_LS_20 <- read.table("W_LS_29_brainweb_20.0.txt")
W_OSL_20 <- read.table("W_OSL_29_brainweb_20.0.txt")
W_AECM_20 <- read.table("W_AECM_29_brainweb_20.0.txt")

for(i in 1:3){
	temp_test <- W_LS_20[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_OSL_20[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
for(i in 1:3){
	temp_test <- W_AECM_20[,i]
	dim(temp_test) <- c(n_x, n_y, n_z)
	image(temp_test[91,seq_y,seq_z], zlim=range(temp_test[91,seq_y,seq_z]), xaxt='n', yaxt='n', ann=FALSE,
	  		col=colorRampPalette(c("#FFFFFF", (brewer.pal(n=9, name="Purples")), "#000000"))(16))
}
dev.off()

