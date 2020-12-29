LS   <- read.table("Least_sq_time.txt")
EM   <- read.table("EM_time.txt")
OSL  <- read.table("OSL_time.txt")
AECM <- read.table("AECM_time.txt")

dat <- data.frame(Method = rep(c("LS", "EM", "OSL", "AECM"), each=length(LS[,1])), 
				  time_spent = c(LS[,1], EM[,1], OSL[,1], AECM[,1]))

dat <- data.frame(Method = rep(c("OSL", "AECM"), each=length(LS[,1])), 
				  time_spent = c(OSL[,1], AECM[,1]))

library(ggplot2)
pdf("times.pdf")
ggplot(dat, aes(x=Method, y=time_spent, fill=Method)) + 
  geom_violin() + ylab("Time (microseconds)") + theme_bw() + theme(legend.position = "bottom") + 
  scale_y_log10() + coord_trans(y = "log10")

ggplot(dat, aes(x=Method, y=time_spent, fill=Method)) + 
  geom_violin() + ylab("Time (microseconds)") + theme_bw() + theme(legend.position = "bottom") 
  
ggplot(dat, aes(x=Method, y=time_spent, fill=Method)) + 
  geom_violin() + ylab("Time (microseconds)") + theme_bw() + theme(legend.position = "bottom") + 
  scale_y_log10()

ggplot(dat, aes(x=Method, y=time_spent, fill=Method)) + 
  geom_violin() + ylab("Time (microseconds)") + theme_bw() + theme(legend.position = "bottom") + 
  coord_trans(y = "log10")

dev.off()
