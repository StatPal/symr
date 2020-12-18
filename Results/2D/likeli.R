likeli_26 <- read.table("likeli_val_26_refined.txt")[,1]
likeli_29 <- read.table("likeli_val_29_refined.txt")[,1]

## Remove first cases:
likeli_26 <- likeli_26[-1]
likeli_29 <- likeli_29[-1]

pdf("likeli_plot.pdf")
plot(likeli_26[1:100], type="l", ylab="log-likelihood", xlab=print("iteration+1"), ylim=range(c(likeli_26, likeli_29)), lty=2)
lines(likeli_29, col=2, lty = 1)

legend(80, -707158, legend=c("OSL-EM", "AECM"),
       col=1:2, lty=2:1, cex=0.8)


library(ggplot2)
likeli_26 <- read.table("likeli_val_26_refined.txt")
likeli_29 <- read.table("likeli_val_29_refined.txt")
## Remove first cases:
likeli_26 <- data.frame(V1 = likeli_26[-1,])
likeli_26 <- data.frame(V1 = likeli_26[1:100,])
likeli_29 <- data.frame(V1 = likeli_29[-1,])


cols <- c("OSL" = "black", "AECM" = "red")
ltypes <- c("OSL-EM" = "twodash", "AECM" = "solid")

ggplot() + 
  geom_line(data = likeli_26, mapping=aes(x = 1:100, y = V1, color = "OSL", linetype = "OSL-EM")) +
  geom_line(data = likeli_29, mapping=aes(x = 1:3, y = V1, color = "AECM", linetype="AECM")) +  
  scale_colour_manual(name = 'Method', 
                      values = cols, 
                      labels = c('OSL-EM','AECM')) + 
  scale_linetype_manual(values=ltypes, name="Method")+
  guides(linetype = guide_legend(override.aes=list(color=c(2,1))), 
         color=FALSE) +
  xlab("iteration+1") + 
  ylab("log-likelihood") + 
  theme_bw() + theme(legend.position = "bottom")


## Does not work
## though stackoverflow ans says this:
## google "scale_colour_manual and scale_linetype_manual"
## Oh, name, labels both should be same
ggplot() +
  geom_line(data = likeli_26, mapping=aes(x = 1:100, y = V1, color = "black", linetype = "dashed")) + ## Weird
  geom_line(data = likeli_29, mapping=aes(x = 1:3, y = V1, color = "red", linetype="solid")) +
  scale_colour_manual(name = 'Method',
                      values = c("black", "red"),
                      labels = c('OSL-EM','AECM')) +
  scale_linetype_manual(name = 'Method', 
                        values=c("twodash", "solid"), 
                        labels = c('OSL-EM','AECM')) +
  xlab("iteration+1") +
  ylab("log-likelihood") +
  theme_bw() + theme(legend.position = "bottom")


dev.off()
