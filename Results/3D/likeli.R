likeli_26 <- read.table("likeli_val_26_refined.txt")[,1]
likeli_29 <- read.table("likeli_val_29_refined.txt")[,1]


pdf("likeli_plot.pdf", height=3*1.25, width=6)

plot(-log(-likeli_26[1:30]), type="l", ylab="-log(-log-likelihood)", xlab=print("iteration+1"), ylim=range(-log(-c(likeli_26, likeli_29))), lty=2, col="blue")
lines(-log(-likeli_29), col="red", lty = 1)
legend(21, -21, legend=c("OSL-EM", "AECM"), col=c("blue", "red"), lty=2:1, cex=0.8)



## Remove first cases:
likeli_26 <- likeli_26[-1]
likeli_29 <- likeli_29[-1]

plot(likeli_26[1:30], type="l", ylab="Penalized log-likelihood", xlab=print("iteration+1"), ylim=range(c(likeli_26, likeli_29)), lty=2, col="blue")
lines(likeli_29, col="red", lty = 1)

legend(21, -12111964, legend=c("OSL-EM", "AECM"), col=c("blue", "red"), lty=2:1, cex=0.8)






library(ggplot2)
likeli_26 <- read.table("likeli_val_26_refined.txt")
likeli_29 <- read.table("likeli_val_29_refined.txt")
## Remove first cases:
likeli_26 <- data.frame(V1 = likeli_26[-1,])
likeli_26 <- data.frame(V1 = likeli_26[1:30,])
likeli_29 <- data.frame(V1 = likeli_29[-1,])


cols <- c("OSL" = "blue", "AECM" = "#009E73")
ltypes <- c("OSL-EM" = "twodash", "AECM" = "solid")

ggplot() + 
  geom_line(data = likeli_26, mapping=aes(x = 1:30, y = V1, color = "OSL", linetype = "OSL-EM")) +
  geom_line(data = likeli_29, mapping=aes(x = 1:nrow(likeli_29), y = V1, color = "AECM", linetype="AECM")) +  
  scale_colour_manual(name = 'Method', 
                      values = cols, 
                      labels = c('OSL-EM','AECM')) + 
  scale_linetype_manual(values=ltypes, name="Method")+
  guides(linetype = guide_legend(override.aes=list(color=c("#009E73", "blue"))), 
         color=FALSE) +
  xlab("iteration+1") + 
  ylab("Penalized log-likelihood") + 
  theme_bw() + theme(legend.position = c(0.8, 0.8), text = element_text(size=15))


## Does not work
## though stackoverflow ans says this:
## google "scale_colour_manual and scale_linetype_manual"
## Oh, name, labels both should be same
ggplot() +
  geom_line(data = likeli_29, mapping=aes(x = 1:nrow(likeli_29), y = V1, color = "#009E73", linetype="solid")) +
  geom_line(data = likeli_26, mapping=aes(x = 1:30, y = V1, color = "blue", linetype = "dashed")) + ## Weird
  scale_colour_manual(name = 'Method',
                      values = c("#009E73", "blue"),
                      labels = c('OSL-EM','AECM')) +
  scale_linetype_manual(name = 'Method', 
                        values=c("twodash", "solid"), 
                        labels = c('OSL-EM','AECM')) +
  xlab("iteration+1") +
  ylab("Penalized log-likelihood") +
  theme_bw() + theme(legend.position = c(0.8, 0.8), text = element_text(size=15))


dev.off()
