##
## Stage 0: pre-processing images
##
library(png)
library(SpatialPack)

# loading images
data(texmos2) # 'texmos2' is available in SpatialPack
lena   <- readPNG("lena_color.png")
baboon <- readPNG("baboon_color.png")

# converting RGB image to grayscale ('texmos2' is already in grayscale)
lena   <- RGB2gray(lena, method = "ITU")
baboon <- RGB2gray(baboon, method = "ITU")

# normalizing images to [0,1]
lena    <- normalize(lena)
baboon  <- normalize(baboon)
texmos2 <- normalize(texmos2)

##
## Stage 1: processing images
##
dyn.load("SSIM.so")

# reading R sources
source("../code/SSIM_FIT.R")
source("../code/SSIM_TEST.R")
source("gamma.simul.R")

# processing 'texmos2'
texmos2.01 <- gamma.simul(texmos2, Nsize = 100, looks = 1)
texmos2.02 <- gamma.simul(texmos2, Nsize = 100, looks = 2)
texmos2.04 <- gamma.simul(texmos2, Nsize = 100, looks = 4)
texmos2.08 <- gamma.simul(texmos2, Nsize = 100, looks = 8)
texmos2.16 <- gamma.simul(texmos2, Nsize = 100, looks = 16)
texmos2.32 <- gamma.simul(texmos2, Nsize = 100, looks = 32)

##
## Stage 2: output
##

## Table 1: no filter
rbind(apply(texmos2.01$no.filter$pars, 2, mean),
apply(texmos2.02$no.filter$pars, 2, mean),
apply(texmos2.04$no.filter$pars, 2, mean),
apply(texmos2.08$no.filter$pars, 2, mean),
apply(texmos2.16$no.filter$pars, 2, mean),
apply(texmos2.32$no.filter$pars, 2, mean))

#looks    alpha     beta    gamma
# 1    1.000000 1.000000 1.000000
# 2    1.000000 1.000000 1.000000
# 4    1.036826 1.050887 1.036529
# 8    1.081879 1.112977 1.081548
# 16   1.159801 1.214050 1.159511
# 32   1.297262 1.385889 1.296982

## Table 1: Lee filter
rbind(apply(texmos2.01$lee$pars, 2, mean),
apply(texmos2.02$lee$pars, 2, mean),
apply(texmos2.04$lee$pars, 2, mean),
apply(texmos2.08$lee$pars, 2, mean),
apply(texmos2.16$lee$pars, 2, mean),
apply(texmos2.32$lee$pars, 2, mean))

#looks    alpha     beta    gamma
# 1    1.178049 1.232842 1.165766
# 2    1.247468 1.325862 1.241308
# 4    1.308159 1.400332 1.304093
# 8    1.430952 1.540320 1.429447
# 16   1.613290 1.748392 1.612049
# 32   1.824969 1.957541 1.824176

## Table 1: Enhanced Lee filter
rbind(apply(texmos2.01$enhanced$pars, 2, mean),
apply(texmos2.02$enhanced$pars, 2, mean),
apply(texmos2.04$enhanced$pars, 2, mean),
apply(texmos2.08$enhanced$pars, 2, mean),
apply(texmos2.16$enhanced$pars, 2, mean),
apply(texmos2.32$enhanced$pars, 2, mean))

#looks    alpha     beta    gamma
# 1    1.294993 1.375068 1.278478
# 2    1.472828 1.597601 1.460297
# 4    1.573338 1.711984 1.566482
# 8    1.763058 1.919510 1.758788
# 16   1.945709 2.121587 1.942399
# 32   2.122518 2.329744 2.120488

## Table 1: Kuan filter
rbind(apply(texmos2.01$kuan$pars, 2, mean),
apply(texmos2.02$kuan$pars, 2, mean),
apply(texmos2.04$kuan$pars, 2, mean),
apply(texmos2.08$kuan$pars, 2, mean),
apply(texmos2.16$kuan$pars, 2, mean),
apply(texmos2.32$kuan$pars, 2, mean))

#looks    alpha     beta    gamma
# 1    1.254216 1.132388 1.101939
# 2    1.496419 1.409050 1.316629
# 4    1.502808 1.610939 1.488109
# 8    1.745729 1.913125 1.742593
# 16   1.905577 2.065849 1.905297
# 32   2.190247 2.355337 2.191109

## Table 4: no filter
rbind(apply(cbind(texmos2.01$no.filter$H0[,1], texmos2.01$no.filter$H1[,1]), 2, mean),
apply(cbind(texmos2.02$no.filter$H0[,1], texmos2.02$no.filter$H1[,1]), 2, mean),
apply(cbind(texmos2.04$no.filter$H0[,1], texmos2.04$no.filter$H1[,1]), 2, mean),
apply(cbind(texmos2.08$no.filter$H0[,1], texmos2.08$no.filter$H1[,1]), 2, mean),
apply(cbind(texmos2.16$no.filter$H0[,1], texmos2.16$no.filter$H1[,1]), 2, mean),
apply(cbind(texmos2.32$no.filter$H0[,1], texmos2.32$no.filter$H1[,1]), 2, mean))

#looks        H0        H1
# 1    0.5508976 0.5508976
# 2    0.6898613 0.6898613
# 4    0.8024050 0.7959750
# 8    0.8821055 0.8731275
# 16   0.9330943 0.9228411
# 32   0.9633495 0.9527184

## Table 4: Lee filter
rbind(apply(cbind(texmos2.01$lee$H0[,1], texmos2.01$lee$H1[,1]), 2, mean),
apply(cbind(texmos2.02$lee$H0[,1], texmos2.02$lee$H1[,1]), 2, mean),
apply(cbind(texmos2.04$lee$H0[,1], texmos2.04$lee$H1[,1]), 2, mean),
apply(cbind(texmos2.08$lee$H0[,1], texmos2.08$lee$H1[,1]), 2, mean),
apply(cbind(texmos2.16$lee$H0[,1], texmos2.16$lee$H1[,1]), 2, mean),
apply(cbind(texmos2.32$lee$H0[,1], texmos2.32$lee$H1[,1]), 2, mean))

#looks        H0        H1
# 1    0.7582669 0.7220599
# 2    0.8565451 0.8237672
# 4    0.9174231 0.8929375
# 8    0.9533598 0.9336020
# 16   0.9741070 0.9583676
# 32   0.9858325 0.9742058

## Table 4: Enhanced Lee filter
rbind(apply(cbind(texmos2.01$enhanced$H0[,1], texmos2.01$enhanced$H1[,1]), 2, mean),
apply(cbind(texmos2.02$enhanced$H0[,1], texmos2.02$enhanced$H1[,1]), 2, mean),
apply(cbind(texmos2.04$enhanced$H0[,1], texmos2.04$enhanced$H1[,1]), 2, mean),
apply(cbind(texmos2.08$enhanced$H0[,1], texmos2.08$enhanced$H1[,1]), 2, mean),
apply(cbind(texmos2.16$enhanced$H0[,1], texmos2.16$enhanced$H1[,1]), 2, mean),
apply(cbind(texmos2.32$enhanced$H0[,1], texmos2.32$enhanced$H1[,1]), 2, mean))

#looks        H0        H1
# 1    0.7721481 0.7148530
# 2    0.8826419 0.8301850
# 4    0.9404165 0.9065011
# 8    0.9691358 0.9454405
# 16   0.9836372 0.9679935
# 32   0.9910211 0.9808030

## Table 4: Kuan filter
rbind(apply(cbind(texmos2.01$kuan$H0[,1], texmos2.01$kuan$H1[,1]), 2, mean),
apply(cbind(texmos2.02$kuan$H0[,1], texmos2.02$kuan$H1[,1]), 2, mean),
apply(cbind(texmos2.04$kuan$H0[,1], texmos2.04$kuan$H1[,1]), 2, mean),
apply(cbind(texmos2.08$kuan$H0[,1], texmos2.08$kuan$H1[,1]), 2, mean),
apply(cbind(texmos2.16$kuan$H0[,1], texmos2.16$kuan$H1[,1]), 2, mean),
apply(cbind(texmos2.32$kuan$H0[,1], texmos2.32$kuan$H1[,1]), 2, mean))

#looks        H0        H1
# 1    0.7948251 0.7705098
# 2    0.8880360 0.8505996
# 4    0.9404560 0.9110584
# 8    0.9683920 0.9445934
# 16   0.9831806 0.9678051
# 32   0.9910150 0.9802475

##
## REMARK: Modify the above commands to build Tables 2, 3 and 5, 6 with images
##         Lena and Baboon.
##
