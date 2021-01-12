## Id: gamma.simul.R, last updated 2020/07/06, F.Osorio

gamma.simul <- function(img, Nsize = 1000, looks = 1, alpha = 0.05, L = 255)
{ # function to perform to simulation experiment
  ok   <- matrix(0, nrow = Nsize, ncol = 14)
  lee  <- matrix(0, nrow = Nsize, ncol = 14)
  lee2 <- matrix(0, nrow = Nsize, ncol = 14)
  kuan <- matrix(0, nrow = Nsize, ncol = 14)
  now <- proc.time()

  img0 <- round(L * img)

  df <- 3
  cutoff <- qchisq(1 - alpha, df)

  cat("  Progress:\n")
  pb <- txtProgressBar(min = 0, max = Nsize, style = 3)

  # Monte Carlo iterations
  for (i in 1:Nsize) {
    set.seed(101 + i) # setting seeds

    # adding multiplicative noise
    img1 <- imnoise(img, type = "gamma", looks = looks)

    # denoising with Lee filter
    img2 <- denoise(img1, type = "Lee", looks = looks)

    # denoising with enhanced Lee filter
    img3 <- denoise(img1, type = "enhanced", looks = looks)

    # denoising with Kuan filter
    img4 <- denoise(img1, type = "Kuan", looks = looks)

    # SSIM coefficients estimation (no filter) #################################
    img1 <- round(L * img1)
    fm <- SSIM.fit(img0, img1, maxiter = 50)
    cf <- fm$coefficients
    ok[i,1:3] <- cf

    # gradient statistic to test
    z <- SSIM.test(fm)
    ok[i,4] <- ifelse(z$stat > cutoff, 1, 0)
    ok[i,5] <- z$stat
    ok[i,6] <- z$logLik

    # SSIM under H0: alpha = beta = gamma = 1
    o <- SSIM(img0, img1)
    ok[i,7] <- o$SSIM
    ok[i,8:10] <- o$comps

    # SSIM under alternative hypothesis
    o <- SSIM(img0, img1, alpha = cf[1], beta = cf[2], gamma = cf[3])
    ok[i,11] <- o$SSIM
    ok[i,12:14] <- o$comps

    # SSIM coefficients estimation (Lee filter) ################################
    img2 <- round(L * img2)
    fm <- SSIM.fit(img0, img2, maxiter = 50)
    cf <- fm$coefficients
    lee[i,1:3] <- cf

    # gradient statistic to test
    z <- SSIM.test(fm)
    lee[i,4] <- ifelse(z$stat > cutoff, 1, 0)
    lee[i,5] <- z$stat
    lee[i,6] <- z$logLik

    # SSIM under H0: alpha = beta = gamma = 1
    o <- SSIM(img0, img2)
    lee[i,7] <- o$SSIM
    lee[i,8:10] <- o$comps

    # SSIM under alternative hypothesis
    o <- SSIM(img0, img2, alpha = cf[1], beta = cf[2], gamma = cf[3])
    lee[i,11] <- o$SSIM
    lee[i,12:14] <- o$comps

    # SSIM coefficients estimation (enhanced Lee filter) #######################
    img3 <- round(L * img3)
    fm <- SSIM.fit(img0, img3, maxiter = 50)
    cf <- fm$coefficients
    lee2[i,1:3] <- cf

    # gradient statistic to test
    z <- SSIM.test(fm)
    lee2[i,4] <- ifelse(z$stat > cutoff, 1, 0)
    lee2[i,5] <- z$stat
    lee2[i,6] <- z$logLik

    # SSIM under H0: alpha = beta = gamma = 1
    o <- SSIM(img0, img3)
    lee2[i,7] <- o$SSIM
    lee2[i,8:10] <- o$comps

    # SSIM under alternative hypothesis
    o <- SSIM(img0, img3, alpha = cf[1], beta = cf[2], gamma = cf[3])
    lee2[i,11] <- o$SSIM
    lee2[i,12:14] <- o$comps

    # SSIM coefficients estimation (Kuan filter) ###############################
    img4 <- round(L * img4)
    fm <- SSIM.fit(img0, img4, maxiter = 50)
    cf <- fm$coefficients
    kuan[i,1:3] <- cf

    # gradient statistic to test
    z <- SSIM.test(fm)
    kuan[i,4] <- ifelse(z$stat > cutoff, 1, 0)
    kuan[i,5] <- z$stat
    kuan[i,6] <- z$logLik

    # SSIM under H0: alpha = beta = gamma = 1
    o <- SSIM(img0, img4)
    kuan[i,7] <- o$SSIM
    kuan[i,8:10] <- o$comps

    # SSIM under alternative hypothesis
    o <- SSIM(img0, img4, alpha = cf[1], beta = cf[2], gamma = cf[3])
    kuan[i,11] <- o$SSIM
    kuan[i,12:14] <- o$comps

    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  speed <- proc.time() - now

  percentage <- 100 * sum(ok[,4]) / Nsize
  pars <- ok[,1:3]
  colnames(pars) <- c("alpha","beta","gamma")
  stats <- ok[,5:6]
  colnames(stats) <- c("stat","logLik")
  H0 <- ok[,7:10]
  colnames(H0) <- c("SSIM","luminance","contrast","structure")
  H1 <- ok[,11:14]
  colnames(H1) <- c("SSIM","luminance","contrast","structure")
  out0 <- list(pars = pars, stats = stats, H0 = H0, H1 = H1, percentage = percentage)

  percentage <- 100 * sum(ok[,4]) / Nsize
  pars <- lee[,1:3]
  colnames(pars) <- c("alpha","beta","gamma")
  stats <- lee[,5:6]
  colnames(stats) <- c("stat","logLik")
  H0 <- lee[,7:10]
  colnames(H0) <- c("SSIM","luminance","contrast","structure")
  H1 <- lee[,11:14]
  colnames(H1) <- c("SSIM","luminance","contrast","structure")
  out1 <- list(pars = pars, stats = stats, H0 = H0, H1 = H1, percentage = percentage)

  percentage <- 100 * sum(ok[,4]) / Nsize
  pars <- lee2[,1:3]
  colnames(pars) <- c("alpha","beta","gamma")
  stats <- lee2[,5:6]
  colnames(stats) <- c("stat","logLik")
  H0 <- lee2[,7:10]
  colnames(H0) <- c("SSIM","luminance","contrast","structure")
  H1 <- lee2[,11:14]
  colnames(H1) <- c("SSIM","luminance","contrast","structure")
  out2 <- list(pars = pars, stats = stats, H0 = H0, H1 = H1, percentage = percentage)

  percentage <- 100 * sum(ok[,4]) / Nsize
  pars <- kuan[,1:3]
  colnames(pars) <- c("alpha","beta","gamma")
  stats <- kuan[,5:6]
  colnames(stats) <- c("stat","logLik")
  H0 <- kuan[,7:10]
  colnames(H0) <- c("SSIM","luminance","contrast","structure")
  H1 <- kuan[,11:14]
  colnames(H1) <- c("SSIM","luminance","contrast","structure")
  out3 <- list(pars = pars, stats = stats, H0 = H0, H1 = H1, percentage = percentage)

  out <- list(no.filter = out0, lee = out1, enhanced = out2, kuan = out3, cutoff = cutoff, speed = speed)
  out
}
