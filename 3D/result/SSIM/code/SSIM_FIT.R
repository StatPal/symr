## $ID: SSIM_FIT.R, last updated 2018/07/30, F. Osorio and W. Barraza

## SSIM parameter estimation using a heteroscedastic regression model

# arguments:
# x, y        nrow-by-ncol images to be compared
# size        subimages window of size x size
# numIter     maximun number of iterations
# tol         relative tolerance for the iterative algorithm

SSIM.fit <- function(x, y, size = 16, maxiter = 200, tol = 1e-6)
{
  # local functions
  luminance <- function(x, y) {
    c1 <- (255 * 0.01)^2
    xbar <- mean(x)
    ybar <- mean(y)
    (2 * xbar * ybar + c1) / (xbar^2 + ybar^2 + c1)
  }

  contrast <- function(x, y) {
    c2 <- (255 * 0.03)^2
    xvar <- var(x)
    yvar <- var(y)
    (2 * sqrt(xvar * yvar) + c2) / (xvar + yvar + c2)
  }

  structure <- function(x, y) {
    c3 <- .5 * (255 * 0.03)^2
    xvar <- var(x)
    yvar <- var(y)
    (cov(x,y) + c3) / (sqrt(xvar * yvar) + c3)
  }

  MSE <- function(x, y) mean((x - y)^2)

  ssim.model <- function(x, y, n)
  {
    subfig <- function(x, n) {
      i <- (row(x) - 1L) %/% n + 1
      j <- (col(x) - 1L) %/% n + 1
      b <- i + (j - 1L) * max(i)
      fig <- split(x, b)
      names(fig) <- NULL
      fig
    }

    sub.x <- subfig(x, n)
    sub.y <- subfig(y, n)

    lxy  <- mapply(luminance, x = sub.x, y = sub.y)
    cxy  <- mapply(contrast,  x = sub.x, y = sub.y)
    sxy  <- mapply(structure, x = sub.x, y = sub.y)
    mse  <- mapply(MSE, x = sub.x, y = sub.y)

    npos <- which((lxy <= 0) | (cxy <= 0) | (sxy <= 0) | (mse <= 0))
    nsum <- sum((lxy <= 0) | (cxy <= 0) | (sxy <= 0) | (mse <= 0))

    if (nsum) {
      lxy <- lxy[-npos]
      cxy <- cxy[-npos]
      sxy <- sxy[-npos]
      mse <- mse[-npos]
    }
    rmse <- sqrt(mse)

    X <- cbind(lxy, cxy, sxy) # luminance = lxy, contrast = cxy, structure = sxy
    z <- cbind(X, rmse)
    colnames(z) <- c("luminance","contrast","structure","RMSE")
    z <- as.data.frame(z)
    z
  }

  mean.fnc <- function(lumin, contr, struct, cf) {
    prod <- (lumin^cf$alpha) * (contr^cf$beta) * (struct^cf$gamma)
    prod
  }

  # pre-processing stage
  now <- proc.time()
  db <- ssim.model(x, y, n = size)
  stage1 <- proc.time() - now
  dims <- c(nrow(db), ncol(db) - 1)

  # initialization
  now <- proc.time()
  fm <- lm(-log(RMSE) ~ log(luminance) + log(contrast) + log(structure), data = db)
  coef <- fm$coef
  phi <- exp(coef[1])
  if (phi < 1)
    phi <- .5 * (1 + sqrt(5))
  coef <- rep(1, 3)
  start <- c(coef, phi)
  z <- 1 / db$RMSE
  control <- c(maxiter, tol, 0)

  # call C routine
  fit <- .C("SSIM_fit",
            z = as.double(z),
            luminance = as.double(db$luminance),
            contrast  = as.double(db$contrast),
            structure = as.double(db$structure),
            dims = as.integer(dims),
            coefficients = as.double(coef),
            scale = as.double(phi),
            control = as.double(control),
            logLik = double(1))

  stage2 <- proc.time() - now
  speed <- rbind(stage1, stage2)[,1:3]
  colnames(speed) <- c("user","system","elapsed")
  rownames(speed) <- c("preprocessing","estimation")

  o <- list(model = db, coefficients = fit$coefficients, scale = fit$scale,
            logLik = fit$logLik, start = start, iterations = fit$control[3],
            speed = speed)
  names(o$coefficients) <- c("luminance","contrast","structure")
  class(o) <- "SSIM.nls"
  o
}

print.SSIM.nls <-
function(x, digits = 4, ...)
{
  cat("\nCoefficients:\n ")
  print(format(round(x$coefficients, digits = digits)), quote = F, ...)
  cat("\nIterations:", x$iterations, "\n")
  invisible(x)
}
