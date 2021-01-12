## $ID: SSIM_TEST.R, last updated 2018/07/24, F. Osorio

SSIM.test <- function(object)
{
  # copying elements from 'object'
  z <- 1 / object$model$RMSE
  luminance <- object$model$luminance
  contrast  <- object$model$contrast
  structure <- object$model$structure
  dims <- dim(object$model)
  dims[2] <- dims[2] - 1
  coef <- object$coef
  scale <- .5 * (1 + sqrt(5))
  control <- c(1.,1.e-6)

  # call C routine
  test <- .C("gradient_test",
             z = as.double(z),
             luminance = as.double(luminance),
             contrast  = as.double(contrast),
             structure = as.double(structure),
             dims = as.integer(dims),
             coefficients = as.double(coef),
             scale = as.double(scale),
             control = as.double(control),
             score = double(3),
             stat = double(1),
             logLik = double(1))

  df <- 3
  pval <- pchisq(test$stat, df, lower.tail = FALSE)
  o <- list(stat = test$stat, df = 3, p.value = pval, coefficient = coef,
            logLik = test$logLik, score = test$score)
  class(o) <- "SSIM.test"
  return(o)
}

print.SSIM.test <- function(x, digits = 4, ...)
{
  cat("\n")
  cat("Gradient test to assess the SSIM balance\n")
  cat("\n")
  cat("gradient statistic:", format(round(x$stat, digits = digits)), "on",
      format(round(x$df, digits = digits)), "degrees of freedom, p-value:",
      format(round(x$p.value, digits = digits)), "\n")
  cat("alternative hypothesis: SSIM coefficients are not equal to 1\n")
  cat("\nEstimated coefficients:\n ")
  print(format(round(x$coefficient, digits = digits)), quote = F, ...)
  cat("\n")
  invisible(x)
}
