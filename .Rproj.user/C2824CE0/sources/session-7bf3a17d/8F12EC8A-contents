summary.DV <- function(object, ...) {
  a <- list()
  class(a) <- "htest"
  a$statistic <- object$statistic
  a$parameter <- object$parameter
  a$p.value <- object$p.value
  a$null.value <- object$null.value
  a$alternative <- object$alternative
  a$method <- object$method
  a$data.name <- object$data.name
  return(a)
}

print.summary.DV <- function(x, ...) {
  cat("Test for linearity of PNAR(p) versus the non-linear ST-PNAR(p)", "\n")
  cat("\n")
  cat("Test statistic value = ", x$statistic, "df = ", x$parameter, "p-value = ", x$p.value, "\n")
  cat("Alternative hypothesis: At least one coefficient of the non-linear component is not zero")
}

print.DV <- function(x, ...) {
  cat("Results: \n")
  a <- c(x$statistic, x$parameter, x$p.value)
  print(a)
}
