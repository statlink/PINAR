summary.nonlin <- function(object, ...) {
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

print.summary.nonlin <- function(x, ...) {
  cat("Linearity test against non-linear ID-PNAR(p) model", "\n")
  cat("\n")
  cat("Test statistic value = ", x$statistic, "df = ", x$parameter, "p-value = ", x$p.value, "\n")
  cat("Alternative hypothesis: True gamma parameter is greater than zero")
}

print.nonlin <- function(x, ...) {
  cat("Results: \n")
  a <- c(x$statistic, x$parameter, x$p.value)
  print(a)
}
