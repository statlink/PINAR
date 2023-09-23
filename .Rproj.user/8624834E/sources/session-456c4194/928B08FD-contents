summary.PNAR <- function(object, ...) {
  a <- list()
  class(a) <- "summary.PNAR"
  a$coefs <- object$coefs
  a$score <- object$score
  a$loglik <- object$loglik
  a$ic <- object$ic
  return(a)
}

print.summary.PNAR <- function(x, ...) {
  cat("Coefficients: \n")
  print(x$coefs)
  cat("---", "\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  cat("\n")
  cat("Log-likelihood value: ", x$loglik, "\n")
  cat("AIC:", x$ic[1], "BIC:", x$ic[2], "QIC:", x$ic[3])
}

print.PNAR <- function(x, ...) {
  cat("Coefficients: \n")
  a <- x$coefs[, 1]
  names(a) <- rownames(x$coefs)
  print(a)
}
