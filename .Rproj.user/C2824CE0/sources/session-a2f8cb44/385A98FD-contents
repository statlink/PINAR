lin_ic_plot <- function(y, W, p = 1:10, Z = NULL, uncons = FALSE) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  if ( !is.null(Z) ) {
    if ( min(Z) < 0 ) {
      stop('The matrix of covariates Z contains negative values.')
    }
    Z <- model.matrix(~., as.data.frame(Z))
    Z <- Z[1:dim(y)[2], -1, drop = FALSE]
  }

  y <- t(y)
  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]
  z <- W %*% y

  ic <- matrix(nrow = length(p), ncol = 3)

  for ( i in 1:length(p) ) {

    wy <- NULL
    for ( ti in (p[i] + 1):TT )  wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p[i])], y[, (ti - 1):(ti - p[i])], Z) )
    wy <- cbind(1, wy)

    XX <- crossprod(wy)
    Xy <- Rfast::eachcol.apply(wy, as.vector( y[, -c(1:p[i])] ) )
    x0 <- solve(XX, Xy)
    x0[x0 < 0] <- 0.001

    m <- length(x0)
    # Lower and upper bounds (positivity constraints)
    lb <- rep(0, m)
    ub <- rep(Inf, m)

    # algorithm and relative tolerance
    opts <- list( "algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-8 )

    if ( uncons ) {

      s_qmle <- nloptr::nloptr(x0 = x0, eval_f = .logl_linpq, eval_grad_f = .scor_linpq, lb = lb,
                               ub = ub, opts = opts, N = N, TT = TT, y = y, W = W, wy = wy, p = p[i], Z = Z)
    } else {
      # Inequality constraints (parameters searched in the stationary region)
      # b are the parameters to be constrained
      constr <- function(b, N, TT, y, W, wy, p, Z) {
        con <- sum( b[2:(2 * p + 1)] ) - 1
        return(con)
      }
      # Jacobian of constraints
      # b are the parameters to be constrained
      j_constr <- function(b, N, TT, y, W, wy, p, Z) {
        j_con <- rep(1, m)
        j_con[1] <- 0
        return(j_con)
      }

      s_qmle <- nloptr::nloptr(x0 = x0, eval_f = .logl_linpq, eval_grad_f = .scor_linpq,
                               lb = lb, ub = ub, eval_g_ineq = constr, eval_jac_g_ineq = j_constr,
                               opts = opts, N = N, TT = TT, y = y, W = W, wy = wy, p = p[i], Z = Z)
    }

    coeflin <- s_qmle$solution

    ic[i, 1] <- aic_lins <- 2 * m + 2 * s_qmle$objective
    ic[i, 2] <- log(TT) * m + 2 * s_qmle$objective
    ola <- .scor_hess_outer_linpq(coeflin, N, TT, y, wy, p[i], Z)
    S_lins <- ola$scor
    H_lins <- ola$hh
    G_lins <- ola$out
    solveH_lins <- solve(H_lins)
    V_lins <- solveH_lins %*% G_lins %*% solveH_lins
    ic[i, 3] <- 2 * sum( H_lins * V_lins ) + 2 * s_qmle$objective

  }

  plot(p, ic[, 1], type = "b", xlab = "Lag", ylab = "AIC", cex.lab = 1.3, cex.axis = 1.3, xaxt = "n")
  abline(v = p, col = "lightgrey", lty = 2)
  abline(h = seq(min(ic[, 1]), max(ic[, 1]), length.out = length(p) ), lty = 2, col = "lightgrey")
  points(p, ic[, 1], pch = 16, col = "blue")
  mtext(text = p, side = 1, at = p, las = 1, font = 2, line = 0.7 )

  dev.new()

  plot(p, ic[, 2], type = "b", xlab = "Lag", ylab = "BIC", cex.lab = 1.3, cex.axis = 1.3, xaxt = "n")
  abline(v = p, col = "lightgrey", lty = 2)
  abline(h = seq(min(ic[, 2]), max(ic[, 2]), length.out = length(p) ), lty = 2, col = "lightgrey")
  points(p, ic[, 2], pch = 16, col = "blue")
  mtext(text = p, side = 1, at = p, las = 1, font = 2, line = 0.7 )

  dev.new()

  plot(p, ic[, 3], type = "b", xlab = "Lag", ylab = "QIC", cex.lab = 1.3, cex.axis = 1.3, xaxt = "n")
  abline(v = p, col = "lightgrey", lty = 2)
  abline(h = seq(min(ic, 3), max(ic[, 3]), length.out = length(p) ), lty = 2, col = "lightgrey")
  points(p, ic[, 3], pch = 16, col = "blue")
  mtext(text = p, side = 1, at = p, las = 1, font = 2, line = 0.7 )

  colnames(ic) <- c("AIC", "BIC", "QIC")
  rownames(ic) <- paste("Lag=", p, sep = "")
  ic
}
