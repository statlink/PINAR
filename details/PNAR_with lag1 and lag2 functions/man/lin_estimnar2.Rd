\name{Constrained estimation of the linear PNAR(2) model}
\alias{lin_estimnar2}

\title{
Constrained estimation of the linear PNAR(2) model
}

\description{
Constrained estimation of the linear PNAR(2) model.
}

\usage{
lin_estimnar2(y, W)
}

\arguments{
\item{y}{
An \eqn{N \times TT} numerical matrix with the \eqn{N} multivariate
count time series over \eqn{TT} time periods.
}
\item{W}{
The \eqn{N \times N} adjacency matrix describing the network.
This must have non-negative values. It may have 0s and 1s or
any positive number. The function row-normalizes the matrix.
}
}

\details{
This function performs a constrained estimation of the Poisson
linear PNAR(2) model.
}

\value{
A list including:
\item{coeflin}{
The vector of estimated QMLE coefficients.
}
\item{selin}{
The standard errors of the estimated coefficients.
}
\item{tlin}{
The t-test statistics of the estimates coefficients.
}
\item{score}{
The value of the score function at the optimization point.
}
\item{loglik}{
The value of the maximized log-likelihood.
}
\item{aic_lin}{
The Akaike information criterion (AIC).
}
\item{bic_lin}{
The Bayesian information criterion (BIC).
}
\item{qic_lin}{
The Quasi information criterion (QIC).
}
}

\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence
in time series. \emph{Technometrics}, 159, 262-3270.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{lin_estimnar1}, \link{lin_estimnarpq}, \link{log_lin_estimnarpq}
}
}

\examples{
W <- adja_gnp(20, 0.3)
y <- matrix(rpois(20 * 100, 10), ncol = 100)
mod2 <- lin_estimnar2(y, W)
}


