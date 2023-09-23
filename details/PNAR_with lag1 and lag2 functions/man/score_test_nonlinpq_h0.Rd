\name{Quasi score test for testing linearity of the nonlinear (ID) PNAR(p) model}
\alias{score_test_nonlinpq_h0}

\title{
Quasi score test of the nonlinear Intercept Drift (ID) PNAR(p)
model under the null assumption of linearity
}

\description{
Quasi score test of the nonlinear Intercept Drift (ID) PNAR(p)
model under the null assumption of linearity.
}

\usage{
score_test_nonlinpq_h0(b, y, W, p, d, Z = NULL)
}

\arguments{
\item{b}{
The estimated parameters from the linear model.
}
\item{y}{
An \eqn{N \times TT} numerical matrix with the multivariate count time series,
where the \eqn{N} denotes the number of nodes in the network and \eqn{TT} denotes
the temporal size of the series.
}
\item{W}{
The \eqn{N \times N} adjacency matrix describing the network.
This must have non-negative values. It may have 0s and 1s or
any positive number. The function row-normalizes the matrix.
}
\item{p}{
The number of lags in the model.
}
\item{d}{
The lag parameter of nonlinear variable (should be between 1 and p).
}
\item{Z}{
An \eqn{N  \times q} matrix of covariates (one for each column),
where \eqn{q} is the number of covariates in the model.
Note that they must be non-negative.
}
}

\details{
The function computes the quasi score test of the nonlinear Intercept Drift
(ID) PNAR(p) model under the null assumption of linearity. q non time-varying
covariates may also be included.
}

\value{
A list including:
\item{stat}{
The value of the test statistic.
}
\item{pvalue}{
The p-value of the test statistic.
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
\code{\link{score_test_stnar1_j}, \link{score_test_tnar1_j},
\link{log_lin_estimnarpq}, \link{lin_estimnarpq}
}
}

\examples{
data(crime)
data(crime_W)
mod1 <- lin_estimnar1(crime, crime_W)
ca <- mod1$coeflin
score_test_nonlinpq_h0(ca, crime, crime_W, p = 2, d = 1)
}


