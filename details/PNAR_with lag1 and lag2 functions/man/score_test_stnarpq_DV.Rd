\name{Davies bound of sup-type test for testing linearity of PNAR(p) versus the STNAR alternative}
\alias{score_test_stnarpq_DV}

\title{
Davies bound of sup-type test for testing linearity of PNAR(p) versus the
STNAR alternative
}

\description{
Davies bound of sup-type test for testing linearity of PNAR(p) versus the
STNAR alternative.
}

\usage{
score_test_stnarpq_DV(b, y, W, p, d, Z = NULL, gama_L, gama_U, len = 100)
}

\arguments{
\item{b}{
The estimated parameters from the linear model.
}
\item{y}{
An \eqn{N \times TT} numerical matrix with the multivariate count time series,
where the \eqn{N} denotes the number of nodes in the network and \eqn{TT}
denotes the temporal size of the series.
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
\item{gama_L}{
The lower value of the \eqn{\gamma} values to consider.
}
\item{gama_U}{
The upper value of the \eqn{\gamma} values to consider.
}
\item{len}{
The length of the grid of values of \eqn{\gamma} values to consider.
}
}

\details{
The function computes the Davies bound of sup-type test for testing
linearity of PNAR(p) model versus the nonlinear Smooth Transition model
(STNAR) alternative. Non time-varying covariates may also be included.
}

\value{
A list including:
\item{DV}{
The Davies bound of p-values for sup test.
}
\item{supLM}{
The value of the sup test statistic in the sample y.
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


