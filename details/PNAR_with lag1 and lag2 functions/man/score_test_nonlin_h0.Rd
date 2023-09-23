\name{Computation of quasi score test for testing linearity of the PNAR(1) model}
\alias{score_test_nonlin_h0}

\title{
Computation of quasi score test for testing linearity of the PNAR model versus
the nonlinear intercept drift alternative, with 1 lag
}

\description{
Computation of quasi score test for testing linearity of the PNAR model versus
the nonlinear intercept drift alternative, with 1 lag.
}

\usage{
score_test_nonlin_h0(ca, y, W)
}

\arguments{
\item{ca}{
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
}

\details{
This function computes the quasi score test for testing linearity of the PNAR
model versus the nonlinear intercept drift alternative, with 1 lag.
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
\link{poisson.MOD2}, \link{log_lin_estimnarpq}, \link{lin_estimnarpq}
}
}

\examples{
data(crime)
data(crime_W)
mod1 <- lin_estimnar1(crime, crime_W)
ca <- mod1$coeflin
score_test_nonlin_h0(ca, crime, crime_W)
}


