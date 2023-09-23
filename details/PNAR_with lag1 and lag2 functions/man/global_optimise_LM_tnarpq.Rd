\name{Optimization of the score test statistic under null hypothesis for the TNAR model}
\alias{global_optimise_LM_tnarpq}

\title{
Optimization of the score test statistic under null hypothesis for the TNAR model
}

\description{
Optimization of the score test statistic under null hypothesis for the TNAR model.
}

\usage{
global_optimise_LM_tnarpq(gama_L, gama_U, len = 100, b, y, W,
p, d, Z = NULL, tol = 1e-9)
}

\arguments{
\item{gama_L}{
The lower value of the \eqn{\gamma} values to consider.
}
\item{gama_U}{
The upper value of the \eqn{\gamma} values to consider.
}
\item{len}{
The number of increments to consider for the \eqn{\gamma} parameter.
}
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
\item{tol}{
The tolerance level for the otimizer.
}
}

\details{
The function optimizes the score test statistic under the null assumption
with respect to the unknown nuisance parameter (gamma). The optimization
employes the Brent algorithm applied on the interval [gama_L, gama_U].
To be sure that the a global optimum is found, the optimization is performed
at (len-1) consecutive equidistant sub-intervals and then the minimum/maximum
over them is taken as global optimum.
}

\value{
A list including:
\item{gama}{
The optimum value of the \eqn{\gamma} parameter.
}
\item{supLM}{
The value of the objective function at the optimum.
}
\item{int}{
A vector with the extremes points of subintervals.
}
}

\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to Time Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise
dependence in time series. \emph{Technometrics}, 159, 262-3270.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{score_test_tnar1_j}, \link{score_test_nonlin_h0},
\link{log_lin_estimnarpq}
}
}

\examples{
data(crime)
data(crime_W)
mod1 <- lin_estimnarpq(crime, crime_W, p = 2)
b <- mod1$coeflin
score_test_tnarpq_j(5, b, crime, crime_W, p = 2, d = 1, J = 5)
}


