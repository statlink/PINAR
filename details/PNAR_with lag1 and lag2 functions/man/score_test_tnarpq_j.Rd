\name{Bootstrap p-value for the sup-type test for testing linearity of PNAR(p) model II}
\alias{score_test_tnarpq_j}

\title{
Bootstrap p-value for the sup-type test for testing linearity of PNAR(p) model
versus the nonlinear TPNAR alternative
}

\description{
Bootstrap p-value for the sup-type test for testing linearity of PNAR(p) model
versus the nonlinear TPNAR alternative.
}

\usage{
score_test_tnarpq_j(supLM, b, y, W, p, d, Z = NULL, J = 499,
gama_L = 0.01, gama_U = 5, tol = 1e-9, ncores = 1)
}

\arguments{
\item{supLM}{
The sup-test statistic.
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
\item{J}{
The number of bootstrap samples to draw.
}
\item{gama_L}{
The lower value of the \eqn{\gamma} values to consider.
}
\item{gama_U}{
The upper value of the \eqn{\gamma} values to consider.
}
\item{tol}{
The tolerance level for the otimizer.
}
\item{ncores}{
The number of cores to use for parallel computing. By default the number of
cores is set to 1 (no parallel computing). \bold{Note:} If for some reason the
parallel does not work then load the doParallel package yourseleves
(\bold{library(doParallel)}).
}
}

\details{
The function computes a bootstrap p-value for the sup-type test for testing
linearity of PNAR(p) model versus the nonlinear Threshold model
(TPNAR) alternative. Non time-varying covariates may also be included.
}

\value{
A list including:
\item{spj}{
The bootstrap p-value of the sup test.
}
\item{pJ}{
The adjusted version of bootstrap p-value of the sup test.
}
\item{gamaj}{
The optimal values of gamma parameter for score test boostrap replications.
}
\item{supLMj}{
The values of perturbed test statistic at the optimum point gamaj.
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


