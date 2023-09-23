\name{Simulation of counts from a Poisson log-linear PNAR(p) model}
\alias{poisson.MODpq.log}

\title{
Simulation of counts from a Poisson log-linear PNAR(p) model with q non
time-varying covariates
}

\description{
Simulation of counts from a Poisson log-linear PNAR(p) model with q non
time-varying covariates.
}

\usage{
poisson.MODpq.log(b, W, p, Z = NULL, TT, N, copula = "gaussian",
corrtype = "equicorrelation", rho, dof = 1)
}

\arguments{
\item{b}{
The coefficients of the model in the following order:
(intercept, p network effects, p AR effects, covariates).
}
\item{W}{
The \eqn{N \times N} adjacency matrix describing the network.
This must have non-negative values. It may have 0s and 1s or
any positive number. The function row-normalizes the matrix.
}
\item{p}{
The number of lags in the model.
}
\item{Z}{
An \eqn{N  \times q} matrix of covariates (one for each column),
where \eqn{q} is the number of covariates in the model.
Note that they must be non-negative.
}
\item{TT}{
The temporal sample size.
}
\item{N}{
The number of nodes on the network.
}
\item{copula}{
Which copula function to use? The "normal", "t", or "clayton"?
}
\item{corrtype}{
The covariance matrix will either be the "compund" or "toeplitz".
}
\item{rho}{
The the value of the correlation of the copula.
}
\item{dof}{
The degrees of freedom for Student's t copula.
}
}

\details{
This function generates counts from a Poisson log-linear PNAR(p) model, where
q non time-varying covariates are allowed as well.
}

\value{
A list including:
\item{p2R}{
The Toeplitz correlation matrix, if employed in the copula or NULL else.
}
\item{lambda}{
An \eqn{N x TT} matrix of generated Poisson means for N time series over TT.
}
\item{y}{
An \eqn{N x TT} matrix of generated counts for N time series over TT.
}
}

\references{
Edelmann, D, K. Fokianos. and M. Pitsillou. (2019). An Updated Literature
Review of Distance Correlation and Its Applications to TT Series.
\emph{International Statistical Review}, 87, 237-262.

Fokianos K. and M. Pitsillou (2017). Consistent testing for pairwise dependence
in time series. \emph{Technometrics}, 159, 262-3270.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{poisson.MODpq}, \link{poisson.MOD2.stnar}, \link{log_lin_estimnarpq}
}
}

\examples{

}


