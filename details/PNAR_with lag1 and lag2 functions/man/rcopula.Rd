\name{Random number generation of copula functions}
\alias{rcopula}

\title{
Random number generation of copula functions
}

\description{
Random number generation of copula functions.
}

\usage{
rcopula(n, N, copula = "gaussian", corrtype = "equicorrelation",
rho, dof = 1, cholR = NULL)
}

\arguments{
\item{n}{
The number of random values to generate.
}
\item{N}{
The number of variables for which random values will be generated.
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
\item{cholR}{
If you have the Cholesky decomposition of the correlation matrix R pass
it here, otherwise leave it NULL.
}
}

\details{
This function generates random copula values.
}

\value{
An \eqn{n x N} matrix with the simulated copula values.
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
\code{\link{getN}, \link{poisson.MOD2}, \link{poisson.MOD2.log},
\link{poisson.MODpq}, \link{poisson.MODpq.log}
}
}

\examples{
u <- rcopula(n = 100, N = 50, rho = 0.3)
}


