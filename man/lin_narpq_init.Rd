\name{lin_narpq_init}
\alias{lin_narpq_init}

\title{
Starting values for the linear Poisson NAR(p) model model with p lags and q
covariates (PNAR(p))
}

\description{
Starting values for the linear Poisson Network Autoregressive model of order
\eqn{p} with \eqn{q} covariates (PNAR(\eqn{p})).
}

\usage{
lin_narpq_init(y, W, p, Z = NULL)
}

\arguments{
\item{y}{
A \eqn{TT} x \eqn{N} time series object or a \eqn{TT} x \eqn{N} numerical matrix with the \eqn{N} multivariate
count time series over \eqn{TT} time periods.
}
\item{W}{
The \eqn{N} x \eqn{N} row-normalized non-negative adjacency matrix describing
the network. The main diagonal entries of the matrix should be zeros, all the
other entries should be non-negative and the maximum sum of elements over the
rows should equal one. The function row-normalizes the matrix if a non-normalized
adjacency matrix is provided.
}
\item{p}{
The number of lags in the model.
}
\item{Z}{
An \eqn{N} x \eqn{q} matrix of covariates (one for each column), where \eqn{q}
is the number of covariates in the model. Note that they must be non-negative.
}
}

\details{
The function computes starting values to be used in the function \code{\link{lin_estimnarpq}}.
These are simply the ordinary least squares estimators with a correction.
If any of the the resulting coefficients is negative they become equal to 0.001
}

\value{
A vector with the initial values.
}

\references{
Armillotta, M. and K. Fokianos (2024). Count network autoregression. \emph{Journal of Time Series Analysis}, 45(4): 584--612.

Armillotta, M. and K. Fokianos (2023). Nonlinear network autoregression. \emph{Annals of Statistics}, 51(6): 2526--2552.

Armillotta, M., Tsagris, M. and Fokianos, K. (2024). Inference for Network Count Time Series with the R Package PNAR. \emph{The R Journal}, 15/4: 255--269.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{lin_estimnarpq}
}
}

\examples{
data(crime)
data(crime_W)
x0 <- lin_narpq_init(crime, crime_W, p = 2)
}


