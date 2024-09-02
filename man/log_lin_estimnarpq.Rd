\name{log_lin_estimnarpq}
\alias{log_lin_estimnarpq}

\title{
Estimation of the log-linear Poisson NAR(p) model with p lags and q
covariates (log-PNAR(p))
}

\description{
Estimation of the log-linear Poisson Network Autoregressive model of order
\eqn{p} with \eqn{q} covariates (log-PNAR(\eqn{p})).
}

\usage{
log_lin_estimnarpq(y, W, p, Z = NULL, uncons = FALSE, init = NULL,
xtol_rel = 1e-8, maxeval = 100)
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
is the number of covariates in the model.
}
\item{uncons}{
logical, if TRUE an unconstrained optimization is performed (default is FALSE).
}
\item{init}{
A vector of starting values for the optimization algorithm. If this is NULL,
the function computes them internally.
}
\item{xtol_rel}{
The stopping tolerance of the optimization algorithm.
}
\item{maxeval}{
The maximum number of evalutions the optimization algorithm will perform.
}
}

\details{
This function performs a constrained estimation of the linear Poisson NAR(\eqn{p}) model with \eqn{q} non-negative valued covariates, for each node of the network \eqn{i=1,...,N} over the time sample \eqn{t=1,...,TT}, defined as
\deqn{
\nu_{i,t}=\beta_{0}+\sum_{h=1}^{p}(\beta_{1h}X_{i,t-h}+\beta_{2h}Y_{i,t-h})+\sum_{l=1}^{q}\delta_{l}Z_{i,l},
}
where \eqn{X_{i,t}=\sum_{j=1}^{N}W_{ij}Y_{j,t}} is the network effect, i.e. the weighted average impact of node \eqn{i} connections, with the weights of the mean being \eqn{W_{ij}}, the single element of the network matrix \eqn{W}. The sequence \eqn{\nu_{i,t}} is the log of the expectation of \eqn{Y_{i,t}}, conditional to its past values. The parameter \eqn{\beta_{0}} is the intercept of the model, \eqn{\beta_{1h}} are the network coefficients, \eqn{\beta_{2h}} are the autoregressive parameters, and \eqn{\delta_{l}} are the coefficients assocciated to the covariates \eqn{Z_{i,l}}.

The estimation of the parameters of the model is performed by Quasi Maximum Likelihood Estimation (QMLE), maximizing the following quasi log-likelihood
\deqn{
l(\theta)=\sum_{t=1}^{TT}\sum_{i=1}^{N}\left[Y_{i,t}\nu_{i,t}(\theta)-e^{\nu_{i,t}(\theta)}\right]
}
with respect to the vector of unknown parameters \eqn{\theta} described above.

By default, the optimization is constrained in the stationary region where \eqn{\sum_{h=1}^{p}(|\beta_{1h}|+|\beta_{2h}|)<1}; this can be removed by setting \code{uncons = TRUE}. However, the model estimates might be inconsistent if the estimated parameters lie outside the stationary region.

The ordinary least squares estimates are employed as starting values of the optimization procedure. Robust standard errors and z-tests are also returned.
}

\value{
A list with attribute class "PNAR" including:
\item{coefs}{
A matrix with the estimated QMLE coefficients, their standard errors, their Z-test statistics and the relevant p-values computed via the standard normal approximation.
}
\item{score}{
The value of the quasi score function at the optimization point.
It should be close to 0 if the optimization is successful.
}
\item{loglik}{
The value of the maximized quasi log-likelihood.
}
\item{ic}{
A vector with the Akaike information criterion (AIC), the Bayesian information criterion (BIC) and the Quasi information criterion (QIC).
}
Alternatively, these can be printed via the function \code{\link{summary.PNAR}}.
}

\references{
Armillotta, M. and K. Fokianos (2024a). Count network autoregression. \emph{Journal of Time Series Analysis}, 45(4): 584--612.

Armillotta, M. and K. Fokianos (2023). Nonlinear network autoregression. \emph{Annals of Statistics}, 51(6): 2526--2552.

Armillotta, M., Tsagris, M. and Fokianos, K. (2024b). Inference for Network Count Time Series with the R Package PNAR. \emph{The R Journal}, 15/4: 255--269.
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
mod1 <- log_lin_estimnarpq(crime, crime_W, p = 2)
summary(mod1)
}


