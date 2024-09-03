\name{log_lin_ic_plot}
\alias{log_lin_ic_plot}

\title{
Scatter plot of information criteria versus the number of lags in the log-linear Poisson NAR(p) model with p lags and q
covariates (log-PNAR(p))
}

\description{
Scatter plot of information criteria versus the number of lags in log-linear Poisson Network Autoregressive model of order
\eqn{p} with \eqn{q} covariates (log-PNAR(\eqn{p})).
}

\usage{
log_lin_ic_plot(y, W, p = 1:10, Z = NULL, uncons = FALSE, ic = "QIC")
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
A vector with integer numbers, the range of lags in the model, for which the AIC,
BIC and QIC will be computed.
}
\item{Z}{
An \eqn{N} x \eqn{q} matrix of covariates (one for each column), where \eqn{q}
is the number of covariates in the model. Note that they must be non-negative.
}
\item{uncons}{
Logical, if TRUE an unconstrained optimization without stationarity
constraints is performed (default is FALSE).
}
\item{ic}{
The information criterion you want to plot, "QIC" (default value), "AIC" or "BIC".
}
}

\details{
The function computes the AIC, BIC or QIC for a range of lag orders of the
log-linear Poisson Network Autoregressive model of order \eqn{p} with \eqn{q}
covariates (PNAR(\eqn{p})).
}

\value{
A scatter plot with the lag order versus either QIC (default), AIC or BIC, and a vector
with their values, for each lag order.
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
\code{\link{log_lin_estimnarpq}, \link{lin_ic_plot}
}
}

\examples{
data(crime)
data(crime_W)
log_lin_ic_plot(crime, crime_W, p = 1:3)
}


