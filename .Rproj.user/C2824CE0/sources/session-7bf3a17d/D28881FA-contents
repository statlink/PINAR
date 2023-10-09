\name{global_optimise_LM_stnarpq}
\alias{global_optimise_LM_stnarpq}

\title{
Optimization of the score test statistic for the ST-PNAR(p) model
}

\description{
Global optimization of the linearity test statistic for the Smooth Transition
Poisson Network Autoregressive model of order \eqn{p} with \eqn{q} covariates
(ST-PNAR(\eqn{p})) with respect to the nuisance scale parameter \eqn{\gamma}.
}

\usage{
global_optimise_LM_stnarpq(gama_L = NULL, gama_U = NULL, len = 10, b, y, W,
p, d, Z = NULL, tol = 1e-9)
}

\arguments{
\item{gama_L}{
The lower value of the \eqn{\gamma} values to consider.
Use NULL if there is not information about its value.
See the details for default computation.
}
\item{gama_U}{
The upper value of the \eqn{\gamma} values to consider.
Use NULL if there is not information about its value.
See the details for default computation.
}
\item{len}{
The number of increments to consider for the \eqn{\gamma} parameter.
}
\item{b}{
The estimated parameters from the linear model, in the following order:
(intercept, network parameters, autoregressive parameters, covariates).
The length of the vector should be \eqn{2p + 1 + q}, where \eqn{q} denotes the number of covariates.
}
\item{y}{
A \eqn{TT \times N} time series object or a \eqn{TT \times N} numerical matrix with the \eqn{N} multivariate count time series over \eqn{TT} time periods.
}
\item{W}{
The \eqn{N \times N} row-normalized non-negative adjacency matrix describing
the network. The main diagonal entries of the matrix should be zeros, all the
other entries should be non-negative and the maximum sum of elements over the
rows should equal one. The function row-normalizes the matrix if a non-normalized
adjacency matrix is provided.
}
\item{p}{
The number of lags in the model.
}
\item{d}{
The lag parameter of non-linear variable (should be between 1 and \eqn{p}).
}
\item{Z}{
An \eqn{N \times q} matrix of covariates (one for each column), where \eqn{q}
is the number of covariates in the model. Note that they must be non-negative.
}
\item{tol}{
Tolerance level for the optimizer.
}
}

\details{
The function optimizes the quasi score test statistic, under the null assumption of linearity, for testing linearity of Poisson Network Autoregressive model of order \eqn{p} against the following ST-PNAR(\eqn{p}) model, with respect to the unknown nuisance parameter (\eqn{\gamma}). For each node of the network \eqn{i=1,...,N} over the time sample \eqn{t=1,...,TT}
\deqn{
\lambda_{i,t}=\beta_{0}+\sum_{h=1}^{p}(\beta_{1h}X_{i,t-h}+\beta_{2h}Y_{i,t-h}+\alpha_{h}e^{-\gamma X_{i,t-d}^{2}}X_{i,t-h})+\sum_{l=1}^{q}\delta_{l}Z_{i,l}
}
where \eqn{X_{i,t}=\sum_{j=1}^{N}W_{ij}Y_{j,t}} is the network effect, i.e. the weighted average impact of node \eqn{i} connections, with the weights of the mean being \eqn{W_{ij}}, the single element of the network matrix \eqn{W}. The sequence \eqn{\lambda_{i,t}} is the expectation of \eqn{Y_{i,t}}, conditional to its past values.

The null hypothesis of the test is defined as \eqn{H_{0}: \alpha_{1}=...=\alpha_{p}=0}, versus the alternative that at least one among \eqn{\alpha_{h}} is not 0. The test statistic has the form
\deqn{
LM(\gamma)=S^{'}(\hat{\theta},\gamma)\Sigma^{-1}(\hat{\theta},\gamma)S(\hat{\theta},\gamma)
}
where
\deqn{
S(\hat{\theta},\gamma)=\sum_{t=1}^{TT}\sum_{i=1}^{N}\left(\frac{Y_{i,t}}{\lambda_{i,t}(\hat{\theta},\gamma)}-1\right) \frac{\partial\lambda_{i,t}(\hat{\theta},\gamma)}{\partial\alpha}
}
is the partition of the quasi score related to the vector of non-linear parameters \eqn{\alpha=(\alpha_{1},...,\alpha_{p})}, evaluated at the estimated parameters \eqn{\hat{\theta}} under the null assumption \eqn{H_{0}} (linear model) and \eqn{\Sigma(\hat{\theta},\gamma)} is the variance of \eqn{S(\hat{\theta},\gamma)}.

The optimization employes the Brent algorithm (Brent, 1973) applied in the interval from \code{gama_L} to \code{gama_U}. To be sure that the global optimum is found, the optimization is performed at (\code{len}-1) consecutive equidistant sub-intervals and then the maximum over them is taken as global optimum.

The values of \code{gama_L} and \code{gama_U} are computed internally as \code{gama_L} \eqn{=-\log(0.9)/X^{2}} and \code{gama_U} \eqn{=-\log(0.1)/X^{2}}, where \eqn{X} is the overall mean of \eqn{X_{i,t}} over the nodes \eqn{i=1,...,N} and times \eqn{t=1,...,TT}. Since the non-linear function \eqn{e^{-\gamma X_{i,t-d}^{2}}} ranges between 0 and 1, by considering \eqn{X} to be a representative value for the network mean, \code{gama_U} and \code{gama_L} would be the values of \eqn{\gamma} leading the non-linear switching function to be 0.1 and 0.9, respectively, so that in the optimization procedure the extremes of the function domain are excluded. Alternatively, their value can be supplied by the user. For details see Armillotta and Fokianos (2022b, Sec. 4-5).
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
A vector with the extremes points of sub-intervals.
}
}

\references{
Armillotta, M. and K. Fokianos (2022a). Poisson network autoregression.
\url{https://arxiv.org/abs/2104.06296}

Armillotta, M. and K. Fokianos (2022b). Testing linearity for network
autoregressive models.
\url{https://arxiv.org/abs/2202.03852}

Armillotta, M., Tsagris, M. and Fokianos, K. (2022c). The R-package PNAR
for modelling count network time series.
\url{https://arxiv.org/abs/2211.02582}

Brent, R. (1973) Algorithms for Minimization without Derivatives.
Prentice-Hall, Englewood Cliffs N.J.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{score_test_stnarpq_j}, \link{global_optimise_LM_tnarpq},
\link{score_test_tnarpq_j}
}
}

\examples{
data(crime)
data(crime_W)
mod1 <- lin_estimnarpq(crime, crime_W, p = 1)
b <- mod1$coefs[, 1]
global_optimise_LM_stnarpq(b = b, y = crime, W = crime_W, p = 1, d = 1)
}


