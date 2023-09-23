\name{Bootstrap test for testing linearity of PNAR(1) model}
\alias{score_test_stnar1_j}
\alias{score_test_tnar1_j}

\title{
Bootstrap sup-type test for testing linearity of PNAR model,
with 1 lag, versus the STNAR alternative or versus the TNAR alternative
}

\description{
Bootstrap sup-type test for testing linearity of PNAR model,
with 1 lag, versus the STNAR alternative or versus the TNAR alternative.
}

\usage{
score_test_stnar1_j(ca, y, W, J = 499, gama_L, gama_U, len = 10)
score_test_tnar1_j(ca, y, W, J = 499, gama_L, gama_U, len = 10)
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
\item{J}{
The number of bootstrap samples to draw.
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
The score_test_stnar1_j function performs a bootstrap sup-type test for testing
linearity of PNAR model, with 1 lag, versus the nonlinear Smooth Transition
model (STNAR) alternative.

The score_test_tnar1_j performs a bootstrap sup-type test for testing
linearity of PNAR model, with 1 lag, versus the nonlinear Threshold model
(TNAR) alternative.
}

\value{
A list including:
\item{supLM}{
The value value of the sup test statistic.
}
\item{aveLM}{
The value of the mean test statistic.
}
\item{pJ}{
The bootstrap p-value of the sup test.
}
\item{apJ}{
The bootstrap p-value for the mean test statistic.
}
\item{cpJ}{
The adjusted bootstrap p-value of the sup test statistic.
}
\item{capJ}{
The adjusted bootstrap corrected p-value for the mean test statistic.
}
\item{LM}{
The values of test statistic for different values in the grid of the
\eqn{\gamma} parameter.
}
\item{supLMj}{
The values of perturbed sup test statistic for different values in the grid of
the \eqn{\gamma} parameter.
}
\item{aveLMj}{
The values of perturbed mean test statistic for different values
in the grid of the \eqn{\gamma} values.
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
mod1 <- lin_estimnar1(crime, crime_W)
ca <- mod1$coeflin
score_test_stnar1_j(ca, crime, crime_W, J = 100, len = 10,
gama_L = 1, gama_U = 1.5)
}


