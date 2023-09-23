\name{Count the number of events before a specified time}
\alias{getN}

\title{
Count the number of events before a specified time
}

\description{
Count the number of events before a specified time.
}

\usage{
getN(x, tt = 1)
}

\arguments{
\item{x}{
A matrix of (positive) inter-event times.
}
\item{tt}{
A positive time.
}
}

\details{
}

\value{
The number of events before time "tt"" (possibly 0), for each column of x.
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
\code{\link{rcopula}, \link{poisson.MOD2}, \link{poisson.MOD2.log}
}
}

\examples{
x <- rcopula(n = 100, N = 50, rho = 0.3, dof = 1)
getN(x)
}


