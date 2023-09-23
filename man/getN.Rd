\name{getN}
\alias{getN}

\title{
Count the number of events within a specified time
}

\description{
This function counts the number of events within a specified time.
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

\value{
The number of events within time \code{tt} (possibly 0), for each column of
\code{x}.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{rcopula}, \link{poisson.MODpq}, \link{poisson.MODpq.log}
}
}

\examples{
x <- rcopula(n = 100, N = 50, rho = 0.3)
getN(x)
}



