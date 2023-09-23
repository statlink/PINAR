\name{Generation of a network from the Stochastic Block Model (SBM)}
\alias{adja}

\title{
Generation of a network from the Stochastic Block Model (SBM)
}

\description{
Generation of a network from the Stochastic Block Model (SBM).
}

\usage{
adja(N, K, alpha, directed = FALSE)
}

\arguments{
\item{N}{
The number of nodes on the network.
}
\item{K}{
The number of blocks.
}
\item{alpha}{
The network density.
}
\item{directed}{
This generates by default an undirected network. If you want a
directed network set this equal to TRUE.
}
}

\details{
This function generates a network from the Stochastic Block Model (SBM).
To do so, for each couple of nodes it performs a Bernoulli trial with values
1 "draw an edge", 0 "otherwise". The probabilities of these trials are bigger
if the two nodes are in the same block, lower otherwise, and they are specified
based on the number of nodes on the network (N), the number of blocks (K) and
network density (alpha).
}

\value{
A row-normalized weighted adjacency matrix describing the network.
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
\code{\link{adja_gnp}, \link{Gmatr}
}
}

\examples{
W <- adja(20, K = 2, alpha = 0.3)
}


