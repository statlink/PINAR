\name{Generation of a network from the Erdos-Renyi model (ER)}
\alias{adja_gnp}

\title{
Generation of a network from the Erdos-Renyi model (ER)
}

\description{
Generation of a network from the Erdos-Renyi model (ER).
}

\usage{
adja_gnp(N, alpha, directed = FALSE)
}

\arguments{
\item{N}{
The number of nodes on the network.
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
This function generates a network from the Erdos-Renyi model (ER). To do so,
for each couple of nodes it performs a Bernoulli trial with values 1
"draw an edge", 0 "otherwise". Each trial has the same probability of having
an edge specified based on the number of nodes on the network (N) and the
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
\code{\link{adja}, \link{Gmatr}
}
}

\examples{
W <- adja_gnp(20, 0.3)
}


