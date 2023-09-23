\name{Generation of the matrix of coefficients of the linear PNAR model}
\alias{Gmatr}

\title{
Generation of the matrix of coefficients of the linear PNAR model
}

\description{
Generation of the matrix of coefficients of the linear PNAR model.
}

\usage{
Gmatr(b0, b1, b2, W, N)
}

\arguments{
\item{b0}{
The intercept coefficient.
}
\item{b1}{
The network coefficient.
}
\item{b2}{
The autoregressive coefficient.
}
\item{W}{
The row-normalized weighted adjacency matrix describing the network.
}
\item{N}{
The number of nodes on the network.
}
}

\details{
The function generates the matrix of coefficients of the linear PNAR
model, called G matrix.
}

\value{
A list including:
\item{b0}{
A vector with the intercept coefficient.
}
\item{G}{
The matrix of coefficients of the linear PNAR model.
}
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
\code{\link{adja}, \link{adja_gnp}
}
}

\examples{

}


