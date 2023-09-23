\name{Chicago crime dataset}
\alias{crime}
\docType{data}
\title{
Chicago crime dataset
}
\description{
Chicago crime for 72 months...
}
\usage{crime}
\format{
A matrix with 552 rows and 72 columns.
\describe{
\item{\code{date}}{a numeric vector}
\item{\code{ibm}}{a numeric vector}
\item{\code{sp}}{a numeric vector}
}
}

\source{
Where di you get the data from?
}

\references{
Tsay, R. S. (2010). \emph{Analysis of Financial Time Series}.
Hoboken, NJ: Wiley. Third edition.
}

\seealso{
\code{\link{crime_W}, \link{adja}, \link{lin_estimnar1}, \link{lin_estimnar2}
}
}

\examples{
data(crime)
data(crime_W)

mod1 <- lin_estimnar1(crime, crime_W)
mod2 <- log_lin_estimnar1(crime, crime_W)
}
\keyword{datasets}
