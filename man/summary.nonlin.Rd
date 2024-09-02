\name{summary.nonlin}
\alias{summary.nonlin}
\alias{print.summary.nonlin}
\alias{print.nonlin}

\title{S3 methods for extracting the results of the non-linear hypothesis test
}

\description{
S3 methods for extracting the results of the non-linear hypothesis test.
}

\usage{
\method{summary}{nonlin}(object, ...)
\method{print}{summary.nonlin}(x, ...)
\method{print}{nonlin}(x, ...)
}

\arguments{
\item{object}{
An object containing the results of the function \code{\link{score_test_nonlinpq_h0}}.
}
\item{x}{
An object containing the results of the function \code{\link{score_test_nonlinpq_h0}}.
}
\item{...}{
Extra arguments the user can pass.
}
}

\details{
The functions print the output of the non-linear hypothesis test.
}

\value{
The functions print the results of the function \code{\link{score_test_nonlinpq_h0}}.
}

\references{
Armillotta, M., Tsagris, M. and Fokianos, K. (2024b). Inference for Network Count Time Series with the R Package PNAR. \emph{The R Journal}, 15/4: 255--269.
}

\author{
Mirko Armillotta, Michail Tsagris and Konstantinos Fokianos.
}

\seealso{
\code{\link{score_test_nonlinpq_h0}
}
}

\examples{
data(crime)
data(crime_W)
mod1 <- lin_estimnarpq(crime, crime_W, p = 2)
ca <- mod1$coefs[, 1]
a <- score_test_nonlinpq_h0(ca, crime, crime_W, p = 2, d = 1)
print(a)
summary(a)
}

