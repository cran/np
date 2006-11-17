\name{w}
\title{where}
\usage{
w(x)
}
\alias{w}
\arguments{
  \item{x}{an array of logical values}
}
\value{
  an array of indices where \code{x == TRUE}.
}
\description{
  \code{w} computes indices where \code{x == TRUE}.
}
\details{
  \code{w} accepts an array of type \code{logical} and returns
  an array of indices where \code{x == TRUE}. If \code{x} contains no
  \code{TRUE} entries, \code{NULL} is returned.
}
\seealso{
  \code{\link{match}}
}
\author{
  Tristen Hayfield \email{hayfietj@mcmaster.ca}, Jeffrey S. Racine \email{racinej@mcmaster.ca}
}

\examples{
x <- rnorm(100)
w(x > 1.0)
}

\keyword{manip}
