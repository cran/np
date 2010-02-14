% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npunitest}
\alias{npunitest}
\title{ Kernel Consistent Density Equality Test with Mixed Data Types}
\description{
  \code{npunitest} implements the consistent metric entropy test of
  Maasoumi and Racine (2002) for two arbitrary, stationary
  univariate nonparametric densities on common support.
}
\usage{
npunitest(data.x = NULL,
          data.y = NULL,
          bootstrap = TRUE,
          boot.num = 399,
          bw.x = NULL,
          bw.y = NULL,
          random.seed = 42,
          \dots)
}

\arguments{
  \item{data.x, data.y}{
    univariate vectors containing the variables.
  }
 
  \item{bootstrap}{
    a logical value which specifies whether to conduct the bootstrap
    test or not. If set to \code{FALSE}, only the statistic will be
    computed. Defaults to \code{TRUE}.
  }
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
 
  \item{bw.x, bw.y}{
    numeric (scalar) bandwidths. Defaults to plug-in (see details below).
  }
 
  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }

  \item{\dots}{ additional arguments supplied to specify the bandwidth
    type, kernel types, and so on.  This is used since we specify bw as
    a numeric scalar and not a \code{bandwidth} object, and is of
    interest if you do not desire the default behaviours. To change the
    defaults, you may specify any of \code{bwscaling}, \code{bwtype},
    \code{ckertype}, \code{ckerorder}, \code{ukertype},
    \code{okertype}.}
 
}

\value{
  \code{npunitest} returns an object of type \code{unitest} with the
  following components
 
  \item{Srho}{ the statistic \code{Srho} }
  \item{Srho.bootstrap}{ contains the bootstrap replications of \code{Srho} }
  \item{P}{ the P-value of the statistic }
  \item{boot.num}{ number of bootstrap replications }
  \item{bw.x, bw.y}{ scalar bandwidths for \code{data.x, data.y} }
   
  \code{\link{summary}} supports object of type \code{unitest}.
}
\references{

  Granger, C.W. and E. Maasoumi and J.S. Racine (2004), \dQuote{A
  dependence metric for possibly nonlinear processes}, Journal of Time
  Series Analysis, 25, 649-669.
  
  Maasoumi, E. and J.S. Racine (2002), \dQuote{Entropy and
    Predictability of Stock Market Returns,} Journal of Econometrics,
    107, 2, pp 291-312.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{

  \code{npunitest} computes the nonparametric metric entropy (normalized
   Hellinger of Granger, Maasoumi and Racine (2004)) for testing
   equality of two univariate density/probability functions,
   \eqn{D[f(x), f(y)]}{D[f(x), f(y)]}. See Maasoumi and Racine (2002)
   for details. Default bandwidths are of the plug-in variety
   (\code{\link{bw.SJ}} for continuous variables and direct plug-in for
   discrete variables). The bootstrap is conducted via simple resampling
   with replacement from the joint sample.

}

\section{Usage Issues}{
 
  See the example below for proper usage.
 
}

\seealso{
  \code{\link{npdeneqtest},\link{npdeptest},\link{npsdeptest}, \link{npsymtest}}
}

\examples{
set.seed(1234)
n <- 1000

## Compute the statistic only, different densities

x <- rchisq(n,df=10)
y <- rnorm(n,sd=10000,mean=-10000)

npunitest(x,y,bootstrap=FALSE)

## Compute the statistic only, same densities

x <- rnorm(n)
y <- rnorm(n)

npunitest(x,y,bootstrap=FALSE)

# Data drawn from same density function

x <- rnorm(n)
y <- rnorm(n)
npunitest(x,y,boot.num=99)

Sys.sleep(5)

# Data drawn from different density functions

x <- rnorm(n)
y <- rchisq(n,df=5)
npunitest(x,y,boot.num=99)

# Data drawn from same probability function

x <- factor(rbinom(n,2,.5))
y <- factor(rbinom(n,2,.5))
npunitest(x,y,boot.num=99)

Sys.sleep(5)

# Data drawn from same probability functions

x <- factor(rbinom(n,2,.5))
y <- factor(rbinom(n,2,.1))
npunitest(x,y,boot.num=99)
}
\keyword{ univar }
\keyword{ nonparametric }
