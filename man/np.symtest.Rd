% $Id: np.cmstest.Rd,v 1.58 2006/11/03 21:17:20 tristen Exp $
\name{npsymtest}
\alias{npsymtest}
\title{ Kernel Consistent Asymmetry Test with Mixed Data Types }
\description{
  \code{npsymtest} implements a consistent entropy-based test for
  asymmetry as described in Maasoumi and Racine (2009).
}
\usage{
npsymtest(data = NULL,
          boot.num = 399,
          bw = NULL,
          random.seed = 42)
}

\arguments{
  \item{data}{
    a vector containing the variable.
  }
  
  \item{boot.num}{
    an integer value specifying the number of bootstrap
    replications to use. Defaults to \code{399}.
  }
  
  \item{bw}{
    a numeric (scalar) bandwidth. Defaults to plug-in (see details below).
  }
  
  \item{random.seed}{
    an integer used to seed R's random number generator. This is to
    ensure replicability. Defaults to 42.
  }
  
}

\value{
  \code{npsymtest} returns an object of type \code{symtest} with the
  following components
  
  \item{Srho}{ the statistic \code{Srho} }
  \item{Srho.bootstrap}{ contains the bootstrap replications of \code{Srho} }
  \item{P}{ the P-value of the statistic }
  \item{boot.num}{ number of bootstrap replications }
  \item{ data.rotate }{ the rotated data series }
  \item{ bw }{ the numeric (scalar) bandwidth }
    
  \code{\link{summary}} supports object of type \code{symtest}.
}
\references{

  Maasoumi, E. and J. S. Racine(2009), \dQuote{A robust entropy-based
  test of asymmetry for discrete and continuous processes}, Econometric
  Reviews, 28, 246-261.

}
\author{
    Tristen Hayfield \email{hayfield@phys.ethz.ch}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\details{
  
  \code{npsymtest} computes the Bhattacharya-Matusita-Hellinger metric
  entropy measure of divergence between the densities/probabilities of
  the data and the rotated data. See Maasoumi and Racine (2009) for
  details.
 
}

\section{Usage Issues}{
  
  When using data of type \code{\link{factor}} it is crucial that the
  variable not be an alphabetic character string (i.e. the factor must
  be integer-valued). The rotation is conducted about the median after
  conversion to type \code{\link{numeric}} which is then converted back
  to type \code{\link{factor}}. Failure to do so will have unpredictable
  results. See the example below for proper usage.
  
}

\examples{
set.seed(1234)

# Asymmetric discrete probability function

x <- factor(rbinom(50,2,.8))
npsymtest(x,boot.num=99)

Sys.sleep(5)

# Symmetric discrete probability function

x <- factor(rbinom(50,2,.5))
npsymtest(x,boot.num=99)

Sys.sleep(5)

# Asymmetric continuous density function

y <- rchisq(50,df=2)
npsymtest(y,boot.num=99)

Sys.sleep(5)

# Symmetric continuous density function

y <- rnorm(50)
npsymtest(y,boot.num=99)

}
\keyword{ univar }
\keyword{ nonparametric }
