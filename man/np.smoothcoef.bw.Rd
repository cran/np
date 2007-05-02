% $ID$
\name{npscoefbw}
\alias{npscoefbw}
\alias{npscoefbw.NULL}
\alias{npscoefbw.default}
\alias{npscoefbw.formula}
\alias{npscoefbw.scbandwidth}
\title{Smooth Coefficient Kernel Regression Bandwidth Selection}

\description{
  \code{npscoefbw} computes a bandwidth object for a smooth
  coefficient kernel regression estimate of a one (1) dimensional
  dependent variable on 
  \eqn{p+q}-variate explanatory data, using the model
  \eqn{Y_i = W_{i}^{\prime} \gamma (Z_i) + u_i}{Y_i = t(W_i) * gamma(Z_i) +
    u_i} where \eqn{W_i'=(1,X_i')}{t(W_i) = (1,t(X_i))}
  given training points (consisting of explanatory data and dependent
  data), and a bandwidth specification, which can be a \code{rbandwidth}
  object, or a bandwidth vector, bandwidth type and kernel type. 
}

\usage{
npscoefbw(\dots)

\method{npscoefbw}{formula}(formula, data, subset, na.action,\dots)

\method{npscoefbw}{NULL}(xdat = stop("invoked without data 'xdat'"),
          ydat = stop("invoked without data 'ydat'"),
          zdat = NULL,
          bws,
          \dots)

\method{npscoefbw}{default}(xdat = stop("invoked without data 'xdat'"),
          ydat = stop("invoked without data 'ydat'"),
          zdat = NULL,
          bws, 
          nmulti,
          cv.iterate,
          cv.num.iterations,
          backfit.iterate,
          backfit.maxiter,
          backfit.tol,
          ftol,
          tol,
          bandwidth.compute = TRUE,
          bwmethod,
          bwscaling,
          bwtype,
          ckertype,
          ckerorder,
          \dots)

\method{npscoefbw}{scbandwidth}(xdat = stop("invoked without data 'xdat'"),
          ydat = stop("invoked without data 'ydat'"),
          zdat = NULL,
          bws,
          nmulti,
          cv.iterate = FALSE,
          cv.num.iterations = 1,
          backfit.iterate = FALSE,
          backfit.maxiter = 100,
          backfit.tol = .Machine$double.eps,
          ftol = .Machine$double.eps,
          tol = sqrt(.Machine$double.eps),
          bandwidth.compute = TRUE,
          \dots)
}
\arguments{
  \item{formula}{
    a symbolic description of variables on which bandwidth selection is
    to be performed. The details of constructing a formula are
    described below.
  }

  \item{data}{
    an optional data frame, list or environment (or object
    coercible to a data frame by \code{\link{as.data.frame}}) containing the variables
    in the model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which the
    function is called.
  }

  \item{subset}{
    an optional vector specifying a subset of observations to be used in
    the fitting process. 
  }

  \item{na.action}{
    a function which indicates what should happen when the data contain
   \code{NA}s. The default is set by the \code{\link{na.action}} setting of options, and is
    \code{\link{na.fail}} if that is unset. The (recommended) default is
    \code{\link{na.omit}}.  
  }

  \item{xdat}{
    a \eqn{p}-variate data frame of explanatory data (training data),
    which, by default, populates the columns \eqn{2} through \eqn{p+1}
    of \eqn{W} in the model equation, and in the 
    absence of \code{zdat}, will also correspond to 
    \eqn{Z} from the model equation.
  }

  \item{ydat}{
    a one (1) dimensional numeric or integer vector of dependent data, each
    element \eqn{i} corresponding to each observation (row) \eqn{i} of
    \code{xdat}. 
  }

  \item{zdat}{
    an optionally specified \eqn{q}-variate data frame of explanatory
    data (training data), which corresponds to \eqn{Z}
    in the model equation. Defaults to be the same as \code{xdat}.
  }


  \item{bws}{
    a bandwidth specification. This can be set as a \code{scbandwidth}
    object returned from a previous invocation, or as a vector of
    bandwidths, with each element \eqn{i} corresponding to the bandwidth
    for column \eqn{i} in \code{xdat}. In either case, the bandwidth
    supplied will serve as a starting point in the numerical search for
    optimal bandwidths. If specified as a vector, then additional
    arguments will need to be supplied as necessary to specify the
    bandwidth type, kernel types, selection methods, and so on. This can
    be left unset.
  }

  \item{\dots}{
  additional arguments supplied to specify the regression type,
  bandwidth type, kernel types, selection methods, and so on, detailed
  below. 
  }
  

  \item{bandwidth.compute}{
    a logical value which specifies whether to do a numerical search for
    bandwidths or not. If set to \code{FALSE}, a \code{scbandwidth} object
    will be returned with bandwidths set to those specified
    in \code{bws}. Defaults to \code{TRUE}.
  }

  \item{bwmethod}{
    which method was used to select bandwidths. \code{cv.ls}
    specifies least-squares cross-validation, which is all that is
    currently supported. Defaults to \code{cv.ls}.
  }

  \item{bwscaling}{
    a logical value that when set to \code{TRUE} the supplied
    bandwidths are interpreted as `scale factors' (\eqn{c_j}{c[j]}),
    otherwise when the value is \code{FALSE} they are interpreted as `raw
    bandwidths' (\eqn{h_j}{h[j]} for continuous datatypes,
		\eqn{\lambda_j}{lambda[j]} for discrete datatypes). For continuous
		datatypes, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are related by the
		formula \eqn{h_j = c_j \sigma_j n^{-1/(2P+l)}}{h[j] =
			c[j]*sigma[j]*n^(-1/(2*P+l))}, where \eqn{\sigma_j}{sigma[j]}
    is the standard deviation of continuous variable \eqn{j},
    \eqn{n} the number of observations, \eqn{P}{P} the order of the kernel,
    and \eqn{l}{l} the number of continuous variables. For discrete
    datatypes, \eqn{c_j}{c[j]} and \eqn{h_j}{h[j]} are related by the
		formula \eqn{h_j = c_jn^{-2/(2P+l)}}{h[j] = c[j]*n^(-2/(2*P+l))},
		where here \eqn{j}{[j]} denotes discrete variable \eqn{j}{j}.
    Defaults to \code{FALSE}.
  }

  \item{bwtype}{
    character string used for the continuous variable bandwidth type,
    specifying the type of bandwidth provided. Defaults to
    \code{fixed}. Option summary:\cr
    \code{fixed}: fixed bandwidths or scale factors \cr
    \code{generalized_nn}: generalized nearest neighbors \cr
    \code{adaptive_nn}:  adaptive nearest neighbors
  }
  
  \item{ckertype}{
    character string used to specify the continuous kernel type.
    Can be set as \code{gaussian}, \code{epanechnikov}, or
    \code{uniform}. Defaults to \code{gaussian}.
  }

  \item{ckerorder}{
    numeric value specifying kernel order (one of
    \code{(2,4,6,8)}). Kernel order specified along with a
  \code{uniform} continuous kernel type will be ignored. Defaults to
  \code{2}.
  }

  \item{nmulti}{
    integer number of times to restart the process of finding extrema of
    the cross-validation function from different (random) initial
    points. Defaults to \code{ncol(xdat)}.
  }

  \item{ftol}{
    tolerance on the value of the cross-validation function 
    evaluated at located minima. Defaults to \code{1.19e-07
      (FLT_EPSILON)}
  }

  \item{tol}{
    tolerance on the position of located minima of the 
    cross-validation function. Defaults to \code{1.49e-08
      (sqrt(DBL_EPSILON))}
  }

  \item{cv.iterate}{
    boolean value specifying whether or not to perform iterative,
    cross-validated backfitting on the data. See details for limitations
    of the backfitting procedure. Defaults to \code{FALSE}.
  }

  \item{cv.num.iterations}{
    integer specifying the number of times to iterate the backfitting
    process over all covariates. Defaults to \code{1}.
  }

  \item{backfit.iterate}{
    boolean value specifying whether or not to iterate evaluations of
    the smooth coefficient estimator, for extra accuracy, during the
    cross-validated backfitting procedure. Defaults to \code{FALSE}.
  }

  \item{backfit.maxiter}{
    integer specifying the maximum number of times to iterate the
    evaluation of the smooth coefficient estimator in the attempt to
    obtain the desired accuracy. Defaults to \code{100}.
  }

  \item{backfit.tol}{
    tolerance to determine convergence of iterated evaluations of the
    smooth coefficient estimator. Defaults to \code{.Machine$double.eps}.
  }
}

\details{
 \code{npscoefbw} implements a variety of methods for semiparametric regression on
 multivariate (\eqn{p+q}-variate) explanatory data defined over a set of possibly
 continuous data. The approach is
 based on Li and Racine (2003) who employ \sQuote{generalized product kernels}
 that admit a mix of continuous and discrete datatypes.

 Three classes of kernel estimators for the continuous datatypes are
 available: fixed, adaptive nearest-neighbor, and generalized
 nearest-neighbor. Adaptive nearest-neighbor bandwidths change with
 each sample realization in the set, \eqn{x_i}{x[i]}, when estimating the
 density at the point \eqn{x}. Generalized nearest-neighbor bandwidths change
 with the point at which the density is estimated, \eqn{x}. Fixed bandwidths
 are constant over the support of \eqn{x}.

 \code{npscoefbw} may be invoked \emph{either} with a formula-like
 symbolic 
 description of variables on which bandwidth selection is to be
 performed \emph{or} through a simpler interface whereby data is passed
 directly to the function via the \code{xdat}, \code{ydat}, and
 \code{zdat} 
 parameters. Use of these two interfaces is \bold{mutually exclusive}.

 Data contained in the data frame \code{xdat} may be continuous and
 in \code{zdat} may be of mixed type. Data can be entered in an
 arbitrary order
 and data
 types will be detected automatically by the routine (see \code{\link{np}} for
 details).

 Data for which bandwidths are to be estimated may be specified
 symbolically. A typical description has the form \code{dependent}
 \code{data} \code{~} \code{parametric} \code{explanatory} \code{data}
 \code{|} \code{nonparametric} \code{explanatory} \code{data},
 where \code{dependent} \code{data} is a univariate response, and
 \code{parametric} \code{explanatory} \code{data} and
\code{nonparametric} \code{explanatory}
 \code{data} are both series of variables specified by name, separated by
 the separation character '+'. For example, \code{ y1 ~ x1 + x2 | z1 }
 specifies that the bandwidth object for the smooth coefficient model
 with response \code{y1}, linear parametric regressors \code{x1} and
 \code{x2}, and nonparametric regressor (that is, the slope-changing
 variable) \code{z1} is to be estimated. See below for further examples.
 In the case where the nonparametric (slope-changing) variable is not
 specified, it is assumed to be the same as the parametric variable.

 A variety of kernels may be specified by the user. Kernels implemented
 for continuous datatypes include the second, fourth, sixth, and eighth
 order Gaussian and Epanechnikov kernels, and the uniform
 kernel. Unordered discrete datatypes use a variation on Aitchison and
 Aitken's (1976) kernel, while ordered datatypes use a variation of the
 Wang and van Ryzin (1981) kernel.
}
\value{
  if \code{bwtype} is set to \code{fixed}, an object containing bandwidths
  (or scale factors if \code{bwscaling = TRUE}) is returned. If it is set to
  \code{generalized_nn} or \code{adaptive_nn}, then instead the \eqn{k}th nearest
  neighbors are returned for the continuous variables while the discrete
  kernel bandwidths are returned for the discrete variables. Bandwidths
  are stored in a vector under the component name \code{bw}. Backfitted
  bandwidths are stored under the component name \code{bw.fitted}.

  The functions \code{\link{predict}}, \code{\link{summary}}, and
  \code{\link{plot}} support 
  objects of this class.

}



\references{

 Aitchison, J. and C.G.G. Aitken (1976), \dQuote{Multivariate binary
   discrimination by the kernel method,} Biometrika, 63, 413-420.

 Cai Z. (2007), \dQuote{Trending time-varying coefficient time series
   models with serially correlated errors,} Journal of Econometrics,
 136, 163-188.

 Hastie, T. and R. Tibshirani (1993), \dQuote{Varying-coefficient
   models,} Journal of the Royal Statistical Society, B 55, 757-796.

 Li, Q. and J.S. Racine (2007), \emph{Nonparametric Econometrics: Theory
 and Practice,} Princeton University Press.

  Li, Q. and J.S. Racine (2007), \dQuote{Smooth varying-coefficient
    nonparametric models for qualitative and quantitative data,}
  manuscript.

  Pagan, A. and A. Ullah (1999), \emph{Nonparametric Econometrics,}
 Cambridge University Press.
 
 Racine, J.S. and D. Ouyang and Q. Li (2007), \dQuote{Nonparametric
   multilevel models: a smoothing approach,} manuscript.

 Wang, M.C. and J. van Ryzin (1981), \dQuote{A class of smooth estimators
 for discrete distributions,}  Biometrika, 68, 301-309.

}
\author{
    Tristen Hayfield \email{hayfietj@mcmaster.ca}, Jeffrey S. Racine
    \email{racinej@mcmaster.ca}
}

\section{Usage Issues}{
  If you are using data of mixed types, then it is advisable to use the
  \code{\link{data.frame}} function to construct your input data and not
  \code{\link{cbind}}, since \code{\link{cbind}} will typically not work as
  intended on mixed data types and will coerce the data to the same
  type.

  Caution: multivariate data-driven bandwidth selection methods are, by
  their nature, \emph{computationally intensive}. Virtually all methods
  require dropping the \eqn{i}th observation from the data set, computing an
  object, repeating this for all observations in the sample, then
  averaging each of these leave-one-out estimates for a \emph{given}
  value of the bandwidth vector, and only then repeating this a large
  number of times in order to conduct multivariate numerical
  minimization/maximization. Furthermore, due to the potential for local
  minima/maxima, \emph{restarting this procedure a large number of times may
  often be necessary}. This can be frustrating for users possessing
  large datasets. For exploratory purposes, you may wish to override the
  default search tolerances, say, setting ftol=.01 and tol=.01 and
  conduct multistarting (the default is to restart ncol(zdat)
  times) as is done for a number of examples. Once the procedure
  terminates, you can restart search with default tolerances using those
  bandwidths obtained from the less rigorous search (i.e., set
  \code{bws=bw} on subsequent calls to this routine where \code{bw} is
  the initial bandwidth object).  A version of this package using the
  \code{Rmpi} wrapper is under development that allows one to deploy
  this software in a clustered computing environment to facilitate
  computation involving large datasets.

  Support for backfitted bandwidths is experimental and is limited in
  functionality. The code does not support asymptotic standard errors
  or out of sample estimates with backfitting.

}

\seealso{
  \code{\link{npregbw}}, \code{\link{npreg}} }
\examples{
# EXAMPLE 1 (INTERFACE=FORMULA): For this example, we compute a
# smooth coefficient model with differing X and Z data

n <- 100
x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
bw <- npscoefbw(formula=y~x|z)
summary(bw)

# EXAMPLE 1 (INTERFACE=DATAFRAME): For this example, we compute a
# smooth coefficient model with differing X and Z data

n <- 100
x <- runif(n)
z <- runif(n, min=-2, max=2)
y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
bw <- npscoefbw(xdat=x, ydat=y, zdat=z)
summary(bw)
}
\keyword{nonparametric}
