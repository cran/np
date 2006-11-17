smoothcoefficient <- 
  function(bws, eval, mean, merr = NA, beta = NA,
           grad = NA, gerr = NA, resid = NA,
           ntrain, trainiseval = FALSE, residuals = FALSE,
           betas = FALSE,
           xtra = rep(NA, 6)){

    if (missing(bws) | missing(eval) | missing(ntrain))
      stop("improper invocation of smoothcoefficient constructor")

    d = list(
      bw = bws$bw,
      bws = bws,
      xnames = bws$xnames,
      ynames = bws$ynames,
      nobs = nrow(eval),
      ndim = bws$ndim,
      nord = bws$nord,
      nuno = bws$nuno,
      ncon = bws$ncon,
      pscaling = bws$pscaling,
      ptype = bws$ptype,
      pckertype = bws$pckertype,
      pukertype = bws$pukertype,
      pokertype = bws$pokertype,
      eval = eval,
      mean = mean,
      merr = merr,
      ntrain = ntrain,
      trainiseval = trainiseval,
      residuals = residuals,
      betas = betas,
      beta = beta,
      grad = grad,
      gerr = gerr,
      resid = resid,
      R2 = xtra[1],
      MSE = xtra[2],
      MAE = xtra[3],
      MAPE = xtra[4],
      CORR = xtra[5],
      SIGN = xtra[6]
      )

    class(d) = "smoothcoefficient"

    d
  }

print.smoothcoefficient <- function(x, digits=NULL, ...){
  cat("\nSmooth Coefficient Model",
      "\nRegression data: ", x$ntrain, " training points,",
      ifelse(x$trainiseval, "",
             paste(" and ", x$nobs," evaluation points,", sep="")),
      " in ",x$ndim," variable(s)\n",sep="")

  print(matrix(x$bw,ncol=x$ndim,dimnames=list(paste(x$pscaling,":",sep=""),x$xnames)))

  cat(genRegEstStr(x))

  cat(genBwKerStrs(x$bws))
  cat('\n\n')  

  if(!missing(...))
    print(...,digits=digits)
  invisible(x)
}

coef.smoothcoefficient <- function(object, ...) {
 object$beta 
}
fitted.smoothcoefficient <- function(object, ...){
 object$mean 
}
plot.smoothcoefficient <- function(x, ...) { npplot(bws = x$bws, ...) }
residuals.smoothcoefficient <- function(object, ...) {
 if(object$residuals) { return(object$resid) } else { return(npscoef(bws = object$bws, residuals =TRUE)$resid) } 
}
se.smoothcoefficient <- function(x){ x$merr }
predict.smoothcoefficient <- function(object, ...) {
  fitted(eval(npscoef(bws = object$bws, ...), env = parent.frame()) )
}

summary.smoothcoefficient <- function(object, ...){
  cat("\nSmooth Coefficient Model",
      "\nRegression data: ", object$ntrain, " training points,",
      ifelse(object$trainiseval, "",
             paste(" and ", object$nobs," evaluation points,", sep="")),
      " in ",object$ndim," variable(s)\n",sep="")

  cat(genOmitStr(object))
  cat("\n")

  print(matrix(object$bw,ncol=object$ndim,dimnames=list(paste(object$pscaling,":",sep=""),object$xnames)))

  cat(genRegEstStr(object))
  cat("\n")
  cat(genGofStr(object))

  cat(genBwKerStrs(object$bws))
  cat('\n\n')  
}
