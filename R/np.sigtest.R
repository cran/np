# $Id: np.sigtest.R,v 1.44 2006/11/03 20:41:14 jracine Exp $

# This function implements a test of significance for both discrete
# and continuous variables. It accepts a data frame for explanatory
# data (mixed datatypes allowed), a vector for y for a regression
# model, an npregbw object, and a set of indices for the
# columns of X for which the test is to be run (default = all).

# Note - this conducts _individual_ tests of significance only. It
# uses a wild bootstrap to handle potential heteroskedasticity (though
# it perhaps could be readily changed to resample (y.star, X) pairs
# and perhaps this is desirable).

# Tristen XXX - this one is trivial to clean up - no plotting etc. All
# that is needed is a simple pretty print method like that for
# npcmstest. The difference here is that by default it conducts one
# test for each column of xdat and therefore returns a vector of In
# test statistics and a vector of P values (along with the actual
# indices used).

# Tristen XXX - we need to pass args for things like bandwidth options
# etc. for bootstrap method II, and for the npreg stuff as
# well...

npsigtest <- function(xdat,
                      ydat,
                      bws=NULL,
                      data=NULL,
                      boot.num=399,
                      boot.method=c("iid","wild","wild-rademacher"),
                      boot.type=c("I","II"),
                      index=seq(1,ncol(xdat)),
                      random.seed = 42,
                      ...) {
  miss.xy <- c(missing(xdat),missing(ydat))
  miss.bw <- is.null(bws)

  if (all(miss.xy) && !is.null(bws$formula)) {
    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf)

    ydat <- model.response(tmf)
    xdat <- tmf[, attr(attr(tmf, "terms"),"term.labels"), drop = FALSE]
  } else {
    if(all(miss.xy) && !is.null(bws$call)){
      xdat = eval(bws$call[["xdat"]], environment(bws$call))
      ydat = eval(bws$call[["ydat"]], environment(bws$call))
    }
    xdat = toFrame(xdat)
    
    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]
  }

  if (is.factor(ydat))
    stop("dependent variable must be continuous.")
  
  set.seed(random.seed)
  

  boot.type <- match.arg(boot.type)
  boot.method <- match.arg(boot.method)  

  if(boot.type=="II") {
    ## Store a copy of the bandwidths passed in
    bws.original <- bws
  }

  num.obs <- nrow(xdat)
  
  In <- numeric(length(index))
  P <- numeric(length(index))

  ## Some constants for the wild bootstrap

  a <- -0.6180339887499  # (1-sqrt(5))/2
  b <- 1.6180339887499   # (1+sqrt(5))/2
  P.a <-0.72360679774998 # (1+sqrt(5))/(2*sqrt(5))

  ## A vector for storing the resampled statistics

  In.vec <- numeric(boot.num)

  ## ii is the counter for successive elements of In and P...

  In.mat = matrix(data = 0, ncol = length(index), nrow = boot.num)

  ii <- 0

  console <- newLineConsole()

  for(i in index) {

    ## Increment counter...

    ii <- ii + 1

    if(boot.type=="II") {

      ## Reset bw vals to original as the ith component of bws gets
      ## overwritten when index changes so needs to be seet to its
      ## original value

      bws <- bws.original

    }

    ## Note - xdat must be a data frame

    ## Construct In, the average value of the squared derivatives of
    ## the jth element, discrete or continuous

    In[ii] <- mean((npreg(txdat = xdat, tydat = ydat,
                          bws=bws, gradients=TRUE,...)$grad[,i])^2)

    ## We now construct mhat.xi holding constant the variable whose
    ## significance is being tested at its median. First, make a copy
    ## of the data frame xdat

    xdat.eval <- xdat

    ## Check whether variable being tested is a factor (unordered or
    ## ordered) or not. If it is, cast the median as a factor. Note
    ## that ## all we do is to replace the ith column of xdat with its
    ## median, then ## evaluate the mean holding xdat[,i] constant at
    ## this value.

    ##if(is.factor(xdat[,i])) {
    ##      xdat.eval[,i] <- cast(levels(xdat[,i])[median(as.integer(xdat[,i]))], xdat[,i])
    ##} else {
    ##  xdat.eval[,i] <- median(xdat[,i])
    ##}
    ##xdat.eval[,i] <- uocquantile(xdat[,i], 0.5)

    xdat.eval[,i] <- sample(xdat[,i], replace=TRUE)

    mhat.xi <-  npreg(txdat = xdat,
                      tydat = ydat,
                      exdat=xdat.eval,
                      bws=bws,...)$mean
    
    delta.bar <- mean(ydat-mhat.xi)
  
    ei <- ydat - mhat.xi - delta.bar

    for(i.star in 1:boot.num) {

      if(boot.type=="I") {
        msg <- paste("Bootstrap replication ",
                     i.star,
                     "/",
                     boot.num,
                     " for variable ",
                     i,
                     " of (",
                     paste(index,collapse=","),
                     ")... ",
                     sep="")
      } else {
        msg <- paste("Bootstrap rep. ",
                     i.star,
                     "/",
                     boot.num,
                     " for variable ",
                     i,
                     " of (",
                     paste(index,collapse=","),
                     ")... ",
                     sep="")
      }

      console <- printPush(msg = msg, console)

      if(boot.method == "iid") {

        ydat.star <- mhat.xi + sample(ei, replace=TRUE)

      } else if(boot.method == "wild") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated holding
        ## the variable tested at its median, and add to that a wild
        ## bootstrap draw from the original disturbance vector
        
##        ydat.star <- mhat.xi + (ei-mean(ei))*
##          ifelse(rbinom(num.obs, 1, P.a) == 1, a, b)+
##            mean(ei)

        ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, a, b)

      } else if(boot.method == "wild-rademacher") {

        ## Conduct a wild bootstrap. We generate a sample for ydat
        ## (ydat.star) drawn from the conditional mean evaluated holding
        ## the variable tested at its median, and add to that a wild
        ## bootstrap draw from the original disturbance vector
        
##        ydat.star <- mhat.xi + (ei-mean(ei))*
##          ifelse(rbinom(num.obs, 1,0.5) == 1, -1, 1)+
##            mean(ei)

        ydat.star <- mhat.xi + ei*ifelse(rbinom(num.obs, 1, P.a) == 1, a, b)
        

      } 

      if(boot.type=="II") {

        ## For Bootstrap II method, starting values are taken from
        ## bandwidths passed in (bws.original). We then conduct
        ## cross-validation for the bootstrap sample and use only the
        ## new bw for variable i along with the original bandwidths
        ## for the remaining variables

        bws.boot <- npregbw(xdat = xdat, ydat = ydat.star,
                            bws=bws.original,...)        

        ## Copy the new cross-validated bandwidth for variable i into
        ## bw.original and use this below.
        
        bws <- bws.original

        bws$bw[i] <- bws.boot$bw[i]

      }

      In.vec[i.star] <- mean((npreg(txdat = xdat,
                                    tydat = ydat.star,
                                    bws=bws,
                                    gradients=TRUE,...)$grad[,i])^2)

      console <- printPop(console)
    }

    ## Compute the P-value

    P[ii] <- mean(ifelse(In.vec>In[ii],1,0))

    In.mat[,ii] = In.vec

  }

  ## Return a list containing the statistic and its P-value
  ##cat("\n")
  
  ## Tristen XXX - I would also like, if/when you implement a pretty
  ## print method, to pass back a matrix each column containing the
  ## bootstrapped In.vec for each variable...

  sigtest(In=In, In.mat, P=P,
          bws = if(miss.bw) NA else bws,
          ixvar = index,
          boot.method, boot.type, boot.num)

}
