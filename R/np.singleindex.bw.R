npindexbw <-
  function(...){
    args = list(...)
    if (is(args[[1]],"formula"))
      UseMethod("npindexbw",args[[1]])
    else if (!is.null(args$formula))
      UseMethod("npindexbw",args$formula)
    else
      UseMethod("npindexbw",args[[w(names(args)=="bws")[1]]])
  }

npindexbw.formula <-
  function(formula, data, subset, na.action, ...){
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(mf), nomatch = 0)
    mf <- mf[c(1,m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    
    ydat <- model.response(mf)
    xdat <- mf[, attr(attr(mf, "terms"),"term.labels"), drop = FALSE]
    
    tbw = npindexbw(xdat = xdat, ydat = ydat, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    tbw$call <- match.call(expand.dots = FALSE)
    environment(tbw$call) <- parent.frame()
    tbw$formula <- formula
    tbw$rows.omit <- as.vector(attr(mf,"na.action"))
    tbw$nobs.omit <- length(tbw$rows.omit)
    tbw$terms <- attr(mf,"terms")

    tbw <-
      updateBwNameMetadata(nameList =
                           list(ynames =
                                attr(mf, "names")[attr(tbw$terms, "response")]),
                           bws = tbw)

    tbw
  }


npindexbw.NULL <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, bandwidth.compute = TRUE, ...){

    xdat <- toFrame(xdat)

    bws <- double(ncol(xdat)+1)

    tbw <- npindexbw.default(xdat = xdat,
                             ydat = ydat,
                             bws = bws, ...)

    ## clean up (possible) inconsistencies due to recursion ...
    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    tbw <- updateBwNameMetadata(nameList =
                                list(ynames = deparse(substitute(ydat))),
                                bws = tbw)

    tbw
  }
  

npindexbw.default <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, bandwidth.compute = TRUE,
           nmulti, ...){

    xdat <- toFrame(xdat)

    if(!(is.vector(ydat) | is.factor(ydat)))
      stop("'ydat' must be a vector")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")
      
      warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    if (coarseclass(bws) != "numeric" | length(bws) != ncol(xdat)+1)
      stop(paste("manually specified 'bws' must be a numeric vector of length ncol(xdat)+1.",
                 "See documentation for details."))

    tbw <- sibandwidth(beta = bws[1:ncol(xdat)],
                       h = bws[ncol(xdat)+1], ...,
                       nobs = dim(xdat)[1],
                       xdati = untangle(xdat),
                       ydati = untangle(data.frame(ydat)),
                       xnames = names(xdat),
                       ynames = deparse(substitute(ydat)),
                       bandwidth = bws[ncol(xdat)+1],
                       bandwidth.compute = bandwidth.compute)


    if (tbw$method == "kleinspady" & !setequal(ydat,c(0,1))) 
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("nmulti")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)
    
    if (bandwidth.compute)
      tbw <- eval(parse(text=paste("npindexbw.sibandwidth(xdat=xdat, ydat=ydat, bws=tbw",
                          ifelse(any.m, ",",""),
                          paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                          ")")))

    mc <- match.call(expand.dots = FALSE)
    environment(mc) <- parent.frame()
    tbw$call <- mc

    return(tbw)
  }

  
npindexbw.sibandwidth <-
  function(xdat = stop("training data xdat missing"),
           ydat = stop("training data ydat missing"),
           bws, nmulti, random.seed = 42, ...){

    set.seed(random.seed)
    xdat = toFrame(xdat)

    if (missing(nmulti))
      nmulti = ncol(xdat)
    
    if (bws$method == "kleinspady" & !setequal(ydat,c(0,1))) 
      stop("Klein and Spady's estimator requires binary ydat with 0/1 values only")

    if (ncol(xdat) < 2) {
      if (coarseclass(xdat[,1]) != "numeric")
        stop("xdat must contain at least one continuous variable")
      
      warning(paste("xdat has one dimension. Using a single index model to reduce",
                    "dimensionality is unnecessary."))
    }

    ## catch and destroy NA's
    goodrows = 1:dim(xdat)[1]
    rows.omit = attr(na.omit(data.frame(xdat,ydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    xdat = xdat[goodrows,,drop = FALSE]
    ydat = ydat[goodrows]
    
    ## convert to numeric
    if (is.factor(ydat))
      ydat <- dlev(ydat)[as.integer(ydat)]
    else
      ydat <- as.double(ydat)
    
    xdat = toMatrix(xdat)


    ## Note - there are two methods currently implemented, Ichimura's
    ## least squares approach and Klein and Spady's likelihood approach.

    ##We define ichimura's objective function. Since we normalize beta_1
    ##to equal 1, we only worry about beta_2...beta_k in the index
    ##function for these internals, and when computing the leave-one-out
    ##objective function use c(1,beta). However, we do indeed return
    ##c(1,beta) which can be used in the index.model function above.

    ichimura <- function(param) {

      ##Define the leave-one-out objective function, sum (y - \hat
      ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
      ## parameters in `param') and then let h denote the kth column.
      if (ncol(xdat) == 1)
        beta = numeric(0)
      else 
        beta <- param[1:(ncol(xdat)-1)]
      
      h <- param[ncol(xdat)]

      ## Next we define the sum of squared leave-one-out residuals

      sum.squares.leave.one.out <- function(xdat,ydat,beta,h) {

        ## Normalize beta_1 = 1 hence multiply X by c(1,beta)
        
        index <- as.matrix(xdat) %*% c(1,beta)

        fit.loo <-  npksum(txdat=index,tydat=ydat,leave.one.out=TRUE,bandwidth.divide=TRUE,bws=c(h))$ksum/npksum(txdat=index,leave.one.out=TRUE,bandwidth.divide=TRUE,bws=c(h))$ksum

        return(sum((ydat-fit.loo)^2))

      }

      ## For the objective function, we require a positive bandwidth, so
      ## return an infinite penalty for negative h

      if(h > 0) {
        return(sum.squares.leave.one.out(xdat,ydat,beta,h))
      } else {
        return(.Machine$double.xmax)
      }

    }

    ## We define ichimura's objective function. Since we normalize beta_1
    ## to equal 1, we only worry about beta_2...beta_k in the index
    ## function for these internals, and when computing the leave-one-out
    ## objective function use c(1,beta). However, we do indeed return
    ## c(1,beta) which can be used in the index.model function above.

    kleinspady <- function(param) {

      ## Define the leave-one-out objective function, sum (y - \hat
      ## G(X\hat\beta))^2. We let beta denote beta_2...beta_k (first k-1
      ## parameters in `param') and then let h denote the kth column.
      if (ncol(xdat) == 1)
        beta = numeric(0)
      else 
        beta <- param[1:(ncol(xdat)-1)]

      h <- param[ncol(xdat)]

      ## Next we define the sum of squared leave-one-out residuals

      sum.log.leave.one.out <- function(xdat,ydat,beta,h) {

        ## Normalize beta_1 = 1 hence multiply X by c(1,beta)
        
        index <- as.matrix(xdat) %*% c(1,beta)

        fit.loo <-  npksum(txdat=index,
                                 tydat=ydat,
                                 leave.one.out=TRUE,
                                 bandwidth.divide=TRUE,
                                 bws=c(h))$ksum / npksum(txdat=index,
                                   leave.one.out=TRUE,
                                   bandwidth.divide=TRUE,
                                   bws=c(h))$ksum

        ## Maximize the log likelihood, therefore minimize minus...
        
        return( - sum(ifelse(ydat==1,log(fit.loo),log(1-fit.loo))) )

      }

      ## For the objective function, we require a positive bandwidth, so
      ## return an infinite penalty for negative h

      if(h > 0) {
        return(sum.log.leave.one.out(xdat,ydat,beta,h))
      } else {
        return(.Machine$double.xmax)
      }

    }

    ## Now we implement multistarting

    min <- .Machine$double.xmax
    numimp <- 0
    value <- numeric(nmulti)

    console <- newLineConsole()

    for(i in 1:nmulti) {

      console <- printPush(paste(sep="", "Multistart ", i, " of ", nmulti, "..."), console)
      ##cv.console <- newLineConsole(console)
      
      n <- nrow(xdat)

      ## We use the nlm command to minimize the objective function using
      ## starting values. Note that since we normalize beta_1=1 here beta
      ## is the k-1 vector containing beta_2...beta_k

      if(i == 1) {

        ## Initial values taken from OLS fit with a constant used for
        ## multistart 1
        ols.fit <- lm(ydat~xdat,x=TRUE)
        fit <- fitted(ols.fit)
        
        if (ncol(xdat) != 1){
          if (setequal(bws$beta[2:ncol(xdat)],c(0)))
            beta <- coef(ols.fit)[3:ncol(ols.fit$x)]
          else
            beta = bws$beta[2:ncol(xdat)]
        } else { beta = numeric(0) }

        if (bws$bw == 0)
          h <- 1.06*(min(sd(fit),IQR(fit)/1.34))*n^{-1/5}
        else
          h = bws$bw
          
        
        if(bws$method == "ichimura") {

          suppressWarnings(nlm.return <- nlm(ichimura, c(beta,h)))
          
        } else if(bws$method == "kleinspady") {

          suppressWarnings(nlm.return <- nlm(kleinspady, c(beta,h)))
          
        }

        value[i] <- nlm.return$minimum

        param <- nlm.return$estimate
        min <- nlm.return$minimum
        numimp <- numimp + 1
        best <- i

      } else {

        ## Random initialization used for remaining multistarts

        beta <- rnorm(ncol(xdat)-1, sd = i)
        h <- runif(1,min=0,max=2)*(min(sd(fit),IQR(fit)/1.34))*n^{-1/5}

        if(bws$method == "ichimura") {

          suppressWarnings(nlm.return <- nlm(ichimura, c(beta,h)))
          
        } else if(bws$method == "kleinspady") {

          suppressWarnings(nlm.return <- nlm(kleinspady, c(beta,h)))
          
        }

        value[i] <- nlm.return$minimum

        if(nlm.return$minimum < min) {
          param <- nlm.return$estimate
          min <- nlm.return$minimum
          numimp <- numimp + 1
          best <- i
        }

      }
      console <- printPop(console)
    }
    ##cat('\n')

    ## Return a list with beta (we append the restricted value of
    ## beta_1=1), the bandwidth h, the value of the objective function at
    ## its minimum, the number of restarts that resulted in an improved
    ## value of the objective function, the restart that resulted in the
    ## smallest value, and the vector of objective function values.

    ## TRISTEN XXX - kindly make sure I have named return vals consistent
    ## with other functions, e.g., ifval, fval etc? Also, kindly return a
    ## `pretty print' type of object

    bws <- sibandwidth(beta =
                       if(ncol(xdat) == 1) 1.0
                       else c(1,param[1:(ncol(xdat)-1)]),
                       h = param[ncol(xdat)],
                       method = bws$method,
                       ckertype = bws$ckertype,
                       ckerorder = bws$ckerorder,
                       fval = min/bws$nobs,
                       ifval = best,
                       numimp = numimp,
                       fval.vector = value,
                       nobs = bws$nobs,
                       xdati = bws$xdati,
                       ydati = bws$ydati,
                       xnames = bws$xnames,
                       ynames = bws$ynames,
                       bandwidth = param[ncol(xdat)],
                       rows.omit = rows.omit)


    bws

  }
