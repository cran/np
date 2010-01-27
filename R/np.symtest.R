## Test for asymmetry applicable to both real-values and categorical
## univariate data.

npsymtest <- function(data = NULL,
                      boot.num = 399,
                      bw = NULL,
                      random.seed = 42) {

  if(is.data.frame(data)) stop("you must enter a data vector (not data frame)")
  if(is.null(data)) stop("you must enter a data vector")
  if(ncol(data.frame(data)) != 1) stop("data must have one dimension only")
  if(boot.num < 9) stop("number of bootstrap replications must be >= 9")

  set.seed(random.seed)

  ## Remove NA

  data <- na.omit(data)

  ## Note - this function accepts numeric and factors/ordered,
  ## however, factors must be integers.

  if(is.numeric(data)) {
    ## Numerical stability - integrate does best centered around zero,
    ## so standardize numeric data - this affects nothing else.
    data <- (data-mean(data))/sd(data)
    data.rotate <- -(data-mean(data))+mean(data)
    if(is.null(bw)) bw <- bw.SJ(data)
  } else {
    if(is.ordered(data)) {
      ## Rotate around the median, ordered case.
      data.levels <- levels(data)
      tmp <- as.numeric(data.matrix(data))
      location <- median(sort(unique(tmp)))
      data.rotate <- ordered(-(tmp-location)+location,levels=data.levels)
    } else {
      ## unordered case.
      data.levels <- levels(data)
      tmp <- as.numeric(data.matrix(data))
      location <- median(sort(unique(tmp)))
      data.rotate <- factor(-(tmp-location)+location,levels=data.levels)
    }
    ## Optimal bandwidth for factor (unordered) used for both
    ## (plug-in).
    if(is.null(bw)) {
      n <- length(data)
      c <- length(unique(data))
      xeval <- unique(data)
      p <- fitted(npudens(tdat=data,edat=xeval,bws=0))
      sum.Lambda3 <- c/(c-1)*sum(p*(1-p))
      sum.Lambda2.minus.Lambda1.sq <- c^2/((c-1)^2)*sum(p*(1-p))
      n.sum.Lambda1.sq <- n*sum(((1-c*p)/(c-1))^2)
      bw <- sum.Lambda3/(sum.Lambda2.minus.Lambda1.sq+n.sum.Lambda1.sq)
    }
  }
	
  ##  Impose the null and generate a data of resampled statistics

  if(is.numeric(data)) {
    data.null <- c(data,data.rotate)
  } else {
    if(is.ordered(data)) {
      data.null <- ordered(c(as.character(data),as.character(data.rotate)))
    } else {
      data.null <- factor(c(as.character(data),as.character(data.rotate)))
    }
  }

  ## Hellinger distance (Srho) between the densities for the data and
  ## the rotated data. This function accepts numeric and
  ## factor/ordered.
  
  Srho.sym <- function(data,data.rotate,bw) {
    if(ncol(data.frame(data)) != 1)  stop("data must have one dimension only") 
    if(is.numeric(data)) {
      h <- function(x,data,data.rotate) {
        f.data <- fitted(npudens(tdat=data,edat=x,bws=bw))
        f.data.rotate <- fitted(npudens(tdat=data.rotate,edat=x,bws=bw))
        return(0.5*(sqrt(f.data)-sqrt(f.data.rotate))**2)
      }
      return(integrate(h,-Inf,Inf,data=data,data.rotate=data.rotate)$value)
    } else {
      xeval <- unique(data)
      p.data <- fitted(npudens(tdat=data,edat=xeval,bws=bw))
      p.data.rotate <- fitted(npudens(tdat=data.rotate,edat=xeval,bws=bw))   
      ## Sum the information function over the unique probabilities and
      ## return.
      return(sum(0.5*(sqrt(p.data)-sqrt(p.data.rotate))**2))
    }
  }
  
  ## Compute Srho using kernel estimates, same bw is appropriate for
  ## both densities (they are the same data). Use plug-in estimators
  ## for both. 

  test.stat <- Srho.sym(data,data.rotate,bw)

  ## Function to be fed to tsboot - accepts a vector of integers
  ## corresponding to all observations in the sample (1,2,...) that
  ## get permuted/rearranged to define resampled data.

	boot.fun <- function(ii,data.null,bw)
	{
	  null.sample1 <- data.null[ii]
    if(is.numeric(data.null)) {
      null.sample2 <- -(null.sample1-mean(null.sample1))+mean(null.sample1)
    } else {
      if(is.ordered(data)) {
        ## Rotate around the median, ordered case
        null.sample1.levels <- levels(null.sample1)
        tmp <- as.numeric(data.matrix(null.sample1))
        location <- median(sort(unique(tmp)))
        null.sample2 <- ordered(-(tmp-location)+location,levels=data.levels)
      } else {
        ## unordered case
        null.sample1.levels <- levels(null.sample1)
        tmp <- as.numeric(data.matrix(null.sample1))
        location <- median(sort(unique(tmp)))
        null.sample2 <- factor(-(tmp-location)+location,levels=data.levels)
      }
    }
    return(Srho.sym(null.sample1,null.sample2,bw))
	}

  if(is.numeric(data)) {
    ## Optimal block length is based upon `data' (not data.null).
    boot.blocklen <- b.star(data,round=TRUE)[1,1]
  } else {
    boot.blocklen <- b.star(as.numeric(data.matrix(data)),round=TRUE)[1,1]
  }

  ## Need to bootstrap integers for data.null to accommodate both
  ## numeric and factor/ordered.

  ii <- 1:(2*length(data))

  console <- newLineConsole()
  console <- printPush(paste(sep="", "Working..."), console)

  resampled.stat <- tsboot(ii, boot.fun,
                           R = boot.num,
                           n.sim = length(data),
                           l = boot.blocklen,
                           sim = "geom",
                           data.null = data.null,
                           bw=bw)$t

  console <- printClear(console)
  console <- printPop(console)  

  p.value <- mean(ifelse(resampled.stat > test.stat, 1, 0))

  symtest(Srho = test.stat,
          Srho.bootstrap = resampled.stat,
          P = p.value,
          boot.num = boot.num,
          data.rotate = data.rotate,
          bw = bw)

}
