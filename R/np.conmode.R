npconmode <-
  function(bws = stop(paste("bandwidths are required to perform the estimate!",
             "please set 'bws'")), ...){
    args = list(...)

    if (is.recursive(bws)){
      if (!is.null(bws$formula) && is.null(args$txdat))
        UseMethod("npconmode",bws$formula)
      else if (!is.null(args$data) || !is.null(args$newdata))
        stop("data and newdata specified, but bws has no formula")
      else if (!is.null(bws$call) && is.null(args$txdat))
        UseMethod("npconmode",bws$call)
      else
        UseMethod("npconmode",bws)
    } else {
      UseMethod("npconmode",bws)
    }

  }

npconmode.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    umf <- tmf <- eval(tmf)

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]

    if ((has.eval <- !is.null(newdata))) {
      has.ey <- succeedWithResponse(tt, newdata)

      if (has.ey){
        umf <- emf <- model.frame(tt, data = newdata)
        eydat <- emf[, bws$variableNames[["response"]], drop = FALSE]
      } else {
        umf <- emf <- model.frame(formula(bws)[-2], data = newdata)
      }

      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }
    
    ev <-
      eval(parse(text=paste("npconmode(txdat = txdat, tydat = tydat,",
                   ifelse(has.eval,paste("exdat = exdat,",ifelse(has.ey,"eydat = eydat,","")),""),
                   "bws = bws, ...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npconmode.call <-
  function(bws, ...) {
    npconmode(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
              tydat = eval(bws$call[["ydat"]], environment(bws$call)),
              bws = bws, ...)
  }


npconmode.conbandwidth <-
  function (bws,
            txdat = stop("invoked without training data 'txdat'"),
            tydat = stop("invoked without training data 'tydat'"),
            exdat, eydat,
            ...){

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    no.ex = missing(exdat)
    no.ey = missing(eydat)

    if (!no.ex)
      exdat = toFrame(exdat)

    if (!no.ey)
      eydat = toFrame(eydat)

    ## catch and destroy NA's
    goodrows = 1:dim(txdat)[1]
    rows.omit = attr(na.omit(data.frame(txdat,tydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Tranining data has no rows without NAs")

    txdat = txdat[goodrows,,drop = FALSE]
    tydat = tydat[goodrows,,drop = FALSE]

    if (!no.ex){
      goodrows = 1:dim(exdat)[1]
      rows.omit = eval(parse(text=paste('attr(na.omit(data.frame(exdat',
                               ifelse(no.ey,"",",eydat"),')), "na.action")')))

      goodrows[rows.omit] = 0

      exdat = exdat[goodrows,,drop = FALSE]

      if (!no.ey)
        eydat = eydat[goodrows,,drop = FALSE]

      if (all(goodrows==0))
        stop("Evaluation data has no rows without NAs")
    }


    tnrow = dim(txdat)[1]
    enrow = ifelse(no.ex, tnrow, dim(exdat)[1])

    if (!no.ey & no.ex)
      stop("npconmode: invalid invocation: 'eydat' provided but not 'exdat'")

    if (bws$yndim != 1 | bws$yncon > 0)
      stop("'tydat' must consist of one (1) discrete variable")

    mdens = double(enrow)
    tdens = double(enrow)
    tf = logical(enrow)
    indices = integer(enrow)

    tdensE <- parse(text = paste("npcdens(txdat = txdat, tydat = tydat,",
                      "exdat = ", ifelse(no.ex, "txdat", "exdat"), ",",
                      "eydat = rep((bws$ydati$all.ulev[[1]])[i], enrow),",
                      "bws = bws)$condens"))


    for(i in 1:nlevels(tydat[,1])){
      tdens <- eval(tdensE)

      tf = tdens > mdens
      indices[tf] = i
      mdens[tf] = tdens[tf]
    }

    con.mode <- eval(parse(text = paste("conmode(bws = bws,",
                             "xeval = ", ifelse(no.ex, "txdat", "exdat"), ",",
                             ifelse(no.ey & !no.ex, "",
                                    paste("yeval = ", ifelse(no.ey, "tydat",
                                                             "eydat"), ",")),
                             "conmode = (bws$ydati$all.ulev[[1]])[indices],",
                             "condens = mdens, ntrain = nrow(txdat),",
                             "trainiseval = no.ex)")))
    
    if (!(no.ey & !no.ex)){
      confusion.matrix <- 
        table(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL),
              factor(con.mode$conmode,exclude = NULL), dnn=c("Actual", "Predicted"))

      cj <- match(levels(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL)),
                  levels(factor(con.mode$conmode,exclude = NULL)), nomatch = 0)
      rj <- cj > 0

      t.diag <- cj
      t.diag[rj] <-  diag(confusion.matrix[rj,cj,drop=FALSE])
      
      CCR.overall <- sum(t.diag)/enrow
      
      CCR.byoutcome <- t.diag/rowSums(confusion.matrix)
      names(CCR.byoutcome) <- levels(factor(if (no.ex) tydat[,1] else eydat[,1], exclude = NULL))

      con.mode$confusion.matrix <- confusion.matrix
      con.mode$CCR.overall <- CCR.overall
      con.mode$CCR.byoutcome <- CCR.byoutcome

      confusion.matrix <- confusion.matrix/enrow
      t.diag <- t.diag/enrow

      fit.mcfadden <- sum(t.diag) - (sum(confusion.matrix^2)-sum(t.diag^2))
      con.mode$fit.mcfadden <- fit.mcfadden
    }
    con.mode
  }

npconmode.default <-
  function (bws,
            txdat = stop("invoked without training data 'txdat'"),
            tydat = stop("invoked without training data 'tydat'"),
            exdat, eydat,
            ...){

    ## maintain x names and 'toFrame'
    txdat <- toFrame(txdat)

    ## maintain y names and 'toFrame'
    tydat <- toFrame(tydat)

    tbw <- npcdensbw(bws = bws,
                     xdat = txdat,
                     ydat = tydat,
                     bandwidth.compute = FALSE,
                     ...)

    mc.names <- names(match.call(expand.dots = FALSE))
    margs <- c("exdat", "eydat")
    m <- match(margs, mc.names, nomatch = 0)
    any.m <- any(m != 0)

    eval(parse(text=paste("npconmode.conbandwidth(txdat=txdat, tydat=tydat, bws=tbw",
                 ifelse(any.m, ",",""),
                 paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
                 ")")))
  }
