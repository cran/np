npcdens <-
  function(bws = stop(paste("bandwidths are required to perform the estimate!",
             "please set 'bws'")), ...){
    args = list(...)

    if (is.recursive(bws)){
      if (!is.null(bws$formula) && is.null(args$txdat))
        UseMethod("npcdens",bws$formula)
      else if (!is.null(args$data) || !is.null(args$newdata))
        stop("data and newdata specified, but bws has no formula")
      else if (!is.null(bws$call) && is.null(args$txdat))
        UseMethod("npcdens",bws$call)
      else
        UseMethod("npcdens",bws)
    } else {
      UseMethod("npcdens",bws)
    }

  }

npcdens.formula <-
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
      umf <- emf <- model.frame(tt, data = newdata)

      eydat <- emf[, bws$variableNames[["response"]], drop = FALSE]
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    ev <-
      eval(parse(text=paste("npcdens(txdat = txdat, tydat = tydat,",
                   ifelse(has.eval,"exdat = exdat, eydat = eydat,",""), "bws = bws, ...)")))
    ev$rows.omit <- as.vector(attr(umf,"na.action"))
    ev$nobs.omit <- length(ev$rows.omit)
    return(ev)
  }

npcdens.call <-
  function(bws, ...) {
    npcdens(txdat = eval(bws$call[["xdat"]], environment(bws$call)),
            tydat = eval(bws$call[["ydat"]], environment(bws$call)),
            bws = bws, ...)
  }


npcdens.conbandwidth <- function(bws,
                                 txdat = stop("invoked without training data 'txdat'"),
                                 tydat = stop("invoked without training data 'tydat'"),
                                 exdat, eydat, gradients = FALSE, ...){

  if (xor(missing(exdat),missing(eydat)))
    stop("evaluation data must be supplied for both 'exdat' and 'eydat'")

  no.exy = missing(exdat)

  txdat = toFrame(txdat)
  tydat = toFrame(tydat)

  if (!no.exy){
    exdat = toFrame(exdat)
    eydat = toFrame(eydat)

    if (! txdat %~% exdat )
      stop("'txdat' and 'exdat' are not similar data frames!")

    if (! tydat %~% eydat )
      stop("'tydat' and 'eydat' are not similar data frames!")

  }

  if (length(bws$xbw) != length(txdat))
    stop("length of bandwidth vector does not match number of columns of 'txdat'")

  if (length(bws$ybw) != length(tydat))
    stop("length of bandwidth vector does not match number of columns of 'tydat'")

  xccon = unlist(lapply(txdat[,bws$ixcon, drop = FALSE],class))
  if ((any(bws$ixcon) && !all((xccon == class(integer(0))) | (xccon == class(numeric(0))))) ||
      (any(bws$ixord) && !all(unlist(lapply(txdat[,bws$ixord, drop = FALSE],class)) ==
                             class(ordered(0)))) ||
      (any(bws$ixuno) && !all(unlist(lapply(txdat[,bws$ixuno, drop = FALSE],class)) ==
                             class(factor(0)))))
    stop("supplied bandwidths do not match 'txdat' in type")

  yccon = unlist(lapply(tydat[,bws$iycon, drop = FALSE],class))
  if ((any(bws$iycon) && !all((yccon == class(integer(0))) | (yccon == class(numeric(0))))) ||
      (any(bws$iyord) && !all(unlist(lapply(tydat[,bws$iyord, drop = FALSE],class)) ==
                             class(ordered(0)))) ||
      (any(bws$iyuno) && !all(unlist(lapply(tydat[,bws$iyuno, drop = FALSE],class)) ==
                             class(factor(0)))))
    stop("supplied bandwidths do not match 'tydat' in type")
  
  ## catch and destroy NA's
  goodrows = 1:dim(txdat)[1]
  rows.omit = attr(na.omit(data.frame(txdat,tydat)), "na.action")
  goodrows[rows.omit] = 0

  if (all(goodrows==0))
    stop("Data has no rows without NAs")

  txdat = txdat[goodrows,,drop = FALSE]
  tydat = tydat[goodrows,,drop = FALSE]

  if (!no.exy){
    goodrows = 1:dim(exdat)[1]
    rows.omit = attr(na.omit(data.frame(exdat,eydat)), "na.action")
    goodrows[rows.omit] = 0

    if (all(goodrows==0))
      stop("Data has no rows without NAs")

    exdat = exdat[goodrows,,drop = FALSE]
    eydat = eydat[goodrows,,drop = FALSE]
  }


  tnrow = nrow(txdat)
  enrow = ifelse(no.exy,tnrow,nrow(exdat))

  ## re-assign levels in training and evaluation data to ensure correct
  ## conversion to numeric type.
  
  txdat <- adjustLevels(txdat, bws$xdati)
  tydat <- adjustLevels(tydat, bws$ydati)
  
  if (!no.exy){
    exdat <- adjustLevels(exdat, bws$xdati, allowNewCells = TRUE)
    eydat <- adjustLevels(eydat, bws$ydati, allowNewCells = TRUE)
  }

  ## grab the evaluation data before it is converted to numeric
  if(no.exy){
    txeval <- txdat
    tyeval <- tydat
  } else {
    txeval <- exdat
    tyeval <- eydat
  }

  ## at this stage, data to be sent to the c routines must be converted to
  ## numeric type.
  
  tydat = toMatrix(tydat)

  tyuno = tydat[, bws$iyuno, drop = FALSE]
  tycon = tydat[, bws$iycon, drop = FALSE]
  tyord = tydat[, bws$iyord, drop = FALSE]


  txdat = toMatrix(txdat)

  txuno = txdat[, bws$ixuno, drop = FALSE]
  txcon = txdat[, bws$ixcon, drop = FALSE]
  txord = txdat[, bws$ixord, drop = FALSE]

  if (!no.exy){
    eydat = toMatrix(eydat)

    eyuno = eydat[, bws$iyuno, drop = FALSE]
    eycon = eydat[, bws$iycon, drop = FALSE]
    eyord = eydat[, bws$iyord, drop = FALSE]


    exdat = toMatrix(exdat)

    exuno = exdat[, bws$ixuno, drop = FALSE]
    excon = exdat[, bws$ixcon, drop = FALSE]
    exord = exdat[, bws$ixord, drop = FALSE]
  } else {
    eyuno = data.frame()
    eycon = data.frame()
    eyord = data.frame()

    exuno = data.frame()
    excon = data.frame()
    exord = data.frame()
  }

  myopti = list(
    num_obs_train = tnrow,
    num_obs_eval = enrow,
    int_LARGE_SF = ifelse(bws$scaling, SF_NORMAL, SF_ARB),
    BANDWIDTH_den_extern = switch(bws$type,
      fixed = BW_FIXED,
      generalized_nn = BW_GEN_NN,
      adaptive_nn = BW_ADAP_NN),
    int_MINIMIZE_IO=ifelse(options('np.messages'), IO_MIN_FALSE, IO_MIN_TRUE),
    xkerneval = switch(bws$cxkertype,
      gaussian = CKER_GAUSS + bws$cxkerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cxkerorder/2 - 1,
      uniform = CKER_UNI),
    ykerneval = switch(bws$cykertype,
      gaussian = CKER_GAUSS + bws$cykerorder/2 - 1,
      epanechnikov = CKER_EPAN + bws$cykerorder/2 - 1,
      uniform = CKER_UNI),
    uxkerneval = switch(bws$uxkertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    uykerneval = switch(bws$uykertype,
      aitchisonaitken = UKER_AIT,
      liracine = UKER_LR),
    oxkerneval = switch(bws$oxkertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR),
    oykerneval = switch(bws$oykertype,
      wangvanryzin = OKER_WANG,
      liracine = OKER_LR),
    num_yuno = bws$ynuno,
    num_yord = bws$ynord,
    num_ycon = bws$yncon,
    num_xuno = bws$xnuno,
    num_xord = bws$xnord,
    num_xcon = bws$xncon,
    no.exy = no.exy,
    gradients = gradients,
    ymcv.numRow = attr(bws$ymcv, "num.row"),
    xmcv.numRow = attr(bws$xmcv, "num.row"),
    densOrDist = NP_DO_DENS)
  
  myout=
    .C("np_density_conditional",
       as.double(tyuno), as.double(tyord), as.double(tycon),
       as.double(txuno), as.double(txord), as.double(txcon),
       as.double(eyuno), as.double(eyord), as.double(eycon),
       as.double(exuno), as.double(exord), as.double(excon),
       as.double(c(bws$xbw[bws$ixcon],bws$ybw[bws$iycon],
                   bws$ybw[bws$iyuno],bws$ybw[bws$iyord],
                   bws$xbw[bws$ixuno],bws$xbw[bws$ixord])),
       as.double(bws$ymcv), as.double(attr(bws$ymcv, "pad.num")),
       as.double(bws$xmcv), as.double(attr(bws$xmcv, "pad.num")),
       as.integer(myopti),
       condens = double(enrow),
       conderr = double(enrow),
       congrad = double(enrow*bws$xndim),
       congerr = double(enrow*bws$xndim),
       log_likelihood = double(1),
       PACKAGE="np" )[c("condens", "conderr", "congrad", "congerr", "log_likelihood")]

  if(gradients){
    myout$congrad = matrix(data=myout$congrad, nrow = enrow, ncol = bws$xndim, byrow = FALSE) 
    rorder = numeric(bws$xndim)
    rorder[c((1:bws$xndim)[bws$ixcon], (1:bws$xndim)[bws$ixuno], (1:bws$xndim)[bws$ixord])]=1:bws$xndim
    myout$congrad = myout$congrad[, rorder, drop = FALSE]

    myout$congerr = matrix(data=myout$congerr, nrow = enrow, ncol = bws$xndim, byrow = FALSE)
    myout$congerr = myout$congerr[, rorder, drop = FALSE]
  } else {
    myout$congrad = NA
    myout$congerr = NA
  }


  return( condensity(bws = bws,
                     xeval = txeval,
                     yeval = tyeval,
                     condens = myout$condens, conderr = myout$conderr,
                     congrad = myout$congrad, congerr = myout$congerr,
                     ll = myout$log_likelihood,
                     ntrain = tnrow, trainiseval = no.exy,
                     rows.omit = rows.omit) )

}

npcdens.default <- function(bws,
                            txdat = stop("training data 'txdat' missing"),
                            tydat = stop("training data 'tydat' missing"),
                            exdat, eydat, gradients, ...){

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
  margs <- c("exdat", "eydat", "gradients")
  m <- match(margs, mc.names, nomatch = 0)
  any.m <- any(m != 0)
  

  eval(parse(text=paste("npcdens.conbandwidth(txdat=txdat, tydat=tydat, bws=tbw",
               ifelse(any.m, ",",""),
               paste(mc.names[m], ifelse(any.m,"=",""), mc.names[m], collapse=", "),
               ")")))
}

