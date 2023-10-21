

# Produce Index ----------------------------------------------------------------


#' Calculate survey indices by age.
#'
#' This is based on the methods described in
#' Berg et al. (2014): "Evaluation of alternative age-based methods for estimating relative abundance from survey data in relation to assessment models",
#' Fisheries Research 151(2014) 91-99.
#' @title Calculate survey indices by age.
#' @param x DATRASraw object
#' @param ages vector of ages
#' @param myids haul.ids for grid
#' @param kvecP vector with spatial smoother max. basis dimension for each age group, strictly positive part of model
#' @param kvecZ vector with spatial smoother max. basis dimension for each age group, presence/absence part of model (ignored for Tweedie models)
#' @param gamma model degress of freedom inflation factor (see 'gamma' argument to mgcv::gam() )
#' @param cutOff treat observations below this value as zero
#' @param fam distribution, either "Gamma","LogNormal", or "Tweedie".
#' @param useBIC use BIC for smoothness selection (overrides 'gamma' argument)
#' @param nBoot number of bootstrap samples used for calculating index confidence intervals
#' @param mc.cores number of cores for parallel processing
#' @param method smoothness selection method used by 'gam'
#' @param predD optional DATRASraw object or data.frame (or named list with such objects, one for each year with names(predD) being the years) , defaults to NULL. If not null this is used as grid.
#' @param modelZ vector of model formulae for presence/absence part, one pr. age group (ignored for Tweedie models)
#' @param modelP vector of model formulae for strictly positive repsonses, one pr. age group
#' @param knotsP optional list of knots to gam, strictly positive repsonses
#' @param knotsZ optional list of knots to gam, presence/absence
#' @param predfix optional named list of extra variables (besides Gear, HaulDur, Ship, and TimeShotHour),  that should be fixed during prediction step (standardized)
#' @param linkZ link function for the grDevices::dev.new part of the model, default: "logit" (not used for Tweedie models).
#' @param CIlevel Confidence interval level, defaults to 0.95.
#' @param ... Optional extra arguments to "gam"
#' @return A survey index (list)
#' @examples
#' \dontrun{
#' library(DATRAS) # example data
#' library(sdmgamindex)
#' library(tidyverse)
#' load("data/ebs_example_datras.rda") #object: dat_wrangled. Source: FOSS; see https://emilymarkowitz-noaa.github.io/sdmgamindex/articles/A-data-prep.html for how these data were wrangled in preparation for this example.
#' # Megsie todo: add example here using EBS data!
#' }
#' @importFrom MASS mvrnorm
#' @export
get_surveyidx <- function(x,
                          ages,
                          myids,
                          kvecP = rep(12 * 12, length(ages)),
                          kvecZ = rep(8 * 8, length(ages)),
                          gamma = 1.4,
                          cutOff = 1,
                          fam = "Gamma",
                          useBIC = FALSE,
                          nBoot = 1000,
                          mc.cores = 1,
                          method = "ML",
                          predD = NULL,
                          modelZ = NULL,
                          modelP = NULL,
                          knotsP = NULL,
                          knotsZ = NULL,
                          predfix = NULL,
                          linkZ = "logit",
                          CIlevel = 0.95,
                          ...) {
  if (is.null(modelZ) & is.null(modelP) & fam != "Tweedie") {
    modelZ <-
      modelP <-
      rep("Year+s(lon,lat,k=kvecZ[a],bs='ts')+s(Ship,bs='re',by=dum)+s(Depth,bs='ts')+s(TimeShotHour,bs='cc')",
          length(ages) )
  }
  
  if (is.null(x$Nage)) {
    stop("No age matrix 'Nage' found.")
  }
  if (is.null(colnames(x$Nage))) {
    stop("No colnames found on 'Nage' matrix.")
  }
  if (length(modelP) < length(ages)) {
    stop(" length(modelP) < length(ages)")
  }
  if (length(kvecP) < length(ages)) {
    stop(" length(kvecP) < length(ages)")
  }
  if (fam[1] != "Tweedie") {
    if (length(modelZ) < length(ages)) {
      stop(" length(modelZ) < length(ages)")
    }
    if (length(kvecZ) < length(ages)) {
      stop(" length(kvecZ) < length(ages)")
    }
  }
  ## check for valid family names
  stopifnot(fam[1] %in% c("Gamma", "LogNormal", "Tweedie", "negbin"))
  if (length(fam) < length(ages)) {
    famVec = rep(fam[1], length(ages))
  } else {
    famVec <- fam
  }
  
  dataAges <- as.numeric(gsub("[+]", "", colnames(x$Nage)))
  if (!all(ages %in% dataAges)) {
    stop(paste0(
      "age(s) ",
      setdiff(ages, dataAges),
      " not found in 'Nage' matrix"
    ))
  }
  x[[1]]$Year <- as.factor(x[[1]]$Year)
  x[[2]]$Year <- as.factor(x[[2]]$Year)
  pModels <- zModels <-
    gPreds <-  ##last data year's predictions
    gPreds2 <-  ## all years predictions
    gPreds2.CV <- ## coefficient of variation of previous
    allobs <- ## response vector (zeroes and positive)
    resid <- list() ## residuals
  predDc <- predD ##copy of predD
  
  if (exists(".Random.seed")) {
    oldseed <- get(".Random.seed", .GlobalEnv)
    oldRNGkind <- RNGkind()
    on.exit({
      do.call("RNGkind", as.list(oldRNGkind))
      assign(".Random.seed", oldseed, .GlobalEnv)
    })
  }
  set.seed(314159265)
  
  yearNum <- as.numeric(as.character(x$Year))
  yearRange <- min(yearNum):max(yearNum)
  
  ## Choose most frequent gear as reference gear
  gearNames <- names(stats::xtabs( ~ Gear, data = x[[2]]))
  myGear <-
    names(stats::xtabs( ~ Gear, data = x[[2]]))[which.max(stats::xtabs( ~ Gear, data = x[[2]]))]
  
  if (!is.null(predfix$Gear)) {
    myGear <- predfix$Gear
  }
  
  resMat <-
    matrix(data = NA,
           nrow = length(yearRange),
           ncol = length(ages))
  upMat <- loMat <- resMat
  
  do_one_a <- function(a, cutoff) {
    age = which(dataAges == ages[a])
    ddd <- x[[2]]
    ddd$dum <- 1.0
    ddd$A1 <- ddd$Nage[, age]
    gammaPos <- gamma
    gammaZ <- gamma
    
    if (useBIC) {
      nZ <- nrow(ddd)
      nPos <-
        nrow(ddd[ddd$A1 > cutoff, ]) # nrow(subset(ddd,A1>cutOff))
      gammaPos <- log(nPos) / 2
      gammaZ <- log(nZ) / 2
      cat("gammaPos: ", gammaPos, " gammaZ: ", gammaZ, "\n")
    }
    pd <- ddd[ddd$A1 > cutOff, ] # subset(ddd,A1>cutOff)
    if (famVec[a] == "LogNormal") {
      f.pos <- stats::as.formula(paste("log(A1) ~", modelP[a]))
      
      f.0 <-
        stats::as.formula(paste("A1>", cutOff, " ~", modelZ[a]))
      
      
      print(system.time(m_pos <-
                          DATRAS:::tryCatch.W.E(
                            mgcv::gam(
                              f.pos,
                              data = ddd[ddd$A1 > cutoff, ],
                              #subset(ddd,A1>cutOff),
                              gamma = gammaPos,
                              method = method,
                              knots = knotsP,
                              na.action = stats::na.fail,
                              ...
                            )
                          )$value))
      
      
      if (class(m_pos)[2] == "error") {
        print(m_pos)
        stop(
          "Error occured for age ",
          a,
          " in the positive part of the model\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
      
      print(system.time(m0 <-
                          DATRAS:::tryCatch.W.E(
                            mgcv::gam(
                              f.0,
                              gamma = gammaZ,
                              data = ddd,
                              family = stats::binomial(link = linkZ),
                              method = method,
                              knots = knotsZ,
                              na.action = stats::na.fail,
                              ...
                            )
                          )$value))
      
      
      if (class(m0)[2] == "error") {
        print(m0)
        stop(
          "Error occured for age ",
          a,
          " in the binomial part of the model\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
      
    } else if (famVec[a] == "Gamma") {
      f.pos <- stats::as.formula(paste("A1 ~", modelP[a]))
      
      f.0 <-
        stats::as.formula(paste("A1>", cutOff, " ~", modelZ[a]))
      
      
      print(system.time(m_pos <- DATRAS:::tryCatch.W.E(
        mgcv::gam(
          formula = f.pos,
          data = ddd[ddd$A1 >
                       cutoff, ],
          #subset(ddd,A1>cutOff),
          family = stats::Gamma(link =
                                  "log"),
          gamma = gammaPos,
          method = method,
          knots = knotsP,
          na.action = stats::na.fail,
          ...
        )
      )$value))
      
      
      if (class(m_pos)[2] == "error") {
        print(m_pos)
        stop(
          "Error occured for age ",
          a,
          " in the positive part of the model\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
      
      print(system.time(m0 <-
                          DATRAS:::tryCatch.W.E(
                            mgcv::gam(
                              formula = f.0,
                              gamma = gammaZ,
                              data = ddd,
                              family = stats::binomial(link =
                                                         linkZ),
                              method = method,
                              knots = knotsZ,
                              na.action = stats::na.fail,
                              ...
                            )
                          )$value))
      
      
      if (class(m0)[2] == "error") {
        print(m0)
        stop(
          "Error occured for age ",
          a,
          " in the binomial part of the model\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
    } else if (famVec[a] == "Tweedie") {
      ddd$A1[ddd$A1 < cutOff] = 0
      pd <- ddd
      f.pos <- stats::as.formula(paste("A1 ~", modelP[a]))
      
      print(system.time(m_pos <-
                          DATRAS:::tryCatch.W.E(
                            mgcv::gam(
                              formula = f.pos,
                              data = ddd,
                              family = mgcv::tw,
                              gamma = gammaPos,
                              method = method,
                              knots = knotsP,
                              na.action =
                                stats::na.fail,
                              ...
                            )
                          )$value))
      if (class(m_pos)[2] == "error") {
        print(m_pos)
        stop(
          "Error occured for age ",
          a,
          ".\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
      m0 <- NULL
    } else if (famVec[a] == "negbin") {
      pd <- ddd
      f.pos <- stats::as.formula(paste("A1 ~", modelP[a]))
      
      print(system.time(m_pos <-
                          DATRAS:::tryCatch.W.E(
                            mgcv::gam(
                              formula = f.pos,
                              data = ddd,
                              family = mgcv::nb,
                              gamma = gammaPos,
                              method = method,
                              knots = knotsP,
                              na.action = stats::na.fail,
                              ...
                            )
                          )$value))
      if (class(m_pos)[2] == "error") {
        print(m_pos)
        stop(
          "Error occured for age ",
          a,
          ".\n",
          "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n"
        )
      }
      m0 = NULL
      
    }
    ## Calculate total log-likelihood
    if (famVec[a] == "Tweedie" || famVec[a] == "negbin") {
      totll <- stats::logLik(m_pos)[1]
    } else {
      p0p <- (1 - stats::predict(m0, type = "response"))
      ppos <- p0p[ddd$A1 > cutOff]
      p0m1 <- p0p[ddd$A1 <= cutOff]
      if (famVec[a] == "Gamma")  {
        totll <- sum(log(p0m1)) + sum(log(1 - ppos)) + stats::logLik (m_pos)[1]
      }
      ## if logNormal model, we must transform til log-likelihood to be able to use AIC
      ## L(y) = prod( dnorm( log y_i, mu_i, sigma^2) * ( 1 / y_i ) ) => stats::logLik (y) = sum( log[dnorm(log y_i, mu_i, sigma^2)]  - log( y_i ) )
      if (famVec[a] == "LogNormal") {
        totll <-
          sum(log(p0m1)) + sum(log(1 - ppos)) + stats::logLik (m_pos)[1] - sum(m_pos$y)
      }
    }
    
    if (is.null(predD)) {
      predD <-
        ddd[ddd$haul.id %in% myids, ] # = subset(ddd,haul.id %in% myids)
    }
    res <- numeric(length(yearRange))
    lores <- upres <- res
    gp2 <- gp2_cv <- list()
    
    for (y in levels(ddd$Year)) {
      ## take care of years with all zeroes
      if (!any(ddd$A1[ddd$Year == y] > cutOff)) {
        res[which(as.character(yearRange) == y)] <-
          upres[which(as.character(yearRange) == y)] <-
          lores[which(as.character(yearRange) == y)] <- 0
        next
        
      }
      if (is.list(predDc) &&
          !class(predDc) %in% c("data.frame", "DATRASraw")) {
        predD <- predDc[[as.character(y)]]
      }
      if (is.null(predD)) {
        stop(paste("Year", y, " not found in predD"))
      }
      ## OBS: effects that should be removed should be included here
      predD$Year <- y
      predD$dum <- 0
      predD$ctime <- as.numeric(as.character(y))
      predD$TimeShotHour <- mean(ddd$TimeShotHour)
      predD$Ship <- names(which.max(summary(ddd$Ship)))
      predD$timeOfYear <- mean(ddd$timeOfYear)
      predD$HaulDur <- 30.0
      predD$Gear <- myGear
      
      if (!is.null(predfix)) {
        ##optional extra variables for standardization
        stopifnot(is.list(predfix))
        for (n in names(predfix)) {
          predD[, n] <- predfix[[n]]
        }
      }
      p.1 <- p.0 <- NULL
      try({
        Xp.1 <- stats::predict(m_pos, newdata = predD, type = "lpmatrix")
        
        OS_pos <- numeric(nrow(predD))
        
        terms_pos <- stats::terms(m_pos)
        if (!is.null(m_pos$offset)) {
          off.num.pos <- attr(terms_pos, "offset")
          for (i in off.num.pos) {
            OS_pos <-
              OS_pos + eval(attr(terms_pos, "variables")[[i + 1]], predD)
          }
        }
        p.1 <- Xp.1 %*% stats::coef(m_pos) + OS_pos
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          Xp_0 <- stats::predict(m0, newdata = predD, type = "lpmatrix")
          
          brp_0 <- stats::coef(m0)
          
          OS0 <- numeric(nrow(predD))
          terms_0 <- stats::terms(m0)
          if (!is.null(m0$offset)) {
            off_num_0 <- attr(terms_0, "offset")
            for (i in off_num_0)
              OS0 <- OS0 + eval(attr(terms_0,
                                     "variables")[[i + 1]], predD)
          }
          p.0 <- m0$family$linkinv(Xp_0 %*% brp_0 + OS0)
          
        }
      })
      
      ## take care of failing predictions
      if (!is.numeric(p.1) |
          (!famVec[a] %in% c("Tweedie", "negbin") &&
           !is.numeric(p.0))) {
        res[which(as.character(yearRange) == y)] <-
          upres[which(as.character(yearRange) == y)] <-
          lores[which(as.character(yearRange) == y)] <- 0
        next
        
      }
      sig2 <- m_pos$sig2
      
      if (famVec[a] == "Gamma") {
        res[which(as.character(yearRange) == y)] <-
          sum(p.0 * exp(p.1))
        gPred = p.0 * exp(p.1)
      }
      if (famVec[a] == "LogNormal")  {
        res[which(as.character(yearRange) == y)] <-
          sum(p.0 * exp(p.1 + sig2 / 2))
        gPred = p.0 * exp(p.1 + sig2 / 2)
      }
      if (famVec[a] %in% c("Tweedie", "negbin")) {
        res[which(as.character(yearRange) == y)] <-
          sum(exp(p.1))
        gPred = exp(p.1)
      }
      
      gp2[[y]] <- gPred
      
      if (nBoot > 10) {
        brp_1 <- MASS::mvrnorm(n = nBoot, stats::coef(m_pos), m_pos$Vp)
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          brp_0 <- MASS::mvrnorm(n = nBoot, stats::coef(m0), m0$Vp)
          
          OS0 = matrix(0, nrow(predD), nBoot)
          
          terms_0 <- stats::terms(m0)
          if (!is.null(m0$offset)) {
            off_num_0 <- attr(terms_0, "offset")
            for (i in off_num_0) {
              OS0 <- OS0 + eval(attr(terms_0, "variables")[[i + 1]], predD)
            }
          }
          rep0 = m0$family$linkinv(Xp_0 %*% t(brp_0) + OS0)
        }
        OS_pos <- matrix(0, nrow(predD), nBoot)
        
        terms_pos <- stats::terms(m_pos)
        
        if (!is.null(m_pos$offset)) {
          off.num.pos <- attr(terms_pos, "offset")
          for (i in off.num.pos) {
            OS_pos <-
              OS_pos + eval(attr(terms_pos, "variables")[[i + 1]], predD)
          }
        }
        if (famVec[a] == "LogNormal") {
          rep1 <- exp(Xp.1 %*% t(brp_1) + sig2 / 2 + OS_pos)
        } else {
          rep1 <- exp(Xp.1 %*% t(brp_1) + OS_pos)
        }
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          idxSamp <- base::colSums(rep0 * rep1)
          
          gp_sd <- apply(rep0 * rep1, 1, stats::sd)
        } else {
          idxSamp <- base::colSums(rep1)
          
          gp_sd <- apply(rep1, 1, stats::sd)
        }
        gp2_cv[[y]] <- (gp_sd / gPred)
        
        halpha <- (1 - CIlevel) / 2
        upres[which(as.character(yearRange) == y)] <-
          stats::quantile(idxSamp, 1 - halpha, na.rm = TRUE)
        lores[which(as.character(yearRange) == y)] <-
          stats::quantile(idxSamp, halpha, na.rm = TRUE)
      }
    } ## rof years
    return(list(
      res = res,
      m_pos = m_pos,
      m0 = m0,
      lo = lores,
      up = upres,
      gp = gPred,
      ll = totll,
      pd = pd,
      gp2 = gp2,
      gp2_cv = gp2_cv
    ))
  }## end do.one
  noAges = length(ages)
  
  rr = parallel::mclapply(1:noAges, do_one_a, mc.cores = mc.cores)
  
  logl = 0
  
  for (a in 1:noAges) {
    resMat[, a] = rr[[a]]$res
    
    zModels[[a]] = rr[[a]]$m0
    
    pModels[[a]] = rr[[a]]$m_pos
    
    loMat[, a] = rr[[a]]$lo
    
    upMat[, a] = rr[[a]]$up
    
    gPreds[[a]] = rr[[a]]$gp
    
    logl = logl + rr[[a]]$ll
    gPreds2[[a]] = rr[[a]]$gp2
    gPreds2.CV[[a]] = rr[[a]]$gp2_cv
    allobs[[a]] = x[[2]]$Nage[, a]
  }
  
  get_edf <- function(m) {
    return(sum(m$edf))
  }
  
  totEdf <-
    sum(unlist(lapply(zModels, get_edf))) + sum(unlist(lapply(pModels, get_edf)))
  rownames(resMat) <- yearRange
  colnames(resMat) <- ages
  
  out <- list(
    idx = resMat,
    zModels = zModels,
    pModels = pModels,
    lo = loMat,
    up = upMat,
    gPreds = gPreds,
    logLik  = logl,
    edfs = totEdf,
    gPreds2 = gPreds2,
    gPreds2.CV  =  gPreds2.CV,
    family = famVec,
    cutOff = cutOff,
    dataAges = dataAges,
    yearNum = yearNum,
    refGear = myGear,
    predfix  =  predfix,
    knotsP = knotsP,
    knotsZ = knotsZ,
    allobs = allobs,
    CIlevel = CIlevel
  )
  
  class(out) <- "surveyIdx"
  set.seed(314159265) ## reset seed here (in case multicore option is used)
  for (a in 1:noAges) {
    resid[[a]] = stats::residuals(out, a)
  }
  out$residuals <- resid
  
  return(out)
}




#' Calculate confidence intervals for a named parameter in a survey index model.
#'
#' @title Calculate confidence intervals for a named parameter in a survey index model.
#' @param x survey index
#' @param dat DATRASraw object
#' @param parName name of the parameter, e.g. "Gear"
#' @param cutOff see getsdmgamindex()
#' @param nboot see getsdmgamindex()
#' @param pOnly only calculate for positive part of model, defaults to FALSE.
#' @return list of estimates + ci bounds for each age group.
#' @importFrom MASS mvrnorm
#' @export
get_effect <- function(x,
                       dat,
                       parName = "Gear",
                       cutOff,
                       nboot = 1000,
                       pOnly = FALSE) {
  noAges <- length(x$pModels)
  res <- list()
  
  for (a in 1:noAges) {
    cat("Age ", a, "\n")
    shipSelP <- grep(parName, names(stats::coef(x$pModels[[a]])))
    shipSelZ <- grep(parName, names(stats::coef(x$zModels[[a]])))
    
    dd <- dat[dat$Nage[, a] > cutOff] # subset(dat,Nage[,a]>cutOff)
    zNam <- levels(dat$Ship)
    pNam <- levels(dd$Ship)
    dif <- setdiff(zNam, pNam)
    remo <- which(zNam %in% dif)
    if (length(remo) > 0) {
      shipSelZ <- shipSelZ[-remo]
    }
    
    if (length(shipSelP) != length(shipSelZ)) {
      print("unequal number of ship effects")
    }
    
    brp_1 <-
      MASS::mvrnorm(n = nboot, stats::coef(x$pModels[[a]]), x$pModels[[a]]$Vp)
    brp_0 <-
      MASS::mvrnorm(n = nboot, stats::coef(x$zModels[[a]]), x$zModels[[a]]$Vp)
    
    shipE <- exp(brp_1[, shipSelP, drop = FALSE])
    if (!pOnly) {
      shipE <-
        shipE * x$zModels[[a]]$family$linkinv(brp_0[, shipSelZ, drop = FALSE])
    }
    
    upres <- apply(shipE, 2, stats::quantile, probs = 0.975)
    lores <- apply(shipE, 2, stats::quantile, probs = 0.025)
    
    tmp <- cbind(colMeans(shipE), upres, lores)
    
    
    res[[a]] <- tmp
  }
  return(res)
}





#' Re-compute standardized survey indices for an alternative grid from a previous fitted "surveyIdx" model.
#'
#' @title Re-compute standardized survey indices for an alternative grid from a previous fitted "surveyIdx" model.
#' @param x DATRASraw dataset
#' @param model object of class "surveyIdx" as created by "get_surveyidx"
#' @param predD optional DATRASraw object, defaults to NULL. If not null this is used as grid.
#' @param myids haul.ids for grid
#' @param nBoot number of bootstrap samples used for calculating index confidence intervals
#' @param predfix optional named list of extra variables (besides Gear, HaulDur, Ship, and TimeShotHour),  that should be fixed during prediction step (standardized)
#' @param mc.cores mc.cores number of cores for parallel processing
#' @return An object of class "surveyIdx"
#' @importFrom MASS mvrnorm
#' @export
redo_surveyidx <- function(x,
                           model,
                           predD = NULL,
                           myids,
                           nBoot = 1000,
                           predfix = list(),
                           mc.cores = 1) {
  ages <- as.numeric(colnames(model$idx))
  dataAges <- model$dataAges
  famVec <- model$family
  cutOff <- model$cutOff
  
  yearNum <- model$yearNum
  yearRange <- min(yearNum):max(yearNum)
  
  
  gPreds <- list() ##last data year's predictions
  gPreds2 <- list() ## all years predictions
  gPreds2.CV <- list()
  predDc <- predD
  
  myGear <- model$refGear
  
  resMat <- matrix(NA, nrow = length(yearRange), ncol = length(ages))
  upMat <- resMat <- resMat
  
  do_one_a <- function(a) {
    age <- which(dataAges == ages[a])
    ddd <- x[[2]]
    ddd$dum <- 1.0
    ddd$A1 <- ddd$Nage[, age]
    m_pos <- model$pModels[[a]]
    m0 = NULL
    if (!famVec[a] %in% c("Tweedie", "negbin")) {
      m0 <- model$zModels[[a]]
    }
    if (is.null(predD)) {
      predD <-
        ddd[ddd$haul.id %in% myids, ] #=subset(ddd,haul.id %in% myids)
    }
    res <- numeric(length(yearRange))
    lores <- upres <- res
    gp2 <- gp2_cv <- list()
    
    do_one_y <- function(y) {
      cat("Doing year ", y, "\n")
      
      if (is.list(predDc) &&
          !class(predDc) %in% c("data.frame", "DATRASraw")) {
        predD <- predDc[[as.character(y)]]
      }
      if (is.null(predD)) {
        stop(paste("Year", y, " not found in predD"))
      }
      
      ## take care of years with all zeroes
      if (!any(ddd$A1[ddd$Year == y] > cutOff)) {
        return(list(
          res = 0,
          upres = 0,
          lores = 0,
          gp2 = NULL,
          gp2_cv = NULL
        ))
      }
      
      ## OBS: effects that should be removed should be included here
      predD$Year <- y
      predD$dum <- 0
      predD$ctime <- as.numeric(as.character(y))
      
      predD$TimeShotHour <- mean(ddd$TimeShotHour)
      predD$Ship <- names(which.max(summary(ddd$Ship)))
      predD$timeOfYear <- mean(ddd$timeOfYear)
      
      predD$HaulDur <- 30.0
      predD$Gear <- myGear
      
      if (!is.null(predfix)) {
        ##optional extra variables for standardization
        stopifnot(is.list(predfix))
        for (n in names(predfix)) {
          predD[, n] <- predfix[[n]]
        }
      }
      p.1 <- p.0 <- NULL
      try({
        Xp.1 = stats::predict(m_pos, newdata = predD, type = "lpmatrix")
        
        OS_pos = numeric(nrow(predD))
        
        terms_pos = stats::terms(m_pos)
        if (!is.null(m_pos$offset)) {
          off.num.pos <- attr(terms_pos, "offset")
          for (i in off.num.pos)
            OS_pos <- OS_pos +
              eval(attr(terms_pos, "variables")[[i + 1]], predD)
        }
        p.1 <- Xp.1 %*% stats::coef(m_pos) + OS_pos
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          Xp_0 <- stats::predict(m0, newdata = predD, type = "lpmatrix")
          brp_0 <- stats::coef(m0)
          OS0 <- numeric(nrow(predD))
          terms_0 <- stats::terms(m0)
          if (!is.null(m0$offset)) {
            off_num_0 <- attr(terms_0, "offset")
            for (i in off_num_0) {
              OS0 <- OS0 + eval(attr(terms_0, "variables")[[i + 1]], predD)
            }
          }
          p.0 <- m0$family$linkinv(Xp_0 %*% brp_0 + OS0)
          
        }
      })
      
      ## take care of failing predictions
      if (!is.numeric(p.1) |
          (!famVec[a] %in% c("Tweedie", "negbin") && !is.numeric(p.0))) {
        return(list(
          res = 0,
          upres = 0,
          lores = 0,
          gp2 = NULL,
          gp2_cv = NULL
        ))
      }
      
      sig2 <- m_pos$sig2
      idx <- NA
      if (famVec[a] == "Gamma") {
        idx <- sum(p.0 * exp(p.1))
        gPred <- p.0 * exp(p.1)
      }
      if (famVec[a] == "LogNormal")  {
        idx <- sum(p.0 * exp(p.1 + sig2 / 2))
        gPred <- p.0 * exp(p.1 + sig2 / 2)
      }
      if (famVec[a] %in% c("Tweedie", "negbin"))  {
        idx <- sum(exp(p.1))
        gPred <- exp(p.1)
      }
      
      if (nBoot > 10) {
        brp_1 <- MASS::mvrnorm(n = nBoot, stats::coef(m_pos), m_pos$Vp)
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          brp_0 <- MASS::mvrnorm(n = nBoot, stats::coef(m0), m0$Vp)
          
          OS0 <- matrix(0, nrow(predD), nBoot)
          
          terms_0 = stats::terms(m0)
          if (!is.null(m0$offset)) {
            off_num_0 <- attr(terms_0, "offset")
            
            for (i in off_num_0) {
              OS0 <- OS0 + eval(attr(terms_0, "variables")[[i + 1]], predD)
            }
          }
          rep0 <- m0$family$linkinv(Xp_0 %*% t(brp_0) + OS0)
          
        }
        OS_pos <- matrix(0, nrow(predD), nBoot)
        
        terms_pos <- stats::terms(m_pos)
        if (!is.null(m_pos$offset)) {
          off.num.pos <- attr(terms_pos, "offset")
          for (i in off.num.pos) {
            OS_pos <- OS_pos + eval(attr(terms_pos, "variables")[[i + 1]], predD)
          }
        }
        
        if (famVec[a] == "LogNormal") {
          rep1 <- exp(Xp.1 %*% t(brp_1) + sig2 / 2 + OS_pos)
          
        } else {
          rep1 <- exp(Xp.1 %*% t(brp_1) + OS_pos)
          
        }
        
        if (!famVec[a] %in% c("Tweedie", "negbin")) {
          idxSamp <- base::colSums(rep0 * rep1)
          
          gp_sd <- apply(rep0 * rep1, 1, stats::sd)
          gp_cv <- gp_sd / gPred
        } else {
          idxSamp <- base::colSums(rep1)
          
          gp_sd <- apply(rep1, 1, stats::sd)
          gp_cv <- gp_sd / gPred
        }
        halpha <- (1 - model$CIlevel) / 2
        
        return(
          list(
            res = idx,
            upres = stats::quantile(idxSamp, 1 - halpha, na.rm = TRUE),
            lores = stats::quantile(idxSamp, halpha, na.rm = TRUE),
            gp2 = gPred,
            gp2_cv = gp_cv
          )
        )
      }
    } ## rof years
    yres <-
      parallel::mclapply(levels(ddd$Year), do_one_y, mc.cores = mc.cores)
    
    for (y in levels(ddd$Year)) {
      ii = which(as.character(yearRange) == y)
      res[ii] = yres[[ii]]$res
      upres[ii] = yres[[ii]]$upres
      lores[ii] = yres[[ii]]$lores
      gp2[[y]] = yres[[ii]]$gp2
      gp2_cv[[y]] = yres[[ii]]$gp2_cv
    }
    
    return(
      list(
        res = res,
        m_pos = m_pos,
        m0 = m0,
        lo = lores,
        up = upres,
        gp = utils::tail(gp2, 1),
        gp2 = gp2,
        gp2_cv = gp2_cv
      )
    )
  }## end do.one
  
  rr <- lapply(1:length(age), do_one_a)
  logl = 0
  
  for (a in 1:length(age)) {
    resMat[, a] <- rr[[a]]$res
    loMat[, a] <- rr[[a]]$lo
    upMat[, a] <- rr[[a]]$up
    gPreds[[a]] <- rr[[a]]$gp
    gPreds2[[a]] <- rr[[a]]$gp2
    gPreds2.CV[[a]] <- rr[[a]]$gp2_cv
  }
  rownames(resMat) <- yearRange
  colnames(resMat) <- ages
  
  out <- list(
    idx = resMat,
    zModels = model$zModels,
    pModels = model$pModels,
    lo = loMat,
    up = upMat,
    gPreds = gPreds,
    logLik = model$logLik,
    edfs = model$edfs,
    pData = model$pData,
    gPreds2 = gPreds2,
    gPreds2.CV = gPreds2.CV,
    family = famVec,
    cutOff = cutOff,
    dataAges = dataAges,
    yearNum = yearNum,
    refGear = myGear,
    predfix = predfix,
    knotsP = model$knotsP,
    knotsZ = model$knotsZ
  )
  
  class(out) <- "surveyIdx"
  
  return(out)
}


#' Survey index using the stratified mean method using ICES statistical rectangles as strata.
#'
#' Previously named get_surveyidxStratMean.
#'
#' @title Survey index using the stratified mean method using ICES statistical rectangles as strata.
#' @param x DATRASraw object. Must contain a matrix: x[[2]]$Nage.
#' @param ageCols which columns of the Nage matrix should be included?
#' @param doLog log-transform?
#' @return a matrix with survey indices
#' @export
get_surveyidx_stratmean <-
  function(x, ageCols, doLog = FALSE) {
    ysplit <- split(x, x$Year)
    
    res <- matrix(NA, nrow = length(ysplit), ncol = length(ageCols))
    for (y in 1:length(ysplit)) {
      if (!doLog) {
        byRec <- stats::aggregate(
          x = ysplit[[y]][[2]]$Nage[, ageCols],
          by = list(ysplit[[y]][[2]]$StatRec),
          FUN = "mean")
      } else {
        byRec <- stats::aggregate(
          x = log(ysplit[[y]][[2]]$Nage[, ageCols] + 1),
          by = list(ysplit[[y]][[2]]$StatRec),
          FUN = "mean")
      }
      
      res[y, ] <- colMeans(byRec[, -1])
    }
    return(res)
  }



# Post-process survey index ----------------------------------------------------


#' Likelihood ratio test for comparing two survey indices.
#'
#' Previously called anova_SI. TOLDEO.
#'
#' @param m1 model 1
#' @param m2 model 2
#'
#' @return A p-value.
#' @export
#'
#' @examples
#' # Simple example
#' dat <- datasets::airquality
#' m1 <- mgcv::gam(Ozone ~ s(Solar.R), data = dat)
#' m1
#' m2 <- mgcv::gam(Ozone ~ Wind + s(Solar.R), data = dat)
#' m2
#' # anova_likelihood_ratio_test(m1 = m1, m2 = m2)
anova_likelihood_ratio_test <- function(m1, m2) {
  ll1 <- m1$logLik
  ll2 <- m2$logLik
  edfs1 <- m1$edfs
  edfs2 <- m2$edfs
  
  if (edfs2 > edfs1) {
    out <- 1 - stats::pchisq(2 * (ll2 - ll1), edfs2 - edfs1)
  } else {
    out <- 1 - stats::pchisq(2 * (ll1 - ll2), edfs1 - edfs2)
  }
  return(out)
}




#' Concentration transform
#'
#' Previously called concTransform. Helper function for plotting survey indices.
#' This analysis makes it easier to transform X values that are concentrations.
#' Because it is common to transform X values to their logarithms, which is
#' required before fitting some models to your data. Since the logarithm of zero
#' is undefined, if you enter X=0 that value will be empty (missing) after
#' transformation. This analysis  lets you substitute some other value
#' (a tiny concentration) for zero before taking the logarithm.
#'
#' @param x a vector of log-responses
#'
#' @return vector of transformed responses
#' @export
#'
#' @examples
#' plot(concentration_transform(x = seq(1,100,5)))
#' plot(concentration_transform(x = stats::rnorm(n = 123, mean = 5, sd = 25)))
concentration_transform <- function (x) {
  i <- order(x)
  ys <- sort(exp(x))
  p <- ys / sum(ys)
  x[i] <- cumsum(p)
  return(x)
}


#' Akaike Information Criterion (or BIC) for survey index models
#'
#' previously named AIC.surveyIdx.
#'
#' @title Akaike Information Criterion (or BIC) for survey index models
#' @param x survey index as return from get_surveyidx
#' @param BIC if TRUE compute BIC instead of AIC
#' @return numeric value
#' @export
get_surveyidx_aic <- function(x, BIC = FALSE) {
  if (!BIC) {
    out <- (2 * x$edfs - 2 * x$logLik)
  }
  
  if (pmatch("Tweedie", x$pModels[[1]]$family$family, nomatch = -1) == 1 ||
      pmatch("negbin", x$pModels[[1]]$family$family, nomatch = -1) == 1) {
    out <-
      log(length(x$pModels[[1]]$y) * (length(x$pModels))) * x$edfs - 2 * x$logLik
  } else {
    out <-
      log(length(x$zModels[[1]]$y) * (length(x$zModels))) * x$edfs - 2 * x$logLik
  }
  
  return(out)
}



#' Helper function to "borrow" missing age groups from other years
#'
#' In years where there are less than 'n' individuals of age 'age',
#' add fake individuals of that age such that there are 'n'.
#' The length of the individuals are set to the mean (or whatever 'fun' specifies)
#' of all other individuals of the same age.
#' For the minimum and maximum age groups fun it is reasonable to replace 'mean' with 'min' and 'max' respectively.
#' Note, that you might need to call 'addSpectrum' on the object again.
#' @title Helper function to "borrow" missing age groups from other years
#' @param x DATRASraw object
#' @param age age to impute
#' @param n at least this many individuals in each year
#' @param fun A function such as 'mean','median','min', or 'max'.
#' @return a DATRASraw object
#' @export
fix_age_group <- function(x,
                          age = 0,
                          n = 3,
                          fun = "mean") {
  
  cm.breaks <- attr(x, "cm.breaks")
  f <- match.fun(fun)
  d <- split(x, x$Year)
  subsLength <- f(x[[3]]$LngtCm, na.rm = TRUE)
  for (y in 1:length(d)) {
    nobs <- sum(d[[y]][[1]]$Age == age, na.rm = TRUE)
    if (nobs < n) {
      sel <- sample(1:nrow(d[[y]][[1]]), n - nobs)
      d[[y]][[1]] <- rbind(d[[y]][[1]][sel, ], d[[y]][[1]])
      
      d[[y]][[1]][1:(n - nobs), "Age"] <- age
      
      d[[y]][[1]][1:(n - nobs), "LngtCm"] <- subsLength
      
      d[[y]][[1]][1:(n - nobs), "NoAtALK"] <- 1
      
    }
  }
  dd <- do.call("c", d)
  if (!is.null(cm.breaks)) {
    dd <- DATRAS::addSpectrum(dd, cm.breaks = cm.breaks)
  }
  return(dd)
}



#' Calculate internal consistency of a survey index.
#'
#' Previously called internalCons.
#'
#' @title Calculate internal consistency of a survey index.
#' @param tt A matrix with survey indices (rows=years, cols=ages)
#' @param print_plot Plot it?
#' @return a vector of consistencies
#' @export
consistency_internal <- function(tt, print_plot = FALSE) {
  tt[tt == 0] <- NA
  sel1 <- 1:(nrow(tt) - 1)
  sel2 <- 2:nrow(tt)
  
  if (print_plot) {
    grDevices::dev.new()
    b <- ceiling((ncol(tt) - 1) / 2)
    graphics::par(mfrow = c(b, b))
  }
  for (a in 1:(ncol(tt) - 1)) {
    cat("Age ", a," vs ", a + 1, " : ",
        stats::cor(log(tt[sel1, a]), log(tt[sel2, a + 1]),
                   use = "pairwise.complete.obs"), "\n")
    if (print_plot) {
      plot(x = log(tt[sel1, a]), y = log(tt[sel2, a + 1]))
      graphics::abline(0, 1)
    }
  }
  return(sapply(1:(ncol(tt) - 1), function(a)
    stats::cor(log(tt[sel1, a]), log(tt[sel2, a + 1]), use = "pairwise.complete.obs")))
  
}


#' Calculate external consistencies between two survey indices.
#'
#' Previously called externalCons.
#'
#' Proper alignment of years and ages must be ensured by the user.
#' @title Calculate external consistencies between two survey indices.
#' @param tt A matrix with survey indices (rows=years, cols=ages)
#' @param tt2 A matrix with survey indices (rows=years, cols=ages)
#' @param print_plot plot it?
#' @return A vector of correlations (consistencies)
#' @export
consistency_external <- function(tt, tt2, print_plot = FALSE) {
  tt[tt == 0] <-
    tt2[tt2 == 0] <- NA
  if (print_plot) {
    grDevices::dev.new()
    b = ceiling((ncol(tt)) / 2)
    graphics::par(mfrow = c(b, b))
  }
  for (a in 1:ncol(tt)) {
    cat( "Survey 1 Age ", a, " vs Survey 2 ", a, " : ",
         stats::cor(log(tt[, a]), log(tt2[, a]),
                    use = "pairwise.complete.obs"), "\n")
    if (print_plot) {
      plot(x = log(tt[, a]), y = log(tt2[, a]))
      graphics::abline(0, 1)
    }
  }
  out <- sapply(1:ncol(tt),
                function(a)
                  stats::cor(
                    x = log(tt[, a]),
                    y = log(tt2[, a]),
                    use = "pairwise.complete.obs"))
  return(out)
}



# Simulate ---------------------------------------------------------------------


#' Simulate data from a surveyIdx model (experimental and subject to change)
#'
#' Previously named simulate surveyIdx.simulate.
#'
#' @title Simulate data from a sdmgamindex model (experimental and subject to change)
#' @param model object of class 'surveyIdx'
#' @param d A dataset (DATRASraw object)
#' @param sampleFit Use a random sample from the gaussian approximation to distribution of the estimated parameter vector. Default: FALSE.
#' @param condSim optional results of previous call to this function. Use this if you want to generate many datasets (much faster, since mean predictions are re-used).
#' @return list with  1) simulated observations with noise 2) mean (no noise) 3) zero probability.
#' @export
get_surveyidx_sim <- function(model,
                              d,
                              sampleFit = FALSE,
                              condSim = NULL) {
  ages <- as.numeric(colnames(model$idx))
  dataAges <- model$dataAges
  famVec <- model$family
  
  out <- d$Nage
  out_mu <- list()
  out_mu0 <- list()
  
  for (a in 1:length(ages)) {
    if (sampleFit && is.null(condSim)) {
      m_pos <- model$pModels[[a]]
      
      Xp.1 <- stats::predict(object = m_pos,
                             newdata = d[[2]],
                             type = "lpmatrix")
      brp_1 <- MASS::mvrnorm(n = 1,
                             mu = stats::coef(m_pos),
                             Sigma = m_pos$Vp)
      OS_pos <- matrix(data = 0,
                       nrow = nrow(d[[2]]),
                       ncol = 1)
      terms_pos <- stats::terms(m_pos)
      
      if (!is.null(m_pos$offset)) {
        off.num.pos <- attr(terms_pos, "offset")
        for (i in off.num.pos) {
          OS_pos <-
            OS_pos + eval(attr(terms_pos, "variables")[[i + 1]], d[[2]])
        }
      }
      mu <- Xp.1 %*% brp_1 + OS_pos
      out_mu[[a]] <- exp(mu) ## obs, assumes log link!
      
      if (!famVec[a] %in% c("Tweedie", "negbin")) {
        m0 <- model$zModels[[a]]
        Xp_0 <- stats::predict(object = m0,
                               newdata = d[[2]],
                               type = "lpmatrix")
        brp_0 <- MASS::mvrnorm(n = 1,
                               mu = stats::coef(m0),
                               Sigma = m0$Vp)
        OS0 <- matrix(0, nrow(d[[2]]), 1)
        
        terms_0 <- stats::terms(m0)
        if (!is.null(m0$offset)) {
          off_num_0 <- attr(terms_0, "offset")
          for (i in off_num_0) {
            OS0 <- OS0 + eval(attr(terms_0, "variables")[[i + 1]], d[[2]])
          }
        }
        mu0 <- m0$family$linkinv(Xp_0 %*% brp_0 + OS0)
        
        out_mu0[[a]] = mu0
      }
      
    }
    if (!is.null(condSim)) {
      out_mu[[a]] <- condSim[[2]][[a]]
    }
    if (!is.null(condSim) &&
        !famVec[a] %in% c("Tweedie", "negbin")) {
      out_mu0[[a]] <- condSim[[3]][[a]]
    }
    
    if (famVec[a] == c("Tweedie")) {
      p <- model$pModels[[a]]$family$getTheta(TRUE)
      phi <- model$pModels[[a]]$scale
      if (!sampleFit && is.null(condSim)) {
        out_mu[[a]] <- stats::predict(
          object = model$pModels[[a]],
          newdata = d[[2]],
          type = "response"
        )
      }
      out[, a] <- mgcv::rTweedie(mu = out_mu[[a]], p = p, phi = phi)
    } else if (famVec[a] == "LogNormal") {
      if (!sampleFit && is.null(condSim)) {
        out_mu0[[a]] <- stats::predict(
          object = model$zModels[[a]],
          newdata = d[[2]],
          type = "response"
        )
      }
      pa <- stats::rbinom(nrow(out), 1, prob = out_mu0[[a]])
      sig2 <- model$pModels[[a]]$sig2
      if (!sampleFit && is.null(condSim)) {
        out_mu[[a]] <- exp(stats::predict(object = model$pModels[[a]],
                                          newdata = d[[2]]))
      }
      pos <- exp(x = stats::rnorm(
        n = nrow(out),
        mean = log(out_mu[[a]]),
        sd = sqrt(sig2)
      ))
      out[, a] <- pa * pos
    } else if (famVec[a] == "Gamma") {
      if (!sampleFit && is.null(condSim)) {
        out_mu0[[a]] <- stats::predict(
          object = model$zModels[[a]],
          newdata = d[[2]],
          type = "response"
        )
      }
      pa <- stats::rbinom(n = nrow(out),
                          size = 1,
                          prob = out_mu0[[a]])
      
      if (!sampleFit && is.null(condSim)) {
        out_mu[[a]] <- stats::predict(
          object = model$pModels[[a]],
          newdata = d[[2]],
          type = "response"
        )
      }
      shape <- (1 / model$pModels[[a]]$scale)
      pos <- stats::rgamma(n = nrow(out),
                           rate = shape / out_mu[[a]],
                           shape = shape)
      out[, a] <- pa * pos
    }
    
  }
  return(list(sim = out,
              mu = out_mu,
              prob0 = out_mu0))
}


##' Plot survey index list (e.g. retrospective analysis)
##'
##' @title Plot survey index list (e.g. retrospective analysis)
##' @param x (named) list of "surveyIdx" objects for example from "retro.surveyIdx" or "leaveout.surveyIdx"
##' @param base Either index of x that should considered the "base run" (integer), OR object of class "surveyIdx". Confidence bounds will be shown for this model only.
##' @param rescale Should indices be rescaled to have mean 1 (over the set of intersecting years)? Default: FALSE
##' @param lwd line width argument to plot
##' @param main if not NULL override main plotting default title of "Age group a"
##' @param allCI show 95\% confidence lines for all indices? Default FALSE.
##' @param includeCI Show confidence intervals? Default TRUE.
##' @param ylim Y axis range. If NULL (default) then determine automatically.
##' @return nothing
##' @export
plot_simulation_list <-
  function(x,
           base = 1,
           rescale = FALSE,
           lwd = 1.5,
           main = NULL,
           allCI = FALSE,
           includeCI = TRUE,
           ylim = NULL) {
    
    if (class(base) == "surveyIdx") {
      x <- c(list(base), x)
      base <- 1
    }
    stopifnot(is.numeric(base))
    nx <- length(x)
    mainwasnull <- is.null(main)
    n <- ncol(x[[base]]$idx)
    if (n > 1) {
      op <- par(mfrow = n2mfrow(n))
      on.exit(par(op))
    }
    
    cols <- rainbow(nx)
    if (nx == 2)
      cols <- 2:3
    cols[base] <- "black"
    allyears <- lapply(x, function(x) rownames(x$idx))
    rsidx <- 1:nrow(x[[nx]]$idx)
    if (rescale) {
      commonyears <- allyears[[1]]
      if (nx > 1) {
        for (i in 2:nx) {
          commonyears <- intersect(commonyears, allyears[[i]])
        }
        if (length(commonyears) == 0)
          stop("rescaling not possible because the set of common years is empty")
      }
    }
    
    ss <- ssbase <- 1
    
    for (aa in 1:n) {
      rangevec <- x[[1]]$idx[, aa]
      for (xx in 2:nx)
        rangevec <- c(rangevec, x[[xx]]$idx[, aa])
      if (includeCI) {
        for (xx in 1:nx)
          rangevec <- c(rangevec, x[[xx]]$lo[, aa], x[[xx]]$up[, aa])
      }
      
      yl = range(rangevec)
      if (rescale) {
        rsidx <- which(rownames(x[[base]]$idx) %in% commonyears)
        ssbase <- mean(x[[base]]$idx[rsidx, aa], na.rm = TRUE)
        yl <- yl / ssbase
      }
      if (!is.null(ylim)){
        yl = ylim
      }
      if (mainwasnull) {
        main <- paste("Age group", colnames(x[[base]]$idx)[aa])
        y <- as.numeric(rownames(x[[base]]$idx))
        plot(
          y,
          x[[base]]$idx[, aa] / ssbase,
          type = "b",
          ylim = yl,
          main = main,
          xlab = "Year",
          ylab = "Index"
        )
      }
      if (includeCI) {
        polygon(c(y, rev(y)),
                c(x[[base]]$lo[, aa], rev(x[[base]]$up[, aa])) / ssbase,
                col = "lightgrey",
                border = NA)
      }
      for (i in 1:length(x)) {
        y <- as.numeric(rownames(x[[i]]$idx))
        if (rescale) {
          rsidx <- which(rownames(x[[i]]$idx) %in% commonyears)
          ss <- mean(x[[i]]$idx[rsidx, aa], na.rm = TRUE)
        }
        lines(y,
              x[[i]]$idx[, aa] / ss,
              col = cols[i],
              type = "b",
              lwd = lwd)
        
        if (includeCI && allCI && i != base) {
          lines(
            y,
            x[[i]]$lo[, aa] / ss,
            col = cols[i],
            lwd = lwd * 0.6,
            lty = 2
          )
          lines(
            y,
            x[[i]]$up[, aa] / ss,
            col = cols[i],
            lwd = lwd * 0.6,
            lty = 2
          )
        }
        
      }
      y <- as.numeric(rownames(x[[base]]$idx))
      lines(y, x[[base]]$idx[, aa] / ssbase, type = "b", lwd = lwd)
      
    }
    if (!is.null(names(x))) {
      legend(
        "topleft",
        legend = names(x),
        col = cols,
        lty = 1,
        lwd = lwd,
        pch = 1
      )
    }
  }

# Randomized quantile residuals ------------------------------------------------

#' Randomized quantile residuals for class 'sdmgamindex'
#'
#' Previously named residuals.surveyIdx.
#'
#' @title Randomized quantile residuals for class 'sdmgamindex'
#' @param x An object of type 'sdmgamindex' as created by 'get_surveyidx'
#' @param a age group
#' @return A vector of residuals, which should be iid standard normal distributed
#' @export
get_surveyidx_resid <- function(x, a = 1) {
  if (pmatch("Tweedie", x$pModels[[a]]$family$family,
             nomatch = -1) == 1) {
    resi <- qres_tweedie(x$pModels[[a]])
  } else if (pmatch("gaussian", x$pModels[[a]]$family$family,
                    nomatch = -1) == 1) {
    ## Delta-Lognormal
    y <- x$allobs[[a]]
    logy <- y
    psel <- (y > x$cutOff)
    logy[psel] <- log(y[psel])
    
    p0 <- stats::predict(x$zModels[[a]], type = "response") ## P( obs > cutOff )
    vari <- x$pModels[[a]]$sig2
    cdfpos <- rep(0, length(y))
    cdfpos[psel] <- stats::pnorm(
      q = logy[psel],
      mean = stats::predict(x$pModels[[a]], type = "link"),
      sd = sqrt(vari)
    )
    u <- (1 - p0) + cdfpos * p0
    u[!psel] <- stats::runif(sum(!psel), min = 0, max = u[!psel])
    resi <- stats::qnorm(u)
    
  } else if (pmatch("Gamma", x$pModels[[a]]$family$family,
                    nomatch = -1) == 1) {
    # requireNamespace("MASS")
    y <- x$allobs[[a]]
    psel <- y > x$cutOff
    
    p0 <- stats::predict(x$zModels[[a]], type = "response")
    means <- stats::predict(x$pModels[[a]], type = "response")
    shape <- MASS::gamma.shape(x$pModels[[a]])$alpha
    cdfpos <- rep(0, length(y))
    cdfpos[psel] <- stats::pgamma(q = y[psel],
                                  shape = shape,
                                  rate = shape / means)
    u <- (1 - p0) + cdfpos * p0
    u[!psel] <- stats::runif(sum(!psel), min = 0, max = u[!psel])
    resi <- stats::qnorm(u)
    
  } else if (pmatch("Negative Binomial", x$pModels[[a]]$family$family,
                    nomatch = -1) == 1) {
    y <- x$pModels[[a]]$y
    size <- x$pModels[[a]]$family$getTheta(TRUE)
    mu <- stats::fitted(x$pModels[[a]])
    p <- size / (mu + size)
    a <- ifelse(y > 0, stats::pbeta(p, size, pmax(y, 1)), 0)
    b <- stats::pbeta(p, size, y + 1)
    u <- stats::runif(n = length(y),
                      min = a,
                      max = b)
    resi <- stats::qnorm(u)
  }
  
  return(resi)
}

#' Randomized quantile residuals for Tweedie models
#'
#' Previously named qres.tweedie.
#'
#' @title Randomized quantile residuals for Tweedie models
#' @param gam.obj A gam object (mgcv package)
#' @return A vector of residuals, which should be iid standard normal distributed
#' @export
#' @import tweedie
qres_tweedie <- function (gam.obj) {
  
  # requireNamespace("tweedie")
  mu <- stats::fitted(gam.obj)
  y <- gam.obj$y
  df <- gam.obj$df.residual
  w <- gam.obj$prior.weights
  if (is.null(w))
    w <- rep(1, length(y))
  p <- gam.obj$family$getTheta(TRUE)
  dispersion <- gam.obj$scale
  u <- numeric(length(y))
  fitt <- stats::fitted(gam.obj)
  for (i in 1:length(u)) {
    u[i] <-
      tweedie::ptweedie(
        q = y[i],
        power = p,
        mu = fitt[i],
        phi = dispersion / w[i]
      )
  }
  
  if (p > 1 && p < 2) {
    u[y == 0] <- stats::runif(sum(y == 0), min = 0, max = u[y == 0])
  }
  
  return(stats::qnorm(u))
}


# Standardized model outputs ---------------------------------------------------




#' Plot a sdmgamindexGrid
#'
#' @param grid a sdmgamindexGrid (as created by the "get_grid" function)
#' @param pch Inherited from base::plot(). plotting ‘character’, i.e., symbol to use.
#' @param gridCol Color of grid plot output.
#'
#' @return nothing
#' @export
#'
#' @examples
#' # TOLEDO
plot_surveyidx_grid <- function(grid,
                                pch = 1,
                                gridCol = "lightgrey") {
  plot(grid[[1]],
       grid[[2]],
       xlab = "Longitude",
       ylab = "Latitude",
       pch = pch)
  lonbps <- grid[[4]]
  latbps <- grid[[5]]
  
  for (i in 1:nrow(lonbps)) {
    graphics::abline(v = lonbps[i, 1], col = gridCol)
  }
  
  graphics::abline(v = utils::tail(lonbps, 1)[2], col = gridCol)
  
  for (i in 1:nrow(latbps)) {
    graphics::abline(h = latbps[i, 1], col = gridCol)
  }
  
  graphics::abline(h = utils::tail(latbps, 1)[2], col = gridCol)
  
  maps::map(
    "worldHires",
    fill = TRUE,
    plot = TRUE,
    add = TRUE,
    col = grDevices::grey(0.5)  )
}

#' Write survey index to file in standard XSA/SAM format
#'
#' Previously named exportSI.
#'
#' @title Write survey index to file in standard XSA/SAM format
#' @param x matrix with survey indices
#' @param ages vector of ages
#' @param years vector of years
#' @param toy fraction of year the survey is conducted (between 0 and 1)
#' @param file filename to write to
#' @param nam file description header
#' @return nothing
#' @export
export_surveyidx <- function(x,
                             ages,
                             years,
                             toy,
                             file,
                             nam = "") {
  cat(nam, "\n", file = file)
  cat(range(as.numeric(as.character(years))),
      "\n",
      file = file,
      append = TRUE)
  cat("1 1 ", rep(toy, 2), "\n", file = file, append = TRUE)
  cat(min(ages), max(ages), "\n",
      file = file,
      append = TRUE)
  utils::write.table(
    round(cbind(1, x[, ]), 4),
    file = file,
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE )
}

#' Visualize results from a survey index model fitted with get_surveyidx().
#'
#' @title Visualize results from a survey index model fitted with get_surveyidx().
#' @param x Survey index as produced by getsdmgamindex()
#' @param dat DATRASraw object
#' @param alt.idx optional matrix with alternative index
#' @param myids vector of haul ids that constitute the grid
#' @param cols which age columns to consider?
#' @param select character vector of chosen plots. Either one of "index","map","absolutemap","CVmap","residuals","fitVsRes",""resVsYear","resVsShip","spatialResiduals", or a number. Numbers refer to smooths in the order they appear in the formula.
#' @param par 'par' settings for plotting (a named list).
#' @param colors colors for spatial effect.
#' @param map.cex size of grid points on maps
#' @param plotByAge boolean (default=TRUE). If true, par(par) is called for each age group.
#' @param legend boolean (default=TRUE). add legends to plot?
#' @param predD DATRASraw object with grid (optional). Overrides 'myids' if supplied.
#' @param year numeric scalar or vector (default=NULL). If 'select' equals 'map' a specific year can be chosen (only meaningful for time-varying spatial effects). If select equals 'absolutemap' or 'CVmap' then year must be a vector.
#' @param main optional main title (overrides default title)
#' @param legend.signif Number of significant digits in map legends
#' @param legend.pos Position of legend (e.g. "bottomleft") see ?legend
#' @param restoreOldPar restore old par() on exit? Default=FALSE
#' @param mapBubbles boolean (default=FALSE) add observation bubbles?
#' @param cutp optional vector of break points for the color scale on maps
#' @param ... Additional parameters for plot()
#'
#' @return series of production plots
#' @import maps mapdata tweedie
#' @export
#'
#' @examples
#' # plot_surveyidx()
plot_surveyidx <- function (x,
                            dat,
                            alt.idx = NULL,
                            myids,
                            cols = 1:length(x$pModels),
                            select = c("index", "map", "residuals", "fitVsRes"),
                            par = list(mfrow = c(3, 3)),
                            colors = rev(grDevices::heat.colors(6)),
                            map.cex = 1,
                            plotByAge = TRUE,
                            legend = TRUE,
                            predD = NULL,
                            year = NULL,
                            main = NULL,
                            legend.signif = 3,
                            legend.pos = "topright",
                            restoreOldPar = FALSE,
                            mapBubbles = FALSE,
                            cutp = NULL,
                            ...) {
  
  if (!plotByAge & !is.null(par)) {
    op <- par(par)
    if (restoreOldPar){
      on.exit(par(op))
    }
  }
  mainwasnull <- is.null(main)
  for (a in cols) {
    if (mainwasnull) {
      main <- paste("Age group", colnames(dat$Nage)[a])
    }
    if (plotByAge & !is.null(par)) {
      op <- par(par)
      if (restoreOldPar) {
        on.exit(par(op))
      }
    }
    if (any(select == "index")) {
      ys <- range(as.numeric(levels(dat$Year)))
      ys <- ys[1]:ys[2]
      yl <- range(c(x$idx[, a], 0, x$lo[, a], x$up[, a]) /
                    mean(x$idx[, a]), na.rm = TRUE)
      if (!is.null(alt.idx) && a <= ncol(alt.idx)) {
        yl <- range(c(alt.idx[, a] / mean(alt.idx[, a]),
                      yl)) * 1.1
        plot(
          x = ys,
          y = alt.idx[, a] / mean(alt.idx[, a], na.rm = TRUE),
          ylim = yl,
          col = 2,
          ylab = "Index",
          xlab = "Year",
          main = main )
      } else {
        plot(
          x = ys,
          y = rep(NA, length(ys)),
          ylim = yl,
          col = 2,
          ylab = "Index",
          xlab = "Year",
          main = main
        )
      }
      idx <- x$idx
      lo <- x$lo
      up <- x$up
      idx[x$idx <= 0] <- lo[x$idx <= 0] <- up[x$idx <= 0] <- NA
      graphics::lines(x = ys,
                      idx[, a] / mean(idx[, a], na.rm = TRUE),
                      lwd = 2)
      graphics::lines(
        x = ys,
        lo[, a] / mean(idx[, a], na.rm = TRUE),
        lwd = 2,
        lty = 2
      )
      graphics::lines(ys,
                      up[, a] / mean(idx[, a], na.rm = TRUE),
                      lwd = 2,
                      lty = 2)
      
      if (legend && !is.null(alt.idx)) {
        legend(
          legend.pos,
          pch = c(1, NA),
          lty = c(NA, 1),
          col = c(2, 1),
          legend = c("alt.idx", "GAM")
        )
      }
    }
    if (any(select == "map")) {
      xlims <- range(dat$lon, na.rm = TRUE)
      ylims <- range(dat$lat, na.rm = TRUE)
      mapvals = NULL
      if (is.null(predD)) {
        tmp <-
          dat[dat$haul.id %in% myids,] #subset(dat, haul.id %in% myids)
      } else {
        tmp <- predD
      }
      if (is.null(year)) {
        concT <- sdmgamindex::concentration_transform(log(x$gPreds[[a]]))
        mapvals <- x$gPreds[[a]]
      } else {
        y <- which(as.numeric(as.character(names(x$gPreds2[[a]]))) == year)
        if (length(y) == 0) {
          stop(paste("Year", year, "age group", a, "not found."))
        }
        concT <- sdmgamindex::concentration_transform(log(x$gPreds2[[a]][[y]]))
        mapvals <- x$gPreds2[[a]][[y]]
      }
      
      if (length(colors) > 1) {
        zFac <- cut(concT, 0:length(colors) / length(colors))
      } else {
        zFac <- 1
      }
      
      if (length(map.cex) > 1) {
        sFac <- cut(log(x$gPreds[[a]]), length(map.cex))
      } else {
        sFac <- 1
      }
      
      myCols <- colors
      plot(
        tmp$lon,
        y = tmp$lat,
        col = 1,
        pch = 1,
        cex = map.cex[sFac],
        xlim = xlims,
        ylim = ylims,
        xlab = "Longitude",
        ylab = "Latitude",
        main = main,
        ... )
      
      graphics::points(
        tmp$lon,
        y = tmp$lat,
        col = myCols[zFac],
        pch = 16,
        cex = map.cex[sFac] )
      
      maps::map(
        database = "worldHires",
        xlim = xlims,
        ylim = ylims,
        fill = TRUE,
        plot = TRUE,
        add = TRUE,
        col = grDevices::grey(0.5) )
      
      if (legend) {
        maxcuts <- stats::aggregate(mapvals ~ zFac, FUN = max)
        mincuts <- stats::aggregate(mapvals ~ zFac, FUN = min)
        mm <- mean(mapvals)
        ml <- signif(mincuts[, 2] / mm, legend.signif)
        ml[1] <- 0
        leg <- paste0("[", ml, ",", signif(maxcuts[, 2] / mm, legend.signif), "]")
        legend(
          legend.pos,
          legend = leg,
          pch = 16,
          col = colors,
          bg = "white" )
      }
    }
    if (any(select == "absolutemap") || any(select == "CVmap")) {
      if (is.null(year) || length(year) < 1) {
        stop("argument 'year' must be vector of length>=1 for type 'absolutemap'")
      }
      if (!all(year %in% levels(dat$Year))) {
        stop("invalid years selected")
      }
      xlims <- range(dat$lon, na.rm = TRUE)
      ylims <- range(dat$lat, na.rm = TRUE)
      
      if (any(select == "absolutemap")) {
        colsel <- "gPreds2"
      } else {
        colsel <- "gPreds2.CV"
      }
      
      ## collect all years as data.frame
      ally <- data.frame(val = x[[colsel]][[a]][[1]],
                         year = as.character(levels(dat$Year)[1]))
      cc <- 0
      
      for (y in levels(dat$Year)) {
        cc <- cc + 1
        ally <- rbind(ally, data.frame(val = x[[colsel]][[a]][[cc]],
                                       year = as.character(levels(dat$Year)[cc])))
      }
      
      ally$conc <- sdmgamindex::concentration_transform(log(ally$val))
      if (is.null(cutp)) {
        ally$zFac <- cut(ally$conc, 0:length(colors) / length(colors))
      } else {
        if (length(cutp) != length(colors) + 1) {
          stop("incompatible number of colors and length of cutp")
        }
        ally$zFac <- cut(ally$val, cutp)
      }
      bubbleScale <- 0.005 * max(dat$Nage[, a])
      for (yy in year) {
        if (is.null(predD)) {
          tmp = dat[dat$haul.id %in% myids,] # subset(dat, haul.id %in% myids)
        } else {
          tmp = predD
          if (is.list(tmp) && !class(tmp) %in% c("data.frame", "DATRASraw")){
            tmp = predD[[as.character(yy)]]
          }
        }
        
        plot(
          x = tmp$lon,
          y = tmp$lat,
          col = 1,
          pch = 1,
          cex = map.cex,
          xlab = "Longitude",
          ylab = "Latitude",
          axes = FALSE
        )
        graphics::box()
        title::graphics(main = yy, line = 1)
        sel = which(ally$year == yy)
        graphics::points(
          x = tmp$lon,
          y = tmp$lat,
          col = colors[as.numeric(ally$zFac[sel])],
          pch = 16,
          cex = map.cex
        )
        maps::map(
          database = 'worldHires',
          xlim = xlims,
          ylim = ylims,
          fill = TRUE,
          plot = TRUE,
          add = TRUE,
          col = grDevices::grey(0.5)
        )
        
        if (mapBubbles) {
          dy <- dat[dat$Year == y,] # subset(dat,Year==yy)
          graphics::points(
            x = dy$lon,
            y = dy$lat,
            cex = sqrt(dy$Nage[, a] / bubbleScale)
          )
        }
        
        if (legend && yy == year[1]) {
          if (is.null(cutp)) {
            maxcuts <- stats::aggregate(val ~ zFac, data = ally, FUN = max)
            mincuts <- stats::aggregate(val ~ zFac, data = ally, FUN = min)
            mm <- base::mean(ally$val)
            ml <- base::signif(mincuts[, 2] / mm, legend.signif)
            ml[1] = 0
            leg <- paste0("[", ml, ",", signif(maxcuts[, 2] / mm, legend.signif), "]")
            legend(
              legend.pos,
              legend = leg,
              pch = 16,
              col = colors,
              bg = "white" )
          }
          legend(
            legend.pos,
            legend = levels(ally$zFac),
            pch = 16,
            col = colors,
            bg = "white" )
        }
      }##rof year
    }## absolutemap
    
    for (k in 1:length(select)) {
      ss <- base::suppressWarnings(as.numeric(select[k]))
      if (!is.na(ss)) {
        mgcv::plot.gam(x$pModels[[a]], select = ss, main = main, ...)
      }
    }
    if (any(select == "residuals") || any(select == "fitVsRes") ||
        any(select == "resVsYear") || any(select == "resVsShip") ||
        any(select == "spatialResiduals")) {
      resi <- x$residuals[[a]] ##residuals(x,a)
    }
    if (any(select == "residuals")) {
      graphics::hist(resi,
                     nclass = 30,
                     main = main,
                     xlab = "Residuals",
                     ...)
    }
    if (any(select == "fitVsRes")) {
      plot(
        x = stats::fitted(x$pModels[[a]]),
        y = stats::residuals(x$pModels[[a]]),
        xlab = "Fitted",
        ylab = "Residuals",
        main = main,
        ...
      ) ## TODO: use quantile residuals here too
    }
    if (any(select == "resVsYear")) {
      plot(
        x = dat$Year,
        y = resi,
        main = main,
        xlab = "Year",
        ylab = "Residuals",
        ...
      )
    }
    if (any(select == "resVsShip")) {
      plot(
        x = dat$Ship,
        y = resi,
        main = main,
        xlab = "Year",
        ylab = "Residuals",
        ...
      )
    }
    if (any(select == "spatialResiduals")) {
      scale <- 3 * map.cex
      if (is.null(year) || length(year) > 1) {
        stop("a single year must be supplied")
      }
      sel <- which(dat[[2]]$Year == as.character(year))
      plot(
        x = dat$lon,
        y = dat$lat,
        type = "n",
        xlab = "Longitude",
        ylab = "Latitude",
        main = main,
        ...
      )
      maps::map(
        database = "worldHires",
        fill = TRUE,
        plot = TRUE,
        add = TRUE,
        col = grDevices::grey(0.5)
      )
      positive <- resi[sel] > 0
      graphics::points(
        dat$lon[sel][positive],
        dat$lat[sel][positive],
        pch = 1,
        cex = scale * sqrt(resi[sel][positive]),
        col = "blue"
      )
      graphics::points(
        dat$lon[sel][!positive],
        dat$lat[sel][!positive],
        pch = 1,
        cex = scale * sqrt(-resi[sel][!positive]),
        col = "red"
      )
    }
  }
  
  return(list())
}


# Covert and process raw data ------------------------------------------------------------------


#' Data frame with one row per haul to DATRASraw alike object
#'
#' @param x Data.frame. Zero filled catch data.
#'
#' @export
#'
#' @examples
#' dat <- data.frame(
#'           species = 1:2,
#'           Year = 2016:2020,
#'           lon = rnorm(n = 10, mean = 170, sd = 10),
#'           lat = rnorm(n = 10, mean = 170, sd = 10),
#'           sx = rnorm(n = 10, mean = 60, sd = 10),
#'           sy = rnorm(n = 10, mean = 170, sd = 10),
#'           covA = rnorm(n = 10, mean = 3, sd = 10),
#'           covB = rnorm(n = 10, mean = 20, sd = 10),
#'           EFFORT = rnorm(n = 10, mean = 20, sd = 10))
#'  ds <- split(dat,dat$species)
#'  lapply(ds, get_datrasraw)
#'
get_datrasraw <- function(x) {
  x$haul.id <- 1:nrow(x)
  dd <- list()
  dd[[1]] <- data.frame()
  dd[[2]] <- x
  dd[[3]] <- data.frame()
  class(dd) <- "DATRASraw"
  return(dd)
}

#' Prediction a grid for a specific year
#'
#' @param year Numeric string. The years the data is available for.
#' @param x Data.frame. The data.frame of the covariate data going into the model.
#' @param subsel Default = NULL.
#' @param varsbyyr Character string. The name of the variables that vary by year.
#' @param vars Character string. The name of the variables that do not vary by year.
#'
#' @return list of data by year
#' @export
#'
#' @examples
#' dat <- data.frame(
#'           lon = rnorm(n = 10, mean = 170, sd = 10),
#'           lat = rnorm(n = 10, mean = 170, sd = 10),
#'           sx = rnorm(n = 10, mean = 60, sd = 10),
#'           sy = rnorm(n = 10, mean = 170, sd = 10),
#'           covA = rnorm(n = 10, mean = 3, sd = 10),
#'           covB2019 = rnorm(n = 10, mean = 20, sd = 10),
#'           covB2020 = rnorm(n = 10, mean = 20, sd = 10),
#'           covB2021 = rnorm(n = 10, mean = 20, sd = 10),
#'           covB2022 = rnorm(n = 10, mean = 20, sd = 10),
#'           EFFORT = rnorm(n = 10, mean = 20, sd = 10))
#'
#' vars <- "covA"
#' varsbyyr <- "covB"
#' YEARS <- 2019:2022
#' pg <- lapply(YEARS,
#'        FUN = get_prediction_grid,
#'        x = dat,
#'        vars = vars,
#'        varsbyyr = varsbyyr)
#' names(pg) <- YEARS
#' pg
get_prediction_grid <- function(year,
                                x,
                                subsel = NULL,
                                varsbyyr,
                                vars) {
  nam <- c("lon", "lat", "sx", "sy", vars)
  nam2 <- tidyr::crossing(varsbyyr, year) %>%
    dplyr::mutate(vars0 = paste0(varsbyyr, year)) %>%
    dplyr::select(vars0) %>%
    unlist()
  # nam2 = paste0(c(vars),year)
  
  pd <- data.frame(x[, c(nam, nam2)])
  colnames(pd) <- c(nam,
                    unique(gsub(
                      pattern = "[0-9]+",
                      replacement = "",
                      x = nam2
                    )))
  
  pd$EFFORT <- 1.0
  ## Note, Better to use median effort here if splines on effort is used.
  ## Otherwise uncertainty is inflated because it is outside normal effort range.
  
  if (!is.null(subsel)) {
    pd = pd[subsel, ]
  }
  return(pd)
}

#' Create a grid of haul positions from a DATRASraw object.
#'
#' @title Create a grid of haul positions from a DATRASraw object.
#' @param dd DATRASraw object
#' @param nLon number of grid cells in the longitude direction.
#' @return A sdmgamindexGrid (a list of coordinates and haul.ids)
#' @export
get_grid <- function(dd,
                     nLon = 20) {
  mlon <- mean(dd$lon)
  mlat <- mean(dd$lat)
  
  kmPerDegLon <- sdmgamindex::calc_distance(
    long1 = convert_deg_rad(mlon),
    lat1 = convert_deg_rad(mlat),
    long2 = convert_deg_rad(mlon + 1),
    lat2 = convert_deg_rad(mlat)
  )
  
  kmPerDegLat <- sdmgamindex::calc_distance(
    long1 = convert_deg_rad(mlon),
    lat1 = convert_deg_rad(mlat),
    long2 = convert_deg_rad(mlon),
    lat2 = convert_deg_rad(mlat + 1)
  )
  
  get_bps <- function(labs) {
    cbind(lower = as.numeric(sub("\\((.+),.*", "\\1", labs)),
          upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", labs)))
  }
  
  lonf <- cut(dd$lon, nLon, dig.lab = 8)
  lonbps <- get_bps(levels(lonf))
  gridsize <- base::diff(lonbps[1, ]) * kmPerDegLon
  nLat <- round(diff(range(dd$lat)) * kmPerDegLat / gridsize)
  latf <- cut(dd$lat, nLat, dig.lab = 8)
  dd$StatRec2 <- as.factor(paste(lonf, latf))
  
  latbps <- get_bps(levels(latf))
  cat("Approximate grid size: ",
      mean(c(base::diff(lonbps[1, ]) * kmPerDegLon,
             diff(latbps[1, ]) * kmPerDegLat )), " km\n")
  
  uRecs <- unique(as.character(dd$StatRec2))
  N <- length(uRecs)
  mylon <- mylat <- myids <- numeric(N)
  k <- 0
  
  for (rec in uRecs) {
    k <- k + 1
    tmp <-
      dd[[2]][dd[[2]]$StatRec2 == rec,] # subset(dd[[2]],StatRec2==rec)
    mlon <- mean(tmp$lon)
    mlat <- mean(tmp$lat)
    dist <- sqrt((mlon - tmp$lon) ^ 2 + (mlat - tmp$lat) ^ 2)
    sel <- which.min(dist)
    mylon[k] <- tmp$lon[sel]
    mylat[k] <- tmp$lat[sel]
    myids[k] <- as.character(tmp$haul.id[sel])
  }
  
  ret <- list(mylon, mylat, myids, lonbps, latbps)
  class(ret) <- "sdmgamindexGrid"
  return(ret)
}


#' Change CRS of coordinates
#'
#' @param x Decimal degrees longititude.
#' @param y Decimal degrees latitude.
#' @param crs_in Default = "+proj=longlat +datum=WGS84".
#' @param crs_out Default = "EPSG:3338".
#'
#' @return data.frame with 3 columns: ID, Longitude, and Latitude
#' @export
#'
#' @examples
#' convert_crs(x = 170, y = 62)
#' dat <- sdmgamindex::noaa_afsc_public_foss[,c("longitude_dd_start", "latitude_dd_start")]
#' head(dat)
#' ll <- sdmgamindex::convert_crs( # project data
#'    x = dat$longitude_dd_start,
#'    y = dat$latitude_dd_start,
#'    crs_in = "+proj=longlat +datum=WGS84",
#'    crs_out = "EPSG:3338")
#' head(ll)
convert_crs <- function(x,
                        y,
                        crs_in = "+proj=longlat +datum=WGS84",
                        crs_out = "EPSG:3338"){ # "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs",) {
  
  xy = data.frame(ID = 1:length(x),
                  X = x,
                  Y = y)
  sp::coordinates(xy) <- c("X", "Y")
  sp::proj4string(xy) <- sp::CRS(crs_in)
  res <- sp::spTransform(xy, sp::CRS(crs_out))
  res <- as.data.frame(res)
  names(res) <- c("ID", "X", "Y")
  
  return(res)
}

#' Convert decimal degrees to radians
#' Previously called deg2rad.
#'
#' @param deg Numeric Decimal degrees of a latitude or longitude.
#'
#' @return The value of degrees in radians.
#' @export
#'
#' @examples
#' convert_deg_rad(170.1)
#' convert_deg_rad(-60.1)
convert_deg_rad <- function(deg) {
  return(deg * pi / 180)
}


#' Get distance (km) between two locations
#' previously called gcd.hf
#'
#' @param long1 Numeric Decimal degrees. Longitude of the first location.
#' @param lat1 Numeric Decimal degrees. Latitude of the first location.
#' @param long2 Numeric Decimal degrees. Longitude of the second location.
#' @param lat2 Numeric Decimal degrees. Latitude of the second location.
#'
#' @return Distance between locations in km.
#' @export
#'
#' @examples
#' calc_distance(long1 = 170, lat1 = 62, long2 = 170.1, lat2 = 62.1)
#' calc_distance(long1 = -170, lat1 = -62, long2 = -170.1, lat2 = -62.1)
#' calc_distance(long1 = 170, lat1 = -62, long2 = 170.1, lat2 = -62.1)
#' calc_distance(long1 = 50, lat1 = 50, long2 = 50.1, lat2 = 50.1)
calc_distance <-
  function(long1, lat1, long2, lat2) {
    R <- 6371 # Earth mean radius [km]
    delta.long <- (long2 - long1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat / 2) ^ 2 + cos(lat1) * cos(lat2) * sin(delta.long / 2)^2
    c <- 2 * asin(min(1, sqrt(a)))
    d <- R * c
    
    return(d) # Distance in km
  }

# Mapping ----------------------------------------------------------------------

#' Get bathymetric prediction grid corresponding to the area for a DATRASraw object using the 'marmap' package
#'
#' @title Get bathymetric prediction grid corresponding to the area for a DATRASraw object using the 'marmap' package
#' @param d DATRASraw object
#' @param minDepth Minimum depth to include
#' @param maxDepth Maximum depth to include
#' @param resolution grid resolution (see marmap::getNOAA.bathy)
#' @param maxDist Do not include grid points farther than maxDist from nearest observation.
#' @param keep Save grid on disk for fast loading next time?
#' @param shapefile extra shapefile information to add (optional)
#' @param select columns to extract from shapefile
#' @return data.frame with depths and geographical coordinates
#' @export
get_bathy_grid <- function(d,
                           minDepth = 10,
                           maxDepth = Inf,
                           resolution = 2,
                           maxDist = Inf,
                           keep = TRUE,
                           shapefile = NULL,
                           select = NULL) {
  if (!requireNamespace("marmap", quietly = TRUE)) {
    stop("Package \"marmap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  bathy <-
    marmap::getNOAA.bathy(
      lon1 = min(d$lon),
      lon2 = max(d$lon),
      lat1 = min(d$lat),
      lat2 = max(d$lat),
      resolution = resolution,
      keep = keep
    )
  xyzbathy <- as.data.frame(marmap::as.xyz(bathy))
  colnames(xyzbathy) <- c("lon", "lat", "Depth")
  xyzbathy$Depth <- -xyzbathy$Depth
  xyzbathy <- xyzbathy[xyzbathy$Depth > minDepth & xyzbathy$Depth < maxDepth, ]#subset(xyzbathy, Depth>minDepth & Depth<maxDepth)
  if (is.finite(maxDist)) {
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop(
        "Package \"RANN\" needed for this function to work when 'maxDist' is used. Please install it.",
        call. = FALSE
      )
    }
    nearest <- RANN::nn2(d[[2]][, c("lon", "lat")], xyzbathy[, 1:2], k =
                           1)
    xyzbathy <- xyzbathy[nearest$nn.dists < maxDist , ]
  }
  xyzbathy <- as.data.frame(xyzbathy)
  if (is.character(shapefile)) {
    if (file.exists(shapefile)) {
      shape <- maptools::readShapeSpatial(shapefile)
    } else {
      stop("shapefile not found")
    }
    tmp <- xyzbathy
    sp::coordinates(tmp) <- ~ lon + lat
    xtra <- sp::over(tmp, shape)
    if (!is.null(select))
      xtra <- xtra[select]
    xyzbathy <- cbind(xyzbathy, xtra)
  }
  
  return(xyzbathy)
}

# Data -------------------------------------------------------------------------

#' @title EBS Prediction Grid
#' @description EBS Prediction Grid
#' @usage data("pred_grid_ebs")
#' @author Emily Markowitz (emily.markowitz AT noaa.gov) and James Thorson (james.thorson AT noaa.gov)
#' @format A data frame with 36690 observations on the following 3 variables.
#' \describe{
#'   \item{\code{lon}}{Numeric; Longitude (one hundred thousandth of a decimal degree).}
#'   \item{\code{lat}}{Numeric; Latitude (one hundred thousandth of a decimal degree).}
#'   \item{\code{Shape_Area}}{Numeric; The area that this location represents. }
#'   }
#' @source https://github.com/James-Thorson-NOAA/VAST
#' @keywords eastern bering sea grid prediciton
#' @examples
#' data(pred_grid_ebs)
#' @details DETAILS
"pred_grid_ebs"

