
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
#' @param gamma model degress of freedom inflation factor (see 'gamma' argument to gam() )
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
#' @param linkZ link function for the binomial part of the model, default: "logit" (not used for Tweedie models).
#' @param CIlevel Confidence interval level, defaults to 0.95.
#' @param ... Optional extra arguments to "gam"
#' @return A survey index (list)
#' @examples
#' \dontrun{
#' library(surveyIndex)
#' ##downloadExchange("NS-IBTS",1994:2014)
#' dAll<-readExchangeDir(".",strict=FALSE)
#' mc.cores<-2; library(parallel)
#' d<-subset(dAll, Species=="Pollachius virens",Quarter==1,HaulVal=="V",StdSpecRecCode==1, Gear=="GOV")
#' dAll<-NULL; gc(); ## lose dAll because it takes up a lot of memory
#' d<-addSpectrum(d,by=1)
#' ## get idea about number of age groups to include
#' agetab<-xtabs(NoAtALK~Year+Age,data=d[[1]])
#' agetab.df<-as.data.frame(agetab)
#' ages<-1:8
#' ## require at least 1 aged individual in each year
#' for(a in ages){
#'     if(any(agetab.df$Freq[agetab.df$Age==a]<1))
#'         d<-fix_age_group(d,age=a,fun=ifelse(a==min(ages),"min","mean"))
#' }
#' d<-subset(d,Age>=min(ages))
#'
#' ###############################
#' ## Convert to numbers-at-age
#' ###############################
#' d.ysplit <- split(d, d$Year)
#' ALK<-mclapply(d.ysplit,fitALK,minAge=min(ages),maxAge=max(ages),autoChooseK=TRUE,useBIC=TRUE,
#'                varCof=FALSE,maxK=50,mc.cores=mc.cores)
#' Nage<-mclapply(ALK,predict,mc.cores=mc.cores)
#' for(i in 1:length(ALK)) d.ysplit[[i]]$Nage=Nage[[i]];
#' dd <- do.call("c",d.ysplit)
#'
#' ## Fit model
#'
#' grid <- get_grid(dd, nLon=40)
#' ## set max basis dim for spatial smooths by age, P=positive and Z=zero/absence.
#' ## These are set relatively low here to speed up the example
#' kvP <- c(50,50,50,40,30,rep(10,length(ages)-5))
#' kvZ <- kvP / 2;
#' mP <- rep("Year+s(lon,lat,k=kvecP[a],bs='ts')+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(ages)  );
#' mZ <- rep("Year+s(lon,lat,k=kvecZ[a],bs='ts')+s(Depth,bs='ts',k=6)+offset(log(HaulDur))",length(ages)  );
#'
#' SIQ1 <- get_surveyidx(dd,ages=ages,myids=grid[[3]],cutOff=0.1,kvecP=kvP,kvecZ=kvZ,
#'          modelZ=mZ,modelP=mP,mc.cores=mc.cores) ## if errors are encountered, debug with mc.cores=1
#'
#' strat.mean<-get_surveyidxStratMean(dd,ages)
#'
#' ## plot indices, distribution map, and estimated depth effects
#' plot_surveyidx(SIQ1,dd,cols=ages,alt.idx=strat.mean,grid[[3]],par=list(mfrow=c(3,3)),legend=FALSE,
#'                select="index",plotByAge=FALSE)
#'
#' plot_surveyidx(SIQ1,dd,cols=ages,alt.idx=NULL,grid[[3]],par=list(mfrow=c(3,3)),legend=FALSE,
#'                 colors=rev(heat.colors(8)),select="map",plotByAge=FALSE)
#'
#' plot_surveyidx(SIQ1,dd,cols=ages,alt.idx=NULL,grid[[3]],par=list(mfrow=c(3,3)),
#'                 legend=FALSE,select="2",plotByAge=FALSE)
#'
#' ## Calculate internal concistency and export to file
#' internalCons(SIQ1$idx)
#' exportSI(SIQ1$idx,ages=ages,years=levels(dd$Year),toy=mean(dd$timeOfYear),file="out.dat",
#'          nam="Survey index demo example")
#' }
#' @importFrom MASS mvrnorm
#' @export
get_surveyidx <- function(
    x,
    ages,
    myids,
    kvecP=rep(12*12,length(ages)),
    kvecZ=rep(8*8,length(ages)),
    gamma=1.4,
    cutOff=1,
    fam="Gamma",
    useBIC=FALSE,
    nBoot=1000,
    mc.cores=1,
    method="ML",
    predD=NULL,
    modelZ=rep("Year+s(lon,lat,k=kvecZ[a],bs='ts')+s(Ship,bs='re',by=dum)+s(Depth,bs='ts')+s(TimeShotHour,bs='cc')",length(ages)  ),
    modelP=rep("Year+s(lon,lat,k=kvecP[a],bs='ts')+s(Ship,bs='re',by=dum)+s(Depth,bs='ts')+s(TimeShotHour,bs='cc')",length(ages)  ),
    knotsP=NULL,
    knotsZ=NULL,
    predfix=NULL,
    linkZ="logit",
    CIlevel=0.95,
    ... ){

  if(is.null(x$Nage)) stop("No age matrix 'Nage' found.");
  if(is.null(colnames(x$Nage))) stop("No colnames found on 'Nage' matrix.");
  if(length(modelP)<length(ages)) stop(" length(modelP) < length(ages)");
  if(length(kvecP)<length(ages)) stop(" length(kvecP) < length(ages)");
  if(fam[1]!="Tweedie"){
    if(length(modelZ)<length(ages)) stop(" length(modelZ) < length(ages)");
    if(length(kvecZ)<length(ages)) stop(" length(kvecZ) < length(ages)");
  }
  ## check for valid family names
  stopifnot(fam[1] %in% c("Gamma","LogNormal","Tweedie","negbin"))
  if(length(fam)<length(ages)) { famVec = rep(fam[1],length(ages)) } else famVec=fam;

  dataAges <- as.numeric(gsub("[+]","",colnames(x$Nage)))
  if(!all(ages%in%dataAges)) stop(paste0("age(s) ",setdiff(ages,dataAges)," not found in 'Nage' matrix"));
  x[[1]]$Year=as.factor(x[[1]]$Year);
  x[[2]]$Year=as.factor(x[[2]]$Year);
  pModels=list()
  zModels=list()
  gPreds=list() ##last data year's predictions
  gPreds2=list() ## all years predictions
  gPreds2.CV=list() ## coefficient of variation of previous
  allobs=list() ## response vector (zeroes and positive)
  resid=list() ## residuals
  predDc = predD ##copy of predD

  if (exists(".Random.seed")) {
    oldseed <- get(".Random.seed", .GlobalEnv)
    oldRNGkind <- RNGkind()
    on.exit( { do.call("RNGkind", as.list(oldRNGkind)); assign(".Random.seed", oldseed, .GlobalEnv) }  )
  }
  set.seed(314159265)

  yearNum=as.numeric(as.character(x$Year));
  yearRange=min(yearNum):max(yearNum);

  ## Choose most frequent gear as reference gear
  gearNames=names(xtabs(~Gear,data=x[[2]]))
  myGear=names(xtabs(~Gear,data=x[[2]]))[which.max(xtabs(~Gear,data=x[[2]]))]
  if(!is.null(predfix$Gear)) {
    myGear = predfix$Gear
  }

  resMat=matrix(NA,nrow=length(yearRange),ncol=length(ages));
  upMat=resMat;
  loMat=resMat;
  do.one.a<-function(a){
    age = which(dataAges==ages[a])
    ddd=x[[2]]; ddd$dum=1.0;
    ddd$A1=ddd$Nage[,age]
    gammaPos=gamma;
    gammaZ=gamma;

    if(useBIC){
      nZ=nrow(ddd);
      nPos=nrow(subset(ddd,A1>cutOff));
      gammaPos=log(nPos)/2;
      gammaZ=log(nZ)/2;
      cat("gammaPos: ",gammaPos," gammaZ: ",gammaZ,"\n");
    }
    pd = subset(ddd,A1>cutOff)
    if(famVec[a]=="LogNormal"){
      f.pos = as.formula( paste( "log(A1) ~",modelP[a]));
      f.0 = as.formula( paste( "A1>",cutOff," ~",modelZ[a]));

      print(system.time(m.pos<-DATRAS:::tryCatch.W.E(gam(f.pos,data=subset(ddd,A1>cutOff),gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));

      if(class(m.pos)[2] == "error") {
        print(m.pos)
        stop("Error occured for age ", a, " in the positive part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }

      print(system.time(m0<-DATRAS:::tryCatch.W.E(gam(f.0,gamma=gammaZ,data=ddd,family=binomial(link=linkZ),method=method,knots=knotsZ,na.action=na.fail,...))$value));

      if(class(m0)[2] == "error") {
        print(m0)
        stop("Error occured for age ", a, " in the binomial part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }

    } else if(famVec[a]=="Gamma"){
      f.pos = as.formula( paste( "A1 ~",modelP[a]));
      f.0 = as.formula( paste( "A1>",cutOff," ~",modelZ[a]));

      print(system.time(m.pos<-DATRAS:::tryCatch.W.E(gam(f.pos,data=subset(ddd,A1>cutOff),family=Gamma(link="log"),gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));

      if(class(m.pos)[2] == "error") {
        print(m.pos)
        stop("Error occured for age ", a, " in the positive part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }

      print(system.time(m0<-DATRAS:::tryCatch.W.E(gam(f.0,gamma=gammaZ,data=ddd,family=binomial(link=linkZ),method=method,knots=knotsZ,na.action=na.fail,...))$value));

      if(class(m0)[2] == "error") {
        print(m0)
        stop("Error occured for age ", a, " in the binomial part of the model\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }
    } else if(famVec[a]=="Tweedie"){
      ddd$A1[ ddd$A1<cutOff ] = 0
      pd = ddd
      f.pos = as.formula( paste( "A1 ~",modelP[a]));
      print(system.time(m.pos<-DATRAS:::tryCatch.W.E(gam(f.pos,data=ddd,family=tw,gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));
      if(class(m.pos)[2] == "error") {
        print(m.pos)
        stop("Error occured for age ", a, ".\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }
      m0=NULL;
    } else if(famVec[a]=="negbin"){
      pd = ddd
      f.pos = as.formula( paste( "A1 ~",modelP[a]));
      print(system.time(m.pos<-DATRAS:::tryCatch.W.E(gam(f.pos,data=ddd,family=nb,gamma=gammaPos,method=method,knots=knotsP,na.action=na.fail,...))$value));
      if(class(m.pos)[2] == "error") {
        print(m.pos)
        stop("Error occured for age ", a, ".\n", "Try reducing the number of age groups or decrease the basis dimension of the smooths, k\n")
      }
      m0=NULL;
    }
    ## Calculate total log-likelihood
    if(famVec[a]=="Tweedie" || famVec[a]=="negbin"){
      totll = logLik(m.pos)[1]
    } else {
      p0p =(1-predict(m0,type="response"))
      ppos=p0p[ddd$A1>cutOff]
      p0m1=p0p[ddd$A1<=cutOff]
      if(famVec[a]=="Gamma")  totll=sum(log(p0m1))+sum(log(1-ppos))+logLik(m.pos)[1];
      ## if logNormal model, we must transform til log-likelihood to be able to use AIC
      ## L(y) = prod( dnorm( log y_i, mu_i, sigma^2) * ( 1 / y_i ) ) => logLik(y) = sum( log[dnorm(log y_i, mu_i, sigma^2)]  - log( y_i ) )
      if(famVec[a]=="LogNormal") totll=sum(log(p0m1))+ sum(log(1-ppos)) + logLik(m.pos)[1] - sum(m.pos$y);
    }

    if(is.null(predD)) predD=subset(ddd,haul.id %in% myids);
    res=numeric(length(yearRange));
    lores=res;
    upres=res;
    gp2=list()
    gp2.cv=list()

    for(y in levels(ddd$Year)){
      ## take care of years with all zeroes
      if(!any(ddd$A1[ddd$Year==y]>cutOff)){
        res[which(as.character(yearRange)==y)]=0;
        upres[which(as.character(yearRange)==y)] = 0;
        lores[which(as.character(yearRange)==y)] = 0;
        next;
      }
      if(is.list(predDc) && !class(predDc)%in%c("data.frame","DATRASraw")) predD = predDc[[as.character(y)]]
      if(is.null(predD)) stop(paste("Year",y," not found in predD"))
      ## OBS: effects that should be removed should be included here
      predD$Year=y; predD$dum=0;
      predD$ctime=as.numeric(as.character(y));
      predD$TimeShotHour=mean(ddd$TimeShotHour)
      predD$Ship=names(which.max(summary(ddd$Ship)))
      predD$timeOfYear=mean(ddd$timeOfYear);
      predD$HaulDur=30.0
      predD$Gear=myGear;
      if(!is.null(predfix)){ ##optional extra variables for standardization
        stopifnot(is.list(predfix))
        for(n in names(predfix)){
          predD[,n] = predfix[[n]]
        }
      }
      p.1 <- p.0 <- NULL
      try({
        Xp.1=predict(m.pos,newdata=predD,type="lpmatrix");
        OS.pos = numeric(nrow(predD));
        terms.pos=terms(m.pos)
        if(!is.null(m.pos$offset)){
          off.num.pos <- attr(terms.pos, "offset")
          for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos,
                                                              "variables")[[i + 1]], predD)
        }
        p.1 =Xp.1%*%coef(m.pos)+OS.pos;
        if(!famVec[a] %in% c("Tweedie","negbin")){
          Xp.0=predict(m0,newdata=predD,type="lpmatrix");
          brp.0=coef(m0);
          OS0 = numeric(nrow(predD))
          terms.0=terms(m0)
          if(!is.null(m0$offset)){
            off.num.0 <- attr(terms.0, "offset")
            for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0,
                                                        "variables")[[i + 1]], predD)
          }
          p.0 = m0$family$linkinv(Xp.0%*%brp.0+OS0);
        }
      });
      ## take care of failing predictions
      if(!is.numeric(p.1) | (!famVec[a] %in% c("Tweedie","negbin") && !is.numeric(p.0))) {
        res[which(as.character(yearRange)==y)]=0;
        upres[which(as.character(yearRange)==y)] = 0;
        lores[which(as.character(yearRange)==y)] = 0;
        next;
      }
      sig2=m.pos$sig2;

      if(famVec[a]=="Gamma") { res[which(as.character(yearRange)==y)] = sum(p.0*exp(p.1)); gPred=p.0*exp(p.1) }
      if(famVec[a]=="LogNormal")  { res[which(as.character(yearRange)==y)] = sum(p.0*exp(p.1+sig2/2)); gPred=p.0*exp(p.1+sig2/2) }
      if(famVec[a] %in% c("Tweedie","negbin"))  { res[which(as.character(yearRange)==y)] = sum(exp(p.1)); gPred=exp(p.1) }
      gp2[[y]]=gPred;
      if(nBoot>10){
        brp.1=mvrnorm(n=nBoot,coef(m.pos),m.pos$Vp);
        if(!famVec[a] %in% c("Tweedie","negbin")){
          brp.0=mvrnorm(n=nBoot,coef(m0),m0$Vp);
          OS0 = matrix(0,nrow(predD),nBoot);
          terms.0=terms(m0)
          if(!is.null(m0$offset)){
            off.num.0 <- attr(terms.0, "offset")
            for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0,
                                                        "variables")[[i + 1]], predD)
          }
          rep0=m0$family$linkinv(Xp.0%*%t(brp.0)+OS0);
        }
        OS.pos = matrix(0,nrow(predD),nBoot);
        terms.pos=terms(m.pos)
        if(!is.null(m.pos$offset)){
          off.num.pos <- attr(terms.pos, "offset")
          for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos,
                                                              "variables")[[i + 1]], predD)
        }
        if(famVec[a]=="LogNormal"){
          rep1=exp(Xp.1%*%t(brp.1)+sig2/2+OS.pos);
        } else {
          rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
        }

        if(!famVec[a] %in% c("Tweedie","negbin")){
          idxSamp = colSums(rep0*rep1);
          gp.sd = apply(rep0*rep1,1,sd)
          gp2.cv[[y]] = gp.sd / gPred
        } else {
          idxSamp = colSums(rep1);
          gp.sd = apply(rep1,1,sd)
          gp2.cv[[y]] = gp.sd / gPred
        }
        halpha = (1-CIlevel)/2
        upres[which(as.character(yearRange)==y)] = quantile(idxSamp,1-halpha,na.rm=TRUE);
        lores[which(as.character(yearRange)==y)] = quantile(idxSamp,halpha,na.rm=TRUE);
      }
    } ## rof years
    list(res=res,m.pos=m.pos,m0=m0,lo=lores,up=upres,gp=gPred,ll=totll,pd=pd,gp2=gp2,gp2.cv=gp2.cv);
  }## end do.one
  noAges=length(ages);
  rr=parallel::mclapply(1:noAges,do.one.a,mc.cores=mc.cores);
  logl=0;
  for(a in 1:noAges){
    resMat[,a]=rr[[a]]$res;
    zModels[[a]]=rr[[a]]$m0;
    pModels[[a]]=rr[[a]]$m.pos;
    loMat[,a]=rr[[a]]$lo;
    upMat[,a]=rr[[a]]$up;
    gPreds[[a]]=rr[[a]]$gp;
    logl=logl+rr[[a]]$ll
    gPreds2[[a]]=rr[[a]]$gp2
    gPreds2.CV[[a]]=rr[[a]]$gp2.cv
    allobs[[a]]=x[[2]]$Nage[,a]
  }
  getEdf<-function(m) sum(m$edf)
  totEdf=sum( unlist( lapply(zModels,getEdf))) + sum( unlist( lapply(pModels,getEdf)));
  rownames(resMat)<-yearRange
  colnames(resMat)<-ages
  out <- list(idx=resMat,zModels=zModels,pModels=pModels,lo=loMat,up=upMat,gPreds=gPreds,logLik=logl,edfs=totEdf,gPreds2=gPreds2,gPreds2.CV = gPreds2.CV, family=famVec, cutOff=cutOff, dataAges=dataAges, yearNum=yearNum, refGear=myGear, predfix = predfix, knotsP=knotsP, knotsZ=knotsZ, allobs=allobs,CIlevel=CIlevel);
  class(out) <- "surveyIdx"
  set.seed(314159265) ## reset seed here (in case multicore option is used)
  for(a in 1:noAges) resid[[a]] = residuals(out,a)
  out$residuals = resid
  out
}




#' Calculate confidence intervals for a named parameter in a survey index model.
#'
#' @title Calculate confidence intervals for a named parameter in a survey index model.
#' @param x survey index
#' @param dat DATRASraw object
#' @param parName name of the parameter, e.g. "Gear"
#' @param cutOff see getSurveyIndex()
#' @param nboot see getSurveyIndex()
#' @param pOnly only calculate for positive part of model, defaults to FALSE.
#' @return list of estimates + ci bounds for each age group.
#' @importFrom MASS mvrnorm
#' @export
get_effect <-
  function(x,dat,parName="Gear",cutOff,nboot=1000,pOnly=FALSE){
    noAges=length(x$pModels);

    res=list()
    for(a in 1:noAges){
      cat("Age ",a,"\n");
      shipSelP = grep(parName,names(coef(x$pModels[[a]])))
      shipSelZ = grep(parName,names(coef(x$zModels[[a]])))

      dd=subset(dat,Nage[,a]>cutOff)
      zNam=levels(dat$Ship)
      pNam=levels(dd$Ship)
      dif=setdiff(zNam,pNam);
      remo = which(zNam %in% dif)
      if(length(remo)>0) shipSelZ = shipSelZ[-remo];

      if(length(shipSelP)!=length(shipSelZ)) { print("unequal number of ship effects"); }


      brp.1=mvrnorm(n=nboot,coef(x$pModels[[a]]),x$pModels[[a]]$Vp);
      brp.0=mvrnorm(n=nboot,coef(x$zModels[[a]]),x$zModels[[a]]$Vp);

      shipE = exp(brp.1[,shipSelP,drop=FALSE]);
      if(!pOnly) shipE = shipE*x$zModels[[a]]$family$linkinv(brp.0[,shipSelZ,drop=FALSE])

      upres  = apply(shipE,2, quantile,probs=0.975);
      lores  = apply(shipE,2, quantile,probs=0.025);

      tmp=cbind(colMeans(shipE),upres,lores);

      res[[a]]=tmp;
    }
    return(res);
  }


#' Create a grid of haul positions from a DATRASraw object.
#'
#' @title Create a grid of haul positions from a DATRASraw object.
#' @param dd DATRASraw object
#' @param nLon number of grid cells in the longitude direction.
#' @return A surveyIndexGrid (a list of coordinates and haul.ids)
#' @export
get_grid <-
  function(dd,nLon=20){
    mlon=mean(dd$lon)
    mlat=mean(dd$lat)
    kmPerDegLon=calc_distance(convert_deg_rad(mlon),convert_deg_rad(mlat),convert_deg_rad(mlon+1),convert_deg_rad(mlat))
    kmPerDegLat=calc_distance(convert_deg_rad(mlon),convert_deg_rad(mlat),convert_deg_rad(mlon),convert_deg_rad(mlat+1))

    get_bps <- function(labs){
      cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
            upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    }

    lonf = cut(dd$lon, nLon,dig.lab=8)
    lonbps = get_bps(levels(lonf))
    gridsize = diff(lonbps[1,])*kmPerDegLon
    nLat = round(diff(range(dd$lat))*kmPerDegLat/gridsize)
    latf = cut(dd$lat, nLat, dig.lab=8)
    dd$StatRec2 = as.factor(paste(lonf,latf))

    latbps=get_bps(levels(latf))
    cat("Approximate grid size: ", mean( c(diff(lonbps[1,])*kmPerDegLon,diff(latbps[1,])*kmPerDegLat))," km\n");

    uRecs=unique(as.character(dd$StatRec2))
    N=length(uRecs);
    mylon=numeric(N);
    mylat=numeric(N);
    myids=numeric(N);
    k=0;
    for(rec in uRecs)
    {
      k=k+1;
      tmp=subset(dd[[2]],StatRec2==rec);
      mlon=mean(tmp$lon);
      mlat=mean(tmp$lat);
      dist=sqrt( (mlon-tmp$lon)^2+(mlat-tmp$lat)^2);
      sel=which.min(dist);
      mylon[k]=tmp$lon[sel];
      mylat[k]=tmp$lat[sel];
      myids[k]=as.character(tmp$haul.id[sel]);
    }
    ret <- list(mylon,mylat,myids,lonbps,latbps)
    class(ret) <- "surveyIndexGrid"
    return(ret);
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
redo_surveyidx <- function(
    x,
    model,
    predD=NULL,
    myids,
    nBoot=1000,
    predfix=list(),
    mc.cores=1){

  ages = as.numeric(colnames(model$idx))
  dataAges <- model$dataAges
  famVec = model$family
  cutOff = model$cutOff

  yearNum=model$yearNum
  yearRange=min(yearNum):max(yearNum);

  gPreds=list() ##last data year's predictions
  gPreds2=list() ## all years predictions
  gPreds2.CV=list()
  predDc = predD

  myGear=model$refGear

  resMat=matrix(NA,nrow=length(yearRange),ncol=length(ages));
  upMat=resMat;
  loMat=resMat;
  do.one.a<-function(a){
    age = which(dataAges==ages[a])
    ddd=x[[2]]; ddd$dum=1.0;
    ddd$A1=ddd$Nage[,age]
    m.pos = model$pModels[[a]]
    m0 = NULL
    if(!famVec[a] %in% c("Tweedie","negbin")) m0 = model$zModels[[a]]
    if(is.null(predD)) predD=subset(ddd,haul.id %in% myids);
    res=numeric(length(yearRange));
    lores=res;
    upres=res;
    gp2=list()
    gp2.cv=list()
    do.one.y<-function(y){

      cat("Doing year ",y,"\n")

      if(is.list(predDc) && !class(predDc)%in%c("data.frame","DATRASraw")) predD = predDc[[as.character(y)]]
      if(is.null(predD)) stop(paste("Year",y," not found in predD"))

      ## take care of years with all zeroes
      if(!any(ddd$A1[ddd$Year==y]>cutOff)){
        return(list(res=0,upres=0,lores=0,gp2=NULL,gp2.cv=NULL))
      }

      ## OBS: effects that should be removed should be included here
      predD$Year=y; predD$dum=0;
      predD$ctime=as.numeric(as.character(y));
      predD$TimeShotHour=mean(ddd$TimeShotHour)
      predD$Ship=names(which.max(summary(ddd$Ship)))
      predD$timeOfYear=mean(ddd$timeOfYear);
      predD$HaulDur=30.0
      predD$Gear=myGear;
      if(!is.null(predfix)){ ##optional extra variables for standardization
        stopifnot(is.list(predfix))
        for(n in names(predfix)){
          predD[,n] = predfix[[n]]
        }
      }
      p.1 <- p.0 <- NULL
      try({
        Xp.1=predict(m.pos,newdata=predD,type="lpmatrix");
        OS.pos = numeric(nrow(predD));
        terms.pos=terms(m.pos)
        if(!is.null(m.pos$offset)){
          off.num.pos <- attr(terms.pos, "offset")
          for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos,
                                                              "variables")[[i + 1]], predD)
        }
        p.1 =Xp.1%*%coef(m.pos)+OS.pos;
        if(!famVec[a] %in% c("Tweedie","negbin")){
          Xp.0=predict(m0,newdata=predD,type="lpmatrix");
          brp.0=coef(m0);
          OS0 = numeric(nrow(predD))
          terms.0=terms(m0)
          if(!is.null(m0$offset)){
            off.num.0 <- attr(terms.0, "offset")
            for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0,
                                                        "variables")[[i + 1]], predD)
          }
          p.0 = m0$family$linkinv(Xp.0%*%brp.0+OS0);
        }
      });
      ## take care of failing predictions
      if(!is.numeric(p.1) | (!famVec[a] %in% c("Tweedie","negbin") && !is.numeric(p.0))) {
        return(list(res=0,upres=0,lores=0,gp2=NULL,gp2.cv=NULL))
      }
      sig2=m.pos$sig2;
      idx = NA
      if(famVec[a]=="Gamma") { idx <- sum(p.0*exp(p.1)); gPred=p.0*exp(p.1) }
      if(famVec[a]=="LogNormal")  { idx <- sum(p.0*exp(p.1+sig2/2)); gPred=p.0*exp(p.1+sig2/2) }
      if(famVec[a] %in% c("Tweedie","negbin"))  { idx <- sum(exp(p.1)); gPred=exp(p.1) }

      if(nBoot>10){
        brp.1=mvrnorm(n=nBoot,coef(m.pos),m.pos$Vp);
        if(!famVec[a] %in% c("Tweedie","negbin")){
          brp.0=mvrnorm(n=nBoot,coef(m0),m0$Vp);
          OS0 = matrix(0,nrow(predD),nBoot);
          terms.0=terms(m0)
          if(!is.null(m0$offset)){
            off.num.0 <- attr(terms.0, "offset")

            for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0,
                                                        "variables")[[i + 1]], predD)
          }
          rep0=m0$family$linkinv(Xp.0%*%t(brp.0)+OS0);
        }
        OS.pos = matrix(0,nrow(predD),nBoot);
        terms.pos=terms(m.pos)
        if(!is.null(m.pos$offset)){
          off.num.pos <- attr(terms.pos, "offset")

          for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos,
                                                              "variables")[[i + 1]], predD)
        }


        if(famVec[a]=="LogNormal"){
          rep1=exp(Xp.1%*%t(brp.1)+sig2/2+OS.pos);
        } else {
          rep1=exp(Xp.1%*%t(brp.1)+OS.pos);
        }

        if(!famVec[a] %in% c("Tweedie","negbin")){
          idxSamp = colSums(rep0*rep1);
          gp.sd = apply(rep0*rep1,1,sd)
          gp.cv = gp.sd / gPred
        } else {
          idxSamp = colSums(rep1);
          gp.sd = apply(rep1,1,sd)
          gp.cv = gp.sd / gPred
        }
        halpha = (1-model$CIlevel)/2
        return(list(res=idx,upres=quantile(idxSamp,1-halpha,na.rm=TRUE),lores=quantile(idxSamp,halpha,na.rm=TRUE),gp2=gPred,gp2.cv=gp.cv))
      }
    } ## rof years
    yres = parallel::mclapply(levels(ddd$Year),do.one.y,mc.cores=mc.cores)
    for(y in levels(ddd$Year)) {
      ii = which(as.character(yearRange)==y)
      res[ii] = yres[[ii]]$res
      upres[ii] = yres[[ii]]$upres
      lores[ii] = yres[[ii]]$lores
      gp2[[y]] = yres[[ii]]$gp2
      gp2.cv[[y]] = yres[[ii]]$gp2.cv
    }
    list(res=res,m.pos=m.pos,m0=m0,lo=lores,up=upres,gp=tail(gp2,1),gp2=gp2,gp2.cv=gp2.cv);
  }## end do.one
  noAges=length(ages);
  rr=lapply(1:noAges,do.one.a);
  logl=0;
  for(a in 1:noAges){
    resMat[,a]=rr[[a]]$res;
    loMat[,a]=rr[[a]]$lo;
    upMat[,a]=rr[[a]]$up;
    gPreds[[a]]=rr[[a]]$gp;
    gPreds2[[a]]=rr[[a]]$gp2
    gPreds2.CV[[a]]=rr[[a]]$gp2.cv
  }
  rownames(resMat)<-yearRange
  colnames(resMat)<-ages


  out <- list(idx=resMat,
              zModels=model$zModels,
              pModels=model$pModels,
              lo=loMat,up=upMat,
              gPreds=gPreds,
              logLik=model$logLik,
              edfs=model$edfs,
              pData=model$pData,
              gPreds2=gPreds2,
              gPreds2.CV=gPreds2.CV,
              family=famVec,
              cutOff=cutOff,
              dataAges=dataAges,
              yearNum=yearNum,
              refGear=myGear,
              predfix = predfix,
              knotsP=model$knotsP,
              knotsZ=model$knotsZ);
  class(out) <- "surveyIdx"

  return(out)
}


#' Survey index using the stratified mean method using ICES statistical rectangles as strata.
#'
#' @title Survey index using the stratified mean method using ICES statistical rectangles as strata.
#' @param x DATRASraw object. Must contain a matrix: x[[2]]$Nage.
#' @param ageCols which columns of the Nage matrix should be included?
#' @param doLog log-transform?
#' @return a matrix with survey indices
#' @export
get_surveyidxStratMean <-
  function(x,ageCols,doLog=FALSE){
    ysplit=split(x,x$Year);
    res=matrix(NA,nrow=length(ysplit),ncol=length(ageCols))
    for(y in 1:length(ysplit)){

      if(!doLog) {
        byRec=aggregate(ysplit[[y]][[2]]$Nage[,ageCols],by=list(ysplit[[y]][[2]]$StatRec),FUN="mean")  } else {
          byRec=aggregate(log(ysplit[[y]][[2]]$Nage[,ageCols]+1),by=list(ysplit[[y]][[2]]$StatRec),FUN="mean")
        }

      res[y,]=colMeans(byRec[,-1])
    }
    res
  }



# Assess Index -----------------------------------------------------------------

#' Calculate internal consistency of a survey index.
#'
#' @title Calculate internal consistency of a survey index.
#' @param tt A matrix with survey indices (rows=years, cols=ages)
#' @param do.plot Plot it?
#' @return a vector of consistencies
#' @export
internalCons <-
  function(tt,do.plot=FALSE){
    tt[tt==0]=NA
    sel1=1:(nrow(tt)-1);
    sel2=2:nrow(tt);
    if(do.plot){ X11(); b=ceiling((ncol(tt)-1)/2); par(mfrow=c(b,b));}
    for(a in 1:(ncol(tt)-1)){
      cat("Age ",a," vs ",a+1," : ",cor(log(tt[sel1,a]),log(tt[sel2,a+1]),use="pairwise.complete.obs"),"\n")
      if(do.plot) {plot(log(tt[sel1,a]),log(tt[sel2,a+1])); abline(0,1);}
    }
    return( sapply(1:(ncol(tt)-1), function(a) cor(log(tt[sel1,a]),log(tt[sel2,a+1]),use="pairwise.complete.obs")));
  }


# Simulate ---------------------------------------------------------------------


#' Simulate data from a surveyIdx model (experimental and subject to change)
#'
#' @title Simulate data from a surveyIndex model (experimental and subject to change)
#' @param model object of class 'surveyIdx'
#' @param d A dataset (DATRASraw object)
#' @param sampleFit Use a random sample from the gaussian approximation to distribution of the estimated parameter vector. Default: FALSE.
#' @param condSim optional results of previous call to this function. Use this if you want to generate many datasets (much faster, since mean predictions are re-used).
#' @return list with  1) simulated observations with noise 2) mean (no noise) 3) zero probability.
#' @export
surveyIdx.simulate<-function(model,d,sampleFit=FALSE,condSim=NULL){
  ages = as.numeric(colnames(model$idx))
  dataAges <- model$dataAges
  famVec = model$family

  out = d$Nage
  out.mu = list()
  out.mu0 = list()

  for(a in 1:length(ages)){


    if(sampleFit && is.null(condSim)){
      m.pos = model$pModels[[a]]

      Xp.1=predict(m.pos,newdata=d[[2]],type="lpmatrix");
      brp.1=MASS::mvrnorm(n=1,coef(m.pos),m.pos$Vp);
      OS.pos = matrix(0,nrow(d[[2]]),1);
      terms.pos=terms(m.pos)
      if(!is.null(m.pos$offset)){
        off.num.pos <- attr(terms.pos, "offset")
        for (i in off.num.pos) OS.pos <- OS.pos + eval(attr(terms.pos,
                                                            "variables")[[i + 1]], d[[2]])
      }
      mu = Xp.1%*%brp.1+OS.pos
      out.mu[[a]] = exp(mu) ## obs, assumes log link!

      if(!famVec[a] %in% c("Tweedie","negbin")){
        m0 = model$zModels[[a]]
        Xp.0=predict(m0,newdata=d[[2]],type="lpmatrix");
        brp.0=MASS::mvrnorm(n=1,coef(m0),m0$Vp);
        OS0 = matrix(0,nrow(d[[2]]),1);
        terms.0=terms(m0)
        if(!is.null(m0$offset)){
          off.num.0 <- attr(terms.0, "offset")
          for (i in off.num.0) OS0 <- OS0 + eval(attr(terms.0,
                                                      "variables")[[i + 1]], d[[2]])
        }
        mu0=m0$family$linkinv(Xp.0%*%brp.0+OS0);
        out.mu0[[a]] = mu0
      }

    }
    if(!is.null(condSim)) out.mu[[a]] = condSim[[2]][[a]]
    if(!is.null(condSim) && !famVec[a] %in% c("Tweedie","negbin")) out.mu0[[a]] = condSim[[3]][[a]]

    if(famVec[a]==c("Tweedie")){
      p = model$pModels[[a]]$family$getTheta(TRUE)
      phi = model$pModels[[a]]$scale
      if(!sampleFit && is.null(condSim)) {  out.mu[[a]] = predict(model$pModels[[a]],newdata=d[[2]],type="response") }
      out[,a] = rTweedie(out.mu[[a]],p,phi)
    } else if(famVec[a]=="LogNormal"){
      if(!sampleFit && is.null(condSim)) out.mu0[[a]] = predict(model$zModels[[a]],newdata=d[[2]],type="response")
      pa = rbinom(nrow(out),1,prob=out.mu0[[a]])
      sig2 = model$pModels[[a]]$sig2
      if(!sampleFit && is.null(condSim)) out.mu[[a]] = exp(predict(model$pModels[[a]],newdata=d[[2]]))
      pos = exp( rnorm(nrow(out),log(out.mu[[a]]),sqrt(sig2)) )
      out[,a] = pa*pos
    } else if(famVec[a]=="Gamma"){
      if(!sampleFit && is.null(condSim)) out.mu0[[a]] = predict(model$zModels[[a]],newdata=d[[2]],type="response")
      pa = rbinom(nrow(out),1,prob=out.mu0[[a]])

      if(!sampleFit && is.null(condSim)){
        out.mu[[a]]= predict(model$pModels[[a]],newdata=d[[2]],type="response")
      }
      shape = 1/model$pModels[[a]]$scale
      pos = rgamma(nrow(out),rate=shape/out.mu[[a]],shape=shape)
      out[,a] = pa*pos
    }

  }
  return(list(sim=out,mu=out.mu,prob0=out.mu0))
}

# Randomized quantile residuals ------------------------------------------------

#' Randomized quantile residuals for class 'surveyIndex'
#'
#' @title Randomized quantile residuals for class 'surveyIndex'
#' @param x An object of type 'surveyIndex' as created by 'get_surveyidx'
#' @param a age group
#' @return A vector of residuals, which should be iid standard normal distributed
#' @export
residuals.surveyIdx <- function(x,a=1){
  if (pmatch("Tweedie", x$pModels[[a]]$family$family,
             nomatch = -1) == 1) {
    resi = qres_tweedie(x$pModels[[a]])
  } else if(pmatch("gaussian", x$pModels[[a]]$family$family,
                   nomatch = -1) == 1 ){## Delta-Lognormal
    y = x$allobs[[a]]
    logy = y
    psel = y>x$cutOff
    logy[psel] = log( y[psel] )

    p0 = predict(x$zModels[[a]],type="response") ## P( obs > cutOff )
    vari = x$pModels[[a]]$sig2
    cdfpos = rep(0,length(y))
    cdfpos[psel] = pnorm( q=logy[psel],mean=predict(x$pModels[[a]],type="link"), sd=sqrt(vari))
    u = (1-p0) + cdfpos*p0
    u[ !psel ] = runif(sum(!psel), min = 0, max = u[!psel])
    resi = qnorm(u)

  } else if(pmatch("Gamma", x$pModels[[a]]$family$family,
                   nomatch = -1) == 1 ){
    requireNamespace("MASS")
    y = x$allobs[[a]]
    psel = y>x$cutOff

    p0 = predict(x$zModels[[a]],type="response")
    means = predict(x$pModels[[a]],type="response")
    shape = MASS::gamma.shape(x$pModels[[a]])$alpha
    cdfpos = rep(0,length(y))
    cdfpos[psel] = pgamma( q=y[psel],shape=shape,rate=shape/means)
    u = (1-p0) + cdfpos*p0
    u[ !psel ] = runif(sum(!psel), min = 0, max = u[!psel])
    resi = qnorm(u)

  } else if(pmatch("Negative Binomial", x$pModels[[a]]$family$family,
                   nomatch = -1) == 1 ){
    y = x$pModels[[a]]$y
    size = x$pModels[[a]]$family$getTheta(TRUE)
    mu = fitted(x$pModels[[a]])
    p = size/(mu + size)
    a = ifelse(y > 0, pbeta(p, size, pmax(y, 1)), 0)
    b = pbeta(p, size, y + 1)
    u = runif(n = length(y), min = a, max = b)
    resi = qnorm(u)
  }

  resi
}

#' Randomized quantile residuals for Tweedie models
#'
#' @title Randomized quantile residuals for Tweedie models
#' @param gam.obj A gam object (mgcv package)
#' @return A vector of residuals, which should be iid standard normal distributed
#' @export
#' @import tweedie
qres_tweedie<-function (gam.obj){
  requireNamespace("tweedie")
  mu <- fitted(gam.obj)
  y <- gam.obj$y
  df <- gam.obj$df.residual
  w <- gam.obj$prior.weights
  if (is.null(w))
    w <- rep(1,length(y))
  p <- gam.obj$family$getTheta(TRUE)
  dispersion <- gam.obj$scale
  u <- numeric(length(y))
  fitt <- fitted(gam.obj)
  for(i in 1:length(u)){
    u[i] <- tweedie::ptweedie(q = y[i], power = p, mu = fitt[i],phi = dispersion/w[i])
  }
  if (p > 1 && p < 2)
    u[y == 0] <- runif(sum(y == 0), min = 0, max = u[y == 0])


  return(qnorm(u))

}





#' Plot a surveyIndexGrid
#'
#' @title Plot a surveyIndexGrid
#' @param grid a surveyIndexGrid (as created by the "get_grid" function)
#' @return nothing
#' @export
plot_surveyidx_grid<-function(grid, pch=1,gridCol="lightgrey"){
  plot(grid[[1]],grid[[2]],xlab="Longitude",ylab="Latitude",pch=pch)
  lonbps = grid[[4]]
  latbps = grid[[5]]
  for(i in 1:nrow(lonbps)) abline(v=lonbps[i,1],col=gridCol)
  abline(v=tail(lonbps,1)[2],col=gridCol)
  for(i in 1:nrow(latbps)) abline(h=latbps[i,1],col=gridCol)
  abline(h=tail(latbps,1)[2],col=gridCol)
  maps::map("worldHires", fill = TRUE, plot = TRUE,
            add = TRUE, col = grey(0.5))
}


# Standardized model outputs ---------------------------------------------------


#' Write survey index to file in standard XSA/SAM format
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
exportSI <- function(
    x,
    ages,
    years,
    toy,
    file,
    nam=""){

  cat(nam,"\n",file=file)
  cat(range(as.numeric(as.character(years))),"\n",file=file,append=TRUE)
  cat("1 1 ",rep(toy,2),"\n",file=file,append=TRUE)
  cat(min(ages),max(ages),"\n",file=file,append=TRUE)
  write.table(round(cbind(1,x[,]),4),file=file,row.names=FALSE,col.names=FALSE,append=TRUE)
}

#' Visualize results from a survey index model fitted with get_surveyidx().
#'
#' @title Visualize results from a survey index model fitted with get_surveyidx().
#' @param x Survey index as produced by getSurveyIndex()
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
plot_surveyidx<-function (x,
                          dat,
                          alt.idx = NULL,
                          myids,
                          cols = 1:length(x$pModels),
                          select = c("index", "map", "residuals", "fitVsRes"),
                          par = list(mfrow = c(3,3)),
                          colors = rev(heat.colors(6)),
                          map.cex = 1,
                          plotByAge = TRUE,
                          legend = TRUE,
                          predD = NULL,
                          year = NULL,
                          main=NULL,
                          legend.signif=3,
                          legend.pos="topright",
                          restoreOldPar=FALSE,
                          mapBubbles=FALSE,
                          cutp = NULL,
                          ...) {

  if (!plotByAge & !is.null(par)){
    op<-par(par)
    if(restoreOldPar) on.exit(par(op))
  }
  mainwasnull <- is.null(main)
  for (a in cols) {
    if(mainwasnull) main <- paste("Age group", colnames(dat$Nage)[a])
    if (plotByAge & !is.null(par)){
      op<-par(par)
      if(restoreOldPar) on.exit(par(op))
    }
    if (any(select == "index")) {
      ys = range(as.numeric(levels(dat$Year)))
      ys = ys[1]:ys[2]
      yl = range(c(x$idx[, a],0,x$lo[,a],x$up[,a])/mean(x$idx[, a]),na.rm=TRUE)
      if (!is.null(alt.idx) && a <= ncol(alt.idx)) {
        yl = range(c(alt.idx[, a]/mean(alt.idx[, a]),
                     yl)) * 1.1
        plot(ys, alt.idx[, a]/mean(alt.idx[, a], na.rm = TRUE),
             ylim = yl, col = 2, ylab = "Index", xlab = "Year",
             main = main)
      }
      else {
        plot(ys, rep(NA, length(ys)), ylim = yl, col = 2,
             ylab = "Index", xlab = "Year", main = main )
      }
      idx = x$idx
      lo = x$lo
      up = x$up
      idx[x$idx <= 0] = NA
      lo[x$idx <= 0] = NA
      up[x$idx <= 0] = NA
      lines(ys, idx[, a]/mean(idx[, a], na.rm = TRUE),
            lwd = 2)
      lines(ys, lo[, a]/mean(idx[, a], na.rm = TRUE), lwd = 2,
            lty = 2)
      lines(ys, up[, a]/mean(idx[, a], na.rm = TRUE), lwd = 2,
            lty = 2)
      if (legend && !is.null(alt.idx))
        legend(legend.pos, pch = c(1, NA), lty = c(NA,
                                                   1), col = c(2, 1), legend = c("alt.idx", "GAM"))
    }
    if (any(select == "map")) {
      xlims = range(dat$lon, na.rm = TRUE)
      ylims = range(dat$lat, na.rm = TRUE)
      mapvals = NULL
      if (is.null(predD)) {
        tmp = subset(dat, haul.id %in% myids)
      }
      else {
        tmp = predD
      }
      if (is.null(year)) {
        concT = surveyIndex:::concentration_transform(log(x$gPreds[[a]]))
        mapvals = x$gPreds[[a]]
      }
      else {
        y = which(as.numeric(as.character(names(x$gPreds2[[a]]))) ==
                    year)
        if (length(y) == 0)
          stop(paste("Year", year, "age group", a, "not found."))
        concT = surveyIndex:::concentration_transform(log(x$gPreds2[[a]][[y]]))
        mapvals = x$gPreds2[[a]][[y]]
      }
      if (length(colors) > 1)
        zFac = cut(concT, 0:length(colors)/length(colors))
      else zFac = 1
      if (length(map.cex) > 1)
        sFac = cut(log(x$gPreds[[a]]), length(map.cex))
      else sFac = 1
      myCols = colors
      plot(tmp$lon, y = tmp$lat, col = 1, pch = 1, cex = map.cex[sFac],
           xlim = xlims, ylim = ylims, xlab = "Longitude",
           ylab = "Latitude", main = main, ...)
      points(tmp$lon, y = tmp$lat, col = myCols[zFac],
             pch = 16, cex = map.cex[sFac])
      maps::map("worldHires", xlim = xlims, ylim = ylims,
                fill = TRUE, plot = TRUE, add = TRUE, col = grey(0.5))
      if (legend){
        maxcuts = aggregate(mapvals ~ zFac, FUN=max)
        mincuts = aggregate(mapvals ~ zFac, FUN=min)
        mm = mean(mapvals)
        ml = signif(mincuts[,2]/mm,legend.signif)
        ml[1] = 0
        leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
        legend(legend.pos, legend = leg, pch = 16, col = colors, bg = "white")
      }
    }
    if (any(select == "absolutemap") || any(select == "CVmap")) {
      if(is.null(year) || length(year)<1) stop("argument 'year' must be vector of length>=1 for type 'absolutemap'")
      if( !all(year %in% levels(dat$Year)) ) stop("invalid years selected")
      xlims = range(dat$lon, na.rm = TRUE)
      ylims = range(dat$lat, na.rm = TRUE)

      if(any(select == "absolutemap")) colsel = "gPreds2" else colsel = "gPreds2.CV"
      ## collect all years as data.frame
      ally = data.frame(val=x[[colsel]][[a]][[1]],year=as.character(levels(dat$Year)[1]))
      cc=0
      for(y in levels(dat$Year)){
        cc=cc+1
        ally = rbind(ally, data.frame(val=x[[colsel]][[a]][[cc]],
                                      year=as.character(levels(dat$Year)[cc])))
      }
      ally$conc = surveyIndex:::concentration_transform(log(ally$val))
      if(is.null(cutp)){
        ally$zFac=cut( ally$conc,0:length(colors)/length(colors))
      } else {
        if(length(cutp) != length(colors) + 1) stop("incompatible number of colors and length of cutp")
        ally$zFac=cut( ally$val,cutp)
      }
      bubbleScale = 0.005*max(dat$Nage[,a])
      for(yy in year){

        if (is.null(predD)) {
          tmp = subset(dat, haul.id %in% myids)
        }
        else {
          tmp = predD
          if(is.list(tmp) && !class(tmp)%in%c("data.frame","DATRASraw")) tmp = predD[[as.character(yy)]]
        }

        plot(tmp$lon,y=tmp$lat,col=1,pch=1,cex=map.cex,xlab="Longitude",ylab="Latitude",axes=FALSE)
        box()
        title(yy,line=1)
        sel = which(ally$year==yy)
        points(tmp$lon,y=tmp$lat,col=colors[as.numeric(ally$zFac[sel])],pch=16,cex=map.cex)
        maps::map('worldHires',xlim=xlims,ylim=ylims,fill=TRUE,plot=TRUE,add=TRUE,col=grey(0.5))
        if(mapBubbles){
          dy = subset(dat,Year==yy)
          points(dy$lon,dy$lat,cex=sqrt(dy$Nage[,a]/bubbleScale))
        }
        if (legend && yy==year[1]){
          if(is.null(cutp)){
            maxcuts = aggregate(val ~ zFac, data=ally, FUN=max)
            mincuts = aggregate(val ~ zFac, data=ally, FUN=min)
            mm = mean(ally$val)
            ml = signif(mincuts[,2]/mm,legend.signif)
            ml[1] = 0
            leg = paste0("[",ml,",",signif(maxcuts[,2]/mm,legend.signif),"]")
            legend(legend.pos, legend = leg, pch = 16, col = colors, bg = "white")
          }
          legend(legend.pos, legend = levels(ally$zFac), pch = 16, col = colors, bg = "white")
        }
      }##rof year
    }## absolutemap

    for (k in 1:length(select)) {
      ss = suppressWarnings(as.numeric(select[k]))
      if (!is.na(ss)) {
        plot.gam(x$pModels[[a]], select = ss, main = main, ...)
      }
    }
    if (any(select == "residuals") || any(select == "fitVsRes") ||
        any(select == "resVsYear") || any(select == "resVsShip") ||
        any(select == "spatialResiduals")) {

      resi <- x$residuals[[a]] ##residuals(x,a)
    }
    if (any(select == "residuals")) {
      hist(resi, nclass = 30, main = main, xlab = "Residuals",...)
    }
    if (any(select == "fitVsRes")) {
      plot(fitted(x$pModels[[a]]), residuals(x$pModels[[a]]), xlab = "Fitted",
           ylab = "Residuals", main = main,...) ## TODO: use quantile residuals here too
    }
    if (any(select == "resVsYear")) {
      plot(dat$Year, resi, main = main, xlab = "Year", ylab = "Residuals",
           ...)
    }
    if (any(select == "resVsShip")) {
      plot(dat$Ship, resi, main = main, xlab = "Year", ylab = "Residuals",
           ...)
    }
    if (any(select == "spatialResiduals")) {
      scale <- 3 * map.cex
      if (is.null(year) || length(year)>1)
        stop("a single year must be supplied")
      sel <- which(dat[[2]]$Year == as.character(year))
      plot(dat$lon, dat$lat, type = "n",
           xlab = "Longitude", ylab = "Latitude", main = main,...)
      maps::map("worldHires", fill = TRUE, plot = TRUE,
                add = TRUE, col = grey(0.5))
      positive = resi[sel] > 0
      points(dat$lon[sel][positive], dat$lat[sel][positive],
             pch = 1, cex = scale * sqrt(resi[sel][positive]),
             col = "blue")
      points(dat$lon[sel][!positive], dat$lat[sel][!positive],
             pch = 1, cex = scale * sqrt(-resi[sel][!positive]),
             col = "red")
    }
  }

  return(list())
}


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



#' Akaike Information Criterion (or BIC) for survey index models
#'
#' @title Akaike Information Criterion (or BIC) for survey index models
#' @param x survey index as return from get_surveyidx
#' @param BIC if TRUE compute BIC instead of AIC
#' @return numeric value
#' @export get_surveyidx_aic
get_surveyidx_aic<-function(x, BIC=FALSE){
  if(!BIC) return(2*x$edfs-2*x$logLik)
  if(pmatch("Tweedie",x$pModels[[1]]$family$family,nomatch=-1)==1 ||
     pmatch("negbin",x$pModels[[1]]$family$family,nomatch=-1)==1 ){
    log(length(x$pModels[[1]]$y)*(length(x$pModels)))*x$edfs-2*x$logLik
  } else {
    log(length(x$zModels[[1]]$y)*(length(x$zModels)))*x$edfs-2*x$logLik
  }
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
#' plot(concentration_transform(x = rnorm(n = 123, mean = 5, sd = 25)))
concentration_transform <- function (x) {
    i <- order(x)
    ys <- sort(exp(x))
    p <- ys/sum(ys)
    x[i] <- cumsum(p)
    return(x)
  }



#' Calculate external consistencies between two survey indices.
#'
#' Proper alignment of years and ages must be ensured by the user.
#' @title Calculate external consistencies between two survey indices.
#' @param tt A matrix with survey indices (rows=years, cols=ages)
#' @param tt2 A matrix with survey indices (rows=years, cols=ages)
#' @param do.plot plot it?
#' @return A vector of correlations (consistencies)
#' @export
externalCons <-
  function(tt,tt2,do.plot=FALSE){
    tt[tt==0]=NA
    tt2[tt2==0]=NA
    if(do.plot){ X11(); b=ceiling((ncol(tt))/2); par(mfrow=c(b,b));}
    for(a in 1:ncol(tt)){
      cat("Survey 1 Age ",a," vs Survey 2 ",a," : ",cor(log(tt[,a]),log(tt2[,a]),use="pairwise.complete.obs"),"\n")
      if(do.plot) {plot(log(tt[,a]),log(tt2[,a])); abline(0,1);}
    }
    return( sapply(1:ncol(tt),function(a) cor(log(tt[,a]),log(tt2[,a]),use="pairwise.complete.obs") ));
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
fix_age_group <-
  function(x,age=0,n=3,fun="mean"){
    cm.breaks<-attr(x,"cm.breaks")
    f <- match.fun(fun)
    d=split(x,x$Year)
    subsLength=f(x[[3]]$LngtCm,na.rm=TRUE)
    for(y in 1:length(d)){
      nobs = sum(d[[y]][[1]]$Age==age,na.rm=TRUE)
      if(nobs<n) {
        sel=sample(1:nrow(d[[y]][[1]]),n-nobs)
        d[[y]][[1]]=rbind(d[[y]][[1]][sel,],d[[y]][[1]]);
        d[[y]][[1]][1:(n-nobs),"Age"]=age;
        d[[y]][[1]][1:(n-nobs),"LngtCm"]=subsLength;
        d[[y]][[1]][1:(n-nobs),"NoAtALK"]=1;
      }
    }
    dd <- do.call("c",d)
    if(!is.null(cm.breaks)) dd<-addSpectrum(dd,cm.breaks=cm.breaks)
    dd
  }


# Covert and process raw data ------------------------------------------------------------------


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
get_bathy_grid <- function(
    d,
    minDepth=10,
    maxDepth=Inf,
    resolution=2,
    maxDist=Inf,
    keep=TRUE,
    shapefile=NULL,
    select=NULL){

  if (!requireNamespace("marmap", quietly = TRUE)) {
    stop("Package \"marmap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  bathy <- marmap::getNOAA.bathy(lon1=min(d$lon),lon2=max(d$lon),lat1=min(d$lat),lat2=max(d$lat),resolution=resolution, keep=keep)
  xyzbathy <- as.data.frame( as.xyz(bathy) )
  colnames(xyzbathy) <- c("lon","lat","Depth")
  xyzbathy$Depth <- -xyzbathy$Depth
  xyzbathy <- subset(xyzbathy, Depth>minDepth & Depth<maxDepth)
  if(is.finite(maxDist)){
    if (!requireNamespace("RANN", quietly = TRUE)) {
      stop("Package \"RANN\" needed for this function to work when 'maxDist' is used. Please install it.",
           call. = FALSE)
    }
    nearest <- RANN::nn2(d[[2]][,c("lon","lat")],xyzbathy[,1:2],k=1)
    xyzbathy <- xyzbathy[ nearest$nn.dists < maxDist ,]
  }
  xyzbathy <- as.data.frame(xyzbathy)
  if (is.character(shapefile)) {
    if (file.exists(shapefile)) {
      shape <- maptools::readShapeSpatial(shapefile)
    } else {
      stop("shapefile not found")
    }
    tmp <- xyzbathy
    sp::coordinates(tmp) <- ~lon + lat
    xtra <- sp::over(tmp, shape)
    if (!is.null(select))
      xtra <- xtra[select]
    xyzbathy <- cbind(xyzbathy, xtra)
  }

  xyzbathy
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
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c <- 2 * asin(min(1,sqrt(a)))
    d = R * c

    return(d) # Distance in km
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
    return(deg*pi/180)
  }



