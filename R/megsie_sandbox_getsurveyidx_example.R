#  library(DATRAS) # example data
#  library(surveyIndex)
#  library(tidyverse)
#  load("data/ebs_example_datras.rda") #object: ds
#  SPECIES <- c("walleye pollock")
#  YEARS <- 2015:2021
#  SRVY <- "EBS"
#  crs_proj <- "EPSG:3338" # NAD83 / Alaska Albers
#  crs_latlon <- "+proj=longlat +datum=WGS84" # decimal degrees
#
#  dat <- surveyIndex::noaa_afsc_public_foss %>%
#    dplyr::filter(srvy == SRVY &
#                    year %in% YEARS &
#                    common_name %in% SPECIES) %>%
#    dplyr::mutate(hauljoin = paste0(stratum, "_", station, "_", date_time)) %>%
#    dplyr::select(
#      year, date_time, latitude_dd, longitude_dd, # spatiotemproal data
#      cpue_kgha, common_name, # catch data
#      bottom_temperature_c, depth_m, # possible covariate data
#      srvy, area_swept_ha, duration_hr, vessel_id, hauljoin # haul/effort data)
#    )
#
#  head(dat)
#
#
#  ll <- surveyIndex::convert_crs(
#    x = dat$longitude_dd,
#    y = dat$latitude_dd,
#    crs_in = crs_latlon,
#    crs_out = crs_proj)
#
#  YEARS <- sort(unique(dat$year))
#
#  # The surveyIndex::get_surveyidx() expects some columns to be named in a specific way
#  dat_wrangled <- dat %>%
#    dplyr::rename(
#      Year = year,
#      wCPUE = cpue_kgha,
#      COMMON_NAME = common_name,
#      GEAR_TEMPERATURE = bottom_temperature_c,
#      BOTTOM_DEPTH = depth_m,
#      HaulDur = duration_hr,
#      EFFORT = area_swept_ha,
#      Ship = vessel_id) %>%
#    dplyr::mutate(
#      # create some other vars
#      Lon = longitude_dd,
#      Lat = latitude_dd,
#      lon = ll$X,
#      lat = ll$Y,
#      sx = ((longitude_dd - mean(longitude_dd, na.rm = TRUE))/1000),
#      sy = ((latitude_dd - mean(latitude_dd, na.rm = TRUE))/1000),
#      ctime = as.numeric(as.character(Year)),
#      date_time = as.Date(x = date_time, format = "%m/%d/%Y %H:%M:%S"),
#      hour = as.numeric(format(date_time,"%H")),
#      minute = as.numeric(format(date_time,"%M")),
#      day = as.numeric(format(date_time,"%d")),
#      month = as.numeric(format(date_time,"%m")),
#      TimeShotHour = hour + minute/60,
#      timeOfYear = (month - 1) * 1/12 + (day - 1)/365,
#
#      # add some dummy vars and create some other vars
#      Country = "USA",
#      Gear = "dummy",
#      Quarter = "2")  %>%
#    dplyr::mutate(across((c("Year", "Ship", "COMMON_NAME")), as.factor)) %>%
#    dplyr::select(wCPUE, GEAR_TEMPERATURE, BOTTOM_DEPTH, COMMON_NAME, EFFORT,
#                  Year, Ship, Lon, Lat, lat, lon, sx, sy,
#                  ctime, TimeShotHour, timeOfYear, Gear, Quarter, HaulDur, hauljoin)
#
#  head(dat_wrangled)
#
#  pred_grid <- surveyIndex::pred_grid_ebs
#
#  ll <- surveyIndex::convert_crs(
#    x = pred_grid$lon,
#    y = pred_grid$lat,
#    crs_in = crs_latlon,
#    crs_out = crs_proj)
#
#  pred_grid <- pred_grid %>%
#    dplyr::mutate(
#      lon = ll$X,
#      lat = ll$Y,
#      sx = ((lon - mean(lon, na.rm = TRUE))/1000),
#      sy = ((lat - mean(lat, na.rm = TRUE))/1000))
#
#  head(pred_grid)
#
#  dat_cov <- surveyIndex::pred_grid_ebs %>%
#    dplyr::select(-Shape_Area) %>%
#    dplyr::mutate(
#      sx = ((lon - mean(lon, na.rm = TRUE))/1000),
#      sy = ((lat - mean(lat, na.rm = TRUE))/1000))
#
#  sp_extrap_raster <- SpatialPoints(
#    coords = coordinates(as.matrix(dat_cov[,c("lon", "lat")])),
#    proj4string = CRS(crs_latlon) )
#
#  x <- dat_wrangled %>%
#    dplyr::select(Lon, Lat, BOTTOM_DEPTH) %>%
#    stats::na.omit()  %>%
#    sf::st_as_sf(x = .,
#                 coords = c(x = "Lon", y = "Lat"),
#                 crs = sf::st_crs(crs_latlon))
#
#  idw_fit <- gstat::gstat(formula = BOTTOM_DEPTH ~ 1,
#                          locations = x,
#                          nmax = 4)
#
#  # stn_predict <- raster::predict(idw_fit, x)
#
#  extrap_data0 <- raster::predict(
#    idw_fit, sp_extrap_raster) %>%
#    # as(sp_extrap_raster, Class = "SpatialPoints")) %>%
#    sf::st_as_sf() %>%
#    sf::st_transform(crs = crs_latlon)  %>%
#    stars::st_rasterize()
#
#  extrap_data <- stars::st_extract(x = extrap_data0,
#                                   at = as.matrix(dat_cov[,c("lon", "lat")]))
#
#  # # to make future runs of this faster:
#  # save(extrap_data0, extrap_data,
#  #      file = paste0("../inst/VigA_bottom_depth_raster_",
#  #                    min(YEARS),"-",max(YEARS), ".rdata"))
#  # EXTRAP_DATA IS IN THE INST/VIGA_BOTTOM_DEPTH_RASTER DIRECTORY
#  dat_cov <- cbind.data.frame(dat_cov,
#                              "BOTTOM_DEPTH" = extrap_data$var1.pred) %>%
#    stats::na.omit()
#
#  head(dat_cov)
#
#  tmp <- c()
#  for (i in 1:length(YEARS)) {
#    tmp <- c(tmp,
#             grep(pattern = YEARS[i], x = names(coldpool::ebs_bottom_temperature)))
#  }
#
#  extrap_data0 <- coldpool::ebs_bottom_temperature[[tmp]] %>%
#    as(., Class = "SpatialPointsDataFrame") %>%
#    sf::st_as_sf() %>%
#    sf::st_transform(crs = crs_latlon)  %>%
#    stars::st_rasterize() %>%
#    stars::st_extract(x = .,
#                      at = as.matrix(dat_cov[,c("lon", "lat")]))
#
#  names(extrap_data0) <- paste0("GEAR_TEMPERATURE", YEARS)
#
#  dat_cov <- dplyr::bind_cols(dat_cov, extrap_data0) %>%
#    na.omit()
#
#  varsbyyr <- unique( # c("GEAR_TEMPERATURE", "cpi")
#    gsub(pattern = "[0-9]+",
#         replacement = "",
#         x = names(dat_cov)[grepl(names(dat_cov),
#                                  pattern = "[0-9]+")]))
#
#  vars <- unique( # c("BOTTOM_DEPTH")
#    names(dat_cov)[!grepl(names(dat_cov),
#                          pattern = "[0-9]+")])
#  vars <- vars[!(vars %in% c("LONG", "LAT", "lon", "lat", "sx", "sy"))]
#
#  data_hauls <- dat_wrangled %>%
#    dplyr::select(Year, sx, sy,
#                  dplyr::all_of(varsbyyr), dplyr::all_of(vars),
#                  Ship, hauljoin,
#                  lat, lon, Lat, Lon,
#                  ctime, TimeShotHour, timeOfYear, Gear, Quarter,
#                  EFFORT, HaulDur)  %>%
#    # dplyr::filter(!is.na(GEAR_TEMPERATURE)) %>%
#    na.omit() %>%
#    dplyr::distinct()
#
#  data_catch <- dat_wrangled %>%
#    dplyr::select(COMMON_NAME, wCPUE, hauljoin)
#
#  dat_catch_haul <- dplyr::left_join(x = data_hauls,
#                                     y = data_catch,
#                                     by = c("hauljoin")) %>%
#    dplyr::mutate(wCPUE = ifelse(is.na(wCPUE), 0, wCPUE))
#
#  head(dat_catch_haul)
#
#  allpd <- lapply(YEARS, FUN = surveyIndex::get_prediction_grid, x = dat_cov,
#                  vars = vars, varsbyyr = varsbyyr)
#  names(allpd) <- as.character(YEARS)
#
#  head(allpd[1][[1]])
#
#  ds <- split(dat_catch_haul,dat_catch_haul$COMMON_NAME)
#  ds <- lapply(ds, surveyIndex::get_datrasraw)
#  ## OBS, response is added here in "Nage" matrix -- use wCPUE
#  ds <- lapply(ds,function(x) { x[[2]]$Nage <- matrix(x$wCPUE,ncol=1); colnames(x[[2]]$Nage)<-1; x } )
#
#  ds
#
#  save(ds, allpd, file = "data/ebspollock_example_getsurveyidx.rda")
#
# # ds and allpd I think are the only data objects you need.

 # HERE IS THE PART WHERE YOU FINALLY FIT THE MODEL
#load("data/ebspollock_example_getsurveyidx.rda")
xxx <- surveyIndex::get_surveyidx(
  x = ds[["walleye pollock"]], # this is a DATRASraw object
  ages = 1, # default
  myids = NULL, # default
  kvecP = rep(12 * 12, length(ages)), # default
  kvecZ = rep(8 * 8, length(ages)), # default
  gamma = 1, # default
  cutOff = 0, # default
  fam = "Tweedie",
  useBIC = FALSE, # default
  nBoot = 1000, # default
  mc.cores = 1, # default
  method = "ML", # default
  predD = allpd,
  modelZ = NULL, # default
  modelP = "Year +    s(sx,sy,bs=c('ts'),k=376) +     s(sx,sy,bs=c('ts'),k=10,by=Year)", # Null model spatial and temporal with an additional year effect
  knotsP = NULL, # default
  knotsZ = NULL, # default
  predfix = NULL, # default
  linkZ = "logit", # default
  CIlevel = 0.95, # default
  control = list(
    trace = TRUE,
    maxit = 20
  )
)


# Megsie todo: add example here using EBS data!
 get_surveyidx <- function(x = ds$`walleye pollock`,
                           ages = 1,
                           myids = NULL,
                           kvecP = rep(12 * 12, length(ages)),
                           kvecZ = rep(8 * 8, length(ages)),
                           gamma = 1,
                           cutOff = 0,
                           fam = "Tweedie",
                           useBIC = FALSE,
                           nBoot = 1000,
                           mc.cores = 1,
                           method = "ML",
                           predD = allpd,
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
   # Turn years into factors - we may be able to take this out
   x[[1]]$Year <- as.factor(x[[1]]$Year)
   x[[2]]$Year <- as.factor(x[[2]]$Year)

   pModels <- zModels <-
     gPreds <-  ##last data year's predictions
     gPreds2 <-  ## all years predictions
     gPreds2.CV <- ## coefficient of variation of previous
     allobs <- ## response vector (zeroes and positive)
     resid <- list() ## residuals
   predDc <- predD ##copy of predD

   if (exists(".Random.seed")) { #Not sure if this is needed either
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

   do_one_a <- function(a, cutoff) { #a is the number of ages, default = 1. cutOff is the threshold below which a value counts as zero.
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

