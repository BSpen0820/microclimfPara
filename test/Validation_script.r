remove.packages("microclimf")
remove.packages("microclimfPara")

pak::pak_cleanup(force = TRUE)


pak::pkg_install("ilyamaclean/microclimf")
library(microclimf)
library(terra)

# data prep 
data(climdata)
data(vegp)
data(soilc)
data(dtmcaerth)

.rast <- function(m,tem) { 
    r<-rast(m) 
    ext(r)<-ext(tem) 
    crs(r)<-crs(tem) 
    r 
    } 
    
.ta<-function(x,dtm,xdim=5,ydim=5) { 
    a<-array(rep(x,each=ydim*xdim),dim=c(ydim,xdim,length(x))) 
    .rast(a,dtm) } 
    

# Create dummy array datasets 

dtm <- rast(dtmcaerth) # unpack raster 

climdata_snow <- climdata
climdata_snow$temp <- climdata$temp - 12 # so cold enough for snow

climarrayr<-list(
    temp = .ta(climdata$temp, dtm), 
    relhum = .ta(climdata$relhum, dtm), 
    pres = .ta(climdata$pres, dtm), 
    swdown = .ta(climdata$swdown, dtm), 
    difrad = .ta(climdata$difrad, dtm), 
    lwdown = .ta(climdata$lwdown, dtm), 
    windspeed = .ta(climdata$windspeed, dtm), 
    winddir = .ta(climdata$winddir, dtm), 
    precip = .ta(climdata$precip, dtm)) 
    

climarrayr_snow<-list(
    temp = .ta(climdata_snow$temp, dtm), 
    relhum = .ta(climdata_snow$relhum, dtm), 
    pres = .ta(climdata_snow$pres, dtm), 
    swdown = .ta(climdata_snow$swdown, dtm), 
    difrad = .ta(climdata_snow$difrad, dtm), 
    lwdown = .ta(climdata_snow$lwdown, dtm), 
    windspeed = .ta(climdata_snow$windspeed, dtm), 
    winddir = .ta(climdata_snow$winddir, dtm), 
    precip = .ta(climdata_snow$precip, dtm)) 


tme <- as.POSIXlt(climdata$obs_time, tz="UTC") 

dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE) 

# Base Run with orginal microclimf package 
dir.create('./Base_Results', showWarnings = FALSE)

## No Snow First

micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "Base_Results/micropointa_nosnow.rds", compress = FALSE, ascii = FALSE)

mout_nosnow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0)
saveRDS(mout_nosnow, "Base_Results/mout_nosnow.rds", compress = FALSE, ascii = FALSE)

## Now with Snow\

micropointa_snow <- runpointmodela(climarrayr_snow, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa_snow, "Base_Results/micropointa_snow.rds", compress = FALSE, ascii = FALSE)

smod_slow <- runsnowmodel(climarrayr_snow, micropointa_snow, vegp, soilc, dtm, dtmc, tme = tme, altcorrect = 0, method = "slow")
saveRDS(smod_slow, "Base_Results/smod_snow.rds", compress = FALSE, ascii = FALSE)

mout_snow <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, snow = T, snowmod = smod_slow)
saveRDS(mout_snow, "Base_Results/mout_snow.rds", compress = FALSE, ascii = FALSE)

rm(micropointa, mout_nosnow, micropointa_snow, smod_slow, mout_snow); gc()

detach("package:microclimf", unload = TRUE)

pak::pkg_install("BSpen0820/microclimfPara")
library(microclimfPara)

# run in serial first 

## No Snow First
dir.create('./MicroPar_Ser', showWarnings = FALSE)

micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "MicroPar_Ser/micropointa_nosnow.rds", compress = FALSE, ascii = FALSE)

mout_nosnow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0)
saveRDS(mout_nosnow, "MicroPar_Ser/mout_nosnow.rds", compress = FALSE, ascii = FALSE)

## Now with Snow

micropointa_snow <- runpointmodela(climarrayr_snow, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa_snow, "MicroPar_Ser/micropointa_snow.rds", compress = FALSE, ascii = FALSE)

smod_slow <- runsnowmodel(climarrayr_snow, micropointa_snow, vegp, soilc, dtm, dtmc, tme = tme, altcorrect = 0, method = "slow")
saveRDS(smod_slow, "MicroPar_Ser/smod_snow.rds", compress = FALSE, ascii = FALSE)

mout_snow <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, snow = T, snowmod = smod_slow)
saveRDS(mout_snow, "MicroPar_Ser/mout_snow.rds", compress = FALSE, ascii = FALSE)

rm(micropointa, mout_nosnow, micropointa_snow, smod_slow, mout_snow); gc()

# run in parallel now 

## No Snow First
dir.create('./MicroPar_Par', showWarnings = FALSE)

micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "MicroPar_Par/micropointa_nosnow.rds", compress = FALSE, ascii = FALSE)

mout_nosnow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, parallel = TRUE, ncores = 4)
saveRDS(mout_nosnow, "MicroPar_Par/mout_nosnow.rds", compress = FALSE, ascii = FALSE)

## Now with Snow

micropointa_snow <- runpointmodela(climarrayr_snow, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa_snow, "MicroPar_Par/micropointa_snow.rds", compress = FALSE, ascii = FALSE)

smod_slow <- runsnowmodel(climarrayr_snow, micropointa_snow, vegp, soilc, dtm, dtmc, tme = tme, altcorrect = 0, method = "slow", parallel = TRUE, ncores = 4)
saveRDS(smod_slow, "MicroPar_Par/smod_snow.rds", compress = FALSE, ascii = FALSE)

mout_snow <- runmicro(micropointa_snow, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, parallel = TRUE, ncores = 4, snow = T, snowmod = smod_slow)
saveRDS(mout_snow, "MicroPar_Par/mout_snow.rds", compress = FALSE, ascii = FALSE)

rm(micropointa, mout_nosnow, micropointa_snow, smod_slow, mout_snow); gc()