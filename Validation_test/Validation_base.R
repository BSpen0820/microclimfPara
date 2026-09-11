devtools::load_all()
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

climdata$temp <- climdata$temp - 12 # so cold enough for snow

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
    
tme <- as.POSIXlt(climdata$obs_time, tz="UTC") 

dtmc <- aggregate(dtm, 10, fun = "mean", na.rm = TRUE) 

# Run and subset point model array (using subset defaults) 
micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "Validation_test/base/micropointa.rds")
str(micropointa[[1]])

smod_slow <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc, tme = tme, altcorrect = 0, method = "slow")
saveRDS(smod_slow, "Validation_test/base/smod_slow.rds")
str(smod_slow)

mout_slow_snow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, snow = T, snowmod = smod_slow)
saveRDS(mout_slow_snow, "Validation_test/base/mout_slow_snow.rds")
str(mout_slow_snow)
rm(smod_slow, mout_slow_snow); gc()

smod_fast <- runsnowmodel(climarrayr, micropointa, vegp, soilc, dtm, dtmc, tme = tme, altcorrect = 0, method = "fast")
saveRDS(smod_fast, "Validation_test/base/smod_fast.rds")
str(smod_fast)

mout_fast_snow <- runmicro(micropointa, reqhgt = 0.05, vegp, soilc, dtm, dtmc, altcorrect = 0, snow = T, snowmod = smod_fast)
saveRDS(mout_fast_snow, "Validation_test/base/mout_fast_snow.rds")
str(mout_fast_snow)
rm(smod_fast, mout_fast_snow); gc()