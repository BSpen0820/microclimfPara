pak::pkg_install("ilyamaclean/microclimf")
library(microclimf)

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

## No Snow First

micropointa <- runpointmodela(climarrayr, tme, reqhgt = 0.05, dtm, vegp, soilc)
saveRDS(micropointa, "Validation_test/base/micropointa.rds")

