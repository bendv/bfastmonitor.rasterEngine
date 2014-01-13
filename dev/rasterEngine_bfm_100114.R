#-----------------------------------------------------------------
# bfastmonitor with rasterEngine
# adapted from code by Jonathan Greenberg
#
# modified by Ben DeVries to work on a Landsat time series brick
# 10 January 2014
#-----------------------------------------------------------------

require("bfast")
#install.packages("spatial.tools", repos="http://cran.r-project.org")
require("spatial.tools")



############################
# {bdv}: get tura (Landsat) dataset
#install.packages("devtools", repos="http://cran.r-project.org")
library(devtools)
# install the tura package from github 
# This contains a Landsat rasterBrick dataset, and a function - getSceneinfo() which will extract scene information, including dates, from Landsat scene ID's)
# the following works on Linux:
install_git("http://github.com/bendv/tura")
# and the following on Windows:
install_github("bendv/tura")
# LazyData = TRUE; so tura dataset will load when library(tura) is called
library(tura)

# tura should now be loaded
class(tura)

# tura is a Landsat-derived NDVI time series taken from a small area in southwestern Ethiopia which has experienced recent and ongoing deforstation (2009 - present)

# getSceneinfo() is a function designed to extract relevant scene information from Landsat sceneID's
# these sceneID's are already provided as layer names in the tura raster Brick
names(tura)

# include only data from ETM+ and assign to new raster brick
s <- getSceneinfo(names(tura))
b <- dropLayer(tura, which(s$sensor=="TM"))
rm(s)

# for convenience, try a crop of b
e <- extent(c(820212, 821509, 830452, 831908))
b <- crop(b, e)
plot(b, 6)
plot(b, 90)
############################



# Default parameters for bfastmonitor:
formula=response ~ trend + harmon
order=3
lag=NULL
slag=NULL
history=c("ROC", "BP", "all") # default
history=c("ROC") # used in example
type="OLS-MOSUM"
h=0.25
end=10
level=0.05
datetype="irregular" # {bdv}: needs to be "irregular" instead of "16-day" for Landsat time series
data_multiplier <- 1/10000 # {bdv}: if original data have been rescaled

# Function to feed to rasterEngine:

bfastmonitor_rasterEngine <- function(rasterTS, start, formula,order, lag, slag, history, type, h, end, level, dates,
                                      datetype, endperiod, data_multiplier)
{
  # If needed, multiply the data:
  if(data_multiplier!=1) 
    rasterTS <- rasterTS*data_multiplier 
  
  # {bdv} NOTE: these lines do not make sense if you are testing the code inside the function line by line
  # In other words, if rasterTS is a rasterBrick object, the line dim(rasterTS) <- c(npixels, ndates) will delete all of the data in the raster!
  # It seems that rasterEngine() internally converts raster objects to arrays, 
  # so the data structure is different when executing this functino inside rasterEngine()
  # So if you want to test this line by line, set rasterTS <- as.array(b) first
  rasterTS_dims <- dim(rasterTS)
  npixels <- prod(dim(rasterTS)[1:2])
  ndates <- dim(rasterTS)[3]
  # "Flatten" the input array to a 2-d matrix:
  dim(rasterTS) <- c(npixels, ndates)
  
  # Run bfastmonitor pixel-by-pixel:
  # {bdv}: technically, these are no longer pixels, but elements in an array
  # as is the result bfm_out. The elements of bfm_out will be reassigned to a rasterBrick via rasterEngine(), but not in this function
  bfm_out <- foreach(i=seq(npixels), .packages=c("bfast"), .combine=rbind) %do% {
    bfts <- bfastts(rasterTS[i,], dates=dates, type=datetype)
    # {bdv}: trim time series using endperiod
    if(!is.null(endperiod))
      bfts <- window(bfts, end=endperiod)
    
    bfm <- bfastmonitor(data=bfts, start=start, formula=formula, order=order, lag=lag, slag=slag, history=history, type=type, h=h, end=end, level=level)
    return(c(bfm$breakpoint, bfm$magnitude)) 
  }
  
  # Coerce the output back to the correct size array:
  dim(bfm_out) <- c(rasterTS_dims[1:2], 2)
  
  # {bdv}:
  # the last line formats bfm_out as a multi-matrix array
  # ie. there are essentially two matrices (bfm_out[,,1] and bfm_out[,,2]) representing bfm$breakpoint and bfm$magnitude respectively
  # rasterEngine() seems to automatically assign these arrays back to a rasterBrick with identical dimensions as input rasterBrick
  
  return(bfm_out)
}

# register parallel backend:
sfQuickInit() # defaults to all available cores
sfQuickInit(cpus=2) # alternatively, specify number of cores to process over

# get dates
dates <- getSceneinfo(names(b))$date

# {bdv}: set list of args for bfastmonitor()
# leave, endperiod optional. TODO: set default endperiod=NULL in the function itself
args_list <- list(dates=dates, datetype=datetype, start=c(2009,1), endperiod=NULL, formula=formula, order=order,
                  lag=lag, slag=slag, history=history, type=type, h=h, end=end, level=level, data_multiplier=data_multiplier)
# optionally set endperiod to limit monitoring period
args_list$endperiod <- c(2010, 1)

# Now use rasterEngine to execute the function on the brick:
system.time(
  bfastmonitor_raster <- rasterEngine(rasterTS=b, args=args_list, fun=bfastmonitor_rasterEngine, debugmode=FALSE)
)
# {bdv}:
# for brick with extent e (defined above), processing takes 132.39 seconds on local Windows machine (1 core)
# for same brick, processing takes 94.66 seconds on local Windows machine with cpus=2 (which is the default for my system)
# on papaya:
# the same brick (extent e) on 4 cores took 28.473 seconds, but with warning message:
# "In .local(x, ...) : min value not known, use setMinMax"
# the full tura brick on 4 cores took 200.91 seconds with the same warning message
# full brick with endperiod=c(2010, 1) [it works nicely!!] on 6 cores took 138.149, again with the same warning message

# To stop parallel engine:
sfQuickStop()

plot(bfastmonitor_raster)

# use only if start=c(2009, 1) and endperiod=NULL:
library(RColorBrewer)
cols <- brewer.pal(4, "Set1")
plot(bfastmonitor_raster, 1, col=cols, legend=FALSE, main="Change Year")
legend("right", fill=cols, legend=seq(2010, 2013, by=1))

######