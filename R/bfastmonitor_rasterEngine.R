bfastmonitor_rasterEngine <-
function(x, startperiod, endperiod=NULL, dates=NULL, cpus="max", sensor="all", formula=response ~ trend + harmon,
                                      order=3, lag=NULL, slag=NULL, history=c("ROC", "BP", "all"),
                                      type="OLS-MOSUM", h=0.25, end=10, level=0.05, datetype=c("16-day", "irregular"),
                                      data_multiplier=1, filename="")
{
  # args:
    # x - rasterBrick to be analyzed
    # startperiod - start of the monitoring period (e.g. c(2009, 1))
    # endperiod - end of the monitoring period (e.g. c(2010, 1), or NULL for full time series)
    # dates - vector of dates corresponding exactly to layers of x (if omitted, then function will attempt to extract dates from layer names, if they are Landsat sceneID's)
    # cpus - number of cores to run on (for parallel processing). Default is "max", which is half of the available cores.
    # sensor - limit analysis to particular sensor (Landsat only - e.g. "TM" or "ETM+"; leave as "all" to use all data, or for MODIS ts)
      # possible values: c("all", "TM", "ETM+", "ETM+ SLC-on", "ETM+ SLC-off")
    # including other arguments for bfastmonitor() ...
    # data_multiplier - data rescaling factor; leave as 1 if not necessary
    # filename - (optional) filename if resulting rasterBrick is to be written to file
  
  # modify Landsat time series if preferred sensor is indicated
  s <- getSceneinfo(names(x))
  if(sensor!="all"){
    if("ETM+" %in% sensor){
      sensor <- unique(c(sensor, "ETM+ SLC-on", "ETM+ SLC-off"))
    }
    x <- dropLayer(x, which(!s$sensor %in% sensor))
    s <- s[which(s$sensor %in% sensor), ]
    names(x) <- row.names(s)
  }
  
  # function to run bfastmonitor() on an array derived from time series brick
  bfastmonitor_array <- function(rasterTS, start, formula, order, lag, slag, history, type, h, end, level, dates,
                                 datetype, endperiod, data_multiplier, ...)
  {
    
    # rescale data (if necessary)
    if(data_multiplier != 1){
      rasterTS <- rasterTS * data_multiplier
    }
    
    # modify dims of input array
    rasterTS_dims <- dim(rasterTS)
    npixels <- prod(dim(rasterTS)[1:2])
    ndates <- dim(rasterTS)[3]
    dim(rasterTS) <- c(npixels, ndates)
    
    bfm_out <- foreach(i=seq(npixels), .packages=c("bfast"), .combine=rbind) %do% {
      bfts <- bfastts(rasterTS[i,], dates=dates, type=datetype)
      # trim time series using endperiod
      if(!is.null(endperiod))
        bfts <- window(bfts, end=endperiod)
      
      bfm <- bfastmonitor(data=bfts, start=start, formula=formula, order=order, lag=lag, slag=slag, history=history, type=type, h=h, end=end, level=level)
      return(c(bfm$breakpoint, bfm$magnitude)) 
    }
    
    # Coerce the output back to the correct size array:
    dim(bfm_out) <- c(rasterTS_dims[1:2], 2)
         
    return(bfm_out)
  }
  
  # get dates from previous getSceneinfo() data.frame
  dates <- s$date
  
  # register parallel backend
  if(cpus == "max"){
    sfQuickInit()
  } else {
    sfQuickInit(cpus=cpus)
  }
  
  # make list of arguments to pass to bfastmonitor() within rasterEngine()
  args_list <- list(dates=dates, datetype=datetype[1], start=startperiod, endperiod=endperiod, formula=formula, order=order,
                    lag=lag, slag=slag, history=history, type=type, h=h, end=end, level=level, data_multiplier=data_multiplier)
  
  # run bfastmonitor() within rasterEngine()
  if(filename=="")
    bfastmonitor_raster <- rasterEngine(rasterTS=x, args=args_list, fun=bfastmonitor_array, debugmode=FALSE)
  else
    bfastmonitor_raster <- rasterEngine(rasterTS=x, args=args_list, fun=bfastmonitor_array, debugmode=FALSE, filename=filename)
  
  
  # stop parallel engine
  sfQuickStop()
  
  return(bfastmonitor_raster)
}
