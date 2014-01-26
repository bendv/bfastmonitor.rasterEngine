nobs_rasterEngine <- function(x, cpus="max", filename="", asPerc=FALSE){
  # function to calculate number of valid observations per pixel
  # args:
    # x - input rasterBrick or Stack
    # cpus - number of cpus to use. Defaults to "max" (half of all available cpus)
    # filename - optional: write to file
    # asPerc - express as a percentage of nlayers(x)?
  
  ### Returns the following error when run on tura:
    ### chunk processing units require array vector outputs.  Please check your function.
    ### Error in spatial.tools:::focal_hpc_test(x, fun, window_center, window_dims, : 21
  
  if(!class(x) %in% c("RasterBrick", "RasterStack"))
    stop("x must be a RasterBrick or RasterStack")
  
  # function to pass to rasterEngine()
  arrayFunction <- function(r, asPerc=FALSE, ...){
    
    # melt the input array to a 2-D matrix (pixel, layers)
    rdims <- dim(r)
    np <- prod(rdims[1:2])
    nl <- rdims[3]
    dim(r) <- c(np, nl)
    
    # output 1-D vector with nobs per pixel
    obs <- foreach(i=seq(np), .combine=rbind) %do% {
      nobs <- length(r[i, !is.na(r[i, ])])
      if(asPerc)
        nobs <- nobs / nl * 100
      return(nobs)
    }
    
    # coerce array to a 2-D array with original spatial dims
    dim(obs) <- rdims[c(1:2)]
    
    return(obs)
  }
  
  # register parallel backend
  if(cpus == "max"){
    sfQuickInit()
  } else {
    sfQuickInit(cpus=cpus)
  }
  
  # run function over raster in rasterEngine()
  args_list <- list(asPerc=asPerc)
  if(filename=="")
    obs_raster <- rasterEngine(r=x, args=args_list, fun=arrayFunction, debugmode=FALSE)
  else
    obs_raster <- rasterEngine(r=x, args=args_list, fun=arrayFunction, debugmode=FALSE, filename=filename)
  
  # stop parallel engine
  sfQuickStop()
  
  return(obs_raster)
}