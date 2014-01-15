# Author: Loic Dutrieux
# June 2013
# loic.dutrieux@wur.nl

## Warning:::: I harcoded the datatype and the bandorder of the intermediary output in this version of the function

mc.calc <- function(x, fun, mc.cores=1, ...) {
    
    if(mc.cores == 1) { # Normal calc
        out <- calc(x=x, fun=fun, ...)
        return(out)
    } else {
        
        s <- blockSize(x, minblocks=mc.cores)
        blocs <- seq(1, s$n)
        
        
        # Create blocks and write to file
        fun2 <- function(i) {
            e <- extent(x, r1=s$row[i], r2=s$row[i]+s$nrows[i]-1)
            # tmp <- rasterTmpFile()
            b <- crop(x, e)
            out <- calc(x=b, fun=fun) # Does this line need an elipsis
            return(out)
        }
        
        listOut <- mclapply(X=blocs, FUN=fun2, mc.cores=mc.cores)
        
        
        # Mosaic and write to filename
        if(hasArg(filename)){
            listOut$filename <- filename
        }
        if(hasArg(overwrite)){
            listOut$overwrite <- overwrite
        }
        listOut$fun <- max
        out <- do.call(mosaic, listOut)
        
        return(out)
    }
    
}