\name{errTest}
\alias{errTest}
\docType{data}
\title{
Small subset of tura with NA's inserted into cell 51 for dates before 2005.
}
\description{
This small rasterBrick can be used for testing error catching. 
}
\usage{data(errTest)}
\format{
}
\details{
bfastmonitor() will return an error for this brick if startperiod=c(2006, 1) (e.g.) since all values for one pixel for dates before 2005 have been removed, making it impossible to construct a stable history period for that. bfastmonitor_rasterEngine_tryCatch() should assign values of -9999 for this pixel and continue processing all other pixels.
}
\source{
}
\references{
}
\examples{
data(errTest)

\dontrun{
# bfastmonitor_rasterEngine() without error catching returns an error without completing the tastk
bfm <- bfastmonitor_rasterEngine(errTest, startperiod=c(2006, 1), cpus=1, datetype="irregular")
}

# with tryCatch() this error is skipped and a -9999 assigned
bfm <- bfastmonitor_rasterEngine_tryCatch(errTest, startperiod=c(2006, 1), cpus=1, datetype="irregular", printErrors=FALSE)

# extract change(breakpoints) and magnitude and remove -9999 values
change <- raster(bfm, 1)
change[change == -9999] <- NA
magn <- raster(bfm, 2)
magn[magn == -9999] <- NA

# produce an error mask from previous results
errMask <- raster(bfm)
errMask[bfm[[1]] == -9999] <- 1

}

\keyword{datasets}
