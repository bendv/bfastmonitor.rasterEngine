\name{bfastmonitor_rasterEngine}
\alias{bfastmonitor_rasterEngine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply bfastmonitor() using rasterEngine()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bfastmonitor_rasterEngine(x, startperiod, endperiod = NULL, dates = NULL, cpus = "max", sensor = "all", formula = response ~ trend + harmon, order = 3, lag = NULL, slag = NULL, history = c("ROC", "BP", "all"), type = "OLS-MOSUM", h = 0.25, end = 10, level = 0.05, datetype = c("16-day", "irregular"), data_multiplier = 1, filename = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{startperiod}{
%%     ~~Describe \code{startperiod} here~~
}
  \item{endperiod}{
%%     ~~Describe \code{endperiod} here~~
}
  \item{dates}{
%%     ~~Describe \code{dates} here~~
}
  \item{cpus}{
%%     ~~Describe \code{cpus} here~~
}
  \item{sensor}{
%%     ~~Describe \code{sensor} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{order}{
%%     ~~Describe \code{order} here~~
}
  \item{lag}{
%%     ~~Describe \code{lag} here~~
}
  \item{slag}{
%%     ~~Describe \code{slag} here~~
}
  \item{history}{
%%     ~~Describe \code{history} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{end}{
%%     ~~Describe \code{end} here~~
}
  \item{level}{
%%     ~~Describe \code{level} here~~
}
  \item{datetype}{
%%     ~~Describe \code{datetype} here~~
}
  \item{data_multiplier}{
%%     ~~Describe \code{data_multiplier} here~~
}
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, startperiod, endperiod = NULL, dates = NULL, cpus = "max", 
    sensor = "all", formula = response ~ trend + harmon, order = 3, 
    lag = NULL, slag = NULL, history = c("ROC", "BP", "all"), 
    type = "OLS-MOSUM", h = 0.25, end = 10, level = 0.05, datetype = c("16-day", 
        "irregular"), data_multiplier = 1, filename = "") 
{
    s <- getSceneinfo(names(x))
    if (sensor != "all") {
        if ("ETM+" \%in\% sensor) {
            sensor <- unique(c(sensor, "ETM+ SLC-on", "ETM+ SLC-off"))
        }
        x <- dropLayer(x, which(!s$sensor \%in\% sensor))
        s <- s[which(s$sensor \%in\% sensor), ]
        names(x) <- row.names(s)
    }
    bfastmonitor_array <- function(rasterTS, start, formula, 
        order, lag, slag, history, type, h, end, level, dates, 
        datetype, endperiod, data_multiplier, ...) {
        if (data_multiplier != 1) {
            rasterTS <- rasterTS * data_multiplier
        }
        rasterTS_dims <- dim(rasterTS)
        npixels <- prod(dim(rasterTS)[1:2])
        ndates <- dim(rasterTS)[3]
        dim(rasterTS) <- c(npixels, ndates)
        bfm_out <- foreach(i = seq(npixels), .packages = c("bfast"), 
            .combine = rbind) \%do\% {
            bfts <- bfastts(rasterTS[i, ], dates = dates, type = datetype)
            if (!is.null(endperiod)) 
                bfts <- window(bfts, end = endperiod)
            bfm <- bfastmonitor(data = bfts, start = start, formula = formula, 
                order = order, lag = lag, slag = slag, history = history, 
                type = type, h = h, end = end, level = level)
            return(c(bfm$breakpoint, bfm$magnitude))
        }
        dim(bfm_out) <- c(rasterTS_dims[1:2], 2)
        return(bfm_out)
    }
    dates <- s$date
    if (cpus == "max") {
        sfQuickInit()
    }
    else {
        sfQuickInit(cpus = cpus)
    }
    args_list <- list(dates = dates, datetype = datetype[1], 
        start = startperiod, endperiod = endperiod, formula = formula, 
        order = order, lag = lag, slag = slag, history = history, 
        type = type, h = h, end = end, level = level, data_multiplier = data_multiplier)
    if (filename == "") 
        bfastmonitor_raster <- rasterEngine(rasterTS = x, args = args_list, 
            fun = bfastmonitor_array, debugmode = FALSE)
    else bfastmonitor_raster <- rasterEngine(rasterTS = x, args = args_list, 
        fun = bfastmonitor_array, debugmode = FALSE, filename = filename)
    sfQuickStop()
    return(bfastmonitor_raster)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
