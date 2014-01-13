\name{getSceneinfo}
\alias{getSceneinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List Landsat Scene Information from Scene IDs
}
\description{
Generate a table (data.frame) of basic scene information for Landsat5/7/8 data given input scene IDs or filenames. This function can be used to extract this information for single scenes or for an entire collection of scenes (supplied as a character vector). In both cases, a single data frame is returned displaying the sceneID, sensor, path, row, and acquisition date. This data.frame can optionally be written to .csv file for future reference.
}
\usage{
getSceneinfo(sourcefile, filename = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sourcefile}{
character vector containing names of scenes to be analyzed. The scene names must adhere to the Landsat convention. Original filenames (with the .tar.gz extension) are allowed
}
  \item{filename}{
optional output filename. If supplied, a .csv is written containing extracted scene information
}
}
\details{
It is important to supply the input character vector without changing the names of the scenes, since the getSceneinfo function looks at certain positions in the input character for the necessary information. The input vector may be comprised of original .tar.gz filenames as downloaded from (e.g.) GloVis, which allows for a quick check of scene information in batch mode directly after downloading the data.
}
\value{
a data.frame with columns representing sensor, path, row, and acquisition date for given scenes.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ben DeVries
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Extract scene information from a known sceneID
sceneID <- "LC80330322013077LGN00"
info <- getSceneinfo(sceneID)
print(info)

# Extract scene information from a scene filename
sceneID <- "LC80330322013077LGN00.tar.gz"
info <- getSceneinfo(sceneID)
print(info)

# Extract scene information from an entire folder containing subfolders with raw data
folders <- c("LE71700551999255AGS00","LE71700551999271EDC00",
"LE71700552000114SGS00","LE71700552000194EDC00","LE71700552000258SGS00",
"LE71700552001036SGS00","LE71700552001180SGS00","LE71700552001228SGS00",
"LE71700552001276SGS00", "LT51700552001060RSA00")
info <- getSceneinfo(folders)
print(info)

# Plot the number of scenes acquired per year
years <- as.numeric(substr(info$date, 1, 4))
plot(table(years), ylab="# of scenes")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~landsat }
