getSceneinfo <-
function (sourcefile, filename = "", ...) 
{
    sourcefile <- sapply(sourcefile, FUN = function(x) {
        if (grepl(".gz", x)) 
            x <- substr(x, 1, nchar(x) - 3)
        if (grepl(".tar", x)) 
            x <- substr(x, 1, nchar(x) - 4)
        return(x)
    })
    dates <- as.Date(substr(sourcefile, 10, 16), format = "%Y%j")
    sensor <- as.character(mapply(substr(sourcefile, 1, 3), dates, 
        FUN = function(x, y) {
            if (x == "LE7" & y <= "2003-03-31") "ETM+ SLC-on" else if (x == 
                "LE7" & y > "2003-03-31") "ETM+ SLC-off" else if (x == 
                "LT5" | x == "LT4") "TM" else if (x == "LC8") "OLI" else stop("Not a recognized Landsat5/7/8 scene ID.")
        }))
    path <- as.numeric(substr(sourcefile, 4, 6))
    row <- as.numeric(substr(sourcefile, 7, 9))
    info <- data.frame(sensor = sensor, path = path, row = row, 
        date = dates)
    row.names(info) <- sourcefile
    if (filename != "") 
        write.csv(info, file = filename, ...)
    return(info)
}
