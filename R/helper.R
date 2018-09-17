# This file contains common helper functions used in this package that
# are not exported

# Check if an input is a boolean (either TRUE or FALSE)
is.bool <- function(x) {
  is.logical(x) && length(x) == 1 && !is.na(x)
}

# Import a narrowPeak file

peak2Granges <- function(peakfile, broad = FALSE) {
  extraCols.narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  extraCols.broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                           qValue = "numeric")
  broad <- ifelse(grepl(".broadPeak", peakfile), TRUE, FALSE)
  
  if (isTRUE(broad)) {
    grg <- rtracklayer::import(peakfile, format = "BED",
                               extraCols = extraCols.broadPeak)
  } else {
    grg <- rtracklayer::import(peakfile, format = "BED",
                               extraCols = extraCols.narrowPeak)
  }
  return(grg)
}

# move files

file.move <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}


