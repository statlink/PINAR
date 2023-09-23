################################################################################
################################################################################
##  getN() --- function to count the number of events before a specified time
##  arguments:
##    x = a vector of inter-event times (assumed to be positive)
##    tt = a positive time
##  return value:
##    the number of events before time tt (possibly 0)
################################################################################

getN <- function(x, tt = 1) {
  a <- Rfast::colsums( Rfast::colCumSums(x) <= tt)  
  a[ a == dim(x)[1] ] <- NA
  a
}
