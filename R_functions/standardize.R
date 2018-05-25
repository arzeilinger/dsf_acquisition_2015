#### Functions to standardize covariates on the same scale

standardize <- function(covar){
  mean.covar <- mean(covar, na.rm = TRUE)
  sd.covar <- sd(covar[!is.na(covar)])
  return((covar - mean.covar)/sd.covar)
}

standardize2 <- function(covar){
  std.covar <- covar/max(covar, na.rm = TRUE)
  return(std.covar)
}

# standardize() is preferrable; equivalent to calculating z-scores