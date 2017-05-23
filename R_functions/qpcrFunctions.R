#### Functions to manipulate qPCR Ct data for Xylella populations in vectors (and plants)

#### Convert a factor variable into a numeric variable
factor2numeric <- function(factorx){
  numx <- as.numeric(levels(factorx))[factorx]
  return(numx)
}


#### Function to transform a plate setup Excel sheet into a "long" data frame and create a variable with Well identifiers
# Assumes that the plate setup has letter rows and number columns, 1 - 12, and that the letter rows are in first column
# Also assumes that the sample names are idenfitied in each well
# R automatically adds an "X" in front of the number column names when imported, need to strip that out
# Returns "long" data frame with columns: "Well" that idenfies each unique well, and "sample" that has the sample names
# Can then easily be joined to the data frame exported from the qPCR software program
transformPlateSetup <- function(plateSetup){
  require(tidyr); require(dplyr)
  plateSetup2 <- plateSetup %>% gather(., key = column, value = sample, X1:X12)
  plateSetup2$Well <- plateSetup2$column %>% gsub(pattern = "X", x = ., replacement = "") %>% paste(plateSetup2[,1], ., sep="")
  return(plateSetup2)
}


#### FUNCTION TO CALCULATE CFU FROM CT VALUES FROM QPCR
#### Using standard curve data from HL5/HL6 primers for Xylella fastidiosa spiked with BGSS extract
calculateCFU <- function(ctdata, serial_dilution, getModel = FALSE){
  require(tidyr); require(dplyr)
  # Function assumes that "ctdata" contains Ct values for experimental samples and standards in a column named "Ct", 
  # and a colmn specifying dilution names "D1" - "D6" named "sample"
  # Also assumes a serial dilution data.frame has been loaded with columns: "dilution" with values D1 - D6 and "cfu"
  # Function also removes rows without Ct values
  # Appends estimated CFU values to the ctdata object, in a column names "cfu"
  ctdata$Ct <- factor2numeric(ctdata$Ct)
  #ctdata <- ctdata[!is.na(ctdata$Ct),]
  scurve <- right_join(serial_dilution, ctdata[,c("sample", "Ct")], by = c("dilution" = "sample"))
  # Standard curve linear regression
  scMod <- lm(Ct ~ log10(cfu), data = scurve)
  modelResults <- summary(scMod)
  # Back-calculate CFU from standard curve
  getCFU <- function(ct, slope, intercept){
    cfu <- 10^((ct-intercept)/slope)
    return(cfu)
  }
  ctdata$cfu <- getCFU(ctdata$Ct, slope = coef(scMod)[2], intercept = coef(scMod)[1])
  ctdata$cfu[is.na(ctdata$cfu)] <- 0
  # If getModel = TRUE, function returns a list containing the model summary and the data set including CFU
  # If getModel = FALSE, function returns just the data set including CFU
  if(getModel == TRUE){
    output <- list(modelResults, ctdata)
  } else {output <- ctdata}
  return(output)
}