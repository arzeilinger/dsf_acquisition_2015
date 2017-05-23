#### Functions to manipulate qPCR Ct data for Xylella populations in vectors (and plants)

#### Convert a factor variable into a numeric variable
factor2numeric <- function(factorx){
  numx <- as.numeric(levels(factorx))[factorx]
  return(numx)
}


#### Function to transform a plate setup Excel sheet into a "long" data frame and create a variable with Well identifiers
# Assumes that the plate setup has letter rows and number columns, 1 - 12, and that the letter rows are in first column
# Also assumes that the sample names are idenfitied in each well
# Returns "long" data frame with columns: "Well" that idenfies each unique well, and "sample" that has the sample names
# Can then easily be joined to the data frame exported from the qPCR software program
transformPlateSetup <- function(plateSetup){
  require(tidyr); require(dplyr)
  plateSetup2 <- plateSetup %>% gather(., key = column, value = sample, 2:13) %>% arrange_(., names(.)[1])
  plateSetup2$wellNumber <- plateSetup2 %>% row.names() %>% as.integer()
  plateSetup2$well <- with(plateSetup2, paste(row, column, sep=""))
  return(plateSetup2)
}


#### Function to read in Cq data output from LinRegPCR software program
#### Only imports the data table in the "compact" sheet
#### Also strips the fluorophore from the wellNumber
#### File name must include ".xlsx" at end
readLinReg <- function(file, dir = "data/qpcr_data/"){
  filePath <- paste(dir, file, sep = "")
  # Specifications to read in only the data table
  cqdata <- read.xlsx(filePath, sheet = 3, cols = c(1:7), startRow = 4, colNames = TRUE) %>%
    # Get well number from "name" column and make it a separate column "wellNumber"
    mutate(., wellNumber = tstrsplit(name, "_", keep = 1, type.convert = TRUE)[[1]])
  return(cqdata)
}

#### FUNCTION TO CALCULATE CFU FROM N0 VALUES FROM QPCR AND LinRegPCR PROGRAM
calculateCFU <- function(qpcrdata, serial_dilution, getModel = FALSE){
  require(tidyr); require(dplyr)
  # Function assumes that "qpcrdata" contains N0 values for experimental samples and standards in a column named "N0".... 
  # and a column specifying dilution names "D1" - "D6" named "sample"
  # N0 = threshold/(mean_efficiency^Cq), from LinRegPCR
  # Also assumes a serial dilution data.frame has been loaded with columns: "dilution" with values D1 - D6 and "cfu"
  # Appends estimated CFU values to the ctdata object, in a column names "cfu"
  # Merge 'serial_dilution' and 'qpcrdata' sets, where 'dilution' and 'sample' columns contain the dilution names D1 - D6
  scurve <- right_join(serial_dilution, qpcrdata[,c("sample", "N0")], by = c("dilution" = "sample"))
  # Standard curve linear regression
  scMod <- lm(log10(N0) ~ log10(cfu), data = scurve)
  modelResults <- summary(scMod)
  # Back-calculate CFU from standard curve
  getCFU <- function(n0, slope, intercept){
    cfu <- 10^((log10(n0)-intercept)/slope)
    return(cfu)
  }
  qpcrdata$cfu <- getCFU(qpcrdata$N0, slope = coef(scMod)[2], intercept = coef(scMod)[1])
  qpcrdata$cfu[is.na(qpcrdata$cfu)] <- 0
  # If getModel = TRUE, function returns a list containing the model summary and the data set including CFU
  # If getModel = FALSE, function returns just the data set including CFU
  if(getModel == TRUE){
    output <- list(modelResults, qpcrdata)
  } else {output <- qpcrdata}
  return(output)
}