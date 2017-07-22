#### Calculation of Xf populations (CFUs) in vectors from qPCR runs
#### Cq values are calculated from LinRegPCR software program....
#### based on the analysis of Ruijter et al. (2009) Nucleic Acids Research

rm(list = ls())
# libraries
# replaced package xlsx with openxlsx because the latter doesn't rely on rJava which was causing problems
my.packages <- c("openxlsx", "tidyr", "dplyr", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")

# Directory for qPCR data
qpcrDir <- "data/qpcr_data/"

#### Import serial dilution
serial_dilution <- read.csv("data/qpcr_data/serial_dilution_cfu.csv", header = TRUE)

#### First qpcr run
#### Import plate setup and transform plate setups
exp1ps <- read.xlsx("data/qpcr_data/101416-BGSS-expe1_plate_setup.xlsx", sheet = 1) 
exp1ps <- transformPlateSetup(exp1ps)

#### Import qpcr data from LinRegPCR
exp1ct <- readLinReg(file = "101416_BGSS_expe1_linreg.xlsx", dir = qpcrDir)

#### Merge plate setup and Ct data
exp1ct <- exp1ct %>% left_join(., exp1ps[,c("wellNumber", "sample")], by = "wellNumber") %>% dplyr::filter(., !is.na(sample))
str(exp1ct)

# Calculate CFUs from standard curve
# Using standard curve data from HL5/HL6 primers for Xylella fastidiosa spiked with BGSS extract

exp1Results <- calculateCFU(qpcrdata = exp1ct, serial_dilution = serial_dilution, getModel = TRUE)
exp1scurve <- exp1Results[[1]]
exp1Mod <- exp1Results[[2]]
exp1cfu <- exp1Results[[3]]
# plot standard curve
plot(x = log10(exp1scurve$cfu), y = log10(exp1scurve$N0))
summary(exp1Mod)
exp1cfu
# Checking results and standard curve coefficients

write.csv(exp1cfu, file = "output/experiment_1_cfu_linreg.csv", row.names = FALSE)

#### Second qpcr run
#### Import plate setup and transform plate setups
exp2ps <- read.xlsx("data/qpcr_data/101916-BGSS-expe_plate_setup.xlsx", sheet = 1) %>% transformPlateSetup()

#### Import data
exp2ct <- readLinReg(file = "101916_BGSS_exp2_linreg.xlsx", dir = qpcrDir)

#### Merge plate setup and Ct data
exp2ct <- exp2ct %>% left_join(., exp2ps[,c("wellNumber", "sample")], by = "wellNumber") %>% dplyr::filter(., !is.na(sample))
str(exp2ct)

# Calculate CFUs from standard curve
exp2Results <- calculateCFU(qpcrdata = exp2ct, serial_dilution = serial_dilution, getModel = TRUE)
exp2scurve <- exp2Results[[1]]
exp2Mod <- exp2Results[[2]]
exp2cfu <- exp2Results[[3]]
# Look at results, plot standard curve
plot(x = log10(exp2scurve$cfu), y = log10(exp2scurve$N0))
summary(exp2Mod)
exp2cfu


#### Third qpcr run
#### No standard curve constructed within qpcr program
#### Import plate setup and transform plate setups
exp3ps <- read.xlsx("data/qpcr_data/102016_BGSS_expe_plate_setup.xlsx", sheet = 1) %>% transformPlateSetup()

#### Import data
exp3ct <- readLinReg(file = "102016_BGSS_exp3_linreg.xlsx", dir = qpcrDir)

#### Merge plate setup and Ct data
exp3ct <- exp3ct %>% left_join(., exp3ps[,c("wellNumber", "sample")], by = "wellNumber") %>% dplyr::filter(., !is.na(sample))
str(exp3ct)

# Calculate CFUs from standard curve
exp3Results <- calculateCFU(qpcrdata = exp3ct, serial_dilution = serial_dilution, getModel = TRUE)
exp3scurve <- exp3Results[[1]]
exp3Mod <- exp3Results[[2]]
exp3cfu <- exp3Results[[3]]
# Look at results, plot standard curve
plot(x = log10(exp3scurve$cfu), y = log10(exp3scurve$N0))
summary(exp3Mod)
exp3cfu


#### Combine data from qPCR runs and save data
pcrResults <- rbind(exp1cfu, exp2cfu, exp3cfu) %>% as.data.frame()
saveRDS(pcrResults, file = "output/dsf_acquisition_qpcr_cfu_linreg_results.rds")

# check negative control
pcrResults[pcrResults$sample == "NTC",]
