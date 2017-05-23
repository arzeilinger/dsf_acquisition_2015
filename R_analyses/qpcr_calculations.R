#### Calculation of Xf populations (CFUs) in vectors from qPCR runs

rm(list = ls())
# libraries
# loading dtplyr that replaces dplyr and data.table
my.packages <- c("xlsx", "tidyr", "dplyr", "ggplot2")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/qpcrFunctions.R")


#### Import serial dilution
serial_dilution <- read.csv("data/qpcr_data/serial_dilution_cfu.csv", header = TRUE)

#### First qpcr run
#### Import plate setup and transform plate setups
exp1ps <- read.xlsx("data/qpcr_data/101416-BGSS-expe1_plate_setup.xlsx", sheetIndex = 1) 
exp1ps <- transformPlateSetup(exp1ps)

#### Import data
exp1ct <- read.csv("data/qpcr_data/101416-BGSS-expe1.csv", header = TRUE)

#### Merge plate setup and Ct data
exp1ct <- exp1ct %>% left_join(., exp1ps[,c("Well", "sample")], by = "Well") %>% dplyr::filter(., !is.na(sample))
str(exp1ct)

# Calculate CFUs from standard curve
exp1cfu <- calculateCFU(ctdata = exp1ct, serial_dilution = serial_dilution, getModel = TRUE)
exp1Mod <- exp1cfu[[1]]
exp1cfu <- exp1cfu[[2]]
exp1cfu
# Checking results and standard curve coefficients
# Compare estimated CFU
plot(x = exp1cfu$Qty, y = exp1cfu$cfu)
abline(a = 0, b = 1)
# Standard curve estimates: slope = -4.179; intercept = 49.11; R2 = 0.976
diff <- exp1cfu$Qty - exp1cfu$cfu
hist(diff)

write.csv(exp1cfu, file = "output/experiment_1_cfu.csv", row.names = FALSE)

#### Second qpcr run
#### Import plate setup and transform plate setups
exp2ps <- read.xlsx("data/qpcr_data/101916-BGSS-expe_plate_setup.xlsx", sheetIndex = 1) %>% transformPlateSetup()

#### Import data
exp2ct <- read.csv("data/qpcr_data/101916-BGSS-exp2.csv", header = TRUE)

#### Merge plate setup and Ct data
exp2ct <- exp2ct %>% left_join(., exp2ps[,c("Well", "sample")], by = "Well") %>% dplyr::filter(., !is.na(sample))
str(exp2ct)

# Calculate CFUs from standard curve
exp2cfu <- calculateCFU(ctdata = exp2ct, serial_dilution = serial_dilution, getModel = TRUE)
exp2Mod <- exp2cfu[[1]]
exp2cfu <- exp2cfu[[2]]
exp2cfu
# Checking results and standard curve coefficients
# Compare estimated CFU
plot(x = exp2cfu$Qty, y = exp2cfu$cfu)
abline(a = 0, b = 1)
# Standard curve estimates: slope = -4.179; intercept = 49.11; R2 = 0.976
diff <- exp2cfu$Qty - exp2cfu$cfu
hist(diff)


#### Third qpcr run
#### No standard curve constructed within qpcr program
#### Import plate setup and transform plate setups
exp3ps <- read.xlsx("data/qpcr_data/102016_BGSS_expe_plate_setup.xlsx", sheetIndex = 1) %>% transformPlateSetup()

#### Import data
exp3ct <- read.csv("data/qpcr_data/102016-BGSS-exp3.csv", header = TRUE)

#### Merge plate setup and Ct data
exp3ct <- exp3ct %>% left_join(., exp3ps[,c("Well", "sample")], by = "Well") %>% dplyr::filter(., !is.na(sample))
str(exp3ct)

# Calculate CFUs from standard curve
exp3cfu <- calculateCFU(ctdata = exp3ct, serial_dilution = serial_dilution, getModel = TRUE)
exp3Mod <- exp3cfu[[1]]
exp3cfu <- exp3cfu[[2]]
exp3cfu


#### Combine data from qPCR runs and save data
pcrResults <- rbind(exp1cfu, exp2cfu, exp3cfu) %>% as.data.frame()
saveRDS(pcrResults, file = "output/dsf_acquisition_qpcr_cfu_results.rds")

# check negative control
pcrResults[pcrResults$sample == "NTC",]
