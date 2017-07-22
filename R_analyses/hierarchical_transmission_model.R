#### Hierarchical Bayesian model for transmission analysis

#### Attempt in nimble first
rm(list = ls())

my_packages<-c('tidyr', 'dplyr', 'data.table', 'nimble', 'coda')
lapply(my_packages, require, character.only=T)

source("R_functions/standardize.R")
source("R_functions/nimble_definitions.R")

acqdat <- read.csv("data/DSF_acquisition_data.csv")
summary(acqdat)
# convert to centimeters
acqdat$distance <- acqdat$distance*2.54
acqdat$sample <- acqdat %>% with(., paste(genotype, rep, cage, sep=""))

#### qPCR data set 
pcrResults <- readRDS("output/dsf_acquisition_qpcr_cfu_linreg_results.rds")
pcrResults$sample <- gsub(pattern = " ", replacement = "", x = pcrResults$sample)
# Remove non-plant samples and select only sample and CFU columns
pcrResults <- pcrResults[grep("F", pcrResults$sample),] %>% dplyr::select(., sample, cfu)
table(pcrResults$sample)
# Create replicate column
pcrResults$rep <- NA
for(i in 1:length(unique(pcrResults$sample))){
  sample.i <- unique(pcrResults$sample)[i]
  data.i <- which(pcrResults$sample == sample.i)
  pcrResults$rep[data.i] <- c(1:length(data.i))
}
# Spread out data set
pcrSpread <- pcrResults %>% spread(., key = rep, value = cfu)
names(pcrSpread) <- c("sample", "cfu1", "cfu2")
# Average the replicates for each sample
pcrSpread$meancfu <- rowMeans(pcrSpread[,c("cfu1", "cfu2")], na.rm = TRUE)

##############################################################################################
#### Comparing linReg results to previous results for CFUs
oldResults <- readRDS("output/dsf_acquisition_qpcr_cfu_results.rds")
oldResults$sample <- gsub(pattern = " ", replacement = "", x = oldResults$sample)
# Remove non-plant samples and select only sample and CFU columns
oldResults <- oldResults[grep("F", oldResults$sample),] %>% dplyr::select(., sample, cfu)
# Create replicate column
oldResults$rep <- NA
for(i in 1:length(unique(oldResults$sample))){
  sample.i <- unique(oldResults$sample)[i]
  data.i <- which(oldResults$sample == sample.i)
  oldResults$rep[data.i] <- c(1:length(data.i))
}

qpcrDiff <- full_join(pcrResults, oldResults, by = c("sample", "rep"))
# Compare CFUs
plot(x = qpcrDiff$cfu.x, y = qpcrDiff$cfu.y, xlab = "CFUs from LinRegPCR", ylab = "Uncorrected CFUs")
abline(a = 0, b = 1, lty = 2)
#### Note: LinRegPCR puts more samples as zero than my CFU calculations. My numbers are generally biased upward too.
############################################################################################


#### combine acquisition data and qpcr data by sample (genotype, cage, rep combination)
acqdat <- acqdat %>% full_join(., pcrSpread, by = "sample")

#### Set up factor variables
acqdat$genotypeFactor <- acqdat$genotype %>% factor(.) %>% as.numeric()
# Double check that factor levels are correct
# Correct alignment is FT = 1, FW = 2
data.frame(acqdat$genotype, acqdat$genotypeFactor)
# Make interaction variable
acqdat <- mutate(acqdat, genotype_distance = genotypeFactor*distance)

saveRDS(acqdat, file = "output/dsf_acquisition_pcr_full_data_set.rds")

#### standardize continuous covariates
# I don't think I need to log transform Xf populations
# acqdat$log.source.plant.pop <- log10(acqdat$source.plant.pop + 1)
# acqdat$log.meancfu <- log10(acqdat$meancfu + 1)
# acqdat$log.cfu1 <- log10(acqdat$cfu1 + 1)
# acqdat$log.cfu2 <- log10(acqdat$cfu2 + 1)
covars <- c("distance", "genotype_distance", "source.plant.pop", "meancfu")
covars.i <- as.numeric(sapply(covars, function(x) which(names(acqdat) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(acqdat)[i]
  stdname.i <- paste("std", var.i, sep = ".")
  stdvar.i <- standardize(acqdat[,var.i])
  acqdat[,stdname.i] <- stdvar.i
}
str(acqdat)

#### More cleaning of data
# include a plant term for the random effect
acqdat$plant <- paste(acqdat$genotype, acqdat$rep, sep = "") %>% factor()
# Remove cages placed below the point of inoculation and an outlier above 100 cm
acqdat <- acqdat %>% dplyr::filter(., distance >= 0 & distance < 100)


N <- nrow(acqdat)
acqdat$plantID <- acqdat$plant %>% factor(.) %>% as.numeric()
nplant <- acqdat$plantID %>% unique() %>% length()

nimbleTransData <- with(acqdat,
                        list(N = N,
                             nplant = nplant,
                             genotype = genotypeFactor,
                             distance = std.distance,
                             genotype_distance = std.genotype_distance,
                             xf_source_plant = source.plant.pop, # response variables don't need to be standardized
                             xf_vector = meancfu, # response variables don't need to be standardized
                             infected = test.plant.infection,
                             plantID = acqdat$plantID))

saveRDS(nimbleTransData, "output/hierarchical_transmission_nimble_data.rds")

nimbleTransData <- readRDS("output/hierarchical_transmission_nimble_data.rds")
str(nimbleTransData)

#### Exporting the data set for nimble and metadata describing the data set
metadata <- c("N = total sample size",
              "nplant = number of distinct source plants. Two insects were caged on each source plant, so source plant ID needs to be included as a random effect",
              "genotype = factor variable on genotype of source plant: 1 = DSF, 2 = WT",
              "distance = distance along stem of source plant from point of inoculation to where insect was caged (in cm)",
              "genotype_distance = interaction between genotype and distance",
              "xf_source_plant = estimated population of Xylella (the bacteria) in the source plant where the insect was caged",
              "xf_vector = estimated population of Xylella in the insect vector",
              "infected = binary variable indicating if the test plant was infected (1) or not infected (0)",
              "plantID = ID of source plant, for random effect",
              "Note: distance and genotype_distance are standardized around the mean and scaled by their standard deviation. Other variables are not standardized.")

dsfData <- list(metadata = metadata, dsfData = nimbleTransData)
saveRDS(dsfData, file = "output/dsf_data_for_nimble.rds")

# Save to DT's 'zeilinger_grapevine' repo
saveRDS(dsfData, file = "../../zeilinger_grapevine/data/dsf_data_for_nimble.rds")

data3 <- data.frame(genotype = nimbleTransData$genotype,
                    infected = nimbleTransData$infected)

####################################################################################################
#### Define model in BUGS/NIMBLE language

code <- nimbleCode({
  mu_alpha ~ dnorm(0, 0.001)
  sigma_alpha ~ dunif(0, 1000)
  for(i in 1:nplant) { 
    alpha[i] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
  }
  for(i in 1:7) {
    beta[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:2) {
    betagenotype[i] ~ dnorm(0, 0.001)
  }
  for(i in 1:N) {
    # Source plant sub-model
    log(lambda_plant[i]) <- alpha[plantID[i]] + betagenotype[genotype[i]] + beta[1]*distance[i]
    xf_source_plant[i] ~ dpois(lambda_plant[i])
    # Vector sub-model with detection probabilty
    log(lambda_vector[i]) <- beta[2] + beta[3]*lamba_plant[i] # Biological process
    logit(p_vector[i]) <- beta[4] + beta[5]*p_trans[i] # Detection probability
    N_vector[i] ~ dpois(lambda_vector[i]) # True Xf population in vectors
    xf_vector[i] ~ dbin(p_vector[i], N_vector[i]) # Data on Xf population in vectors
    # Transmission probability sub-model
    logit(p_trans[i]) <- beta[6] + beta[7]*N_vector[i]
    infected[i] ~ dbern(p_trans[i])
  }
})

constants <- with(nimbleTransData,
                  list(N=N, 
                       nplant=nplant, 
                       genotype = genotype, 
                       distance = distance,  
                       #genotype_distance = genotype_distance, 
                       # xf_source_plant = xf_source_plant,
                       # xf_vector = xf_vector,
                       # infected = infected,
                       plantID = plantID))

data <- with(nimbleTransData, 
             list(xf_source_plant = xf_source_plant,
                  xf_vector = xf_vector,
                  infected=infected))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,nimbleTransData$nplant), beta=rep(0,7), betagenotype=rep(0,2))

modelInfo <- list(code=code, constants=constants, data=data, inits=inits, name='month')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo$code,
                      modelInfo$constants,
                      modelInfo$data,
                      modelInfo$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:9]')
spec$addSampler('beta[1:3]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[4:9]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors(c('p_occ')) # add a monitor to get p_occ in output
#spec$addMonitors(c('p_obs')) # add a monitor to get p_obs in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 1500
burnin <- 500

ti <- Sys.time()
samplesList <- lapply(1, mcmcClusterFunction)
tf <- Sys.time()

# The time it took to run MCMC
tf-ti

save(samplesList, file = 'output/MCMC_list_climate_transmission.RData')
