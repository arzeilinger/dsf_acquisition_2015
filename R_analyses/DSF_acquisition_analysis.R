#### Preliminary look at DSF-BGSS acquisition experiment data

rm(list = ls())
my.packages <- c("lattice", "tidyr", "ggplot2", "gplots", "lme4", 
                 "dplyr", "multcomp", "optimx", "bbmle", "lmerTest")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/standardize.R")

#### Acquisition data set
acqdat <- read.csv("data/DSF_acquisition_data.csv")
summary(acqdat)
# convert to centimeters
acqdat$distance <- acqdat$distance*2.54
acqdat$sample <- acqdat %>% with(., paste(genotype, rep, cage, sep=""))

#### qPCR data set 
pcrResults <- readRDS("output/dsf_acquisition_qpcr_cfu_linreg_results.rds")
pcrResults$sample <- gsub(pattern = " ", replacement = "", x = pcrResults$sample)
table(pcrResults$sample)
# Average the replicates for each sample
pcrMean <- pcrResults %>% group_by(sample) %>% summarise(meancfu = mean(cfu), sdcfu = sd(cfu)) %>% as.data.frame()

#### combine acquisition data and qpcr data by sample (genotype, cage, rep combination)
acqdat <- acqdat %>% full_join(., pcrMean, by = "sample")

#### standardize continuous covariates
acqdat$log.source.plant.pop <- log10(acqdat$source.plant.pop + 1)
acqdat$log.meancfu <- log10(acqdat$meancfu + 1)
covars <- c("distance", "aap", "log.source.plant.pop", "log.meancfu")
covars.i <- as.numeric(sapply(covars, function(x) which(names(acqdat) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(acqdat)[i]
  stdname.i <- paste("std", var.i, sep = ".")
  stdvar.i <- standardize(acqdat[,var.i])
  acqdat[,stdname.i] <- stdvar.i
}

######################################################################
#### Transmission (test plant infection) results
acqdat <- acqdat %>% dplyr::filter(., !is.na(test.plant.infection) & !is.na(distance) & !is.na(meancfu))
# include a plant term for the random effect
acqdat$plant <- paste(acqdat$genotype, acqdat$rep, sep = "") %>% factor()
# Remove cages placed below the point of inoculation
acqdat <- acqdat[acqdat$distance >= 0,]
str(acqdat)
# Remove distance outlier
acqdat <- acqdat[acqdat$distance < 100,]
#acqdat$test.plant.infection <- as.numeric(levels(acqdat$test.plant.infection))[acqdat$test.plant.infection]

#### Model selection to determine if models with cage or distance fit better
dsfModCage <- glmer(test.plant.infection ~ inoculation.date + genotype*cage*std.log.meancfu + (1|plant), 
                    data = acqdat, family = "binomial",
                    control = glmerControl(optimizer = "bobyqa"))
dsfModDistance <- glmer(test.plant.infection ~ inoculation.date + genotype*std.distance*std.log.meancfu + (1|plant), 
                        data = acqdat, family = "binomial",
                        control = glmerControl(optimizer = "bobyqa"))
AICctab(dsfModCage, dsfModDistance, base = TRUE)
# No difference between models with cage and distance effects

#### GLMM with cage nested in plant 
dsfMod2 <- glmer(test.plant.infection ~ genotype*std.distance + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "bobyqa"))
dsfMod3 <- glmer(test.plant.infection ~ genotype*std.log.meancfu + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "bobyqa"))
dsfMod.G <- glmer(test.plant.infection ~ genotype + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "bobyqa"))

# Compare reduced models to both full model with distance and with cage
AICctab(dsfModDistance, dsfModCage, dsfMod2, dsfMod3, dsfMod.G, base = TRUE)
# dsfMod3 seems best -- no distance/cage effect
plot(dsfMod3)
summary(dsfMod3)

#### Best model with and without interactions
dsfMod3.noInterxn <- glmer(test.plant.infection ~ genotype + std.log.meancfu + (1|plant), 
                                      data = acqdat, family = "binomial",
                                      control = glmerControl(optimizer = "bobyqa"))
AICctab(dsfMod3, dsfMod3.noInterxn, base = TRUE)
summary(dsfMod3.noInterxn)


#### Model with combined genotype and cage, for contrast. But contrast doesn't work
acqdat$genotype.cage <- with(acqdat, factor(paste(genotype, cage, sep = ".")))
dsfMod3.cage <- glmer(test.plant.infection ~ genotype.cage*std.log.meancfu + (1|plant), 
                      data = acqdat, family = "binomial",
                      control = glmerControl(optimizer = "bobyqa"))
summary(dsfMod3.cage)


#### Contrasts for across cage main effect
cageContrast <- rbind("FTP - FWP" = c(0,1,0,-1),
                      "FTD - FWD" = c(1,0,-1,0))
cageContrTest <- glht(dsfMod3.cage, linfct = cageContrast)
summary(cageContrTest)
# Doesn't work. Need to figure out how to use glht with ANCOVA


#############################################################################################################
#### Include source plant xf pop
# Have to drop inoculation.date because xf pops weren't measured on first date

xfpopdat <- acqdat %>% dplyr::filter(., !is.na(source.plant.pop))

# Opimizer engine to use
optimizer <- "bobyqa"

xfpopModFull <- glmer(test.plant.infection ~ genotype*std.distance*std.log.source.plant.pop*std.log.meancfu + (1|plant), 
                      data = xfpopdat, family = "binomial",
                      control = glmerControl(optimizer = "Nelder_Mead"))
                      # control = glmerControl(optimizer = "optimx",
                      #                        optCtrl = list(method = optimizer)))
# Model without distance
xfpopMod2 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop*std.log.meancfu + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
                   # control = glmerControl(optimizer = "optimx",
                   #                        optCtrl = list(method = optimizer)))
# Model without vector xf pops
xfpopMod3 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop*std.distance + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
                  # control = glmerControl(optimizer = "optimx",
                  #                        optCtrl = list(method = optimizer)))
# Model without distance and cfu
xfpopMod4 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
                   # control = glmerControl(optimizer = "optimx",
                   #                        optCtrl = list(method = optimizer)))
# Model with only genotype
xfpopMod.G <- glmer(test.plant.infection ~ genotype + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
                   # control = glmerControl(optimizer = "optimx",
                   #                        optCtrl = list(method = optimizer)))

AICctab(xfpopModFull, xfpopMod2, xfpopMod3, xfpopMod4, xfpopMod.G, base = TRUE)
# xfpopMod2 and xfpopMod4 are equivalent but xfpopMod2 has convergence issues, go with xfpopMod4
xfpopMod4 %>% summary()

#### Model selection when dropping interactions
xfpopMod2.noIntrxn <- glmer(test.plant.infection ~ genotype + std.log.source.plant.pop + std.log.meancfu + (1|plant), 
                            data = xfpopdat, family = "binomial",
                            control = glmerControl(optimizer = "Nelder_Mead"))
                            # control = glmerControl(optimizer = "optimx",
                            #                        optCtrl = list(method = optimizer)))
AICctab(xfpopMod2, xfpopMod2.noIntrxn, base = TRUE)

#### Assessing sample size
plants <- as.data.frame(table(acqdat$rep, acqdat$genotype))
plants.incl <- plants[plants$Freq > 0,]
plants.incl %>% group_by(Var2) %>% summarise(n = length(Freq)) # 25 FT plants and 23 FW plants
table(acqdat$cage, acqdat$genotype) # sample size for each cage*genotype

#### Chi square test
table(acqdat$genotype, acqdat$test.plant.infection)
chisq.test(acqdat$genotype, acqdat$test.plant.infection)


#### Predicted (marginal) means by genotype
# model with cage
acqdat$predTrans <- predict(dsfMod3, type = "response", re.form = NA)
acqdat %>% group_by(genotype) %>% summarise(mean = mean(predTrans), se = sd(predTrans)/sqrt(length(predTrans)))
# model with distance
# acqdat$predTrans <- predict(dsfMod2, type = "response", re.form = NA)
# acqdat %>% group_by(genotype) %>% summarise(mean = mean(predTrans), se = sd(predTrans)/sqrt(length(predTrans)))

########################################################################################
#### Figures for transmission results
## Pecent transmission between genotypes
transPerc <- acqdat %>% group_by(genotype, cage) %>% summarise(perc = (sum(test.plant.infection)/length(test.plant.infection))*100)
transPerc$trt <- paste(transPerc$genotype, transPerc$cage, sep="-")
transPerc$trt.full <- c("DSF-Distal", "DSF-Proximal", "WT-Distal", "WT-Proximal")

transBarplot <- ggplot(transPerc, aes(x=trt.full,y=perc)) +
  geom_bar(position=position_dodge(), stat = "identity",
           fill = "black",
           size = 0.3) +
  ylab("% Transmission") + 
  ylim(c(0,75)) +
  xlab("") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) 
ggsave("results/transmission_barplot.jpg", plot = transBarplot,
       width = 5, height = 4, units = "in")  


## Distance and genotype
## Since distance is not significant, don't need any model fit lines
#xdist <- runif(nrow(acqdat), min(acqdat$distance, na.rm=TRUE), max(acqdat$distance, na.rm=TRUE))
acqdat$trans.dummy <- ifelse(acqdat$genotype == "FW", acqdat$test.plant.infection,
                             ifelse(acqdat$test.plant.infection == 0, 0.05, 0.95))

tiff("results/DSF_transmission_distance_plot_monochrome.tif")
  plot(acqdat[acqdat$genotype == "FW",]$distance, jitter(acqdat[acqdat$genotype == "FW",]$trans.dummy, amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), xlim = c(0,100), pch = 1, col = "black",
       ylab = "Probability of transmission", xlab = "Distance from inoculation point (cm)")
  points(acqdat[acqdat$genotype == "FT",]$distance, jitter(acqdat[acqdat$genotype == "FT",]$trans.dummy, amount = 0), 
         pch = 2, col = "black")
  # lines(smooth.spline(acqdat[acqdat$genotype == "FW",]$distance, acqdat[acqdat$genotype == "FW",]$predTrans, tol = 1e-6), 
  #       lty = 1, lwd = 2, col = "black")
  # lines(smooth.spline(acqdat[acqdat$genotype == "FT",]$distance, acqdat[acqdat$genotype == "FT",]$predTrans, tol = 1e-6), 
  #       lty = 2, lwd = 2, col = "black")
dev.off()

#plot(acqdat$distance, acqdat$predTrans, pch = 2)


#### Plot of vector xf pops vs transmission
tiff("results/DSF_transmission_vector_xfpops_plot_monochrome.tif")
  plot(acqdat[acqdat$genotype == "FW",]$log.meancfu, jitter(acqdat[acqdat$genotype == "FW",]$trans.dummy, amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), 
       #xlim = c(0,120), 
       pch = 1, col = "black",
       ylab = "Probability of transmission", xlab = "Xylella populations in vectors (log10 transformed)")
  points(acqdat[acqdat$genotype == "FT",]$log.meancfu, jitter(acqdat[acqdat$genotype == "FT",]$trans.dummy, amount = 0), 
         pch = 2, col = "black")
  lines(smooth.spline(acqdat[acqdat$genotype == "FW",]$log.meancfu, acqdat[acqdat$genotype == "FW",]$predTrans, tol = 1e-6, nknots = 4), 
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(acqdat[acqdat$genotype == "FT",]$log.meancfu, acqdat[acqdat$genotype == "FT",]$predTrans, tol = 1e-6, nknots = 4), 
        lty = 2, lwd = 2, col = "black")
dev.off()


#### Plot of Xf source plant pops vs transmission
xfpopdat$predTrans2 <- predict(xfpopMod2, type = "response", re.form = NA)

tiff("results/DSF_transmission_sourcepop_plot_monochrome.tif")
  plot(jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, amount = 0), jitter(xfpopdat[xfpopdat$genotype == "FW",]$trans.dummy, amount = 0),
       cex.axis = 1.3, cex.lab = 1, cex = 1.3,
       ylim = c(0,1), xlim = c(0,10), pch = 1, col = "black",
       ylab = "Probability of transmission", xlab = "X. fastidiosa populations in source plant (CFU/g, log10 transformed)")
  points(jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, amount = 0), jitter(xfpopdat[xfpopdat$genotype == "FT",]$trans.dummy, amount = 0), 
         pch = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FW",]$predTrans2, nknots = 4, tol = 1e-10), 
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FT",]$predTrans2, nknots = 4, tol = 1e-10), 
        lty = 2, lwd = 2, col = "black")
dev.off()



#### Tri-variate plot with distance, Xf pop, and P(trans)
library(akima)
# WT contourplot
fwdat <- xfpopdat[xfpopdat$genotype == "FW",]
zzFW <- interp(x=fwdat$distance, y=log(fwdat$source.plant.pop+1), 
             fwdat$predTransPop, duplicate = "mean")
tiff("contourplot_WT.tif")
  filled.contour(zzFW, col = topo.colors(24), 
                 xlab = "Distance (cm)", ylab = "Xf pop in source plant (log)")
dev.off()
# DSF contourplot
ftdat <- xfpopdat[xfpopdat$genotype == "FT",]
zzft <- interp(x=ftdat$distance, y=log(ftdat$source.plant.pop+1), 
             ftdat$predTransPop, duplicate = "mean")
tiff("contourplot_DSF.tif")
  filled.contour(zzft, col = topo.colors(24), 
                 xlab = "Distance (cm)", ylab = "Xf pop in source plant (log)")
dev.off()


#### Tri-variate plot with vector Xf pops, source Xf pop, and P(trans)
library(akima)
xfpopdat$predTransPop <- predict(xfpopMod2, type = "response", re.form = NA)

# WT contourplot
fwdat <- xfpopdat[xfpopdat$genotype == "FW",]
zzFW <- interp(x=log10(fwdat$meancfu+1), y=log10(fwdat$source.plant.pop+1), 
               fwdat$predTransPop, duplicate = "mean")
tiff("results/contourplot_WT.tif")
filled.contour(zzFW, col = topo.colors(24), 
               xlab = "Xf pop in vector (log10)", ylab = "Xf pop in source plant (log10)")
dev.off()
# DSF contourplot
ftdat <- xfpopdat[xfpopdat$genotype == "FT",]
zzft <- interp(x=log10(ftdat$meancfu+1), y=log10(ftdat$source.plant.pop+1), 
               ftdat$predTransPop, duplicate = "mean")
tiff("results/contourplot_DSF.tif")
filled.contour(zzft, col = topo.colors(24), 
               xlab = "Xf pop in vector (log10)", ylab = "Xf pop in source plant (log10)")
dev.off()


# Plot scatterplot of source xf pop and vector xf pop for WT plants
tiff("results/source_vs_vector_scatterplots.tif")
  plot(x = log10(ftdat$meancfu+1), y=log10(ftdat$source.plant.pop+1))
dev.off()



###########################################################################
#### Analysis of Source plant Xf populations

# # Number of infections, by genotype and plant position
# infdat <- ppdat[ppdat$source.infected == 1,]
# table(infdat$Genotype, infdat$Cage)
# chisq.test(infdat$Genotype, infdat$Cage, simulate.p.value = TRUE)

# Xf population size in source plant by genotype and distance of all plants
sourcepopMod <- lmer(log.source.plant.pop ~ genotype*std.distance + (1|plant), 
                      data = xfpopdat,
                      control = lmerControl(optimizer = "bobyqa"))
                      # control = glmerControl(optimizer = "optimx",
                      #                        optCtrl = list(method = "bobyqa")))

plot(sourcepopMod)
summary(sourcepopMod)

xfpopdat$predSourcePop <- predict(sourcepopMod, type = "response", re.form = NA)
# source plant populations and distance from inoculation
tiff("results/source_plant_pop_distance_plot.tif")
  plot(x = jitter(xfpopdat[xfpopdat$genotype == "FW",]$distance, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, amount = 0),
       cex.axis = 1.3, cex.lab = 1, cex = 1.3,
       #ylim = c(0,1), xlim = c(0,10), 
       pch = 1, col = "black",
       xlab = "Distance from inoculation (cm)", ylab = "X. fastidiosa populations in source plant (CFU/g, log10 transformed)")
  points(x = jitter(xfpopdat[xfpopdat$genotype == "FT",]$distance, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, amount = 0), 
         pch = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FW",]$distance, xfpopdat[xfpopdat$genotype == "FW",]$predSourcePop, nknots = 4, tol = 1e-10), 
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FT",]$distance, xfpopdat[xfpopdat$genotype == "FT",]$predSourcePop, nknots = 4, tol = 1e-10), 
        lty = 2, lwd = 2, col = "black")
dev.off()


############################################################################################################
#### Analysis of Xf populations in vectors

xfvectorMod1 <- lmer(log.meancfu ~ genotype*distance*log.source.plant.pop + (1|plant),
                    data = xfpopdat,
                    control = lmerControl(optimizer = "bobyqa"))
xfvectorMod2 <- lmer(log.meancfu ~ genotype*log.source.plant.pop + (1|plant),
                     data = xfpopdat,
                     control = lmerControl(optimizer = "bobyqa"))
AICctab(xfvectorMod1, xfvectorMod2, base = TRUE)
# Dropping distance is clearly better
plot(xfvectorMod2)
summary(xfvectorMod2)

# Plotting
xfpopdat$predVectorPop <- predict(xfvectorMod, type = "response", re.form = NA)
# Xf populations in vectors and in source plants from inoculation
tiff("results/xf_vector_source_plant_plot.tif")
  plot(x = jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.meancfu, amount = 0),
       cex.axis = 1.3, cex.lab = 1, cex = 1.3,
       #ylim = c(0,1), xlim = c(0,10), 
       pch = 1, col = "black",
       ylab = "X. fastidiosa populations in vectors (CFU, log10 transformed)", xlab = "X. fastidiosa populations in source plant (CFU/g, log10 transformed)")
  points(x = jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.meancfu, amount = 0), 
         pch = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FW",]$predVectorPop, nknots = 4, tol = 1e-10), 
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FT",]$predVectorPop, nknots = 4, tol = 1e-10), 
        lty = 2, lwd = 2, col = "black")
dev.off()



#####################################################################################################################
#### False negative vectors

falseNegatives <- acqdat %>% dplyr::filter(., meancfu == 0 & test.plant.infection == 1)
infected <- acqdat %>% dplyr::filter(., meancfu > 0 & test.plant.infection == 1)
outsamples <- rbind(falseNegatives[1:6,c("meancfu", "sample")], infected[1:6,c("meancfu", "sample")])
write.csv(outsamples, file = "output/vector_xf_false_negatives.csv", row.names = FALSE)



#####################################################################################################################
#### Analyses without zeros
#####################################################################################################################
acqdatNZ <- acqdat %>% dplyr::filter(., meancfu > 0)

#### Transmission analysis without source plant populations
#### GLMM with cage nested in plant 
dsfModDistance <- glmer(test.plant.infection ~ inoculation.date + genotype*std.distance*std.log.meancfu + (1|plant), 
                        data = acqdat, family = "binomial",
                        control = glmerControl(optimizer = "bobyqa"))
dsfMod2 <- glmer(test.plant.infection ~ genotype*std.distance + (1|plant), 
                 data = acqdat, family = "binomial",
                 control = glmerControl(optimizer = "bobyqa"))
dsfMod3 <- glmer(test.plant.infection ~ genotype*std.log.meancfu + (1|plant), 
                 data = acqdat, family = "binomial",
                 control = glmerControl(optimizer = "bobyqa"))
dsfMod.G <- glmer(test.plant.infection ~ genotype + (1|plant), 
                  data = acqdat, family = "binomial",
                  control = glmerControl(optimizer = "bobyqa"))

# Compare reduced models to both full model with distance and with cage
AICctab(dsfModDistance, dsfModCage, dsfMod2, dsfMod3, dsfMod.G, base = TRUE)
# dsfMod3 seems best -- no distance/cage effect
plot(dsfMod3)
summary(dsfMod3)

#### Plot of vector xf pops vs transmission excluding zeros
tiff("results/DSF_transmission_vector_xfpops_plot_non-zero.tif")
  plot(acqdatNZ[acqdatNZ$genotype == "FW",]$log.meancfu, jitter(acqdatNZ[acqdatNZ$genotype == "FW",]$trans.dummy, amount = 0),
       cex.axis = 1.3, cex.lab = 1.3, cex = 1.3,
       ylim = c(0,1), 
       #xlim = c(0,120), 
       pch = 1, col = "black",
       ylab = "Probability of transmission", xlab = "Xylella populations in vectors (log10 transformed)")
  points(acqdatNZ[acqdatNZ$genotype == "FT",]$log.meancfu, jitter(acqdatNZ[acqdatNZ$genotype == "FT",]$trans.dummy, amount = 0), 
       pch = 2, col = "black")
  lines(smooth.spline(acqdatNZ[acqdatNZ$genotype == "FW",]$log.meancfu, acqdatNZ[acqdatNZ$genotype == "FW",]$predTrans, tol = 1e-6, nknots = 4), 
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(acqdatNZ[acqdatNZ$genotype == "FT",]$log.meancfu, acqdatNZ[acqdatNZ$genotype == "FT",]$predTrans, tol = 1e-6, nknots = 4), 
        lty = 2, lwd = 2, col = "black")
dev.off()



#### Analysis of Xf populations in vectors
xfpopdatNZ <- acqdatNZ %>% dplyr::filter(., !is.na(source.plant.pop) & source.plant.pop > 0)

xfvectorMod1 <- lmer(log.meancfu ~ genotype*distance*log.source.plant.pop + (1|plant),
                     data = xfpopdatNZ,
                     control = lmerControl(optimizer = "bobyqa"))
xfvectorMod2 <- lmer(log.meancfu ~ genotype*log.source.plant.pop + (1|plant),
                     data = xfpopdatNZ,
                     control = lmerControl(optimizer = "bobyqa"))
AICctab(xfvectorMod1, xfvectorMod2, base = TRUE)
# Dropping distance is clearly better
plot(xfvectorMod2)
summary(xfvectorMod2)

# Plotting -- model fit lines not necessary as relationship is not significant
#xfpopdatNZ$predVectorPop <- predict(xfvectorMod, type = "response", re.form = NA)
# Xf populations in vectors and in source plants from inoculation
tiff("results/xf_vector_source_plant_plot_non-zero.tif")
  plot(x = jitter(xfpopdatNZ[xfpopdatNZ$genotype == "FW",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdatNZ[xfpopdatNZ$genotype == "FW",]$log.meancfu, amount = 0),
       cex.axis = 1.3, cex.lab = 1, cex = 1.3,
       #ylim = c(0,1), xlim = c(0,10), 
       pch = 1, col = "black",
       ylab = "X. fastidiosa populations in vectors (CFU, log10 transformed)", xlab = "X. fastidiosa populations in source plant (CFU/g, log10 transformed)")
  points(x = jitter(xfpopdatNZ[xfpopdatNZ$genotype == "FT",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdatNZ[xfpopdatNZ$genotype == "FT",]$log.meancfu, amount = 0), 
         pch = 2, col = "black")
  # lines(smooth.spline(xfpopdatNZ[xfpopdatNZ$genotype == "FW",]$log.source.plant.pop, xfpopdatNZ[xfpopdatNZ$genotype == "FW",]$predVectorPop, nknots = 4, tol = 1e-10), 
#       lty = 1, lwd = 2, col = "black")
# lines(smooth.spline(xfpopdatNZ[xfpopdatNZ$genotype == "FT",]$log.source.plant.pop, xfpopdatNZ[xfpopdatNZ$genotype == "FT",]$predVectorPop, nknots = 4, tol = 1e-10), 
#       lty = 2, lwd = 2, col = "black")
dev.off()
