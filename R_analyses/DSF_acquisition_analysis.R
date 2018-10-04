#### Preliminary look at DSF-BGSS acquisition experiment data

rm(list = ls())
my.packages <- c("lattice", "tidyr", "ggplot2", "gplots", "lme4", "pROC",
                 "dplyr", "multcomp", "optimx", "bbmle", "lmerTest")
lapply(my.packages, require, character.only = TRUE)

source("R_functions/standardize.R")

#### Load data sets that have been cleaned and with model predictions
acqdat <- readRDS("output/dsf_acquisition_data_set.rds")
xfpopdat <- readRDS("output/dsf_acquisition_data_set_complete.rds")



#################################################################################################
#### Code for cleaning up raw data set and combining with qPCR data
#### Skip this chunk if you just loaded "acqdat" and "xfpopdat" above

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
covars <- c("distance", "aap", "source.plant.pop", "log.source.plant.pop", "log.meancfu")
covars.i <- as.numeric(sapply(covars, function(x) which(names(acqdat) == x), simplify = TRUE))
for(i in covars.i){
  var.i <- names(acqdat)[i]
  stdname.i <- paste("std", var.i, sep = ".")
  stdvar.i <- standardize(acqdat[,var.i])
  acqdat[,stdname.i] <- stdvar.i
}

# Have to drop inoculation.date because xf pops weren't measured on first date
xfpopdat <- acqdat %>% dplyr::filter(., !is.na(source.plant.pop))

###########################################################################################################
#### Transmission (test plant infection) results
###########################################################################################################

acqdat <- acqdat %>% dplyr::filter(., !is.na(test.plant.infection) & !is.na(distance) & !is.na(meancfu))
# include a plant term for the random effect
acqdat$plant <- paste(acqdat$genotype, acqdat$rep, sep = "") %>% factor()
# Remove cages placed below the point of inoculation
acqdat <- acqdat[acqdat$distance >= 0,]
str(acqdat)
# Remove distance outlier
acqdat <- acqdat[acqdat$distance < 100,]

#### Model selection to determine if models with cage or distance fit better
dsfModCage <- glmer(test.plant.infection ~ inoculation.date + genotype*cage*std.log.meancfu + (1|plant), 
                    data = acqdat, family = "binomial",
                    control = glmerControl(optimizer = "Nelder_Mead"))
dsfModDistance <- glmer(test.plant.infection ~ inoculation.date + genotype*std.distance*std.log.meancfu + (1|plant), 
                        data = acqdat, family = "binomial",
                        control = glmerControl(optimizer = "Nelder_Mead"))
AICctab(dsfModCage, dsfModDistance, base = TRUE)
# No difference between models with cage and distance effects

#### GLMM with plant as random effect and dropping inoculation.date
dsfModFull <- glmer(test.plant.infection ~ genotype*std.distance*std.log.meancfu + (1|plant), 
                    data = acqdat, family = "binomial",
                    control = glmerControl(optimizer = "Nelder_Mead"))
dsfMod2 <- glmer(test.plant.infection ~ genotype*std.distance + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "Nelder_Mead"))
dsfMod3 <- glmer(test.plant.infection ~ genotype*std.log.meancfu + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "Nelder_Mead"))
dsfMod.G <- glmer(test.plant.infection ~ genotype + (1|plant), 
                data = acqdat, family = "binomial",
                control = glmerControl(optimizer = "Nelder_Mead"))

# Compare reduced models to both full model with distance and with cage
AICctab(dsfModFull, dsfMod2, dsfMod3, dsfMod.G, base = TRUE)
# dsfMod3 seems best -- no distance/cage effect
plot(dsfMod3)
summary(dsfMod3)

#### Best model with and without interactions
dsfMod3.noInterxn <- glmer(test.plant.infection ~ genotype + std.log.meancfu + (1|plant), 
                                      data = acqdat, family = "binomial",
                                      control = glmerControl(optimizer = "Nelder_Mead"))
AICctab(dsfModFull, dsfMod2, dsfMod3, dsfMod.G, dsfMod3.noInterxn, base = TRUE)
summary(dsfMod3.noInterxn) # Best model is without interaction


# ROC curve
p <- as.numeric(predict(dsfMod3.noInterxn, type = "response")>0.5)
mean(p==acqdat$test.plant.infection)
table(p,acqdat$test.plant.infection)

transVectorROC <- roc(response = acqdat$test.plant.infection, predictor = p)
tiff("results/trans_vector_roc_curve.tiff", compression = "lzw")
  plot(transVectorROC, ylim = c(0,1))
dev.off()



#############################################################################################################
#### Analysis of transmission inluding source plant xf pop
#############################################################################################################

# Note on optimizers: tested bobyqa extensively and Nelder_Mead had fewer convergence issues. Stick with Nelder_Mead

xfpopModFull <- glmer(test.plant.infection ~ genotype*std.distance*std.log.source.plant.pop*std.log.meancfu + (1|plant), 
                      data = xfpopdat, family = "binomial",
                      control = glmerControl(optimizer = "Nelder_Mead"))
# Model without distance
xfpopMod2 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop*std.log.meancfu + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
# Model without vector xf pops
xfpopMod3 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop*std.distance + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
# Model without distance and cfu
xfpopMod4 <- glmer(test.plant.infection ~ genotype*std.log.source.plant.pop + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))
# Model with only genotype
xfpopMod.G <- glmer(test.plant.infection ~ genotype + (1|plant), 
                   data = xfpopdat, family = "binomial",
                   control = glmerControl(optimizer = "Nelder_Mead"))

AICctab(xfpopModFull, xfpopMod2, xfpopMod3, xfpopMod4, xfpopMod.G, base = TRUE)
# xfpopMod2 and xfpopMod4 are equivalent but xfpopMod2 has convergence issues, go with xfpopMod4
xfpopMod4 %>% summary()

## ROC curve
p <- as.numeric(predict(xfpopMod4, type = "response")>0.5)
mean(p==xfpopdat$test.plant.infection)
table(p,xfpopdat$test.plant.infection)

transSourceROC <- roc(response = xfpopdat$test.plant.infection, predictor = p)
tiff("results/trans_source_roc_curve.tiff", compression = "lzw")
  plot(transSourceROC, ylim = c(0,1))
dev.off()


#### Model selection when dropping interactions
xfpopMod4.noIntrxn <- glmer(test.plant.infection ~ genotype + std.log.source.plant.pop + (1|plant), 
                            data = xfpopdat, family = "binomial",
                            control = glmerControl(optimizer = "Nelder_Mead"))
AICctab(xfpopModFull, xfpopMod2, xfpopMod3, xfpopMod4, xfpopMod.G, xfpopMod4.noIntrxn, base = TRUE)
# Best model includes genotype-source plant interaction

#### Assessing sample size
plants <- as.data.frame(table(acqdat$rep, acqdat$genotype))
plants.incl <- plants[plants$Freq > 0,]
plants.incl %>% group_by(Var2) %>% summarise(n = length(Freq)) # 25 FT plants and 23 FW plants
table(acqdat$cage, acqdat$genotype) # sample size for each cage*genotype


#### Predicted (marginal) means by genotype
acqdat$predTrans <- predict(dsfMod3.noInterxn, type = "response", re.form = NA)
acqdat %>% group_by(genotype) %>% summarise(mean = mean(predTrans), se = sd(predTrans)/sqrt(length(predTrans)))

#### Looking at false negatives in xf vector and xf source plant pops
# In vectors
xfpopdat$vector.infection <- xfpopdat %>% with(., ifelse(meancfu > 0, 1, 0))
# In table(), first element is on the side, second element is on the top
table(xfpopdat$test.plant.infection, xfpopdat$vector.infection)

# In source plants
xfpopdat$source.plant.infection <- xfpopdat %>% with(., ifelse(source.plant.pop > 0, 1, 0))
table(xfpopdat$test.plant.infection, xfpopdat$source.plant.infection)

## Summary of transmission for reduced data set
transPercReduced <- xfpopdat %>% group_by(genotype) %>% summarise(perc = (sum(test.plant.infection)/length(test.plant.infection))*100,
                                                                  n = sum(!is.na(test.plant.infection)))


########################################################################################
#### Figures for transmission results
## Pecent transmission between genotypes
# transPerc <- acqdat %>% group_by(genotype, cage) %>% summarise(perc = (sum(test.plant.infection)/length(test.plant.infection))*100)
# transPerc$trt <- paste(transPerc$genotype, transPerc$cage, sep="-")
# transPerc$trt.full <- c("DSF-Distal", "DSF-Proximal", "WT-Distal", "WT-Proximal")
# 
# transBarplot <- ggplot(transPerc, aes(x=trt.full,y=perc)) +
#   geom_bar(position=position_dodge(), stat = "identity",
#            fill = "black",
#            size = 0.3) +
#   ylab("% Transmission") + 
#   ylim(c(0,75)) +
#   xlab("") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank()) 
# 
# transBarplot
# ggsave("results/transmission_barplot.jpg", plot = transBarplot,
#        width = 5, height = 4, units = "in")  
# 
# ## Plot transmission and just genotype
# transPerc2 <- acqdat %>% group_by(genotype) %>% summarise(perc = (sum(test.plant.infection)/length(test.plant.infection))*100,
#                                                           n = sum(!is.na(test.plant.infection)))
# transPerc2$trt <- ifelse(transPerc2$genotype == "FT", "DSF", "WT")
# 
# transBarplot2 <- ggplot(transPerc2, aes(x=trt,y=perc)) +
#   geom_bar(position=position_dodge(), stat = "identity",
#            fill = "black",
#            size = 1.4) +
#   ylab("% test plants infected") + 
#   ylim(c(0,75)) +
#   xlab("Source plant genotype") +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_text(size=20),
#         axis.title=element_text(size=22)) 
# 
# transBarplot2
# ggsave("results/transmission_barplot_just_genotype.jpg", plot = transBarplot2,
#        width = 7, height = 7, units = "in")  
# 
# 
# 
# ## Distance and genotype
# ## Since distance is not significant, don't need any model fit lines
# #xdist <- runif(nrow(acqdat), min(acqdat$distance, na.rm=TRUE), max(acqdat$distance, na.rm=TRUE))
# acqdat$trans.dummy <- ifelse(acqdat$genotype == "FW", acqdat$test.plant.infection,
#                              ifelse(acqdat$test.plant.infection == 0, 0.04, 0.96))
# 
# tiff("results/DSF_transmission_distance_plot_monochrome.tif")
#   plot(acqdat[acqdat$genotype == "FW",]$distance, jitter(acqdat[acqdat$genotype == "FW",]$trans.dummy, amount = 0),
#        cex.axis = 1.4, cex.lab = 1.6, cex = 1.6,
#        ylim = c(0,1), xlim = c(0,100), pch = 1, col = "black",
#        ylab = "Probability of transmission", xlab = "Distance from inoculation point (cm)")
#   points(acqdat[acqdat$genotype == "FT",]$distance, jitter(acqdat[acqdat$genotype == "FT",]$trans.dummy, amount = 0), 
#          pch = 2, col = "black")
#   # lines(smooth.spline(acqdat[acqdat$genotype == "FW",]$distance, acqdat[acqdat$genotype == "FW",]$predTrans, tol = 1e-6), 
#   #       lty = 1, lwd = 2, col = "black")
#   # lines(smooth.spline(acqdat[acqdat$genotype == "FT",]$distance, acqdat[acqdat$genotype == "FT",]$predTrans, tol = 1e-6), 
#   #       lty = 2, lwd = 2, col = "black")
# dev.off()
# 
# plot(acqdat$distance, acqdat$predTrans, pch = 2)
# 
# 
# #### Plot of vector xf pops vs transmission
# tiff("results/DSF_transmission_vector_xfpops_plot_monochrome.tif")
#   plot(acqdat[acqdat$genotype == "FW",]$log.meancfu, jitter(acqdat[acqdat$genotype == "FW",]$trans.dummy, amount = 0),
#        cex.axis = 1.4, cex.lab = 1.6, cex = 1.6,
#        ylim = c(0,1), 
#        #xlim = c(0,120), 
#        pch = 1, col = "black",
#        ylab = "Probability of transmission", xlab = "Xylella populations in vectors (CFU, log10 transformed)")
#   points(acqdat[acqdat$genotype == "FT",]$log.meancfu, jitter(acqdat[acqdat$genotype == "FT",]$trans.dummy, amount = 0), 
#          pch = 2, col = "black")
#   lines(smooth.spline(acqdat[acqdat$genotype == "FW",]$log.meancfu, acqdat[acqdat$genotype == "FW",]$predTrans, tol = 1e-6, nknots = 4), 
#         lty = 1, lwd = 2, col = "black")
#   lines(smooth.spline(acqdat[acqdat$genotype == "FT",]$log.meancfu, acqdat[acqdat$genotype == "FT",]$predTrans, tol = 1e-6, nknots = 4), 
#         lty = 2, lwd = 2, col = "black")
# dev.off()

## In ggplot
vectorTransPlot <- ggplot(data = acqdat) +
  geom_point(data = acqdat[acqdat$genotype == "FW",], aes(y = jitter(test.plant.infection, amount = 0.05), x = log.meancfu), size = 4, pch = 16) +
  geom_point(data = acqdat[acqdat$genotype == "FT",], aes(y = jitter(test.plant.infection, amount = 0.05), x = log.meancfu), size = 4, pch = 1) +
  stat_smooth(data = acqdat[acqdat$genotype == "FW",], aes(y = predTrans, x = log.meancfu), method = "lm", se = FALSE, col = "black", lty = 1) +
  stat_smooth(data = acqdat[acqdat$genotype == "FT",], aes(y = predTrans, x = log.meancfu), method = "lm", se = FALSE, col = "black", lty = 2) +
  xlab("Xylella populations in vectors (CFU, log10)") + ylab("Probability of transmission") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
vectorTransPlot

ggsave(filename = "results/DSF_transmission_vector_xfpops_ggplot.tiff",
       plot = vectorTransPlot,
       width = 7, height = 7, units = "in")


#### Plot of Xf source plant pops vs transmission
xfpopdat$trans.dummy <- ifelse(xfpopdat$genotype == "FW", xfpopdat$test.plant.infection,
                               ifelse(xfpopdat$test.plant.infection == 0, 0.04, 0.96))
xfpopdat$predTrans2 <- predict(xfpopMod4, type = "response", re.form = NA)

## In ggplot
plantTransPlot <- ggplot(data = xfpopdat) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FW",], aes(y = jitter(test.plant.infection, amount = 0.05), x = log.source.plant.pop), size = 4, pch = 16) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FT",], aes(y = jitter(test.plant.infection, amount = 0.05), x = log.source.plant.pop), size = 4, pch = 1) +
  stat_smooth(data = xfpopdat[xfpopdat$genotype == "FW",], aes(y = predTrans2, x = log.source.plant.pop), method = "lm", se = FALSE, col = "black", lty = 1) +
  stat_smooth(data = xfpopdat[xfpopdat$genotype == "FT",], aes(y = predTrans2, x = log.source.plant.pop), method = "lm", se = FALSE, col = "black", lty = 2) +
  xlab("Xylella populations in source plants (CFU/g, log10)") + ylab("Probability of transmission") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
plantTransPlot

ggsave(filename = "results/DSF_transmission_sourcepop_ggplot.tiff",
       plot = plantTransPlot,
       width = 7, height = 7, units = "in")


##############################################################################################################
#### Analysis of Source plant Xf populations
##############################################################################################################

#### Without source.plant.pop zeros
xfpopdatNZ <- xfpopdat %>% dplyr::filter(., source.plant.pop > 0)

# Xf population size in source plant by genotype and distance of all plants
sourcepopModnz1 <- lmer(log.source.plant.pop ~ genotype*std.distance + (1|plant), 
                      data = xfpopdatNZ,
                      control = lmerControl(optimizer = "Nelder_Mead"))
# control = glmerControl(optimizer = "optimx",
#                        optCtrl = list(method = "bobyqa")))
sourcepopModnz2 <- lmer(log.source.plant.pop ~ genotype + std.distance + (1|plant), 
                      data = xfpopdatNZ,
                      control = lmerControl(optimizer = "Nelder_Mead"))
# control = glmerControl(optimizer = "optimx",
#                        optCtrl = list(method = "bobyqa")))
AICctab(sourcepopModnz1, sourcepopModnz2, base = TRUE)
# Best model is sourcepopModnz2

# Residuals vs. predicted plot for Appendix
sourcePlantPlotResid <-   plot(xfpopMod4, pch = 16, col = "black", xlab = "Model predictions", ylab = "Residuals")

tiff("results/source_plant_pop_distance_residuals.tiff", compression = "lzw")
  plot(sourcepopModnz2, pch = 16, col = "black", ylab = "Residuals", xlab = "Model predictions")
dev.off()

# Model results
summary(sourcepopModnz2)
xfpopdatNZ$predSourcePop <- predict(sourcepopModnz2, type = "response", re.form = NA)


## In ggplot
sourceDistancePlot <- ggplot(data = xfpopdat) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FW",], aes(y = log.source.plant.pop, x = distance), size = 4, pch = 16) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FT",], aes(y = log.source.plant.pop, x = distance), size = 4, pch = 1) +
  stat_smooth(data = xfpopdatNZ[xfpopdatNZ$genotype == "FW",], aes(y = predSourcePop, x = distance), method = "lm", se = FALSE, col = "black", lty = 1) +
  stat_smooth(data = xfpopdatNZ[xfpopdatNZ$genotype == "FT",], aes(y = predSourcePop, x = distance), method = "lm", se = FALSE, col = "black", lty = 2) +
  xlab("Distance from inoculation point (cm)") + ylab("Xylella populations in source plants (CFU/g, log10)") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
sourceDistancePlot

ggsave(filename = "results/source_plant_pop_distance_nonzero_ggplot.tiff",
       plot = sourceDistancePlot,
       width = 7, height = 7, units = "in")


## Summary of source.plant.pop between genotypes
summaryPlantPop <- xfpopdatNZ %>% group_by(genotype) %>% summarise(mean = mean(source.plant.pop, na.rm = TRUE),
                                                                 n = sum(!is.na(source.plant.pop)),
                                                                 se = sd(source.plant.pop, na.rm = TRUE)/sqrt(n))

############################################################################################################
#### Analysis of Xf populations in vectors
############################################################################################################

xfvectorMod1 <- lmer(log.meancfu ~ genotype*std.distance*std.log.source.plant.pop + (1|plant),
                    data = xfpopdat,
                    control = lmerControl(optimizer = "Nelder_Mead"))
xfvectorMod2 <- lmer(log.meancfu ~ genotype*std.log.source.plant.pop + (1|plant),
                     data = xfpopdat,
                     control = lmerControl(optimizer = "Nelder_Mead"))
xfvectorMod3 <- lmer(log.meancfu ~ genotype + std.log.source.plant.pop + (1|plant),
                     data = xfpopdat,
                     control = lmerControl(optimizer = "Nelder_Mead"))
AICctab(xfvectorMod1, xfvectorMod2, xfvectorMod3, base = TRUE)
# Dropping distance is clearly better, dropping interaction too

# Plot residuals and predictions
xfvectorPlotResid <-   plot(xfvectorMod3, pch = 16, col = "black", ylab = "Residuals", xlab = "Model predictions")

tiff("results/xf_vector_residuals.tif", compression = "lzw")
  plot(xfvectorMod3, pch = 16, col = "black", ylab = "Residuals", xlab = "Model predictions")
dev.off()

summary(xfvectorMod3)


# Plotting
xfpopdat$predVectorPop <- predict(xfvectorMod3, type = "response", re.form = NA)
# Xf populations in vectors and in source plants from inoculation
tiff("results/xf_vector_source_plant_plot.tif")
  plot(x = jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FW",]$log.meancfu, amount = 0),
       cex.axis = 1.4, cex.lab = 1.2, cex = 1.8,
       #ylim = c(0,1), xlim = c(0,10), 
       pch = 1, col = "black",
       ylab = "Xylella populations in vectors (CFU, log10)", xlab = "Xylella populations in source plant (CFU/g, log10)")
  points(x = jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, amount = 0), y = jitter(xfpopdat[xfpopdat$genotype == "FT",]$log.meancfu, amount = 0), 
         pch = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FW",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FW",]$predVectorPop, nknots = 4, tol = 1e-10),
        lty = 1, lwd = 2, col = "black")
  lines(smooth.spline(xfpopdat[xfpopdat$genotype == "FT",]$log.source.plant.pop, xfpopdat[xfpopdat$genotype == "FT",]$predVectorPop, nknots = 4, tol = 1e-10),
        lty = 2, lwd = 2, col = "black")
dev.off()


## In ggplot
vectorSourcePlot <- ggplot(data = xfpopdat) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FW",], aes(x = log.source.plant.pop, y = log.meancfu), size = 4, pch = 16) +
  geom_point(data = xfpopdat[xfpopdat$genotype == "FT",], aes(x = log.source.plant.pop, y = log.meancfu), size = 4, pch = 1) +
  stat_smooth(data = xfpopdat[xfpopdat$genotype == "FW",], aes(x = log.source.plant.pop, y = log.meancfu), method = "lm", se = FALSE, col = "black", lty = 1) +
  stat_smooth(data = xfpopdat[xfpopdat$genotype == "FT",], aes(x = log.source.plant.pop, y = log.meancfu), method = "lm", se = FALSE, col = "black", lty = 2) +
  ylab("Xylella populations in vectors (CFU, log10)") + xlab("Xylella populations in source plants (CFU/g, log10)") +
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
vectorSourcePlot

ggsave(filename = "results/xf_vector_source_plant_ggplot.tiff",
       plot = vectorSourcePlot,
       width = 7, height = 7, units = "in")


## Summary of vector pop between genotypes
summaryVectorPop <- xfpopdat %>% group_by(genotype) %>% summarise(mean = mean(meancfu, na.rm = TRUE),
                                                                 n = sum(!is.na(meancfu)),
                                                                 se = sd(meancfu, na.rm = TRUE)/sqrt(n))

## Save acqdat and xfpopdat data.frames
saveRDS(acqdat, file = "output/dsf_acquisition_data_set.rds")
saveRDS(xfpopdat, file = "output/dsf_acquisition_data_set_complete.rds")

###########################################################################################################
#### Analysis of Xf populations in vectors without zeros
xfpopdatNZ <- acqdatNZ %>% dplyr::filter(., !is.na(source.plant.pop) & source.plant.pop > 0)

xfvectorMod1 <- lmer(std.log.meancfu ~ genotype*std.distance*std.log.source.plant.pop + (1|plant),
                     data = xfpopdatNZ,
                     control = lmerControl(optimizer = "Nelder_Mead"))
xfvectorMod2 <- lmer(std.log.meancfu ~ genotype*std.log.source.plant.pop + (1|plant),
                     data = xfpopdatNZ,
                     control = lmerControl(optimizer = "Nelder_Mead"))
AICctab(xfvectorMod1, xfvectorMod2, base = TRUE)
# Keeping distance is better
plot(xfvectorMod1)
summary(xfvectorMod1)




######################################################################################################################
#### Combining plots into Supplementary Figure S2 for ms
figs2 <- plot_grid(sourceDistancePlot, vectorSourcePlot, vectorTransPlot, plantTransPlot,
                  align = "v", ncol = 2, nrow = 2, labels = "AUTO", label_size = 22)
figs2
ggsave(filename = "results/figureS2.tiff",
       plot = figs2,
       width = 12, height = 12, units = "in", dpi = 600)


