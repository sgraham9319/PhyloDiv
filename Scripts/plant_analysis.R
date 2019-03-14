#####################
# Plant data analysis
#####################

# Contents of file
# 1.  Formating data for analysis - calculation of PD (Faith's Index)
#     using community data and the plant phylogeny obtained through 
#     phylomatic/phylocom
# 2.  Do distinct organismal groups respond differently? - calculating
#     PD response for each plot pair
# 3.  Are responses context-dependent? - modeling PD responses as a 
#     as a function of environmental variables and type of landuse change
# 4.  Are some lineages more strongly affected than others? - calculating
#     standardized effect sizes of MPD and MNTD for each community
# 5.  What composition changes are not captured by alpha diversity
#     metrics? - decomposing beta phylogenetic diversity into turnover
#     and nestedness components

# Load required packages
library(ape)
library(picante)
library(dplyr)
library(nlme)

# Load functions file
source("R/utils.R")

#===============================
# 1. Formating data for analysis
#===============================

# Load site data
raw_site <- read.csv("Data/plant.csv")

# Load phylogeny
phylo <- read.tree("Data/PlantPhylo")

# Add columns for phylogenetic diversity and species richness
plant <- faith_pd(raw_site, phylo, 4:151)

# Reorder factor levels for landuse so that conserved forms intercept in models
levels(plant$landuse) <- levels(plant$landuse)[c(2, 1, 3, 4)]

#===============
# 2. Modeling PD
#===============

#-----------------
# Data exploration
#-----------------

# Check distribution of dependent variable
hist(plant$PD)
qqnorm(plant$PD)
qqline(plant$PD)

# Check distributions of predictor variables
hist(plant$annual_rainfall)
hist(plant$soil_pc1)

# Check for outliers
dotchart(plant$PD, main = "PD")
dotchart(plant$annual_rainfall, main = "Rain")
dotchart(plant$soil_pc1, main = "Soil")

#-------------------------
# Testing for collinearity
#-------------------------

# Source code from Zuur et al 2009
source("R/HighstatLibV6.R")

# Create matrix of variables
Z <- cbind(plant$PD, plant$annual_rainfall, plant$soil_pc1)
colnames(Z) <- c("PD", "Rain", "Soil")

# Check for correlation among variables
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and soil degradation (-0.5). Check variance 
# inflation factors (VIFs)
corvif(Z)

# No VIFs are greater than 3, suggesting that collinearity is not a problem 
# (Zuur et al. 2009)

#-----------------------------------------
# Selecting model random effects structure
#-----------------------------------------

# This model selection protocol is adapted from the top-down strategy described
# in Zuur et al (2009) pp 121-122. First create models with the same fixed 
# effects structure (most complex structure possible) but with different
# random effects structure and compare them using AIC (not all are nested)

# Create object to be passed to lme to increase number of iterations and 
# achieve convergence (only use for models that fail to converge without this
# additional measure)
lmc <- lmeControl(niterEM = 5200, msMaxIter = 5200)

# Create models
M1 <- gls(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, method = "REML", data = plant)
M2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 | pair_id, method = "REML", data = plant)
M3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 + soil_pc1 | pair_id, method = "REML", data = plant)
M4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 + annual_rainfall | pair_id, 
            control = lmc, method = "REML", data = plant)

# Determine which random effects structure is best according to AIC
AIC(M1, M2, M3, M4) # M2 (random intercepts) is best

# Refit full model with ML then select optimal fixed effects structure
M5 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)

# Does inclusion of landuse:soil_pc1 improve model?
M6 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall, 
          random = ~ 1 | pair_id, method = "ML", data = plant)
anova(M6, M5) # No, new model is M6

# Does inclusion of landuse:annual_rainfall improve model?
M7 <- lme(PD ~ landuse + annual_rainfall + soil_pc1, 
          random = ~ 1 | pair_id, method = "ML", data = plant)
anova(M7, M6) # No, new model is M7

# Does inclusion of soil_pc1 improve model?
M8 <- lme(PD ~ landuse + annual_rainfall, 
          random = ~ 1 | pair_id, method = "ML", data = plant)
anova(M8, M7) # No, new model is M8

# Does inclusion of annual_rainfall improve model?
M9 <- lme(PD ~ landuse, 
          random = ~ 1 | pair_id, method = "ML", data = plant)
anova(M9, M8) # No, new model is M9

# Does inclusion of landuse improve model?
M10 <- lme(PD ~ 1, 
           random = ~ 1 | pair_id, method = "ML", data = plant)
anova(M10, M9) # Yes, final model is M9

# Refit final model with REML
MF <- lme(PD ~ landuse, random = ~ 1 | pair_id, method = "REML", data = plant)

#-----------------
# Model validation
#-----------------

# Check for linearity and equal variance across range of fitted values
plot(fitted(MF), residuals(MF))
abline(0,0)

# Check assumption of normality in residuals
hist(residuals(MF))
qqnorm(residuals(MF))
qqline(residuals(MF))

# Check for relationships between residuals and explanatory variables
plot(plant$landuse, residuals(MF))



#-------------------
# Creating figure 3a
#-------------------

# Make plot
plot(x = byPair$annRain, y = byPair$logPD, pch = 16,
     col = c("red","tan2","gold1")[byPair$landuse],
     ylab = "PD Response", yaxt = "n",
     xlab = "", xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(-1, 0.5, 0.5), labels = seq(-1, 0.5, 0.5), las = 1, 
     cex.axis = 0.8)

# Make global model and obtain confidence interval
mod <- lm(logPD ~ annRain, data = byPair)
modrain <- seq(min(byPair$annRain), max(byPair$annRain), length.out = 1000)
modconf <- predict(mod, newdata = data.frame(annRain = modrain), interval = "confidence")

# Add model line and delimit confidence intervals
lines(modrain, modconf[,1], col = "black", lwd = 1)
lines(modrain, modconf[,2], lty = 2, col = "black", lwd = 1)
lines(modrain, modconf[,3], lty = 2, col = "black", lwd = 1)

# Replot points over lines
points(x = byPair$annRain, y = byPair$logPD, pch = 16,
       col = c("red","tan2","gold1")[byPair$landuse])

#-------------------
# Creating figure S2
#-------------------

# Load function to make transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Create plot
plot(x = byPair$soilDeg, y = byPair$logPD, pch = 16,
     col = c("red","tan2","gold1")[byPair$landuse],
     ylab = "PD Response", yaxt = "n", xlab = "Soil degradation")
axis(side = 2, at = seq(-1, 0.5, 0.5), labels = seq(-1, 0.5, 0.5), las = 1, 
     cex.axis = 0.8)

# Subset data by landuse type and create separate model for each group
AG <- byPair[byPair$landuse == "Agriculture",]
EX <- byPair[byPair$landuse == "Exclosure",]
PA <- byPair[byPair$landuse == "Pastoral",]
AGmod <- lm(logPD ~ soilDeg, data = AG)
EXmod <- lm(logPD ~ soilDeg, data = EX)
PAmod <- lm(logPD ~ soilDeg, data = PA)

# Calculate confidence intervals
AGsoil <- seq(min(AG$soilDeg), max(AG$soilDeg), length.out = 50)
AGconf <- predict(AGmod, newdata = data.frame(soilDeg = AGsoil), interval = "confidence")
EXsoil <- seq(min(EX$soilDeg), max(EX$soilDeg), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(soilDeg = EXsoil), interval = "confidence")
PAsoil <- seq(min(PA$soilDeg), max(PA$soilDeg), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(soilDeg = PAsoil), interval = "confidence")

# Add confidence intervals to plot
polygon(c(rev(AGsoil), AGsoil), c(rev(AGconf[,3]), AGconf[,2]),
        col = makeTransparent("red", alpha = 70), border = NA)
polygon(c(rev(EXsoil), EXsoil), c(rev(EXconf[,3]), EXconf[,2]),
        col = makeTransparent("yellow2"), border = NA)
polygon(c(rev(PAsoil), PAsoil), c(rev(PAconf[,3]), PAconf[,2]),
        col = makeTransparent("tan2"), border = NA)

# Add lines to delimit confidence intervals
lines(AGsoil, AGconf[,1], col = "red", lwd = 2)
lines(AGsoil, AGconf[,2], lty=2, col = "red", lwd = 2)
lines(AGsoil, AGconf[,3], lty=2, col = "red", lwd = 2)

lines(EXsoil, EXconf[,1], col = "gold1", lwd = 2)
lines(EXsoil, EXconf[,2], lty=2, col = "gold1", lwd = 2)
lines(EXsoil, EXconf[,3], lty=2, col = "gold1", lwd = 2)

lines(PAsoil, PAconf[,1], col = "tan2", lwd = 2)
lines(PAsoil, PAconf[,2], lty=2, col = "tan2", lwd = 2)
lines(PAsoil, PAconf[,3], lty=2, col = "tan2", lwd = 2)

# Replot points over lines
points(x = byPair$soilDeg, y = byPair$logPD, pch = 16,
       col = c("red","tan2","gold1")[byPair$landuse])

#=========================================================
# 4. Are some lineages more strongly affected than others?
#=========================================================

# Calculate standardized effect sizes for MPD and MNTD
MPD.ES <- ses.mpd(pres.abs[,-1], cophenetic(phylo), null.model = "independentswap", 
                  abundance.weighted = F, runs = 999, iterations = 1000)
MNTD.ES <- ses.mntd(pres.abs[,-1], cophenetic(phylo), null.model = "independentswap", 
                    abundance.weighted = F, runs = 999, iterations = 1000)

# Combine with PD data (to match with site name)
ES.dat <- cbind(PDData,
                MPD.ES[,c("mpd.obs.p", "mpd.obs.z")],
                MNTD.ES[,c("mntd.obs.p", "mntd.obs.z")])

# Add landuse type data
ESDat <- cbind(byPlot[, c("Site", "Landuse")], ES.dat[match(byPlot$Site, ES.dat$Site), -1])

# Define standard error function
SE <- function(x){sd(x, na.rm = T)/sqrt(length(!is.na(x)))}

# Summarize data for plotting
ESPlot <- ESDat %>%
  group_by(Landuse) %>%
  summarise(MPDz = mean(mpd.obs.z), MPDz.SE = SE(mpd.obs.z),
            MPDcent = mean(mpd.obs.p), MPDcent.SE = SE(mpd.obs.p),
            MNTDz = mean(mntd.obs.z), MNTDz.SE = SE(mntd.obs.z),
            MNTDcent = mean(mntd.obs.p), MNTDcent.SE = SE(mntd.obs.p))

# Reorder treatments for plotting
ESPlot <- ESPlot[c(2, 1, 3, 4),]

# Save summary data to csv for plotting figure 4
write.csv(ESPlot, file = "../Data/Plot data/Plant_MPD_MNTD.csv", row.names = F)

#=========================================================================
# 5. What composition changes are not captured by alpha diversity metrics?
#=========================================================================

# Load functions from Leprieur et al 2012
source("../Data/Leprieuretal2012.R")

# Create dataframe to store decomposition results
TotalBetaPD <- rep(NA, times = nrow(byPair))
Turnover <- rep(NA, times = nrow(byPair))
Nestedness <- rep(NA, times = nrow(byPair))
results <- data.frame(byPair$PairID, byPair$landuse, TotalBetaPD, Turnover, 
                      Nestedness)
colnames(results)[1:2] <- c("PairID", "Landuse")

# Combine raw community data with pair ID data
dat <- cbind(env[, "PairID"], pres.abs[match(env$Site, pres.abs$Site),])
colnames(dat)[1] <- "PairID"

# Calculate beta PD and its components for each plot pair
for(i in unique(byPair$PairID)){
  com <- dat[dat$PairID == i, 3:ncol(dat)]
  results[results$PairID == i,3:5] <- beta.pd.decompo(com = com, tree = phylo, 
                                                      type = "Unifrac")$betadiv
}

# Summarize results by type of land use change
falseComp <- results %>% group_by(Landuse) %>% 
  summarise(Total = mean(TotalBetaPD),
            TotalSD = sd(TotalBetaPD),
            TurnoverMean = mean(Turnover),
            TurnoverSD = sd(Turnover),
            TurnoverSE = SE(Turnover),
            Nested = mean(Nestedness),
            NestedSD = sd(Nestedness))

# Re-order for plotting
falseComp <- falseComp[c(1, 3, 2),]

# Save summary data to csv for plotting figure 6
write.csv(falseComp, file = "../Data/Plot data/Plant_false_comp.csv", row.names = F)
