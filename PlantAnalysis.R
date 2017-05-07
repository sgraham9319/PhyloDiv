#####################
# Plant data analysis
#####################

library(ape)
library(picante)
library(dplyr)

#========================================
# Calculating Phylogenetic Diversity (PD)
#========================================

# Load community data
rawCom <- read.csv("../Data/Plant community data.csv")

# Load phylogeny
phylo <- read.tree("../Data/PlantPhylo")

# Calculate PD (Evolutionary Heritage) for each site and store in data frame
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

#======================================
# Add PD data to experimental plot data
#======================================

# Load site data
byPlot <- read.csv("../Data/Plant data by plot.csv")

# Check site names match
length(byPlot$Site) == sum(byPlot$Site %in% PDData$Site)

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR")])

#=========================
# Format data for modeling
#=========================

byPair <- byPlot %>% 
  group_by(PairID) %>% 
  summarise(landuse = Landuse[Landuse != "Conserved"], 
            annRain = mean(AnnualRainfall), 
            soilDeg = mean(SoilPC1), 
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"])

#================
# Linear modeling
#================

### Data exploration
# Look at distributions of predictor variables
hist(byPair$annRain)
hist(byPair$soilDeg)

# Predictor variables look okay. Check if the response variable has a 
# normal distribution
hist(byPair$PD)
qqnorm(byPair$PD)
qqline(byPair$PD)

# Response variable a little right-skewed. Try log transformation
byPair$logPD <- log(byPair$PD)
hist(byPair$logPD)
qqnorm(byPair$logPD)
qqline(byPair$logPD)

# This looks much better. Continue with log transformed EH as the response
# variable. Now let's look for outliers in each of these variables
dotchart(byPair$logPD, main = "PD")
dotchart(byPair$annRain, main = "Rain")
dotchart(byPair$soilDeg, main = "Soil")

### Testing for collinearity
source("C:/Users/stuart/Documents/Science books/Zuur et al 2009 Mixed Models/HighstatLibV6.R")

Z <- cbind(byPair$logPD, byPair$annRain, byPair$soilDeg)
colnames(Z) <- c("logPD", "Rain", "Soil")
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and PC1 (-0.6). Check variance inflation 
# factors (VIFs)
corvif(Z)

# No VIFs are greater than 3, suggesting that collinearity is not a problem 
# (Zuur et al. 2007)

### Model selection
M1 <- lm(logPD ~ landuse + annRain + soilDeg + landuse:annRain + 
           landuse:soilDeg, data = byPair)

# Use AIC to select final model
step(M1)
Mfinal <- lm(logPD ~ landuse + annRain + soilDeg + landuse:annRain + 
           landuse:soilDeg, data = byPair)

### Model validation

# Check for linearity and equal variance across range of fitted values
plot(fitted(Mfinal), residuals(Mfinal))
abline(0,0)

# Equal variance and linearity assumptions appear justfieid, now check assumption
# of normality in residuals
hist(residuals(Mfinal))
qqnorm(residuals(Mfinal))
qqline(residuals(Mfinal))

# I think this is close enough to a normal distribution. Finally, check for
# influential data points
Influence <- influence.measures(Mfinal)
plot(Influence$infmat[,"cook.d"])

# No points have a Cook's distance greater than 1.

### Model interpretation
summary(Mfinal)

# Significance testing
nolandrain <- lm(logPD ~ landuse + annRain + soilDeg +  
                   landuse:soilDeg, data = byPair)
nolandsoil <- lm(logPD ~ landuse + annRain + soilDeg + landuse:annRain, 
                   data = byPair)
norain <- lm(logPD ~ landuse + soilDeg + landuse:soilDeg, data = byPair)
nosoil <- lm(logPD ~ landuse + annRain + landuse:annRain, data = byPair)
nointer <- lm(logPD ~ landuse + annRain + soilDeg, data = byPair)
noland <- lm(logPD ~ annRain + soilDeg, data = byPair)
# Effect of landuse/rainfall interaction
anova(nolandrain, Mfinal)
# Effect of landuse/soil interaction
anova(nolandsoil, Mfinal)
# Effect of rainfall
anova(norain, nolandrain)
# Effect of soil
anova(nosoil, nolandsoil)
# Effect of landuse
anova(noland, nointer)

# Model selection using LRT might work better for plants because this leads
# to removal of landuse by rainfall interaction which is not significant (p = 0.14).
# The removal of this interaction increases AIC by < 1.

#==============
# Plotting data
#==============
SE <- function(x) {
  sd(x, na.rm = T)/sqrt(length(!is.na(x)))
}

Dat <- byPair %>% group_by(landuse) %>% summarise(mean = mean(PD), SE = SE(PD))

plot(Dat$landuse, Dat$mean, ylim=range(c(Dat$mean-Dat$SE, 
                                         Dat$mean+Dat$SE)),
     pch=19)
arrows(c(1,3,4), Dat$mean-Dat$SE, c(1,3,4), 
       Dat$mean+Dat$SE, code=0)


