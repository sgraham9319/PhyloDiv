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

#===============================
# 1. Formating data for analysis
#===============================

# Load community data
rawCom <- read.csv("../Data/PlantRawCom.csv")

# Convert to presence/absence community composition matrix
pres.abs <- rawCom[, - c(2, 3)]
for(i in 2:ncol(pres.abs)){
  pres.abs[,i][which(pres.abs[, i] > 0)] <- 1
}

# Load phylogeny
phylo <- read.tree("../Data/PlantPhylo")

# Calculate PD for each site and store in data frame
PDData <- data.frame(pres.abs$Site, pd(pres.abs[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

# Calculate PD of a community containing all taxa in regional phylogeny
allTaxaCom <- pres.abs[1, -1]
allTaxaCom[1,] <- 1
pd(allTaxaCom, phylo, include.root = F)

# Load pair ID and environmental data
env <- read.csv("../Data/Plant data by plot.csv")

# Check site names match between PD and environmental data
length(env$Site) == sum(env$Site %in% PDData$Site)

# Combine PD and environmental data
byPlot <- cbind(env, PDData[match(env$Site, PDData$Site), 
                               c("PD", "SR")])

#======================================================
# 2. Do distinct organismal groups respond differently?
#======================================================

# Calculate PD responses
byPair <- byPlot %>% 
  group_by(PairID) %>% 
  summarise(landuse = droplevels(Landuse[Landuse != "Conserved"]), 
            annRain = mean(AnnualRainfall), 
            soilDeg = mean(SoilPC1), 
            PDrawChg = PD[Landuse != "Conserved"] - PD[Landuse == "Conserved"],
            logPD = log(PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"]),
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"],
            SRrawChg = SR[Landuse != "Conserved"] - SR[Landuse == "Conserved"],
            SR = SR[Landuse != "Conserved"] / SR[Landuse == "Conserved"])

# Write PD response data to csv for plotting figures 1 and 2
write.csv(byPair[, c("PairID", "landuse", "PD", "logPD")],
          file = "../Data/Plot data/Plant_PD_response.csv", row.names = F)

#====================================
# 3. Are responses context-dependent?
#====================================

#-----------------
# Data exploration
#-----------------

# Check distributions of predictor variables
hist(byPair$annRain)
hist(byPair$soilDeg)

# Check if the untransformed response variable has a normal distribution
hist(byPair$PD)
qqnorm(byPair$PD)
qqline(byPair$PD)

# Response variable a little right-skewed. Used log transformed values
hist(byPair$logPD)
qqnorm(byPair$logPD)
qqline(byPair$logPD)

# This looks much better. Continue with log transformed PD as the response
# variable. Now look for outliers in each of these variables
dotchart(byPair$logPD, main = "PD")
dotchart(byPair$annRain, main = "Rain")
dotchart(byPair$soilDeg, main = "Soil")

#-------------------------
# Testing for collinearity
#-------------------------

# Source code from Zuur et al 2009
source("HighstatLibV6.R")

# Create matrix of variables
Z <- cbind(byPair$logPD, byPair$annRain, byPair$soilDeg)
colnames(Z) <- c("logPD", "Rain", "Soil")

# Check for correlation among variables
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and soil degradation (-0.6). Check variance 
# inflation factors (VIFs)
corvif(Z)

# No VIFs are greater than 3, suggesting that collinearity is not a problem 
# (Zuur et al. 2007)

#----------------
# Model selection
#----------------

# Create full model
M1 <- lm(logPD ~ landuse + annRain + soilDeg + landuse:annRain + 
           landuse:soilDeg, data = byPair)

# Use AIC to select final model
step(M1)

# Create final model according to AIC (same as full model)
Mfinal <- lm(logPD ~ landuse + annRain + soilDeg + landuse:annRain + 
           landuse:soilDeg, data = byPair)

#-----------------
# Model validation
#-----------------

# Check for linearity and equal variance across range of fitted values
plot(fitted(Mfinal), residuals(Mfinal))
abline(0,0)

# Equal variance and linearity assumptions appear justfied. Check assumption
# of normality in residuals
hist(residuals(Mfinal))
qqnorm(residuals(Mfinal))
qqline(residuals(Mfinal))

# This looks close enough to a normal distribution. Finally, check for
# influential data points
Influence <- influence.measures(Mfinal)
plot(Influence$infmat[,"cook.d"])

# No points have a Cook's distance greater than 1, suggesting influential
# data points are not a problem

#---------------------
# Model interpretation
#---------------------

# Check model summary
summary(Mfinal)

# Create reduced models for significance testing
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
