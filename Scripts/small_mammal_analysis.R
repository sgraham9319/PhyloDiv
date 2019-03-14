############################
# Small mammal data analysis
############################

# Contents of file
# 1.  Calculation of PD (Faith's Index) and SR using community data 
#     and the mammal supertree (Bininda-Emonds et al. 2007)
# 2.  Calculation of MPD and MNTD
# 3.  Combining site environmental data to PD, SR, MPD and MNTD data
# 4.  Calculation of PD, SR, MPD, and MNTD response for each plot pair
# 5.  Linear modeling of PD
# 6.  Determining whether PD, MPD, and MNTD are greater or less than 
#     expected based on SR
# 7.  Summary plots

library(ape)
library(picante)
library(geiger)
library(dplyr)
library(nlme)

# Load functions file
source("R/utils.R")

#===============================
# 1. Formating data for analysis
#===============================

# Load community data
raw_site <- read.csv("Data/small_mammal.csv")

# Load mammal supertree
supertree <- read.nexus("Data/Mammal.supertree.nexus.txt")

# Create tree of sampled small mammal taxa
small_mammal_tree <- subset_supertree(raw_site, supertree, 4:25)

# Add columns for phylogenetic diversity and species richness - warning message
# is saying that PD could not be calculated for communities containing a
# single species
small_mammal <- faith_pd(raw_site, small_mammal_tree, 4:25)

# Reorder factor levels for landuse so that conserved forms intercept in models
levels(small_mammal$landuse) <- levels(small_mammal$landuse)[c(2, 1, 3, 4)]

# Round abundances to whole numbers - MAY NEED THIS LATER!
# rawCom[, 4:ncol(rawCom)] <- round(rawCom[, 4:ncol(rawCom)])

#===============
# 2. Modeling PD
#===============

#-----------------
# Data exploration
#-----------------

# Check distributions of predictor variables
hist(small_mammal$annual_rainfall)
hist(small_mammal$soil_pc1)

# Check if the untransformed response variable has a normal distribution
hist(small_mammal$PD)
qqnorm(small_mammal$PD)
qqline(small_mammal$PD)

# Look for outliers in each of these variables
dotchart(small_mammal$PD, main = "PD")
dotchart(small_mammal$annual_rainfall, main = "Rain")
dotchart(small_mammal$soil_pc1, main = "Soil")

#-------------------------
# Testing for collinearity
#-------------------------

# Source code from Zuur et al 2009
source("R/HighstatLibV6.R")

# Create matrix of variables
Z <- cbind(small_mammal$PD, small_mammal$annual_rainfall, 
           small_mammal$soil_pc1)
colnames(Z) <- c("PD", "Rain", "Soil")

# Check for correlation among variables
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and soil degradation (-0.5). Check variance 
# inflation factors (VIFs)
corvif(Z)

# No VIFs are greater than 3, suggesting that collinearity is not a problem 
# (Zuur et al. 2007)

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
            landuse:soil_pc1, method = "REML", data = small_mammal)
M2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 | pair_id, method = "REML", data = small_mammal)
M3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 + soil_pc1 | pair_id, method = "REML", data = small_mammal)
M4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 + annual_rainfall | pair_id, 
          control = lmc, method = "REML", data = small_mammal)

# Determine which random effects structure is best according to AIC
AIC(M1, M2, M3, M4) # M1 no random effects is best

# Refit full model with ML then select optimal fixed effects structure
M5 <- gls(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, method = "ML", data = small_mammal)

# Does inclusion of landuse:soil_pc1 improve model?
M6 <- gls(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall,
          method = "ML", data = small_mammal)
anova(M6, M5) # No, new model is M6

# Does inclusion of landuse:annual_rainfall improve model?
M7 <- gls(PD ~ landuse + annual_rainfall + soil_pc1,
          method = "ML", data = small_mammal)
anova(M7, M6) # Yes, model is still M6

# Does inclusion of soil_pc1 improve model?
M8 <- gls(PD ~ landuse + annual_rainfall + landuse:annual_rainfall,
          method = "ML", data = small_mammal)
anova(M8, M6) # No, final model is M8

# Refit final model with REML
MF <- gls(PD ~ landuse + annual_rainfall + landuse:annual_rainfall,
          method = "REML", data = small_mammal)

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
plot(small_mammal$landuse, residuals(MF))
plot(small_mammal$annual_rainfall, residuals(MF))






#-------------------
# Creating figure 3b
#-------------------

# Make plot
plot(x = byPair$annRain, y = byPair$logPD, pch = 16,
     col = c("red","tan2","gold1")[byPair$landuse],
     ylab = "PD Response", yaxt = "n",
     xlab = "", xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(-1.5, 1.5, 0.5), labels = seq(-1.5, 1.5, 0.5), las = 1, 
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
ESDat <- cbind(byPlot[, c("Site", "Landuse", "PairID")], 
               ES.dat[match(byPlot$Site, ES.dat$Site), -1])

# Remove pairs not included in analysis
ESDat <- ESDat[ESDat$PairID %in% byPair$PairID,]

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
write.csv(ESPlot, file = "../Data/Plot data/Small_mammal_MPD_MNTD.csv", row.names = F)

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

# Remove pairs not included in analysis
dat <- dat[dat$PairID %in% byPair$PairID,]

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
write.csv(falseComp, file = "../Data/Plot data/Small_mammal_false_comp.csv", row.names = F)
