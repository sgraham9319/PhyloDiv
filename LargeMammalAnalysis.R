############################
# Small mammal data analysis
############################

library(ape)
library(picante)
library(geiger)
library(ncf)

#========================================
# Calculating Phylogenetic Diversity (PD)
#========================================

# Load mammal supertree
mammal.supertree.phylos <- read.nexus("../Data/Mammal.supertree.nexus.txt")

# Select the "bestDates" tree
mammal.tree <- mammal.supertree.phylos$mammalST_bestDates

# Load large mammal community data
rawCom <- read.csv("../Data/Large mammal com data.csv")

# Remove sites where no large mammals were observed
rawCom <- rawCom[-which(apply(rawCom[,-1], MARGIN = 1, FUN = sum) == 0),]

# Remove species that did not appear in any of sampled sites
spsRm <- names(which(apply(rawCom[,-1], MARGIN = 2, FUN = sum) == 0))
rawCom <- rawCom[,-which(names(rawCom) == spsRm)]

# Prune supertree to only the sampled taxa
lmammal <- as.matrix(t(rawCom))
phylo <- treedata(mammal.tree, lmammal)$phy

# Warning returned from above command just says that not all taxa in the mammal 
# supertree were found in our community sample
plot(phylo, type = "fan", cex = 0.6)

# Calculate PD (Evolutionary Heritage) for each site and store in data frame
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = T))
names(PDData)[1] <- "Site"

#==================================
# Add PD data to sampling site data
#==================================

# Load site data
byPlot <- read.csv("../Data/Large mammal data.csv")

# Subset to plots for which PD could be calculated
byPlot <- byPlot[which(byPlot$Site %in% PDData$Site),]

# Check site names match
length(PDData$Site) == sum(PDData$Site %in% byPlot$Site)

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR")])

# Remove the single site that has no rainfall data
byPlot <- byPlot[-which(is.na(byPlot$AnnualRainfall)),]

#===========================================
# Checking for spatial autocorrelation (SAC)
#===========================================

# Plot spatial correlogram of phylogenetic diversity measurements
fit1 <- correlog(x = byPlot$Longitude, y = byPlot$Latitude, 
                 z = byPlot$PD, increment = 1, resamp = 500,
                 latlon = TRUE)
plot(fit1$mean.of.class, fit1$correlation, xlab="Distance class", ylab="Moran's I",
     main="SAC in PD data")

# There appears to be some SAC among sites that are within 10 km of each other, but it
# is not extremely strong. Even among sites within 1 km, Moran's I is just under 0.5.

# Now create full model and check residuals for SAC
M1 <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)
summary(M1)

fit2 <- correlog(x = byPlot$Longitude, y = byPlot$Latitude, z = residuals(M1),
                 increment = 1, resamp = 500, latlon = TRUE)
plot(fit2$mean.of.class, fit2$correlation, xlab = "Distance class", 
     ylab = "Moran's I", main = "SAC in full model residuals")

# SAC appears to be accounted for by model covariates

#================
# Linear modeling
#================

# Data exploration
# Look at distributions of predictor variables
hist(byPlot$AnnualRainfall)
hist(byPlot$SoilPC1)

# Check if the response variable has a normal distribution
hist(byPlot$PD)
qqnorm(byPlot$PD)
qqline(byPlot$PD)

# Look for any outliers
dotchart(byPlot$PD, main = "PD")
dotchart(byPlot$AnnualRainfall, main = "Rain")
dotchart(byPlot$SoilPC1, main = "Soil")

### Testing for collinearity
source("C:/Users/stuart/Documents/Science books/Zuur et al 2009 Mixed Models/HighstatLibV6.R")

Z <- cbind(byPlot$PD, byPlot$AnnualRainfall, byPlot$SoilPC1)
colnames(Z) <- c("PD", "Rain", "Soil")
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and PC1 (-0.5). Check variance inflation factors 
# (VIFs)
corvif(Z)

# No VIFs are greater than 3 (Zuur et al. 2007) which suggests that collinearity is
# not a problem.

### Model selection

# Create full model
M1 <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)

# Use AIC to select final model
step(M1)

Mfinal <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)

### Model validation

# Check for linearity and equal variance across range of fitted values
plot(fitted(Mfinal), residuals(Mfinal))
abline(0,0)

# Equal variance assumption may not be justified, now check assumption
# of normality in residuals
hist(residuals(Mfinal))
qqnorm(residuals(Mfinal))
qqline(residuals(Mfinal))

# I think this is close enough to a normal distribution. Finally, check for
# influential data points
Influence <- influence.measures(Mfinal)
Cooks.distances <- Influence$infmat[,"cook.d"]
plot(Cooks.distances)

# No points have a Cook's distance greater than 1.

### Model interpretation
summary(Mfinal)

# Significance testing
nolandsoil <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall, 
                   data = byPlot)
nolandrain <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 +
                   Landuse:SoilPC1, data = byPlot)
nointer <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1, data = byPlot)
nosoil <- lm(PD ~ Landuse + AnnualRainfall + Landuse:AnnualRainfall, 
             data = byPlot)
norain <- lm(PD ~ Landuse + SoilPC1 + Landuse:SoilPC1, data = byPlot)
noland <- lm(PD ~ AnnualRainfall + SoilPC1, data = byPlot)

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

