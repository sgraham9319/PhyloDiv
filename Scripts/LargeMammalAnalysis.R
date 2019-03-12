############################
# Large mammal data analysis
############################

# Contents of file
# 1.  Calculation of PD (Faith's Index) and SR using community data 
#     and the mammal supertree (Bininda-Emonds et al. 2007)
# 2.  Calculation of MPD and MNTD
# 3.  Combining site environmental data to PD, SR, MPD and MNTD data
# 4.  Check data for spatial autocorrelation
# 5.  Linear modeling of PD
# 6.  Determining whether PD, MPD, and MNTD are greater or less than 
#     expected based on SR


# 4.  Calculation of PD, SR, MPD, and MNTD response for each plot pair
# 7.  Summary plots

library(ape)
library(picante)
library(geiger)
library(ncf)

#===============================
# 1. Formating data for analysis
#===============================

# Load community data
rawCom <- read.csv("../Data/Lg_mammal_com_dom_incl.csv")

# Remove species that did not appear in any of sampled sites
spsRm <- names(which(apply(rawCom[,-1], MARGIN = 2, FUN = sum) == 0))
rawCom <- rawCom[,-which(names(rawCom) == spsRm)]

# Load mammal supertree
mammal.supertree.phylos <- read.nexus("../Data/Mammal.supertree.nexus.txt")

# Select the "bestDates" tree
mammal.tree <- mammal.supertree.phylos$mammalST_bestDates

# Prune supertree to only the sampled taxa
lmammal <- as.matrix(t(rawCom))
phylo <- treedata(mammal.tree, lmammal)$phy # Warning = not all taxa in supertree
                                            # were found in community data

# Calculate PD for each site and store in data frame - warning messages are saying
# that PD could not be calculated for communities containing a single species
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

# Calculate PD of a community containing all taxa in regional phylogeny
allTaxaCom <- rawCom[1, -1]
allTaxaCom[1,] <- 1
pd(allTaxaCom, phylo, include.root = F)

# Manually enter zeroes in the PD column for communities that contained < 2 taxa
PDData$PD[which(PDData$SR < 2)] <- 0

# Load environmental data
env <- read.csv("../Data/Large mammal data.csv")

# Check site names match
length(PDData$Site) == sum(PDData$Site %in% env$Site)

# Combine PD and environment data
byPlot <- cbind(env, PDData[match(env$Site, PDData$Site), 
                               c("PD", "SR")])

#======================================================
# 2. Do distinct organismal groups respond differently?
#======================================================

# Create subsets
Con <- byPlot[byPlot$Landuse == "Conserved" & byPlot$SR >= 2,]
Past <- byPlot[byPlot$Landuse == "Pastoral" & byPlot$SR >= 2,]
Exc <- byPlot[byPlot$Landuse == "Fenced" & byPlot$SR >= 2,]

# Obtain random sample of 1000 sites from each land-use type
resampCon <- Con[sample(x = 1:nrow(Con), size = 1000, replace = T),]
resampPast <- Past[sample(x = 1:nrow(Past), size = 1000, replace = T),]
resampExc <- Exc[sample(x = 1:nrow(Exc), size = 1000, replace = T),]

# Calculate PD response variables for each type of land use change
PDPast <- resampPast$PD / resampCon$PD
PDExc <- resampExc$PD / resampCon$PD

# Calculate log PD responses
logPDPast <- log(PDPast)
logPDExc <- log(PDExc)

# Combine in data frame
PD <- c(PDPast, PDExc)
logPD <- c(logPDPast, logPDExc)
landuse <- rep(c("Pastoral", "Exclosure"), each = 1000)
PairID <- rep(NA, times = 2000)
lMamm <- data.frame(PairID, landuse, PD, logPD)

# Write PD response data to csv for plotting figure 1
write.csv(lMamm, file = "../Data/Plot data/Large_mammal_PD_response.csv",
          row.names = F)

#====================================
# 3. Are responses context-dependent?
#====================================

#-------------------------------------------
# Checking for spatial autocorrelation (SAC)
#-------------------------------------------

# Plot spatial correlogram of phylogenetic diversity measurements
fit1 <- correlog(x = byPlot$Longitude, y = byPlot$Latitude, 
                 z = byPlot$PD, increment = 1, resamp = 500,
                 latlon = TRUE)
plot(fit1$mean.of.class, fit1$correlation, xlab = "Distance class", ylab = "Moran's I",
     main = "SAC in PD data")

# There appears to be some SAC among sites that are within 10 km of each other,
# but it is not extremely strong. Even among sites within 1 km, Moran's I is
# just under 0.5.

# Now create full model and check residuals for SAC
M1 <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)
fit2 <- correlog(x = byPlot$Longitude, y = byPlot$Latitude, z = residuals(M1),
                 increment = 1, resamp = 500, latlon = TRUE)
plot(fit2$mean.of.class, fit2$correlation, xlab = "Distance class", 
     ylab = "Moran's I", main = "SAC in full model residuals")

# SAC appears to be accounted for by model covariates

#-----------------
# Data exploration
#-----------------

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

#-------------------------
# Testing for collinearity
#-------------------------

# Source code from Zuur et al 2009
source("HighstatLibV6.R")

# Create matrix of variables
Z <- cbind(byPlot$PD, byPlot$AnnualRainfall, byPlot$SoilPC1)
colnames(Z) <- c("PD", "Rain", "Soil")

# Check for correlation among variables
pairs(Z, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and PC1 (-0.5). Check variance inflation 
# factors (VIFs)
corvif(Z)

# No VIFs are greater than 3, suggesting that collinearity is not a problem 
# (Zuur et al. 2007)

#----------------
# Model selection
#----------------

# Create full model
M1 <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)

# Use AIC to select final model
step(M1)

# Create final model (same as full model)
Mfinal <- lm(PD ~ Landuse + AnnualRainfall + SoilPC1 + Landuse:AnnualRainfall + 
           Landuse:SoilPC1, data = byPlot)

#-----------------
# Model validation
#-----------------

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

#---------------------
# Model interpretation
#---------------------

# Check model summary
summary(Mfinal)

# Create reduced models for significance testing
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

#-------------------
# Creating figure 3c
#-------------------

# Make plot
plot(x = byPlot$AnnualRainfall, y = byPlot$PD, type = "n",
     ylab = "PD", xlab = "Annual rainfall (mm)",
     ylim = c(0, 500), yaxt = "n", yaxs = "i",
     xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(0,500,100), labels = seq(0,500,100), las = 1,
     cex.axis = 0.8)
axis(side = 1, at = seq(450,750,50), labels = seq(450,750,50),
     cex.axis = 0.8)

# Create land-use type subsets so that separate rainfall models can be made
CN <- byPlot[byPlot$Landuse == "Conserved",]
EX <- byPlot[byPlot$Landuse == "Fenced",]
PA <- byPlot[byPlot$Landuse == "Pastoral",]
CNmod <- lm(PD ~ AnnualRainfall, data = CN)
EXmod <- lm(PD ~ AnnualRainfall, data = EX)
PAmod <- lm(PD ~ AnnualRainfall, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$AnnualRainfall), max(CN$AnnualRainfall), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(AnnualRainfall = CNrain), interval = "confidence")
EXrain <- seq(min(EX$AnnualRainfall), max(EX$AnnualRainfall), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(AnnualRainfall = EXrain), interval = "confidence")
PArain <- seq(min(PA$AnnualRainfall), max(PA$AnnualRainfall), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(AnnualRainfall = PArain), interval = "confidence")

# Load function to make transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Add confidence intervals to plot
polygon(c(rev(CNrain), CNrain), c(rev(CNconf[,3]), CNconf[,2]),
        col = makeTransparent("paleturquoise3", alpha = 70), border = NA)
polygon(c(rev(EXrain), EXrain), c(rev(EXconf[,3]), EXconf[,2]),
        col = makeTransparent("yellow2"), border = NA)
polygon(c(rev(PArain), PArain), c(rev(PAconf[,3]), PAconf[,2]),
        col = makeTransparent("tan2"), border = NA)

# Add lines to delimit confidence intervals
lines(CNrain, CNconf[,1], col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,2], lty=2, col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,3], lty=2, col = "skyblue2", lwd = 2)

lines(EXrain, EXconf[,1], col = "gold1", lwd = 2)
lines(EXrain, EXconf[,2], lty=2, col = "gold1", lwd = 2)
lines(EXrain, EXconf[,3], lty=2, col = "gold1", lwd = 2)

lines(PArain, PAconf[,1], col = "tan2", lwd = 2)
lines(PArain, PAconf[,2], lty=2, col = "tan2", lwd = 2)
lines(PArain, PAconf[,3], lty=2, col = "tan2", lwd = 2)

#-------------------
# Creating figure S3
#-------------------

# Create plot
plot(x = byPlot$SoilPC1, y = byPlot$PD, type = "n",
     ylab = "PD (million years)", xlab = "Soil degradation",
     ylim = c(0, 600), yaxt = "n", yaxs = "i")
axis(side = 2, at = seq(0,600,100), labels = seq(0,600,100), las = 1,
     cex.axis = 0.8)

# Subset data by landuse type and create separate model for each group
CN <- byPlot[byPlot$Landuse == "Conserved",]
EX <- byPlot[byPlot$Landuse == "Fenced",]
PA <- byPlot[byPlot$Landuse == "Pastoral",]
CNmod <- lm(PD ~ SoilPC1, data = CN)
EXmod <- lm(PD ~ SoilPC1, data = EX)
PAmod <- lm(PD ~ SoilPC1, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$SoilPC1), max(CN$SoilPC1), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(SoilPC1 = CNrain), interval = "confidence")
EXrain <- seq(min(EX$SoilPC1), max(EX$SoilPC1), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(SoilPC1 = EXrain), interval = "confidence")
PArain <- seq(min(PA$SoilPC1), max(PA$SoilPC1), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(SoilPC1 = PArain), interval = "confidence")

# Add confidence intervals to plot
polygon(c(rev(CNrain), CNrain), c(rev(CNconf[,3]), CNconf[,2]),
        col = makeTransparent("paleturquoise3", alpha = 70), border = NA)
polygon(c(rev(EXrain), EXrain), c(rev(EXconf[,3]), EXconf[,2]),
        col = makeTransparent("yellow2"), border = NA)
polygon(c(rev(PArain), PArain), c(rev(PAconf[,3]), PAconf[,2]),
        col = makeTransparent("tan2"), border = NA)

# Add lines to delimit confidence intervals
lines(CNrain, CNconf[,1], col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,2], lty=2, col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,3], lty=2, col = "skyblue2", lwd = 2)

lines(EXrain, EXconf[,1], col = "gold1", lwd = 2)
lines(EXrain, EXconf[,2], lty=2, col = "gold1", lwd = 2)
lines(EXrain, EXconf[,3], lty=2, col = "gold1", lwd = 2)

lines(PArain, PAconf[,1], col = "tan2", lwd = 2)
lines(PArain, PAconf[,2], lty=2, col = "tan2", lwd = 2)
lines(PArain, PAconf[,3], lty=2, col = "tan2", lwd = 2)

#=========================================================
# 4. Are some lineages more strongly affected than others?
#=========================================================

# Calculate standardized effect sizes for MPD and MNTD
MPD.ES <- ses.mpd(rawCom[,-1], cophenetic(phylo), null.model = "independentswap", 
                  abundance.weighted = F, runs = 999, iterations = 1000)
MNTD.ES <- ses.mntd(rawCom[,-1], cophenetic(phylo), null.model = "independentswap", 
                    abundance.weighted = F, runs = 999, iterations = 1000)

# Combine with PD data (to match with site name)
ES.dat <- cbind(PDData,
                MPD.ES[,c("mpd.obs.p", "mpd.obs.z")],
                MNTD.ES[,c("mntd.obs.p", "mntd.obs.z")])

# Add landuse type data
ESDat <- cbind(byPlot[, c("Site", "Landuse")], 
               ES.dat[match(byPlot$Site, ES.dat$Site), -1])

# Define standard error function
SE <- function(x){sd(x, na.rm = T)/sqrt(length(!is.na(x)))}

# Summarize data for plotting
ESPlot <- ESDat %>%
  group_by(Landuse) %>%
  summarise(MPDz = mean(mpd.obs.z, na.rm = T), MPDz.SE = SE(mpd.obs.z),
            MPDcent = mean(mpd.obs.p, na.rm = T), MPDcent.SE = SE(mpd.obs.p),
            MNTDz = mean(mntd.obs.z, na.rm = T), MNTDz.SE = SE(mntd.obs.z),
            MNTDcent = mean(mntd.obs.p, na.rm = T), MNTDcent.SE = SE(mntd.obs.p))

# Save summary data to csv for plotting figure 4
write.csv(ESPlot, file = "../Data/Plot data/Large_mammal_MPD_MNTD.csv", row.names = F)