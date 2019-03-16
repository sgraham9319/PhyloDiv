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

library(ape)
library(picante)
library(geiger)
library(ncf)
library(nlme)

# Load functions file
source("R/utils.R")

#===============================
# 1. Formating data for analysis
#===============================

# Load community data
raw_site <- read.csv("Data/large_mammal.csv")

# Load mammal supertree
supertree <- read.nexus("Data/Mammal.supertree.nexus.txt")

# Create tree of sampled large mammal taxa
large_mammal_tree <- subset_supertree(raw_site, supertree, 2:57)

# Add columns for phylogenetic diversity and species richness - warning message
# is saying that PD could not be calculated for communities containing a
# single species
large_mammal <- faith_pd(raw_site, large_mammal_tree, 2:57)

#===============
# 2. Modeling PD
#===============

#-----------------
# Data exploration
#-----------------

# Check distributions of predictor variables
hist(large_mammal$annual_rainfall)
hist(large_mammal$soil_pc1)

# Check if the untransformed response variable has a normal distribution
hist(large_mammal$PD)
qqnorm(large_mammal$PD)
qqline(large_mammal$PD)

# Look for outliers in each of these variables
dotchart(large_mammal$PD, main = "PD")
dotchart(large_mammal$annual_rainfall, main = "Rain")
dotchart(large_mammal$soil_pc1, main = "Soil")

#-------------------------
# Testing for collinearity
#-------------------------

# Source code from Zuur et al 2009
source("R/HighstatLibV6.R")

# Create matrix of variables
var_mat <- cbind(large_mammal$PD, large_mammal$annual_rainfall, 
           large_mammal$soil_pc1)
colnames(var_mat) <- c("PD", "Rain", "Soil")

# Check for correlation among variables
pairs(var_mat, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Some correlation between rainfall and soil degradation (-0.5). Check variance 
# inflation factors (VIFs)
corvif(var_mat)

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
rand1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 | ranch, method = "REML", data = large_mammal)
rand2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
               landuse:soil_pc1, random = ~ 1 + annual_rainfall | ranch, 
             method = "REML", data = large_mammal)
rand3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
            landuse:soil_pc1, random = ~ 1 + soil_pc1 | ranch,
            control = lmc, method = "REML", data = large_mammal)


# Determine which random effects structure is best according to AIC
AICc(rand1, rand2, rand3) # M1 (random intercepts) is best

#----------------------------------
# Check for spatial autocorrelation
#----------------------------------

# Plot spatial correlogram of phylogenetic diversity measurements
fit1 <- correlog(x = large_mammal$longitude, y = large_mammal$latitude, 
                 z = large_mammal$PD, increment = 1, resamp = 500,
                 latlon = TRUE)
plot(fit1$mean.of.class, fit1$correlation, xlab = "Distance (km)", ylab = "Moran's I",
     main = "SAC in PD data")

# There appears to be some SAC among sites that are within 10 km of each other,
# but it is not extremely strong. Even among sites within 1 km, Moran's I is
# just under 0.5.

# Now create full model and check residuals for SAC
full_mod <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
                landuse:soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fit2 <- correlog(x = large_mammal$longitude, y = large_mammal$latitude,
                 z = residuals(full_mod), increment = 1, resamp = 500, latlon = TRUE)
plot(fit2$mean.of.class, fit2$correlation, xlab = "Distance (km)", 
     ylab = "Moran's I", main = "SAC in full model residuals")

# SAC appears to be accounted for by model covariates

#----------------------------------------
# Selecting model fixed effects structure
#----------------------------------------

# Create models off all combinations of fixed effects
fix1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
              landuse:soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall, 
            random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + 
              landuse:soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix5 <- lme(PD ~ landuse + annual_rainfall+ landuse:annual_rainfall,
            random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix6 <- lme(PD ~ landuse + soil_pc1 + landuse:soil_pc1,
            random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix7 <- lme(PD ~ landuse + annual_rainfall, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix8 <- lme(PD ~ landuse + soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix9 <- lme(PD ~ annual_rainfall + soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix10 <- lme(PD ~ landuse, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix11 <- lme(PD ~ annual_rainfall, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix12 <- lme(PD ~ soil_pc1, random = ~ 1 | ranch, method = "ML", data = large_mammal)
fix13 <- lme(PD ~ 1, random = ~ 1 | ranch, method = "ML", data = large_mammal)

# Calculate AICc for each model
AICc(fix1, fix2, fix3, fix4, fix5, fix6, fix7, fix8, fix9, fix10, fix11, fix12, fix13)

# Compare best model (fix1) to model with similar AICc (fix2) using
# likelihood ratio test
anova(fix2, fix1)

# Refit final model with REML
final <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
              landuse:soil_pc1, random = ~ 1 | ranch, method = "REML", data = large_mammal)

#-----------------
# Model validation
#-----------------

# Check for linearity and equal variance across range of fitted values
plot(fitted(final), residuals(final))
abline(0,0)

# Check assumption of normality in residuals
hist(residuals(final))
qqnorm(residuals(final))
qqline(residuals(final))

# Check for relationships between residuals and explanatory variables
plot(large_mammal$landuse, residuals(final))
plot(large_mammal$annual_rainfall, residuals(final))

#---------------------
# Model interpretation
#---------------------

# Extract coefficients from final model
summary(final)


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