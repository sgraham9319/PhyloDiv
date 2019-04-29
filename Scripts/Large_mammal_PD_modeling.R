##########################
# Modeling large mammal PD
##########################

# This script creates and interprets the generalized linear mixed effects model
# of large mammal PD

# Outputs: Figure 2b, Figure S3

# Load required packages
library(ape)
library(dplyr)
library(picante)
library(MuMIn)
library(geiger)
library(ncf)

# Load functions
source("R/utils.R")
source("R/HighstatLibV6.R")

#=====================
# Load and format data
#=====================

# Load community data
l_mamm <- read.csv("Data/large_mammal.csv")

# Change species names to match tip labels in mammal supertree
l_mamm <- match_phylo_names(l_mamm)

# Load mammal supertree
supertree <- read.nexus("Data/mammal_supertree_nexus.txt")

# Create tree of sampled large mammal taxa
large_mammal_tree <- subset_supertree(l_mamm, supertree, 7:62)

# Add columns for phylogenetic diversity and species richness - warning messages
# explain that PD could not be calculated for communities containing a single species
l_mamm <- faith_pd(l_mamm, large_mammal_tree, 7:62)

#=================
# Data exploration
#=================

# Check distribution of dependent variable
hist(l_mamm$PD)
qqnorm(l_mamm$PD)
qqline(l_mamm$PD)

# Check distributions of predictor variables
hist(l_mamm$annual_rainfall)
hist(l_mamm$soil_pc1)

# Check for outliers
dotchart(l_mamm$PD, main = "PD")
dotchart(l_mamm$annual_rainfall, main = "Rain")
dotchart(l_mamm$soil_pc1, main = "Soil")

# Check for absence of strong correlation among variables
var_mat <- cbind(l_mamm$PD, l_mamm$annual_rainfall, l_mamm$soil_pc1)
colnames(var_mat) <- c("PD", "Rain", "Soil")
pairs(var_mat, lower.panel = panel.smooth2,
      upper.panel = panel.cor, diag.panel = panel.hist)

# Check variance inflation factors < 3 (Zuur et al. 2007)
corvif(var_mat)

#=========================================
# Selecting model random effects structure
#=========================================

# This model selection protocol is adapted from the top-down strategy described
# in Zuur et al (2009) pp 121-122. First create models with the same fixed 
# effects structure (most complex structure possible) but with different
# random effects structure and compare them using AIC (not all are nested)

# Create object to be passed to lme to increase number of iterations and 
# achieve convergence (only use for models that fail to converge without this
# additional measure)
lmc <- lmeControl(niterEM = 5200, msMaxIter = 5200)

# Create models
rand1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 | ranch, method = "REML", data = l_mamm)
rand2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + annual_rainfall | ranch, 
             method = "REML", data = l_mamm)
rand3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + soil_pc1 | ranch,
             control = lmc, method = "REML", data = l_mamm)

# Determine which random effects structure is best according to AIC
AICc(rand1, rand2, rand3) # rand1 (random intercepts) is best

#==================================
# Check for spatial autocorrelation
#==================================

# Plot spatial correlogram of phylogenetic diversity measurements
fit1 <- correlog(x = l_mamm$longitude, y = l_mamm$latitude, 
                 z = l_mamm$PD, increment = 1, resamp = 500,
                 latlon = TRUE)
plot(fit1$mean.of.class, fit1$correlation, xlab = "Distance (km)", ylab = "Moran's I",
     main = "SAC in PD data")

# There appears to be some SAC among sites that are within 10 km of each other,
# but it is not extremely strong. Even among sites within 1 km, Moran's I is
# just under 0.5.

# Now create full model and check residuals for SAC
full_mod <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
                  landuse:annual_rainfall + landuse:soil_pc1,
                random = ~ 1 | ranch, method = "ML", data = l_mamm)
fit2 <- correlog(x = l_mamm$longitude, y = l_mamm$latitude,
                 z = residuals(full_mod), increment = 1, resamp = 500, latlon = TRUE)
plot(fit2$mean.of.class, fit2$correlation, xlab = "Distance (km)", 
     ylab = "Moran's I", main = "SAC in full model residuals")

# SAC appears to be accounted for by model covariates

#========================================
# Selecting model fixed effects structure
#========================================

# Create models off all combinations of fixed effects
fix1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
              landuse:soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall, 
            random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + 
              landuse:soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix5 <- lme(PD ~ landuse + annual_rainfall+ landuse:annual_rainfall,
            random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix6 <- lme(PD ~ landuse + soil_pc1 + landuse:soil_pc1,
            random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix7 <- lme(PD ~ landuse + annual_rainfall, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix8 <- lme(PD ~ landuse + soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix9 <- lme(PD ~ annual_rainfall + soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix10 <- lme(PD ~ landuse, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix11 <- lme(PD ~ annual_rainfall, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix12 <- lme(PD ~ soil_pc1, random = ~ 1 | ranch, method = "ML", data = l_mamm)
fix13 <- lme(PD ~ 1, random = ~ 1 | ranch, method = "ML", data = l_mamm)

# Calculate AICc for each model
AICc(fix1, fix2, fix3, fix4, fix5, fix6, fix7, fix8, fix9, fix10, fix11, fix12, fix13)

# Compare best model (fix1) to model with similar AICc (fix2) using
# likelihood ratio test
anova(fix2, fix1) # Significant difference, use more complex model (fix1)

# Refit final model with REML
final <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
               landuse:soil_pc1, random = ~ 1 | ranch, method = "REML", data = l_mamm)

#=================
# Model validation
#=================

# Check for linearity and equal variance across range of fitted values
plot(fitted(final), residuals(final))
abline(0,0)

# Check assumption of normality in residuals
hist(residuals(final))
qqnorm(residuals(final))
qqline(residuals(final))

# Check for relationships between residuals and explanatory variables
plot(l_mamm$landuse, residuals(final))

#=====================
# Model interpretation
#=====================

# Extract coefficients from final model
summary(final)

#=================
# Create Figure 2b
#=================

# Make plot
plot(x = l_mamm$annual_rainfall, y = l_mamm$PD, type = "n",
     ylab = "PD", xlab = "Annual rainfall (mm)",
     ylim = c(0, 500), yaxt = "n", yaxs = "i",
     xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(0,500,100), labels = seq(0,500,100), las = 1,
     cex.axis = 0.8)
axis(side = 1, at = seq(450,750,50), labels = seq(450,750,50),
     cex.axis = 0.8)

# Create land-use type subsets so that separate rainfall models can be made
CN <- l_mamm[l_mamm$landuse == "Conserved",]
EX <- l_mamm[l_mamm$landuse == "Fenced",]
PA <- l_mamm[l_mamm$landuse == "Pastoral",]
CNmod <- lm(PD ~ annual_rainfall, data = CN)
EXmod <- lm(PD ~ annual_rainfall, data = EX)
PAmod <- lm(PD ~ annual_rainfall, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$annual_rainfall), max(CN$annual_rainfall), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(annual_rainfall = CNrain), interval = "confidence")
EXrain <- seq(min(EX$annual_rainfall), max(EX$annual_rainfall), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(annual_rainfall = EXrain), interval = "confidence")
PArain <- seq(min(PA$annual_rainfall), max(PA$annual_rainfall), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(annual_rainfall = PArain), interval = "confidence")

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

#=================
# Create Figure S3
#=================

# Create plot
plot(x = l_mamm$soil_pc1, y = l_mamm$PD, type = "n",
     ylab = "PD (million years)", xlab = "Soil degradation",
     ylim = c(0, 600), yaxt = "n", yaxs = "i")
axis(side = 2, at = seq(0,600,100), labels = seq(0,600,100), las = 1,
     cex.axis = 0.8)

# Subset data by landuse type and create separate model for each group
CN <- l_mamm[l_mamm$landuse == "Conserved",]
EX <- l_mamm[l_mamm$landuse == "Fenced",]
PA <- l_mamm[l_mamm$landuse == "Pastoral",]
CNmod <- lm(PD ~ soil_pc1, data = CN)
EXmod <- lm(PD ~ soil_pc1, data = EX)
PAmod <- lm(PD ~ soil_pc1, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$soil_pc1), max(CN$soil_pc1), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(soil_pc1 = CNrain), interval = "confidence")
EXrain <- seq(min(EX$soil_pc1), max(EX$soil_pc1), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(soil_pc1 = EXrain), interval = "confidence")
PArain <- seq(min(PA$soil_pc1), max(PA$soil_pc1), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(soil_pc1 = PArain), interval = "confidence")

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
