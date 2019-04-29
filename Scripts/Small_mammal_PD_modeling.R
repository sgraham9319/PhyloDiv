##########################
# Modeling small mammal PD
##########################

# This script creates and interprets the generalized linear mixed effects model
# of small mammal PD

# Outputs: Figure 2a

# Load required packages
library(ape)
library(dplyr)
library(picante)
library(MuMIn)
library(geiger)

# Load functions
source("R/utils.R")
source("R/HighstatLibV6.R")

#=====================
# Load and format data
#=====================

# Load community data
s_mamm <- read.csv("Data/small_mammal.csv")

# Load mammal supertree
supertree <- read.nexus("Data/mammal_supertree_nexus.txt")

# Create tree of sampled small mammal taxa
small_mammal_tree <- subset_supertree(s_mamm, supertree, 4:25)

# Add columns for phylogenetic diversity and species richness - warning messages
# explain that PD could not be calculated for communities containing a single species
s_mamm <- faith_pd(s_mamm, small_mammal_tree, 4:25)

# Reorder factor levels for landuse so that conserved forms intercept in models
s_mamm$landuse <- relevel(s_mamm$landuse, "Conserved")

#=================
# Data exploration
#=================

# Check distribution of dependent variable
hist(s_mamm$PD)
qqnorm(s_mamm$PD)
qqline(s_mamm$PD)

# Check distributions of predictor variables
hist(s_mamm$annual_rainfall)
hist(s_mamm$soil_pc1)

# Check for outliers
dotchart(s_mamm$PD, main = "PD")
dotchart(s_mamm$annual_rainfall, main = "Rain")
dotchart(s_mamm$soil_pc1, main = "Soil")

# Check for absence of strong correlation among variables
var_mat <- cbind(s_mamm$PD, s_mamm$annual_rainfall, s_mamm$soil_pc1)
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
             random = ~ 1 | pair_id, method = "REML", data = s_mamm)
rand2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + annual_rainfall | pair_id, 
             control = lmc, method = "REML", data = s_mamm)
rand3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + soil_pc1 | pair_id, method = "REML", data = s_mamm)

# Determine which random effects structure is best according to AIC
AICc(rand1, rand2, rand3) # rand1 (random intercepts) is best

#========================================
# Selecting model fixed effects structure
#========================================

# Create models off all combinations of fixed effects
fix1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
              landuse:soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall, 
            random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + 
              landuse:soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix5 <- lme(PD ~ landuse + annual_rainfall+ landuse:annual_rainfall,
            random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix6 <- lme(PD ~ landuse + soil_pc1 + landuse:soil_pc1,
            random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix7 <- lme(PD ~ landuse + annual_rainfall, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix8 <- lme(PD ~ landuse + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix9 <- lme(PD ~ annual_rainfall + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix10 <- lme(PD ~ landuse, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix11 <- lme(PD ~ annual_rainfall, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix12 <- lme(PD ~ soil_pc1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)
fix13 <- lme(PD ~ 1, random = ~ 1 | pair_id, method = "ML", data = s_mamm)

# Calculate AICc for each model
AICc(fix1, fix2, fix3, fix4, fix5, fix6, fix7, fix8, fix9, fix10, fix11, fix12, fix13)

# Compare best model (fix5) to more complex model with similar AICc (fix7) using
# likelihood ratio test
anova(fix7, fix5) # Significant difference, use more complex model (fix5)

# Refit final model (fix5) with REML
final <- lme(PD ~ landuse + annual_rainfall+ landuse:annual_rainfall,
             random = ~ 1 | pair_id, method = "REML", data = s_mamm)

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
plot(s_mamm$landuse, residuals(final))

#=====================
# Model interpretation
#=====================

# Extract coefficients from final model
summary(final)

#=================
# Create Figure 2a
#=================

# Make plot
plot(x = s_mamm$annual_rainfall, y = s_mamm$PD, type = "n",
     ylab = "PD", xlab = "Annual rainfall (mm)",
     ylim = c(0, 500), yaxt = "n",
     xlim = c(420, 750), xaxt = "n", xaxs = "i",
     main = "Conserved = blue, Ag = red, \nFence = yellow, Pastoral = orange")
axis(side = 2, at = seq(0,500,100), labels = seq(0,500,100), las = 1,
     cex.axis = 0.8)
axis(side = 1, at = seq(450,750,50), labels = seq(450,750,50),
     cex.axis = 0.8)

# Create land-use type subsets so that separate rainfall models can be made
CN <- s_mamm[s_mamm$landuse == "Conserved",]
AG <- s_mamm[s_mamm$landuse == "Agriculture",]
EX <- s_mamm[s_mamm$landuse == "Fenced",]
PA <- s_mamm[s_mamm$landuse == "Pastoral",]
CNmod <- lm(PD ~ annual_rainfall, data = CN)
AGmod <- lm(PD ~ annual_rainfall, data = AG)
EXmod <- lm(PD ~ annual_rainfall, data = EX)
PAmod <- lm(PD ~ annual_rainfall, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$annual_rainfall), max(CN$annual_rainfall), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(annual_rainfall = CNrain), interval = "confidence")
AGrain <- seq(min(AG$annual_rainfall), max(AG$annual_rainfall), length.out = 50)
AGconf <- predict(AGmod, newdata = data.frame(annual_rainfall = CNrain), interval = "confidence")
EXrain <- seq(min(EX$annual_rainfall), max(EX$annual_rainfall), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(annual_rainfall = EXrain), interval = "confidence")
PArain <- seq(min(PA$annual_rainfall), max(PA$annual_rainfall), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(annual_rainfall = PArain), interval = "confidence")

# Add confidence intervals to plot
polygon(c(rev(CNrain), CNrain), c(rev(CNconf[,3]), CNconf[,2]),
        col = makeTransparent("paleturquoise3", alpha = 70), border = NA)
polygon(c(rev(AGrain), AGrain), c(rev(AGconf[,3]), AGconf[,2]),
        col = makeTransparent("red", alpha = 70), border = NA)
polygon(c(rev(EXrain), EXrain), c(rev(EXconf[,3]), EXconf[,2]),
        col = makeTransparent("yellow2"), border = NA)
polygon(c(rev(PArain), PArain), c(rev(PAconf[,3]), PAconf[,2]),
        col = makeTransparent("tan2"), border = NA)

# Add lines to delimit confidence intervals
lines(CNrain, CNconf[,1], col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,2], lty=2, col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,3], lty=2, col = "skyblue2", lwd = 2)

lines(AGrain, AGconf[,1], col = "red", lwd = 2)
lines(AGrain, AGconf[,2], lty=2, col = "red", lwd = 2)
lines(AGrain, AGconf[,3], lty=2, col = "red", lwd = 2)

lines(EXrain, EXconf[,1], col = "gold1", lwd = 2)
lines(EXrain, EXconf[,2], lty=2, col = "gold1", lwd = 2)
lines(EXrain, EXconf[,3], lty=2, col = "gold1", lwd = 2)

lines(PArain, PAconf[,1], col = "tan2", lwd = 2)
lines(PArain, PAconf[,2], lty=2, col = "tan2", lwd = 2)
lines(PArain, PAconf[,3], lty=2, col = "tan2", lwd = 2)