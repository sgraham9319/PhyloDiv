###################
# Modeling plant PD
###################

# This script creates and interprets the generalized linear mixed effects model
# of plant PD

# Load required packages
library(ape)
library(dplyr)
library(picante)
library(MuMIn)

# Load functions
source("R/utils.R")
source("R/HighstatLibV6.R")

#=====================
# Load and format data
#=====================

# Load community data
plant <- read.csv("Data/plant.csv")

# Load phylogeny
phylo <- read.tree("Data/plant_phylo")

# Add columns for phylogenetic diversity and species richness
plant <- faith_pd(plant, phylo, 4:151)

# Reorder factor levels for landuse so that conserved forms intercept in models
plant$landuse <- relevel(plant$landuse, "Conserved")

#=================
# Data exploration
#=================

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

# Check for absence of strong correlation among variables
var_mat <- cbind(plant$PD, plant$annual_rainfall, plant$soil_pc1)
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
             random = ~ 1 | pair_id, method = "REML", data = plant)
rand2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + annual_rainfall | pair_id, 
             control = lmc, method = "REML", data = plant)
rand3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 +
               landuse:annual_rainfall + landuse:soil_pc1,
             random = ~ 1 + soil_pc1 | pair_id, method = "REML", data = plant)

# Determine which random effects structure is best according to AIC
AICc(rand1, rand2, rand3) # rand1 (random intercepts) is best

#========================================
# Selecting model fixed effects structure
#========================================

# Create models off all combinations of fixed effects
fix1 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall + 
              landuse:soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix2 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + landuse:annual_rainfall, 
            random = ~ 1 | pair_id, method = "ML", data = plant)
fix3 <- lme(PD ~ landuse + annual_rainfall + soil_pc1 + 
              landuse:soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix4 <- lme(PD ~ landuse + annual_rainfall + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix5 <- lme(PD ~ landuse + annual_rainfall+ landuse:annual_rainfall,
            random = ~ 1 | pair_id, method = "ML", data = plant)
fix6 <- lme(PD ~ landuse + soil_pc1 + landuse:soil_pc1,
            random = ~ 1 | pair_id, method = "ML", data = plant)
fix7 <- lme(PD ~ landuse + annual_rainfall, random = ~ 1 | pair_id, method = "ML", data = plant)
fix8 <- lme(PD ~ landuse + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix9 <- lme(PD ~ annual_rainfall + soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix10 <- lme(PD ~ landuse, random = ~ 1 | pair_id, method = "ML", data = plant)
fix11 <- lme(PD ~ annual_rainfall, random = ~ 1 | pair_id, method = "ML", data = plant)
fix12 <- lme(PD ~ soil_pc1, random = ~ 1 | pair_id, method = "ML", data = plant)
fix13 <- lme(PD ~ 1, random = ~ 1 | pair_id, method = "ML", data = plant)

# Calculate AICc for each model
AICc(fix1, fix2, fix3, fix4, fix5, fix6, fix7, fix8, fix9, fix10, fix11, fix12, fix13)

# Compare best model (fix10) to models with similar AICc (fix7, fix8) using
# likelihood ratio tests
anova(fix10, fix7) # No difference, use simpler model (fix10)
anova(fix10, fix8) # No difference, use simpler model (fix10)

# Refit final model with REML
final <- lme(PD ~ landuse, random = ~ 1 | pair_id, method = "REML", data = plant)

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
plot(plant$landuse, residuals(final))

#=====================
# Model interpretation
#=====================

# Extract coefficients from final model
summary(final)
