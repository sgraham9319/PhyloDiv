########################
# Prevalence of turnover
########################

# This script addresses the fourth research objective of the manuscript - how 
# prevalent is turnover in community responses to land-use change? It 
# calculates phylogenetic beta-diversity for plant and small mammal communities
# across each plot pair and decomposes this value into its turnover and 
# nestedness components

# Outputs: Figure 6, Table S6

# Load required packages
library(ape)
library(dplyr)
library(geiger)

# Load functions
source("R/utils.R")
source("../Data/Leprieuretal2012.R")

# Load community data
s_mamm <- read.csv("Data/small_mammal.csv", stringsAsFactors = F)
plant <- read.csv("Data/plant.csv", stringsAsFactors = F)

# Load small mammal phylogeny
supertree <- read.nexus("Data/mammal_supertree_nexus.txt")
s_mamm_phylo <- subset_supertree(s_mamm, supertree, 4:25)

# Load plant phylogeny
plant_phylo <- read.tree("Data/plant_phylo")

# Calculate decomposed beta diversity
s_mamm_decomp <- beta_decomp(s_mamm, s_mamm_phylo, 4:25)
plant_decomp <- beta_decomp(plant, plant_phylo, 4:151)

# Combine plant and small mammal data
decomp_dat <- rbind(plant_decomp, s_mamm_decomp)
decomp_dat$group <- rep(c("Plant", "Small Mammal"), each = 3)

# Reorder rows and add x-coordinates for plotting
decomp_dat <- decomp_dat[c(1,4,2,5,3,6), ]
decomp_dat$x <- 1:6

# Create order for colors
cols <- c("green1", "violetred3", "green1", "violetred3", "green1", "violetred3")

# Create plot
plot(decomp_dat$x, decomp_dat$turnover_av, ylim = c(0, 0.7),
     pch = 19, col = cols, xlab = "Type of land use change", yaxs = "i", 
     ylab = "Beta-diversity", xaxt = "n", las = 1, xlim = c(0.5, 6.5), xaxs = "i")

# Add x-axis
axis(side = 1, at = 1:6, labels = rep(c("P", "S"), times = 3), cex.axis = 0.8,
     tck = -0.02, padj = -1)

# Add standard errors
arrows(decomp_dat$x, decomp_dat$turnover_av - decomp_dat$turnover_se, decomp_dat$x, 
       decomp_dat$turnover_av + decomp_dat$turnover_se, code = 0, col = cols, lwd = 1.5)

# Add land-use type labels below x-axis
mtext(text = c("Agriculture", "Fenced", "Pastoral"), side = 1, line = 1.5, at = c(1.5,3.5,5.5))

# Add lines to show mean total beta-diversity for each category
arrows(decomp_dat$x - 0.2, decomp_dat$total_av, decomp_dat$x + 0.2, 
       decomp_dat$total_av, code = 0, col = cols, lwd = 1.5)

# Add lines to separate land-use types
abline(v = 2.5)
abline(v = 4.5)
