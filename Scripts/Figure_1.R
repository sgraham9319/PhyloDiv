
# Load functions file
source("R/utils.R")

#=========================================
# Load community datasets and calculate PD
#=========================================

plant <- read.csv("Data/plant.csv", stringsAsFactors = F)
small_mammal <- read.csv("Data/small_mammal.csv", stringsAsFactors = F)
large_mammal <- read.csv("Data/large_mammal.csv", stringsAsFactors = F)

# Load phylogenies
plant_phylo <- read.tree("Data/plant_phylo")
mammal_supertree <- read.nexus("Data/mammal_supertree_nexus.txt")

# Create phylogenies of sampled small mammal and large mammal taxa
small_mammal_phylo <- subset_supertree(small_mammal, mammal_supertree, 4:25)
large_mammal_phylo <- subset_supertree(large_mammal, mammal_supertree, 7:62)

# Calculate phylogenetic diversity
plant <- faith_pd(plant, plant_phylo, 4:151)
small_mammal <- faith_pd(small_mammal, small_mammal_phylo, 4:25)
large_mammal <- faith_pd(large_mammal, large_mammal_phylo, 7:62)


# 
org_group <- rep(c("plant", "s_mamm", "l_mamm"),
                 times = c(nrow(plant), nrow(small_mammal), nrow(large_mammal)))
landuse <- c(plant$landuse, small_mammal$landuse, large_mammal$landuse)
PD <- c(plant$PD, small_mammal$PD, large_mammal$PD)
comb_dat <- data.frame(org_group, landuse, PD)

comb_dat$treatment <- paste(comb_dat$org_group, comb_dat$landuse, sep = "_")

comb_dat$treatment <- factor(comb_dat$treatment, levels = c("plant_Conserved",
                                                            "s_mamm_Conserved",
                                                            "l_mamm_Conserved",
                                                            "plant_Agriculture",
                                                            "s_mamm_Agriculture", 
                                                            "plant_Fenced",
                                                            "s_mamm_Fenced",
                                                            "l_mamm_Fenced",
                                                            "plant_Pastoral",
                                                            "s_mamm_Pastoral",
                                                            "l_mamm_Pastoral"))

# Create the plot
boxplot(PD ~ treatment, at = c(1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 14), outline = F,
        ylim = c(0,2150),
        col = c("green1", "violetred3", "blue", "green1", "violetred3",
                "green1", "violetred3", "blue", "green1", "violetred3", "blue"),
        xaxt = "n", yaxt = "n",
        ylab = "PD (billion years)",
        xlab = "Land-use type", whisklty = "solid", 
        data = comb_dat)
axis(side = 1, at = c(1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 14), 
     labels = c("P", "S", "L", "P", "S", "P", "S", "L", "P", "S", "L"),
     cex.axis = 0.8, tck = -0.05, padj = -1)
axis(side = 2, at = seq(0, 2000, 500), labels = seq(0, 2, 0.5), cex.axis = 1, las = 1)
mtext(text = c("Conserved", "Agriculture", "Fenced", "Pastoral"),
      side = 1, line = 1.5, at = c(2, 5.5, 9, 13))
abline(v = 4)
abline(v = 7)
abline(v = 11)
