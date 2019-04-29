####################################
# Evaluation of habitat associations
####################################

# This script evaluates habitat associations for all taxa. It determines whether
# each taxon exhibited a significant aversion/affinity to any land-use type and
# also estimates phylogenetic signal in the habitat association trait

# Outputs: Figure 3, Figure 4, Figure 5, Table 2, Table S4, Table S5

# Load required packages
library(dplyr)
library(tidyr)
library(ape)
library(phytools)
library(picante)
library(geiger)

# Load functions
source("R/utils.R")
source("R/habitat_preference_functions.R")

#==========
# 1. Plants
#==========

#-------------------------------------
# Calculating aversions and affinities
#-------------------------------------

# Load community data
plant <- read.csv("Data/plant.csv")

# Summarize abundance data by landuse for each species
abundance <- abund_summ(plant)

# Identify rare species and exclude from habitat preference evaluation
rare_sps <- names(which(apply(plant[, 1:nrow(abundance[[1]]) + 3], 2, sum) < 30))
rare_sps_rows <- which(abundance[[1]]$species %in% rare_sps)
abundance[[1]][, 2:4][rare_sps_rows, ] <- 0
abundance[[2]][, 2:4][rare_sps_rows, ] <- 0
abundance[[3]][, 2:4][rare_sps_rows, ] <- 0

# Estimate habitat preference trait
hab_prefs <- hab_pref_summ(abundance[[1]], abundance[[2]], abundance[[3]])

# Load plant phylogeny
phylo <- read.tree("Data/plant_phylo")

# Calculate phylogenetic signal - add to Table 2
signal_summary <- phy_sig_summ(hab_prefs, phylo)

# Identify significant landuse affinities and aversions
plant_affin <- affin_signif(abundance[[1]], abundance[[2]], abundance[[3]])

# Change affinities of rarely sampled species to "Rare"
plant_affin[which(plant_affin[, "species"] %in% rare_sps), 2:9] <- "RARE"

# Reorder aversion data so they match the order of tip labels in phylogeny
plant_affin <- plant_affin[match(phylo$tip.label, plant_affin[,"species"]),]


#---------------------
# Creating phylogenies
#---------------------

# Plot phylogeny without tip labels - colors will be added in Illustrator
plot(phylo, type = "fan", show.tip.label = F, no.margin = T)

# Create table of affinity data to be consulted when adding colored branches
# in Adobe Illustrator
affin.table <- as.data.frame(plant_affin[, c(1,3,5,7,9)])
affin.table <- affin.table %>% mutate_all(as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  if(sum(affin.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(affin.table[i,])) == 4){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
affin.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(affin.table)){
  affin.table[,i][is.na(affin.table[, i])] <- "-"
  affin.table[,i][which(affin.table[, i] == "RARE")] <- "-"
  affin.table[,i][which(affin.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
affin.table[, 2][which(affin.table[, 2] == "TRUE")] <- "BLUE"
affin.table[, 3][which(affin.table[, 3] == "TRUE")] <- "RED"
affin.table[, 4][which(affin.table[, 4] == "TRUE")] <- "YELLOW"
affin.table[, 5][which(affin.table[, 5] == "TRUE")] <- "ORANGE"

# Repeat to create aversion table
aver.table <- as.data.frame(plant_affin[, c(1,2,4,6,8)])
aver.table <- aver.table %>% mutate_all(as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  if(sum(aver.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(aver.table[i,])) == 4){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
aver.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(aver.table)){
  aver.table[,i][is.na(aver.table[, i])] <- "-"
  aver.table[,i][which(aver.table[, i] == "RARE")] <- "-"
  aver.table[,i][which(aver.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
aver.table[, 2][which(aver.table[, 2] == "TRUE")] <- "BLUE"
aver.table[, 3][which(aver.table[, 3] == "TRUE")] <- "RED"
aver.table[, 4][which(aver.table[, 4] == "TRUE")] <- "YELLOW"
aver.table[, 5][which(aver.table[, 5] == "TRUE")] <- "ORANGE"

# Summarize aversion and affinity data for table S3
num_taxa <- sum(plant_affin[, 2] != "RARE")
raw_nums <- rep(NA, times = ncol(plant_affin))
for(i in 1:ncol(plant_affin)){
  raw_nums[i] <- sum(plant_affin[, i] == "TRUE", na.rm = T)
}
tableS3 <- data.frame(colnames(plant_affin), raw_nums)
tableS3$pct <- (raw_nums/num_taxa) * 100

# Calculate percentage of non-rare taxa showing affinity to at least one
# disturbed land-use
num_dist_affin <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  num_dist_affin[i] <- sum(affin.table[i, 3:5] == "-")
}
100 * (sum(num_dist_affin < 3) / 70) # Add value to table S3

# Repeat for aversion
num_dist_aver <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  num_dist_aver[i] <- sum(aver.table[i, 3:5] == "-")
}
100 * (sum(num_dist_aver < 3) / 70) # Add value to table S3

#=================
# 2. Small mammals
#=================

#-------------------------------------
# Calculating aversions and affinities
#-------------------------------------

# Load community data
s_mamm <- read.csv("Data/small_mammal.csv")

# Round abundances to whole numbers
s_mamm[, 4:25] <- round(s_mamm[, 4:25])

# Summarize abundance data by landuse for each species
abundance <- abund_summ(s_mamm)

# Identify rare species and exclude from habitat preference evaluation
rare_sps <- names(which(apply(s_mamm[, 1:nrow(abundance[[1]]) + 3], 2, sum) < 30))
rare_sps_rows <- which(abundance[[1]]$species %in% rare_sps)
abundance[[1]][, 2:4][rare_sps_rows, ] <- 0
abundance[[2]][, 2:4][rare_sps_rows, ] <- 0
abundance[[3]][, 2:4][rare_sps_rows, ] <- 0

# Estimate habitat preference trait
hab_prefs <- hab_pref_summ(abundance[[1]], abundance[[2]], abundance[[3]])

# Load mammal supertree
supertree <- read.nexus("Data/mammal_supertree_nexus.txt")

# Create tree of sampled small mammal taxa
phylo <- subset_supertree(s_mamm, supertree, 4:25)

# Calculate phylogenetic signal - add to Table 2
signal_summary <- phy_sig_summ(hab_prefs, phylo)

# Identify significant landuse affinities and aversions
s_mamm_affin <- affin_signif(abundance[[1]], abundance[[2]], abundance[[3]])

# Change affinities of rarely sampled species to "Rare"
s_mamm_affin[which(s_mamm_affin[, "species"] %in% rare_sps), 2:9] <- "RARE"

# Reorder aversion data so they match the order of tip labels in phylogeny
s_mamm_affin <- s_mamm_affin[match(phylo$tip.label, s_mamm_affin[, "species"]),]

#---------------------
# Creating phylogenies
#---------------------

# Plot phylogeny - colors will be added in Illustrator
plot(phylo, type = "fan", show.tip.label = T, no.margin = T)

# Create table of affinity data to be consulted when adding colored branches
# in Illustrator
affin.table <- as.data.frame(s_mamm_affin[, c(1,3,5,7,9)])
affin.table <- affin.table %>% mutate_all(as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  if(sum(affin.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(affin.table[i,])) == 4){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
affin.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(affin.table)){
  affin.table[,i][is.na(affin.table[, i])] <- "-"
  affin.table[,i][which(affin.table[, i] == "RARE")] <- "-"
  affin.table[,i][which(affin.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
affin.table[, 2][which(affin.table[, 2] == "TRUE")] <- "BLUE"
affin.table[, 3][which(affin.table[, 3] == "TRUE")] <- "RED"
affin.table[, 4][which(affin.table[, 4] == "TRUE")] <- "YELLOW"
affin.table[, 5][which(affin.table[, 5] == "TRUE")] <- "ORANGE"

# Repeat to create aversion table
aver.table <- as.data.frame(s_mamm_affin[, c(1,2,4,6,8)])
aver.table <- aver.table %>% mutate_all(as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  if(sum(aver.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(aver.table[i,])) == 4){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
aver.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(aver.table)){
  aver.table[,i][is.na(aver.table[, i])] <- "-"
  aver.table[,i][which(aver.table[, i] == "RARE")] <- "-"
  aver.table[,i][which(aver.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
aver.table[, 2][which(aver.table[, 2] == "TRUE")] <- "BLUE"
aver.table[, 3][which(aver.table[, 3] == "TRUE")] <- "RED"
aver.table[, 4][which(aver.table[, 4] == "TRUE")] <- "YELLOW"
aver.table[, 5][which(aver.table[, 5] == "TRUE")] <- "ORANGE"

# Summarize aversion and affinity data for table S3
num_taxa <- sum(s_mamm_affin[, 2] != "RARE")
raw_nums <- rep(NA, times = ncol(s_mamm_affin))
for(i in 1:ncol(s_mamm_affin)){
  raw_nums[i] <- sum(s_mamm_affin[, i] == "TRUE", na.rm = T)
}
tableS3 <- data.frame(colnames(s_mamm_affin), raw_nums)
tableS3$pct <- (raw_nums/num_taxa) * 100

# Calculate percentage of non-rare taxa showing affinity to at least one
# disturbed land-use
num_dist_affin <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  num_dist_affin[i] <- sum(affin.table[i, 3:5] == "-")
}
100 * (sum(num_dist_affin < 3) / num_taxa) # Add value to table S3

# Repeat for aversion
num_dist_aver <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  num_dist_aver[i] <- sum(aver.table[i, 3:5] == "-")
}
100 * (sum(num_dist_aver < 3) / num_taxa) # Add value to table S3

#=================
# 3. Large mammals
#=================

#-------------------------------------
# Calculating aversions and affinities
#-------------------------------------

# Load community data
l_mamm <- read.csv("Data/large_mammal.csv")

# Change species names to match tip labels in mammal supertree
l_mamm <- match_phylo_names(l_mamm)

# Remove outlier for cow abundance
l_mamm$Bos_taurus[which(l_mamm$Bos_taurus > 1000)] <- NA

# Round abundances to whole numbers
l_mamm[, 7:62] <- round(l_mamm[, 7:62])

# Summarize abundance data by landuse for each species
abundance <- abund_summ_lm(l_mamm)

# Identify rare species and exclude from habitat preference evaluation
rare_sps_rows <- which(prop.table(abundance$total) < 0.001)
abundance[, 2:5][rare_sps_rows, ] <- 0

# Estimate habitat preference trait
hab_prefs <- hab_pref_summ_lm(abundance)

# Create tree of sampled large mammal taxa
phylo <- subset_supertree(l_mamm, supertree, 7:62)

# Calculate phylogenetic signal - add to Table 2
signal_summary <- phy_sig_summ(hab_prefs, phylo)

# Identify significant landuse affinities and aversions
l_mamm_affin <- affin_signif_lm(abundance)

# Change affinities of rarely sampled species to "Rare"
l_mamm_affin[rare_sps_rows, 2:7] <- "RARE"

# Reorder aversion data so they match the order of tip labels in phylogeny
l_mamm_affin <- l_mamm_affin[match(phylo$tip.label, l_mamm_affin[, "species"]),]

#---------------------
# Creating phylogenies
#---------------------

# Plot phylogeny - colors will be added in Illustrator
plot(phylo, type = "fan", show.tip.label = F, no.margin = T)

# Create table of affinity data to be consulted when adding colored branches
# in Illustrator
affin.table <- as.data.frame(l_mamm_affin[, c(1,3,5,7)])
affin.table <- affin.table %>% mutate_all (as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  if(sum(affin.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(affin.table[i,])) == 3){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
affin.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(affin.table)){
  affin.table[,i][is.na(affin.table[, i])] <- "-"
  affin.table[,i][which(affin.table[, i] == "RARE")] <- "-"
  affin.table[,i][which(affin.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
affin.table[, 2][which(affin.table[, 2] == "TRUE")] <- "BLUE"
affin.table[, 3][which(affin.table[, 3] == "TRUE")] <- "YELLOW"
affin.table[, 4][which(affin.table[, 4] == "TRUE")] <- "ORANGE"

# Repeat to create aversion table
aver.table <- as.data.frame(l_mamm_affin[, c(1,2,4,6)])
aver.table <- aver.table %>% mutate_all (as.character)

# Add column for no affinity/aversion
total.false.na <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  if(sum(aver.table[i,] == "FALSE", na.rm = T) +
     sum(is.na(aver.table[i,])) == 3){
    total.false.na[i] <- "GREY"
  } else{
    total.false.na[i] <- "-"
  }
}
aver.table$None <- total.false.na

# Replace RARE, NA and FALSE with - (no lines need to be added for these)
for(i in 2:ncol(aver.table)){
  aver.table[,i][is.na(aver.table[, i])] <- "-"
  aver.table[,i][which(aver.table[, i] == "RARE")] <- "-"
  aver.table[,i][which(aver.table[, i] == "FALSE")] <- "-"
}

# Replace "TRUE" with appropriate color
aver.table[, 2][which(aver.table[, 2] == "TRUE")] <- "BLUE"
aver.table[, 3][which(aver.table[, 3] == "TRUE")] <- "YELLOW"
aver.table[, 4][which(aver.table[, 4] == "TRUE")] <- "ORANGE"

# Summarize aversion and affinity data for table S3
num_taxa <- sum(l_mamm_affin[, 2] != "RARE")
raw_nums <- rep(NA, times = ncol(l_mamm_affin))
for(i in 1:ncol(l_mamm_affin)){
  raw_nums[i] <- sum(l_mamm_affin[, i] == "TRUE", na.rm = T)
}
tableS3 <- data.frame(colnames(l_mamm_affin), raw_nums)
tableS3$pct <- (raw_nums/num_taxa) * 100

# Calculate percentage of non-rare taxa showing affinity to at least one
# disturbed land-use
num_dist_affin <- rep(NA, times = nrow(affin.table))
for(i in 1:nrow(affin.table)){
  num_dist_affin[i] <- sum(affin.table[i, 3:4] == "-")
}
100 * (sum(num_dist_affin < 2) / num_taxa) # Add value to table S3

# Repeat for aversion
num_dist_aver <- rep(NA, times = nrow(aver.table))
for(i in 1:nrow(aver.table)){
  num_dist_aver[i] <- sum(aver.table[i, 3:4] == "-")
}
100 * (sum(num_dist_aver < 2) / num_taxa) # Add value to table S3

#=====================
# 4. Creating Table S4
#=====================

# Add ag_aver and ag_pref columns to large mammal data
l_mamm_affin <- cbind(l_mamm_affin, rep(NA, times = nrow(l_mamm_affin)), rep(NA, times = nrow(l_mamm_affin)))

# Reorder columns of large mammal data to match plants and small mammals
l_mamm_affin <- l_mamm_affin[,c(1:3, 8, 9, 4:7)]

# Rename columns
colnames(l_mamm_affin) <- c("species", "con_aver", "con_pref", "ag_aver", "ag_pref", "fen_aver", "fen_pref", "pas_aver", "pas_pref")

# Combine three taxonomic groups
affin <- rbind(plant_affin, s_mamm_affin, l_mamm_affin)
affin <- data.frame(affin)

# Include real name for large mammal taxa as congenerics were used when the sampled
# species did not occur in the tree
l_mamm_names <- read.csv("Data/large_mammal_species_key.csv")
real_name_lm <- as.character(l_mamm_names$binomial[match(l_mamm_affin[, "species"], l_mamm_names$representative_in_supertree)])
real_name_plant <- plant_affin[,"species"]
real_name_sm <- s_mamm_affin[, "species"]
affin$real_name <- c(real_name_plant, real_name_sm, real_name_lm)

# Change name for single small mammal species that had different name in tree
# than community data
affin$real_name[which(affin$real_name == "Crocidura_elgonius")] <- "Crocidura spp."

# Add taxonomic group identifier
affin$taxon <- c(rep("Plant", times = nrow(plant_affin)), 
                 rep("Small mammal", times = nrow(s_mamm_affin)),
                 rep("Large mammal", times = nrow(l_mamm_affin)))

# Remove underscores from species names and capitalize genus names
affin$real_name <- gsub(x = affin$real_name, pattern = "_", replacement = " ")
substr(affin$real_name, 1, 1) <- toupper(substr(affin$real_name, 1, 1))

# Add family information to affinities data
family <- read.csv("Data/family_info_all_taxa.csv")
affin$family <- as.character(family$family[match(affin$real_name, family$species)])

# Reorder columns
tableS2 <- affin[,c(11, 10, 12, 2:9)]

# Rename species column
colnames(tableS2)[2] <- "species"
# Significant preferences and aversions were summarized manually in Excel