#####################################
# Determination of land-use affinites
#####################################

# This script is for the updated (9/16/2018) calculation of land-use affinities
# following reviewer comments pointing out that we should highlight
# any preferences for particular land-uses as well as any aversions. If we
# only look for aversions we are ensuring that our conclusion can only be that
# land-use change is negative, or at best neutral.

# Load required packages
library(dplyr)
library(tidyr)
library(ape)
library(geiger)
library(phytools)
library(picante)
source("R/habitat_preference_functions.R")
source("R/utils.R")

#==========
# 1. Plants
#==========

#-------------------------------------
# Calculating aversions and affinities
#-------------------------------------
# Load site data
plant <- read.csv("Data/plant.csv")

# Create subsets for different disturbed land-use types (and their paired 
# conserved sites)
ag_pairs <- plant$pair_id[which(plant$landuse == "Agriculture")]
pas_pairs <- plant$pair_id[which(plant$landuse == "Pastoral")]
fen_pairs <- plant$pair_id[which(plant$landuse == "Fenced")]
ag <- droplevels(plant[which(plant$pair_id %in% ag_pairs),])
pas <- droplevels(plant[which(plant$pair_id %in% pas_pairs),])
fen <- droplevels(plant[which(plant$pair_id %in% fen_pairs),])

# Summarize abundance data by land-use type for each species
long_ag <- ag %>%
  select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
ag_abund <- long_ag %>% group_by(species) %>%
  summarise(conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            agriculture = tapply(abund, landuse, FUN = sum, na.rm = T)["Agriculture"],
            total = sum(abund, na.rm = T))

long_pas <- pas %>%
  select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
pas_abund <- long_pas %>% group_by(species) %>%
  summarise(conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            pastoral = tapply(abund, landuse, FUN = sum, na.rm = T)["Pastoral"],
            total = sum(abund, na.rm = T))

long_fen <- fen %>%
  select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
fen_abund <- long_fen %>% group_by(species) %>%
  summarise(conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            fenced = tapply(abund, landuse, FUN = sum, na.rm = T)["Fenced"],
            total = sum(abund, na.rm = T))

# Create matrix to store affinity data
plant_affin <- matrix(NA, nrow = nrow(ag_abund), ncol = 9)
colnames(plant_affin) <- c("species", "con_aver", "con_pref", "ag_aver", "ag_pref", "fen_aver", "fen_pref", "pas_aver", "pas_pref")
plant_affin[, "species"] <- ag_abund$species

for(species in 1:nrow(ag_abund)){
  resamples <- matrix(NA, nrow = 999, ncol = 6)
  colnames(resamples) <- c("con_a", "ag", "con_f", "fen", "con_p", "pas")
  for(i in 1:nrow(resamples)) {
    resampA <- table(sample(c("con", "ag"), ag_abund$total[species], replace = T))
    resampF <- table(sample(c("con", "fen"), fen_abund$total[species], replace = T))
    resampP <- table(sample(c("con", "pas"), pas_abund$total[species], replace = T))
    resamples[i, "con_a"] <- resampA["con"]
    resamples[i, "ag"] <- resampA["ag"]
    resamples[i, "con_f"] <- resampF["con"]
    resamples[i, "fen"] <- resampF["fen"]
    resamples[i, "con_p"] <- resampP["con"]
    resamples[i, "pas"] <- resampP["pas"]
  }
  if(sum((ag_abund[species, "conserved"] < quantile(resamples[,"con_a"], probs = 0.025, na.rm = T)), 
         (fen_abund[species, "conserved"] < quantile(resamples[,"con_f"], probs = 0.025, na.rm = T)),
         (pas_abund[species, "conserved"] < quantile(resamples[,"con_p"], probs = 0.025, na.rm = T)), na.rm = T)
     >= 2){plant_affin[species, "con_aver"] <- TRUE}
  else{plant_affin[species, "con_aver"] <- FALSE}
  if(sum((ag_abund[species, "conserved"] > quantile(resamples[,"con_a"], probs = 0.975, na.rm = T)), 
         (fen_abund[species, "conserved"] > quantile(resamples[,"con_f"], probs = 0.975, na.rm = T)),
         (pas_abund[species, "conserved"] > quantile(resamples[,"con_p"], probs = 0.975, na.rm = T)), na.rm = T)
     >= 2){plant_affin[species, "con_pref"] <- TRUE}
  else{plant_affin[species, "con_pref"] <- FALSE}  
  plant_affin[species, "ag_aver"] <- ag_abund[species, "agriculture"] < quantile(resamples[,"ag"], probs = 0.025, na.rm = T)
  plant_affin[species, "ag_pref"] <- ag_abund[species, "agriculture"] > quantile(resamples[,"ag"], probs = 0.975, na.rm = T)
  plant_affin[species, "fen_aver"] <- fen_abund[species, "fenced"] < quantile(resamples[,"fen"], probs = 0.025, na.rm = T)
  plant_affin[species, "fen_pref"] <- fen_abund[species, "fenced"] > quantile(resamples[,"fen"], probs = 0.975, na.rm = T)
  plant_affin[species, "pas_aver"] <- pas_abund[species, "pastoral"] < quantile(resamples[,"pas"], probs = 0.025, na.rm = T)
  plant_affin[species, "pas_pref"] <- pas_abund[species, "pastoral"] > quantile(resamples[,"pas"], probs = 0.975, na.rm = T)
} # Note that NA values occur when a species was never observed in that land-use type

# Change affinities of rarely sampled species to "Rare"
rare_sps <- colnames(plant[4:(nrow(plant_affin) + 3)])[which(apply(plant[, 4:(nrow(plant_affin) + 3)], MARGIN = 2, FUN = sum) < 30)]
plant_affin[which(plant_affin[, "species"] %in% rare_sps), 2:9] <- "RARE"

# Load plant phylogeny
phylo <- read.tree("Data/PlantPhylo")

# Reorder aversion data so they match the order of tip labels in phylogeny
plant_affin <- plant_affin[match(phylo$tip.label, plant_affin[,"species"]),]

plot(phylo, type = "fan", cex = 0.5)

trait1 <- plant_affin[, "con_aver"]
trait1[which(trait1 == "RARE")] <- NA
trait1[which(trait1 == "TRUE")] <- 1
trait1[which(trait1 == "FALSE")] <- 0
trait1 <- as.numeric(trait1)

my.vcv <- vcv.phylo(phylo)
inv.vcv <- solve(my.vcv)

missing <- which(is.na(trait1))
trait1 <- trait1[-missing]
my.vcv1 <- my.vcv[-missing, -missing]
inv.vcv1 <- solve(my.vcv1)

sum(inv.vcv1 %*% as.matrix(trait1)) / sum(inv.vcv1)


#---------------------
# Creating phylogenies
#---------------------

# Plot phylogeny without tip labels - colors will be added in Illustrator
plot(phylo, type = "fan", show.tip.label = F, no.margin = T)

# Create table of affinity data to be consulted when adding colored branches
# in Illustrator
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

# Write affinities and aversions tables to csv
write.csv(affin.table, file = "../Data/Plot data/Plant_affinities.csv", row.names = F)
write.csv(aver.table, file = "../Data/Plot data/Plant_aversions.csv", row.names = F)

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
raw_site <- read.csv("Data/small_mammal.csv")

# Round abundances to whole numbers
raw_site[, 4:25] <- round(raw_site[, 4:25])

# Create subsets for different disturbed land-use types (and their paired 
# conserved sites)
ag_pairs <- raw_site$pair_id[which(raw_site$landuse == "Agriculture")]
pas_pairs <- raw_site$pair_id[which(raw_site$landuse == "Pastoral")]
fen_pairs <- raw_site$pair_id[which(raw_site$landuse == "Fenced")]
ag <- droplevels(raw_site[which(raw_site$pair_id %in% ag_pairs),])
pas <- droplevels(raw_site[which(raw_site$pair_id %in% pas_pairs),])
fen <- droplevels(raw_site[which(raw_site$pair_id %in% fen_pairs),])

# Summarize abundance data by landuse for each species
long_ag <- ag %>%
  select(-longitude, -latitude, -annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
ag_abund <- long_ag %>% group_by(species) %>%
  summarise(Conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            Agriculture = tapply(abund, landuse, FUN = sum, na.rm = T)["Agriculture"],
            total = sum(abund, na.rm = T))

long_pas <- pas %>%
  select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
pas_abund <- long_pas %>% group_by(species) %>%
  summarise(Conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            Pastoral = tapply(abund, landuse, FUN = sum, na.rm = T)["Pastoral"],
            total = sum(abund, na.rm = T))

long_fen <- fen %>%
  select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
         -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
         -soil_pc1, -soil_pc2) %>%
  gather(species, abund, -site, -landuse, -pair_id)
fen_abund <- long_fen %>% group_by(species) %>%
  summarise(Conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
            Fenced = tapply(abund, landuse, FUN = sum, na.rm = T)["Fenced"],
            total = sum(abund, na.rm = T))

# Estimate habitat preference trait
hab_prefs <- hab_pref_summ(ag_abund, fen_abund, pas_abund)

# Load mammal supertree
supertree <- read.nexus("../Data/Mammal.supertree.nexus.txt")

# Create tree of sampled small mammal taxa
phylo <- subset_supertree(raw_site, supertree, 4:25)

# Calculate phylogenetic signal
signal_summary <- phy_sig_summ(hab_prefs, phylo)











# Create matrix to store affinity data
s_mamm_affin <- matrix(NA, nrow = nrow(ag_abund), ncol = 9)
colnames(s_mamm_affin) <- c("species", "con_aver", "con_pref", "ag_aver", "ag_pref", "fen_aver", "fen_pref", "pas_aver", "pas_pref")
s_mamm_affin[, "species"] <- ag_abund$species

for(species in 1:nrow(ag_abund)){
  resamples <- matrix(NA, nrow = 999, ncol = 6)
  colnames(resamples) <- c("con_a", "ag", "con_f", "fen", "con_p", "pas")
  for(i in 1:nrow(resamples)) {
    resampA <- table(sample(c("con", "ag"), ag_abund$total[species], replace = T))
    resampF <- table(sample(c("con", "fen"), fen_abund$total[species], replace = T))
    resampP <- table(sample(c("con", "pas"), pas_abund$total[species], replace = T))
    resamples[i, "con_a"] <- resampA["con"]
    resamples[i, "ag"] <- resampA["ag"]
    resamples[i, "con_f"] <- resampF["con"]
    resamples[i, "fen"] <- resampF["fen"]
    resamples[i, "con_p"] <- resampP["con"]
    resamples[i, "pas"] <- resampP["pas"]
  }
  if(sum((ag_abund[species, "Conserved"] < quantile(resamples[,"con_a"], probs = 0.025, na.rm = T)), 
     (fen_abund[species, "Conserved"] < quantile(resamples[,"con_f"], probs = 0.025, na.rm = T)),
     (pas_abund[species, "Conserved"] < quantile(resamples[,"con_p"], probs = 0.025, na.rm = T)), na.rm = T)
     >= 2){s_mamm_affin[species, "con_aver"] <- TRUE}
  else{s_mamm_affin[species, "con_aver"] <- FALSE}
  if(sum((ag_abund[species, "Conserved"] > quantile(resamples[,"con_a"], probs = 0.975, na.rm = T)), 
     (fen_abund[species, "Conserved"] > quantile(resamples[,"con_f"], probs = 0.975, na.rm = T)),
     (pas_abund[species, "Conserved"] > quantile(resamples[,"con_p"], probs = 0.975, na.rm = T)), na.rm = T)
     >= 2){s_mamm_affin[species, "con_pref"] <- TRUE}
  else{s_mamm_affin[species, "con_pref"] <- FALSE}  
  s_mamm_affin[species, "ag_aver"] <- ag_abund[species, "Agriculture"] < quantile(resamples[,"ag"], probs = 0.025, na.rm = T)
  s_mamm_affin[species, "ag_pref"] <- ag_abund[species, "Agriculture"] > quantile(resamples[,"ag"], probs = 0.975, na.rm = T)
  s_mamm_affin[species, "fen_aver"] <- fen_abund[species, "Fenced"] < quantile(resamples[,"fen"], probs = 0.025, na.rm = T)
  s_mamm_affin[species, "fen_pref"] <- fen_abund[species, "Fenced"] > quantile(resamples[,"fen"], probs = 0.975, na.rm = T)
  s_mamm_affin[species, "pas_aver"] <- pas_abund[species, "Pastoral"] < quantile(resamples[,"pas"], probs = 0.025, na.rm = T)
  s_mamm_affin[species, "pas_pref"] <- pas_abund[species, "Pastoral"] > quantile(resamples[,"pas"], probs = 0.975, na.rm = T)
} # Note that NA values occur when a species was never observed in that land-use type

# Change affinities of rarely sampled species to "Rare"
rare_sps <- colnames(raw_site[4:ncol(raw_site)])[which(apply(raw_site[, 4:ncol(raw_site)], MARGIN = 2, FUN = sum) < 30)]
s_mamm_affin[which(s_mamm_affin[, "species"] %in% rare_sps), 2:9] <- "RARE"

# Load mammal supertree
supertree <- read.nexus("../Data/Mammal.supertree.nexus.txt")

# Create tree of sampled small mammal taxa
phylo <- subset_supertree(raw_site, supertree, 4:25)

# Reorder aversion data so they match the order of tip labels in phylogeny
s_mamm_affin <- s_mamm_affin[match(phylo$tip.label, s_mamm_affin[,"species"]),]

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

# Write affinities and aversions tables to csv
write.csv(affin.table, file = "../Data/Plot data/Small_mammal_affinities.csv", row.names = F)
write.csv(aver.table, file = "../Data/Plot data/Small_mammal_aversions.csv", row.names = F)

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
l_mamm <- read.csv("../Data/Large_mammal_abundance.csv")

# Update landuse column to convert Sanctuary to Conserved and Group to Pastoral
l_mamm$Landuse[which(l_mamm$Landuse == "Sanctuary")] <- "Conservancy"
levels(l_mamm$Landuse)[3] <- "Pastoral"
l_mamm$Landuse <- droplevels(l_mamm$Landuse)

# Remove the species that were not observed at all
no_obs <- which(apply(l_mamm[7:ncol(l_mamm)], MARGIN = 2, FUN = sum) == 0)
l_mamm <- l_mamm[,-(no_obs + 6)]

# Remove outlier for cow abundance
l_mamm$Cow[which(l_mamm$Cow > 1000)] <- NA

# Summarize abundance data by land-use for each species
l_mamm_long <- gather(l_mamm, species, abund, -Landuse, -Site, - Ranch, 
                      -Longitude, -Latitude, -Year)
l_mamm_abund <- l_mamm_long %>% group_by(species) %>%
  summarise(Conserved = tapply(abund, Landuse, FUN = sum, na.rm = T)["Conservancy"],
            Fenced = tapply(abund, Landuse, FUN = sum, na.rm = T)["Fenced"],
            Pastoral = tapply(abund, Landuse, FUN = sum, na.rm = T)["Pastoral"],
            total = sum(abund, na.rm = T))

# Round total abundance values to whole numbers
l_mamm_abund$totalR <- round(l_mamm_abund$total)

# Create matrix to store affinity data
l_mamm_affin <- matrix(NA, nrow = nrow(l_mamm_abund), ncol = 7)
colnames(l_mamm_affin) <- c("species", "con_aver", "con_pref", "fen_aver", "fen_pref", "pas_aver", "pas_pref")
l_mamm_affin[, "species"] <- l_mamm_abund$species

# Calculate affinities
for(species in 1:nrow(l_mamm_abund)){
  resamples <- matrix(NA, nrow = 999, ncol = 3)
  colnames(resamples) <- c("Con", "Fen", "Pas")
  for(i in 1:nrow(resamples)) {
    resamp <- table(sample(c("Con", "Fen", "Pas"), l_mamm_abund$totalR[species], 
                                  replace = T, prob = table(l_mamm$Landuse)/nrow(l_mamm)))
    resamples[i, "Con"] <- resamp["Con"]
    resamples[i, "Fen"] <- resamp["Fen"]
    resamples[i, "Pas"] <- resamp["Pas"]
  }
  l_mamm_affin[species, "con_aver"] <- l_mamm_abund[species, "Conserved"] < quantile(resamples[,"Con"], probs = 0.025, na.rm = T)
  l_mamm_affin[species, "con_pref"] <- l_mamm_abund[species, "Conserved"] > quantile(resamples[,"Con"], probs = 0.975, na.rm = T)
  l_mamm_affin[species, "fen_aver"] <- l_mamm_abund[species, "Fenced"] < quantile(resamples[,"Fen"], probs = 0.025, na.rm = T)
  l_mamm_affin[species, "fen_pref"] <- l_mamm_abund[species, "Fenced"] > quantile(resamples[,"Fen"], probs = 0.975, na.rm = T)
  l_mamm_affin[species, "pas_aver"] <- l_mamm_abund[species, "Pastoral"] < quantile(resamples[,"Pas"], probs = 0.025, na.rm = T)
  l_mamm_affin[species, "pas_pref"] <- l_mamm_abund[species, "Pastoral"] > quantile(resamples[,"Pas"], probs = 0.975, na.rm = T)
}

# Change affinities of rarely sampled species to "Rare"
rare_sps <- l_mamm_abund$species[which(prop.table(l_mamm_abund$totalR) < 0.001)]
l_mamm_affin[which(l_mamm_affin[, "species"] %in% rare_sps), 2:7] <- "RARE"

# Replace common names with latin binomials
taxa_names <- read.csv("../Data/Large_mammal_species_list.csv")
l_mamm_affin[, "species"] <- as.character(taxa_names$Representative_in_supertree[match(l_mamm_affin[,"species"], taxa_names$Spaces_removed_name)])

# Create large mammal phylogeny
raw_com <- l_mamm[, c(3, 7:ncol(l_mamm))]
colnames(raw_com)[-1] <- as.character(taxa_names$Representative_in_supertree[match(colnames(raw_com)[-1], taxa_names$Spaces_removed_name)])
com <- as.matrix(t(raw_com))
phylo <- treedata(mammal.tree, com)$phy

# Reorder aversion data so they match the order of tip labels in phylogeny
l_mamm_affin <- l_mamm_affin[match(phylo$tip.label, l_mamm_affin[,"species"]),]

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

# Write affinities and aversions tables to csv
write.csv(affin.table, file = "../Data/Plot data/Large_mammal_affinities.csv", row.names = F)
write.csv(aver.table, file = "../Data/Plot data/Large_mammal_aversions.csv", row.names = F)

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

#================================
# 4. Creating supplementary table
#================================

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
real_name_lm <- as.character(taxa_names$Binomial[match(l_mamm_affin[, "species"], taxa_names$Representative_in_supertree)])
real_name_plant <- plant_affin[,"species"]
real_name_sm <- s_mamm_affin[, "species"]
affin$real_name <- c(real_name_plant, real_name_sm, real_name_lm)

# Change name for single small mammal species that had different name in tree
# than community data
affin$real_name[which(affin$real_name == "Crocidura_elgonius")] <- "Crocidura spp."

# Add taxonomic group identifier
affin$Taxon <- c(rep("Plant", times = nrow(plant_affin)), 
                 rep("Small mammal", times = nrow(s_mamm_affin)),
                 rep("Large mammal", times = nrow(l_mamm_affin)))

# Remove underscores from species names and capitalize genus names
affin$real_name <- gsub(x = affin$real_name, pattern = "_", replacement = " ")
substr(affin$real_name, 1, 1) <- toupper(substr(affin$real_name, 1, 1))

# Add family information to affinities data
family <- read.csv("../Data/Family_info_all_taxa.csv")
affin$Family <- as.character(family$Family[match(affin$real_name, family$species)])

# Add old aversion data
affin$old_aversion <- as.character(family$Old.aversion[match(affin$real_name, family$species)])

# Format for table S2
tableS2 <- affin[,c(11, 10, 12, 2:9, 13)]
colnames(tableS2)[2] <- "species"

# Write to csv
write.csv(tableS2, file = "../Data/Plot data/Table_S2.csv", row.names = F)
# Significant preferences and aversions will be summarized manually in Excel