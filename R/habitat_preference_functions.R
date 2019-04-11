

#===================================
# Evaluate habitat preference traits
#===================================

hab_pref_summ <- function(ag_dat, fen_dat, pas_dat){
  output <- ag_dat %>%
    select(species) %>%
    mutate(agriculture = ag_dat$Agriculture / ag_dat$total,
           ag_adj = (ag_dat$Agriculture + 0.5) / (ag_dat$total + 1),
           fenced = fen_dat$Fenced / fen_dat$total,
           fen_adj = (fen_dat$Fenced + 0.5) / (fen_dat$total + 1),
           pastoral = pas_dat$Pastoral / pas_dat$total,
           pas_adj = (pas_dat$Pastoral + 0.5) / (pas_dat$total + 1),
           con_ag = ag_dat$Conserved / ag_dat$total,
           con_fen = fen_dat$Conserved / fen_dat$total,
           con_pas = pas_dat$Conserved / pas_dat$total,
           con_ag_adj = (ag_dat$Conserved + 0.5) / (ag_dat$total + 1),
           con_fen_adj = (fen_dat$Conserved + 0.5) / (fen_dat$total + 1),
           con_pas_adj = (pas_dat$Conserved + 0.5) / (pas_dat$total + 1))
  output$conserved <- apply(output[, 8:10], 1, mean, na.rm = T)
  output$con_adj <- apply(output[, 11:13], 1, mean, na.rm = T)
  output <- output %>%
    mutate(agriculture_logit = log(ag_adj / (1 - ag_adj)),
           fenced_logit = log(fen_adj / (1 - fen_adj)),
           pastoral_logit = log(pas_adj / (1 - pas_adj)),
           conserved_logit = log(con_adj / (1 - con_adj))) %>%
    select(-ag_adj, -fen_adj, -pas_adj, -con_ag, -con_fen, -con_pas,
           -con_ag_adj, -con_fen_adj, -con_pas_adj, -con_adj) %>%
    as.data.frame()
  output[which(is.na(output$agriculture)), "agriculture_logit"] <- NaN
  output[which(is.na(output$fenced)), "fenced_logit"] <- NaN
  output[which(is.na(output$pastoral)), "pastoral_logit"] <- NaN
  output[which(is.na(output$conserved)), "conserved_logit"] <- NaN
  output
}

#===========================================
# Calculate p-values for phylogenetic signal
#===========================================

phy_sig_p <- function(x, obs){
  sum(x > obs) / 1000
}

#==============================
# Calculate phylogenetic signal
#==============================

# Calculates phylogenetic signal for each trait provided in a data frame where
# columns 2:ncol(data.frame) each contain data for a different trait, and column
# 1 contains taxa names. The taxa names in the trait data do not need to be in
# the same order as in the tip labels of the provided phylogeny, but the names
# do need to match exactly.

phy_sig_summ <- function(trait_data, phylo){
  
  # Reorder trait data to match phylogeny tip labels
  hab_prefs <- trait_data[match(phylo$tip.label, trait_data$species),]
  
  # Extract trait names
  traits <- colnames(hab_prefs)[2:ncol(hab_prefs)]
  
  # Create output matrix
  output <- matrix(NA, nrow = 4, ncol = length(traits))
  colnames(output) <- traits
  rownames(output) <- c("K", "K0_null_median", "K0_null_sd", "K0_null_p")
  
  # Calculate phylogenetic signal (Blomberg's K) for each trait
  for(trait in 1:length(traits)){
    output[1, trait] <- phylosig(phylo, hab_prefs[, traits[trait]], 
                                 method = "K")
  }
  
  # Create distribution of K values under null model of no phylogenetic signal
  no_signal_null <- matrix(NA, nrow = 1000, ncol = length(traits))
  for(iter in 1:nrow(no_signal_null)){
    new_phylo <- tipShuffle(phylo)
    new_hab_prefs <- hab_prefs[match(new_phylo$tip.label, hab_prefs$species),]
    for(trait in 1:length(traits)){
      k_zero <- phylosig(new_phylo, new_hab_prefs[, traits[trait]], method = "K")
      no_signal_null[iter, trait] <- k_zero
    }
  }
  
  # Summarize null distribution
  output[2, ] <- apply(no_signal_null, 2, median)
  output[3, ] <- apply(no_signal_null, 2, sd)
  
  # Calculate p-value by comparing observed K to null distribution
  output[4, ] <- apply(no_signal_null, 2, phy_sig_p, obs = output[1, ])
  
  # Return output
  output
}

#=============================================
# Convert site x species matrix to long format
#=============================================

long <- function(dat){
  dat %>%
    select(-longitude, - latitude, - annual_rainfall, -frac_veg, -erosion,
           -root_dep_res_50cm, -soil_org_carb, -pH, -sand,
           -soil_pc1, -soil_pc2) %>%
    gather(species, abund, -site, -landuse, -pair_id)
}

#================================================
# Summarize long format abundance data by species
#================================================

abund <- function(long_dat, dist_landuse_type){
  output <- long_dat %>% group_by(species) %>%
    summarise(Conserved = tapply(abund, landuse, FUN = sum, na.rm = T)["Conserved"],
              Dist_land = tapply(abund, landuse, FUN = sum, na.rm = T)[dist_landuse_type],
              total = sum(abund, na.rm = T))
  colnames(output)[3] <- dist_landuse_type
  output
}

#==========================================================
# Summarize abundance data by landuse type for each species
#==========================================================

abund_summ <- function(comm_dat){
  
  # Separate data by landuse type of the plot pair
  ag_pairs <- comm_dat$pair_id[which(comm_dat$landuse == "Agriculture")]
  pas_pairs <- comm_dat$pair_id[which(comm_dat$landuse == "Pastoral")]
  fen_pairs <- comm_dat$pair_id[which(comm_dat$landuse == "Fenced")]
  ag <- droplevels(comm_dat[which(comm_dat$pair_id %in% ag_pairs),])
  pas <- droplevels(comm_dat[which(comm_dat$pair_id %in% pas_pairs),])
  fen <- droplevels(comm_dat[which(comm_dat$pair_id %in% fen_pairs),])
  
  # Summarize abundance for each species
  long_ag <- long(ag)
  ag_abund <- abund(long_ag, "Agriculture")
  long_fen <- long(fen)
  fen_abund <- abund(long_fen, "Fenced")
  long_pas <- long(pas)
  pas_abund <- abund(long_pas, "Pastoral")
  
  # Combine landuse types into single output
  total_abund <- list(ag_abund, fen_abund, pas_abund)
  
  # Return output
  total_abund
}

#======================================================
# Identify significant landuse aversions and affinities
#======================================================

affin_signif <- function(ag_dat, fen_dat, pas_dat){
  
  # Create matrix to store affinity data
  affin <- matrix(NA, nrow = nrow(ag_dat), ncol = 9)
  colnames(affin) <- c("species", "con_aver", "con_pref", "ag_aver", "ag_pref", "fen_aver", "fen_pref", "pas_aver", "pas_pref")
  affin[, "species"] <- ag_dat$species
  
  # Loop through taxa, calculating aversions and affinities for each
  for(species in 1:nrow(ag_dat)){
    resamples <- matrix(NA, nrow = 999, ncol = 6)
    colnames(resamples) <- c("con_a", "ag", "con_f", "fen", "con_p", "pas")
    
    # Loop through landuse types, resampling abundance for each
    for(i in 1:nrow(resamples)) {
      resampA <- table(sample(c("con", "ag"), ag_dat$total[species], replace = T))
      resampF <- table(sample(c("con", "fen"), fen_dat$total[species], replace = T))
      resampP <- table(sample(c("con", "pas"), pas_dat$total[species], replace = T))
      resamples[i, "con_a"] <- resampA["con"]
      resamples[i, "ag"] <- resampA["ag"]
      resamples[i, "con_f"] <- resampF["con"]
      resamples[i, "fen"] <- resampF["fen"]
      resamples[i, "con_p"] <- resampP["con"]
      resamples[i, "pas"] <- resampP["pas"]
    }
    
    # Calculate aversion or affinity for conserved landuse
    if(sum((ag_dat[species, "Conserved"] < quantile(resamples[,"con_a"], probs = 0.025, na.rm = T)), 
           (fen_dat[species, "Conserved"] < quantile(resamples[,"con_f"], probs = 0.025, na.rm = T)),
           (pas_dat[species, "Conserved"] < quantile(resamples[,"con_p"], probs = 0.025, na.rm = T)), na.rm = T)
       >= 2){affin[species, "con_aver"] <- TRUE}
    else{affin[species, "con_aver"] <- FALSE}
    if(sum((ag_dat[species, "Conserved"] > quantile(resamples[,"con_a"], probs = 0.975, na.rm = T)), 
           (fen_dat[species, "Conserved"] > quantile(resamples[,"con_f"], probs = 0.975, na.rm = T)),
           (pas_dat[species, "Conserved"] > quantile(resamples[,"con_p"], probs = 0.975, na.rm = T)), na.rm = T)
       >= 2){affin[species, "con_pref"] <- TRUE}
    else{affin[species, "con_pref"] <- FALSE}
    
    # Record aversions/affinities for disturbed landuses
    affin[species, "ag_aver"] <- ag_dat[species, "Agriculture"] < quantile(resamples[, "ag"], probs = 0.025, na.rm = T)
    affin[species, "ag_pref"] <- ag_dat[species, "Agriculture"] > quantile(resamples[, "ag"], probs = 0.975, na.rm = T)
    affin[species, "fen_aver"] <- fen_dat[species, "Fenced"] < quantile(resamples[, "fen"], probs = 0.025, na.rm = T)
    affin[species, "fen_pref"] <- fen_dat[species, "Fenced"] > quantile(resamples[, "fen"], probs = 0.975, na.rm = T)
    affin[species, "pas_aver"] <- pas_dat[species, "Pastoral"] < quantile(resamples[, "pas"], probs = 0.025, na.rm = T)
    affin[species, "pas_pref"] <- pas_dat[species, "Pastoral"] > quantile(resamples[, "pas"], probs = 0.975, na.rm = T)
  }
  
  # Return output
  affin
}