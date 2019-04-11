

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