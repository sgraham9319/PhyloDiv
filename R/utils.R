

#=========================
# Calculate standard error
#=========================

se <- function(x){
  sd(x, na.rm = T) / sqrt(length(na.omit(x)))
}

#===========================================================================
# Add Faith's PD and species richness as columns to a site by species matrix
#===========================================================================

faith_pd <- function(comm, phylo, abundance_cols){
  comm_mat <- comm[, abundance_cols]
  comm_mat <- comm_mat %>%
    mutate_all(function(x){as.numeric(x > 0)})
  output <- cbind(comm, pd(comm_mat, phylo, include.root = F))
  output[which(is.na(output[, "PD"])), "PD"] <- 0
  output
}

#======================================================
# Subset mammal supertree to taxa included in site data
#======================================================

subset_supertree <- function(comm, supertree, abundance_cols){
  comm_mat <- comm[, abundance_cols]
  comm_mat <- comm_mat %>%
    mutate_all(function(x){as.numeric(x > 0)})
  comm_mat <- as.matrix(t(comm_mat))
  best_dates_tree <- supertree$mammalST_bestDates
  phylo <- suppressWarnings(treedata(best_dates_tree, comm_mat))$phy
  phylo
}

#========================
# Make transparent colors
#========================

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#=========================================================
# Change large mammal taxa names to match mammal supertree
#=========================================================

match_phylo_names <- function(dat){
  species_key <- read.csv("Data/large_mammal_species_key.csv")
  relev_cols <- match(colnames(dat)[7:62], species_key$binomial)
  supertree_names <- as.character(species_key$representative_in_supertree[relev_cols])
  colnames(dat)[7:62] <- supertree_names
  dat
}

#=====================================================
# Summarize decomposed beta diversity by land use type
#=====================================================

beta_decomp <- function(com_mat, phylo, sps_cols){
  
  # Create dataframe to store decomposition results
  results <- as.data.frame(matrix(NA, ncol = 5, nrow = length(unique(s_mamm$pair_id))))
  colnames(results) <- c("pair_id", "landuse", "beta_PD", "turnover", "nestedness")
  results$pair_id <- unique(s_mamm$pair_id)
  results$landuse <- s_mamm[s_mamm$landuse != "Conserved", ]$landuse
  
  # Loop through pairs, decomposing beta pd for each
  for(pair in 1:nrow(results)){
    com <- com_mat[com_mat$pair_id == pair, sps_cols]
    results[results$pair_id == pair, 3:5] <- beta.pd.decompo(com = com, tree = phylo, 
                                                             type = "Unifrac")$betadiv
  }
  
  # Summarize output
  output <- results %>% group_by(landuse) %>% 
    summarise(total_av = mean(beta_PD, na.rm = T),
              total_sd = sd(beta_PD, na.rm = T),
              turnover_av = mean(turnover, na.rm = T),
              turnover_sd = sd(turnover, na.rm = T),
              turnover_se = se(turnover),
              nested_av = mean(nestedness, na.rm = T),
              nested_sd = sd(nestedness, na.rm = T))
  
  # Return output
  output
}
