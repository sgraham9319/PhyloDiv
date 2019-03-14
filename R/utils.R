
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

