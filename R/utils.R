
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