##################################################
# Fig 3. Phylogenetic clustering vs. land use type
##################################################


#===================================
# Evaluating phylogenetic clustering
#===================================

# As expected PD for each site has uncertainty associated with it, when we average
# the observed/expected ratio across sites within a land use type we will need
# to account for error propagation. One way to account for error propagation is
# to simply calculate the average for each land use type using resampling and
# using the variance among these averages as the uncertainty in the group average

# Produce phylogenetic distance matrix from phylogeny
phyDist <- cophenetic(phylo)

# Calculate mean pariwise distance (takes some time to run)
ses.mpd(rawCom[,-1], phyDist, null.model = "taxa.labels", runs = 99)

# Column "mpd.obs.z" in the output is the standardized effect size of MPD - if
# the value is positive the community is overdispersed, but if negative it is
# clustered.

# Instead of using mpd, could also use mntd (mean nearest taxon distance). 
# MPD is generally thought to be more sensitive to tree-wide patterns of 
# phylogenetic clustering and eveness, while MNTD is more sensitive to
# patterns of evenness and clustering closer to the tips of the phylogeny. So
# if community comes from many different orders, but very specific subsets of
# the species within those orders, MPD may be zero but MNTD would be negative.