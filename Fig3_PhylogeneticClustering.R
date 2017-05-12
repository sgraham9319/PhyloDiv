##################################################
# Fig 3. Phylogenetic clustering vs. land use type
##################################################

library(ape)
library(picante)
library(geiger)
library(dplyr)

#======================
# Formatting plant data
#======================

# Load community data
rawCom <- read.csv("../Data/Plant community data.csv")

# Load phylogeny
phylo <- read.tree("../Data/PlantPhylo")

# Calculate PD for each site and store in data frame
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

# Calculate expected PD table
expPD <- expected.pd(phylo)

# Add expected values to PD data
PDData$expPD <- expPD$expected.pd[match(PDData$SR, expPD$n)]

# Load plot data
byPlot <- read.csv("../Data/Plant data by plot.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR", "expPD")])

# Calculate observed / expected PD
byPlot$obsexp <- byPlot$PD / byPlot$expPD

plant <- byPlot %>% 
  group_by(Landuse) %>% 
  summarise(Taxa = "Plant",
            Mean = mean(obsexp, na.rm = T),
            SE = sd(obsexp, na.rm = T) / sqrt(length(!is.na(obsexp))))


#=============================
# Formatting small mammal data
#=============================

# Load mammal supertree
mammal.tree <- read.nexus("../Data/Mammal.supertree.nexus.txt")

# Select the "bestDates" tree
mammal.tree <- mammal.tree$mammalST_bestDates

# Load small mammal community data
rawCom <- read.csv("../Data/Diurnal excluded com data.csv")

# Prune supertree to only the sampled taxa
smammal <- as.matrix(t(rawCom))
phylo <- treedata(mammal.tree, smammal)$phy

# Calculate PD for each site and store in data frame
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

# Calculate expected PD table
expPD <- expected.pd(phylo)

# Add expected values to PD data
PDData$expPD <- expPD$expected.pd[match(PDData$SR, expPD$n)]

# Load site data
byPlot <- read.csv("../Data/Nocturnal only site data.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR", "expPD")])

# Calculate observed / expected PD
byPlot$obsexp <- byPlot$PD / byPlot$expPD

sMamm <- byPlot %>%
  group_by(Landuse) %>% 
  summarise(Taxa = "sMamm",
            Mean = mean(obsexp, na.rm = T),
            SE = sd(obsexp, na.rm = T) / sqrt(length(!is.na(obsexp))))


#=============================
# Formatting large mammal data
#=============================

# Load large mammal community data
rawCom <- read.csv("../Data/Large mammal com data.csv")

# Remove species that did not appear in any of sampled sites
spsRm <- names(which(apply(rawCom[,-1], MARGIN = 2, FUN = sum) == 0))
rawCom <- rawCom[,-which(names(rawCom) == spsRm)]

# Prune supertree to only the sampled taxa
lmammal <- as.matrix(t(rawCom))
phylo <- treedata(mammal.tree, lmammal)$phy

# Calculate PD for each site and store in data frame
PDData <- data.frame(rawCom$Site, pd(rawCom[,-1], phylo, include.root = F))
names(PDData)[1] <- "Site"

# Calculate expected PD table
expPD <- expected.pd(phylo)

# Add expected values to PD data
PDData$expPD <- expPD$expected.pd[match(PDData$SR, expPD$n)]

# Load site data
byPlot <- read.csv("../Data/Large mammal data.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR", "expPD")])

# Calculate observed / expected PD
byPlot$obsexp <- byPlot$PD / byPlot$expPD

lMamm <- byPlot %>%
  group_by(Landuse) %>% 
  summarise(Taxa = "lMamm",
            Mean = mean(obsexp, na.rm = T),
            SE = sd(obsexp, na.rm = T) / sqrt(length(!is.na(obsexp))))


#==================
# Creating the plot
#==================

# Combine data from different organismal groups
plotData <- rbind(plant, sMamm, lMamm)

# Re-order the data to get desired order for plot. Nest organismal group within type
# of land use change, with the order: plant, small mammal, large mammal. Order the
# types of land use as: Agriculture, Exclosure, Pastoral. Also remove NA.
plotData <- plotData[c(2,6,9,1,5,3,7,10,4,8,11),]

# Add xlim values for the plot
plotData$x <- c(1,2,3,5,6,8,9,10,12,13,14)

# Plot the data
cols <- c("green1", "violetred3", "blue", "green1", "violetred3", 
          "green1", "violetred3", "blue", "green1", "violetred3", "blue")
plot(plotData$x, plotData$Mean, ylim=range(c(0.65, 1.1)),
     pch=19, col=cols, xlab="Type of land use change",
     ylab="Observed/Expected PD", xaxt = "n", yaxt = "n", yaxs = "i")
axis(side = 2, at = seq(0.7, 1.1, 0.1), labels = seq(0.7, 1.1, 0.1),
     las = 1)
axis(side = 1, at = c(1,2,3,5,6,8,9,10,12,13,14), 
     labels = c("P", "S", "L", "P", "S", "P", "S", "L", "P", "S", "L"),
     cex.axis = 0.8, tck = -0.02, padj = -1)
arrows(plotData$x, plotData$Mean-plotData$SE, plotData$x, 
       plotData$Mean+plotData$SE, code=0, col=cols)
mtext(text = c("Conserved", "Agriculture", "Exclosure", "Pastoral"), side = 1, 
      line = 1.5, at = c(2,5.5,9,13))
abline(h = 1, lty = "dashed")
abline(v = 4)
abline(v = 7)
abline(v = 11)
