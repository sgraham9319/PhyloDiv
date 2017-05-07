

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

# Load plot data
byPlot <- read.csv("../Data/Plant data by plot.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR")])

# Summarize data by pair
byPair <- byPlot %>% 
  group_by(PairID) %>% 
  summarise(landuse = Landuse[Landuse != "Conserved"], 
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"])

# Add plant data to a data frame for all organismal groups
byPair$Grp <- "plant"
allData <- byPair

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

# Load site data
byPlot <- read.csv("../Data/Nocturnal only site data.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR")])

# Summarize data by pair
byPair <- byPlot %>% 
  group_by(PairID) %>% 
  summarise(landuse = Landuse[Landuse != "Conserved"], 
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"])

# Add to plant data
byPair$Grp <- "sMamm"
allData <- rbind(allData, byPair)

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

# Load site data
byPlot <- read.csv("../Data/Large mammal data.csv")

# Combine PD and plot data
byPlot <- cbind(byPlot, PDData[match(byPlot$Site, PDData$Site), 
                               c("PD", "SR")])

# Remove NAs
byPlot <- byPlot[!is.na(byPlot$PD),]

# Create subsets
Con <- byPlot[byPlot$Landuse == "Conserved",]
Past <- byPlot[byPlot$Landuse == "Pastoral",]
Exc <- byPlot[byPlot$Landuse == "Fenced",]

# Obtain random sample of Conserved sites
resampCon <- Con[sample(x = 1:nrow(Con), size = 1000, replace = T),]

# Obtain random samples of disturbed sites
resampPast <- Past[sample(x = 1:nrow(Past), size = 1000, replace = T),]
resampExc <- Exc[sample(x = 1:nrow(Exc), size = 1000, replace = T),]

# Calculate PD response variables for each type of land use change
PDPast <- resampPast$PD / resampCon$PD
PDExc <- resampExc$PD / resampCon$PD

# Put in data frame
PD <- c(PDPast, PDExc)
Grp <- rep("lMamm", times = 2000)
landuse <- rep(c("Pastoral", "Exclosure"), each = 1000)
PairID <- rep(NA, times = 2000)
lMamm <- data.frame(PairID, landuse, PD, Grp)

# Add to plant and small mammal data
allData <- rbind(allData, lMamm)

#==================
# Creating the plot
#==================

# Create factor with unique value for each organismal group and each
# type of land use change
allData$treat <- paste(allData$Grp, allData$landuse, sep = "")

# Change to factor and manipulate level order for plotting
allData$treat <- as.factor(allData$treat)
allData$treat <- factor(allData$treat, levels = c("plantAgriculture",
                                                  "sMammAgriculture", 
                                                  "plantExclosure",
                                                  "sMammExclosure",
                                                  "lMammExclosure",
                                                  "plantPastoral",
                                                  "sMammPastoral",
                                                  "lMammPastoral"))

# Create the plot
boxplot(PD ~ treat, at = c(1,2,4,5,6,8,9,10), ylim = c(0,2.1), outline = F, 
        col = c("red", "blue", "red", "blue", "gold1", "red", "blue", "gold1"),
        xaxt = "n", yaxt = "n",ylab = "PD response", 
        xlab = "Type of land use change", whisklty = "solid", 
        data = allData)
axis(side = 1, at = c(1,2,4,5,6,8,9,10), labels = c("P", "S", "P", "S", "L", "P", "S", "L"), cex.axis = 0.8,
     tck = -0.05, padj = -1)
axis(side = 2, at = seq(0, 2, 1), labels = seq(0, 2, 1), cex.axis = 1, las = 1)
mtext(text = c("Agriculture", "Exclosure", "Pastoral"), side = 1, line = 1.5, at = c(1.5,5,9))
abline(h = 1, lty = "dashed")
abline(v = 3)
abline(v = 7)
