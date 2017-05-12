# Load function to make transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


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
            rain = mean(AnnualRainfall),
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"])

# Remove conserved from landuse factor levels
byPair$landuse <- droplevels(byPair$landuse)

# Add plant data to combined data frame for later comparison
# of plant and small mammal PD responses
Comb <- byPair
Comb$pPD <- Comb$PD
Comb <- Comb[,-which(names(Comb) == "PD")]

# Make plot
plot(x = byPair$rain, y = byPair$PD, pch = 16,
     col = c("red","gold1","tan2")[byPair$landuse],
     ylab = "PD Change", ylim = c(0.25, 2.3), yaxt = "n", yaxs = "i",
     xlab = "", xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(0.5, 2, 0.5), labels = seq(0.5, 2, 0.5), las = 1, 
     cex.axis = 0.8)

# Make global model and obtain confidence interval
mod <- lm(PD ~ rain, data = byPair)
modrain <- seq(min(byPair$rain), max(byPair$rain), length.out = 1000)
modconf <- predict(mod, newdata = data.frame(rain = modrain), interval = "confidence")

# Add confidence interval to plot
polygon(c(rev(modrain), modrain), c(rev(modconf[,3]), modconf[,2]),
        col = makeTransparent("gray", alpha = 70), border = NA)

# Add model line and delimit confidence intervals
lines(modrain, modconf[,1], col = "black", lwd = 1)
lines(modrain, modconf[,2], lty = 2, col = "black", lwd = 1)
lines(modrain, modconf[,3], lty = 2, col = "black", lwd = 1)

# Replot points over lines
points(x = byPair$rain, y = byPair$PD, pch = 16,
     col = c("red","gold1","tan2")[byPair$landuse])

#======================
# Fig 2b. Small Mammals
#======================

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
            rain = mean(AnnualRainfall),
            PD = PD[Landuse != "Conserved"] / PD[Landuse == "Conserved"])

# Remove conserved from landuse factor levels
byPair$landuse <- droplevels(byPair$landuse)

# Remove outlier
byPair[byPair$PairID == 34, "PD"] <- NA

# Add small mammal data to plant data
Comb$smPD <- byPair$PD

# Make plot
plot(x = byPair$rain, y = byPair$PD, pch = 16,
     col = c("red","gold1","tan2")[byPair$landuse],
     ylab = "PD Change", ylim = c(0, 5), yaxt = "n", yaxs = "i",
     xlab = "", xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(0, 5, 1), labels = seq(0, 5, 1), las = 1, 
     cex.axis = 0.8)

# Make global model and obtain confidence interval
mod <- lm(PD ~ rain, data = byPair)
modrain <- seq(min(byPair$rain), max(byPair$rain), length.out = 1000)
modconf <- predict(mod, newdata = data.frame(rain = modrain), interval = "confidence")

# Add confidence interval to plot
polygon(c(rev(modrain), modrain), c(rev(modconf[,3]), modconf[,2]),
        col = makeTransparent("gray", alpha = 70), border = NA)

# Add model line and delimit confidence intervals
lines(modrain, modconf[,1], col = "black", lwd = 1)
lines(modrain, modconf[,2], lty = 2, col = "black", lwd = 1)
lines(modrain, modconf[,3], lty = 2, col = "black", lwd = 1)

# Replot points over lines
points(x = byPair$rain, y = byPair$PD, pch = 16,
     col = c("red","gold1","tan2")[byPair$landuse])

#======================
# Fig 2c. Large Mammals
#======================

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

# Remove site with no rainfall data
byPlot <- byPlot[!is.na(byPlot$AnnualRainfall),]

# Replace NA with 0
byPlot[is.na(byPlot$PD), "PD"] <- 0

plot(x = byPlot$AnnualRainfall, y = byPlot$PD, type = "n",
     ylab = "PD", xlab = "Annual rainfall (mm)",
     ylim = c(0, 500), yaxt = "n", yaxs = "i",
     xlim = c(420, 750), xaxt = "n", xaxs = "i")
axis(side = 2, at = seq(0,500,100), labels = seq(0,500,100), las = 1,
     cex.axis = 0.8)
axis(side = 1, at = seq(450,750,50), labels = seq(450,750,50),
     cex.axis = 0.8)

CN <- byPlot[byPlot$Landuse == "Conserved",]
EX <- byPlot[byPlot$Landuse == "Fenced",]
PA <- byPlot[byPlot$Landuse == "Pastoral",]
CNmod <- lm(PD ~ AnnualRainfall, data = CN)
EXmod <- lm(PD ~ AnnualRainfall, data = EX)
PAmod <- lm(PD ~ AnnualRainfall, data = PA)

# Calculate confidence intervals
CNrain <- seq(min(CN$AnnualRainfall), max(CN$AnnualRainfall), length.out = 50)
CNconf <- predict(CNmod, newdata = data.frame(AnnualRainfall = CNrain), interval = "confidence")
EXrain <- seq(min(EX$AnnualRainfall), max(EX$AnnualRainfall), length.out = 50)
EXconf <- predict(EXmod, newdata = data.frame(AnnualRainfall = EXrain), interval = "confidence")
PArain <- seq(min(PA$AnnualRainfall), max(PA$AnnualRainfall), length.out = 50)
PAconf <- predict(PAmod, newdata = data.frame(AnnualRainfall = PArain), interval = "confidence")

# Add confidence intervals to plot
polygon(c(rev(CNrain), CNrain), c(rev(CNconf[,3]), CNconf[,2]),
        col = makeTransparent("paleturquoise3", alpha = 70), border = NA)
polygon(c(rev(EXrain), EXrain), c(rev(EXconf[,3]), EXconf[,2]),
        col = makeTransparent("yellow2"), border = NA)
polygon(c(rev(PArain), PArain), c(rev(PAconf[,3]), PAconf[,2]),
        col = makeTransparent("tan2"), border = NA)

# Add lines to delimit confidence intervals
lines(CNrain, CNconf[,1], col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,2], lty=2, col = "skyblue2", lwd = 2)
lines(CNrain, CNconf[,3], lty=2, col = "skyblue2", lwd = 2)

lines(EXrain, EXconf[,1], col = "gold1", lwd = 2)
lines(EXrain, EXconf[,2], lty=2, col = "gold1", lwd = 2)
lines(EXrain, EXconf[,3], lty=2, col = "gold1", lwd = 2)

lines(PArain, PAconf[,1], col = "tan2", lwd = 2)
lines(PArain, PAconf[,2], lty=2, col = "tan2", lwd = 2)
lines(PArain, PAconf[,3], lty=2, col = "tan2", lwd = 2)

#===========================================
# Comparing plant and small mammal responses
#===========================================

# Log transform PD columns
Comb$pPD <- log(Comb$pPD)
Comb$smPD <- log(Comb$smPD)

# Make plot
plot(x = Comb$pPD, y = Comb$smPD, pch = 16,
     #col = c("blue","darkolivegreen2","mediumturquoise")[Comb$landuse],
     col = c("red","gold1","tan2")[Comb$landuse],
     ylab = "ln(Small Mammal PD Change)",
     xlab = "ln(Plant PD Change)")
abline(h = 0, lty = "dashed")
abline(v = 0, lty = "dashed")

# Make global model and obtain confidence interval
mod <- lm(smPD ~ pPD, data = Comb)
modPD <- seq(min(Comb$pPD), max(Comb$pPD), length.out = 1000)
modconf <- predict(mod, newdata = data.frame(pPD = modPD), interval = "confidence")

# Add confidence interval to plot
polygon(c(rev(modPD), modPD), c(rev(modconf[,3]), modconf[,2]),
        col = makeTransparent("gray", alpha = 70), border = NA)

# Add model line and delimit confidence intervals
lines(modPD, modconf[,1], col = "black", lwd = 1)
lines(modPD, modconf[,2], lty = 2, col = "black", lwd = 1)
lines(modPD, modconf[,3], lty = 2, col = "black", lwd = 1)

# Replot points over lines
points(x = Comb$pPD, y = Comb$smPD, pch = 16,
       col = c("red","gold1","tan2")[Comb$landuse])

# Run correlation test
cor.test(x = Comb$pPD, y = Comb$smPD)
summary(mod)$r.squared
