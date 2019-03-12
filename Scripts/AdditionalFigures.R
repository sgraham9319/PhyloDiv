
#====================
# Fig 1. PD responses
#====================

# Load PD response data
plant <- read.csv("../Data/Plot data/Plant_PD_response.csv")
sml.mammal <- read.csv("../Data/Plot data/Small_mammal_PD_response.csv")
lrg.mammal <- read.csv("../Data/Plot data/Large_mammal_PD_response.csv")

# Combine organismal groups into one data frame
PD.resp.plot <- rbind(plant, sml.mammal, lrg.mammal)

# Add organismal group column
PD.resp.plot$org.grp <- c(rep("Plant", times = nrow(plant)),
                          rep("sMamm", times = nrow(sml.mammal)),
                          rep("lMamm", times = nrow(lrg.mammal)))

# Add new column to differentiate all groups for separate boxplots
PD.resp.plot$treat <- paste(PD.resp.plot$org.grp, PD.resp.plot$landuse, sep = "")

# Change new column to factor and manipulate level order for plotting
PD.resp.plot$treat <- as.factor(PD.resp.plot$treat)
PD.resp.plot$treat <- factor(PD.resp.plot$treat, levels = c("PlantAgriculture",
                                                            "sMammAgriculture", 
                                                            "PlantExclosure",
                                                            "sMammExclosure",
                                                            "lMammExclosure",
                                                            "PlantPastoral",
                                                            "sMammPastoral",
                                                            "lMammPastoral"))

# Create the boxplot
boxplot(logPD ~ treat, at = c(1, 2, 4, 5, 6, 8, 9, 10), outline = F,
        ylim = c(-3, 2), 
        col = c("green1", "violetred3", "green1", "violetred3", "blue", "green1", "violetred3", "blue"),
        xaxt = "n", yaxt = "n", ylab = "PD Response",
        xlab = "Type of land use change", whisklty = "solid", 
        data = PD.resp.plot)

# Add labels for organsimal groups on x-axis
axis(side = 1, at = c(1,2,4,5,6,8,9,10), labels = c("P", "S", "P", "S", "L", "P", "S", "L"), cex.axis = 0.8,
     tck = -0.05, padj = -1)

# Add y-axis
axis(side = 2, at = seq(-3, 2, 1), labels = seq(-3, 2, 1), cex.axis = 1, las = 1)

# Add labels for types of land-use change below x-axis
mtext(text = c("Agriculture", "Fenced", "Pastoral"), side = 1, line = 1.5, at = c(1.5,5,9))

# Add line to represent condition of no community change (PD response = 0)
abline(h = 0, lty = "dashed")

# Add lines to separate types on land-use change
abline(v = 3)
abline(v = 7)

#---------------------------------------------------------------
# Determining whether PD responses are statistically significant
#---------------------------------------------------------------

# Create standard error function
SE <- function(x){
  sd(x, na.rm = T) / sqrt(length(!is.na(x)))
}

# Summarize treatments with mean and standard error of PD responses
plotData <- PD.resp.plot %>% 
  group_by(treat) %>%
  summarise(Mean = mean(logPD, na.rm = T),
            St.err = SE(logPD))

# Add column of x-coordinates for the means
plotData$x <- 1:nrow(plotData)

# Plot means and standard errors
plot(x = plotData$x, y = plotData$Mean, ylim = c(min(plotData$Mean-plotData$St.err), 
                                                 max(plotData$Mean+plotData$St.err)),
     ylab = "PD Response", xlab = "Type of land use change", xaxt = "n")
axis(side = 1, at = c(1:8), labels = c("P", "S", "P", "S", "L", "P", "S", "L"), cex.axis = 0.8,
     tck = -0.05, padj = -1)
mtext(text = c("Agriculture", "Exclosure", "Pastoral"), side = 1, line = 1.5, at = c(1.5,4,7), cex = 0.8)
arrows(plotData$x, plotData$Mean-plotData$St.err, plotData$x, 
       plotData$Mean+plotData$St.err, code=0)

# Add dashed line to represent the null hypothesis (zero PD response)
abline(h = 0, lty = "dashed")

# Groups where SEs do not overlap with the dashed line are considered to
# represent statistically significant PD responses to land-use change

#===========================================
# Fig 2. Plant vs. small mammal PD responses
#===========================================

# Extract plant data for sites where small mammal PD responses were analyzed
plant.reduced <- plant[plant$PairID %in% sml.mammal$PairID,]

# Create new data frame containing both plant and small mammal data
Comb <- data.frame(plant.reduced$landuse,
                   plant.reduced$logPD,
                   sml.mammal$logPD)
colnames(Comb) <- c("landuse", "plantPD", "sMammPD")

# Make plot
plot(x = Comb$plantPD, y = Comb$sMammPD, pch = 16,
     col = c("red","gold1","tan2")[Comb$landuse],
     ylab = "Small Mammal PD Response", yaxt = "n",
     xlab = "Plant PD Response", xaxt = "n")
axis(side = 1, at = seq(-1, 0.5, 0.5), labels = seq(-1, 0.5, 0.5), 
     cex.axis = 0.8)
axis(side = 2, at = seq(-1.5, 1.5, 0.5), labels = seq(-1.5, 1.5, 0.5), las = 1, 
     cex.axis = 0.8)

# Make global model and obtain confidence interval
mod <- lm(sMammPD ~ plantPD, data = Comb)
modPD <- seq(min(Comb$plantPD), max(Comb$plantPD), length.out = 1000)
modconf <- predict(mod, newdata = data.frame(plantPD = modPD), interval = "confidence")

# Add model line and delimit confidence intervals
lines(modPD, modconf[,1], col = "black", lwd = 1)
lines(modPD, modconf[,2], lty = 2, col = "black", lwd = 1)
lines(modPD, modconf[,3], lty = 2, col = "black", lwd = 1)

# Replot points over lines
points(x = Comb$plantPD, y = Comb$sMammPD, pch = 16,
       col = c("red","gold1","tan2")[Comb$landuse])

# Run correlation test
cor.test(x = Comb$plantPD, y = Comb$sMammPD)

#====================================
# Fig 4. Mean pairwise distance (MPD)
#====================================

# Load SES data
plant <- read.csv("../Data/Plot data/Plant_MPD_MNTD.csv")
sml.mammal <- read.csv("../Data/Plot data/Small_mammal_MPD_MNTD.csv")
lrg.mammal <- read.csv("../Data/Plot data/Large_mammal_MPD_MNTD.csv")

# Add column to identify taxonomic group
plant$Taxa <- rep("plant", times = nrow(plant))
sml.mammal$Taxa <- rep("sMamm", times = nrow(sml.mammal))
lrg.mammal$Taxa <- rep("lMamm", times = nrow(lrg.mammal))

# Combine taxonomic groups
plotData <- rbind(plant, sml.mammal, lrg.mammal)

# Re-order the data to get desired order for plot
plotData <- plotData[c(1, 5, 9, 2, 6, 3, 7, 10, 4, 8, 11),]

# Add xlim values for the data
plotData$x <- c(1,2,3,5,6,8,9,10,12,13,14)

# Define colors for the data points (one color per taxonomic group)
cols <- c("green1", "violetred3", "blue", "green1", "violetred3", 
          "green1", "violetred3", "blue", "green1", "violetred3", "blue")

# Create plotting function
plot.SES <- function(Means, SEs, y.lab){
  plot(plotData$x, Means, 
       ylim = range(c(Means - SEs, Means + SEs)), yaxt = "n",
       pch = 19, col = cols, xlab = "Type of land use change",
       ylab = y.lab, xaxt = "n")
  axis(side = 1, at = c(1, 2, 3, 5, 6, 8, 9, 10, 12, 13, 14), 
       labels = c("P", "S", "L", "P", "S", "P", "S", "L", "P", "S", "L"),
       cex.axis = 0.8, tck = -0.02, padj = -1)
  arrows(plotData$x, Means - SEs, plotData$x, 
         Means + SEs, code = 0, col = cols)
  mtext(text = c("Conserved", "Agriculture", "Exclosure", "Pastoral"), side = 1, 
        line = 1.5, at = c(2, 5.5, 9, 13))
  abline(h = 0, lty = "dashed")
  abline(v = 4, lwd = 4)
  abline(v = 7)
  abline(v = 11)
}

# Create figure
plot.SES(Means = plotData$MPDz, SEs = plotData$MPDz.SE, y.lab = "MPD (z-score)")

# Add y-axis
axis(side = 2, at = seq(-0.6, 0.6, 0.3), labels = seq(-0.6, 0.6, 0.3), las = 1)

#===========================================
# Fig S4. Mean nearest taxon distance (MNTD)
#===========================================

# Use function defined in previous section to create plot
plot.SES(Means = plotData$MNTDz, SEs = plotData$MNTDz.SE, y.lab = "MNTD (z-score)")

# Add y-axis
axis(side = 2, at = seq(-0.8, 0.4, 0.4), labels = seq(-0.8, 0.4, 0.4), las = 1)

#==========================
# Fig 6. False compensation
#==========================

# Load turnover vs. nestedness data
plant <- read.csv("../Data/Plot data/Plant_false_comp.csv")
sml.mammal <- read.csv("../Data/Plot data/Small_mammal_false_comp.csv")

# Combine plant and small mammal data
plotData <- rbind(plant, sml.mammal)
plotData$Group <- rep(c("Plant", "sMamm"), each = 3)

# Reorder rows and add x-coordinates for plotting
plotData <- plotData[c(1,4,2,5,3,6),]
plotData$x <- 1:6

# Create order for colors
cols <- c("green1", "violetred3", "green1", "violetred3", "green1", "violetred3")

# Create plot
plot(plotData$x, plotData$TurnoverMean, ylim = c(0, 0.7),
     pch = 19, col = cols, xlab = "Type of land use change", yaxs = "i", 
     ylab="Beta-diversity", xaxt = "n", las = 1, xlim = c(0.5, 6.5), xaxs = "i")

# Add x-axis
axis(side = 1, at = 1:6, labels = rep(c("P", "S"), times = 3), cex.axis = 0.8,
     tck = -0.02, padj = -1)

# Add standard errors
arrows(plotData$x, plotData$TurnoverMean - plotData$TurnoverSE, plotData$x, 
       plotData$TurnoverMean + plotData$TurnoverSE, code = 0, col = cols, lwd = 1.5)

# Add land-use type labels below x-axis
mtext(text = c("Agriculture", "Exclosure", "Pastoral"), side = 1, line = 1.5, at = c(1.5,3.5,5.5))

# Add lines to show mean total beta-diversity for each category
arrows(plotData$x - 0.2, plotData$Total, plotData$x + 0.2, 
       plotData$Total, code = 0, col = cols, lwd = 1.5)

# Add lines to separate land-use types
abline(v = 2.5)
abline(v = 4.5)