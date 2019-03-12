#########################################################
# Plotting standardized effect sizes (SES) for PD metrics
#########################################################

# Load SES data
plant <- read.csv("../Data/Plot data/Plant_PD_summary_Ind_Swap.csv")
s.mamm <- read.csv("../Data/Plot data/Sm_Mamm_PD_summary_Ind_Swap.csv")
l.mamm <- read.csv("../Data/Plot data/Lg_Mamm_PD_summary_Ind_Swap.csv")

# Add column to identify taxonomic group
plant$Taxa <- rep("plant", times = nrow(plant))
s.mamm$Taxa <- rep("sMamm", times = nrow(s.mamm))
l.mamm$Taxa <- rep("lMamm", times = nrow(l.mamm))

# Combine taxonomic groups
plotData <- rbind(plant, s.mamm, l.mamm)

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
       ylim = range(c(Means - SEs, Means + SEs)),
       pch = 19, col = cols, xlab = "Type of land use change",
       ylab = y.lab, xaxt = "n")
  # yaxt = "n", yaxs = "i")
  # axis(side = 2, at = seq(0.7, 1.1, 0.1), labels = seq(0.7, 1.1, 0.1), las = 1)
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

# Plot SES figures for PD, MPD, and MNTD
plot.SES(Means = plotData$PDz, SEs = plotData$PDz.SE, y.lab = "PD (z-score)")
plot.SES(Means = plotData$MPDz, SEs = plotData$MPDz.SE, y.lab = "MPD (z-score)")
plot.SES(Means = plotData$MNTDz, SEs = plotData$MNTDz.SE, y.lab = "MNTD (z-score)")

# Plot centile figures for PD, MPD, and MNTD
plot.SES(Means = plotData$PDcent, SEs = plotData$PDcent.SE, y.lab = "PD (centile)")
plot.SES(Means = plotData$MPDcent, SEs = plotData$MPDcent.SE, y.lab = "MPD (centile)")
plot.SES(Means = plotData$MNTDcent, SEs = plotData$MNTDcent.SE, y.lab = "MNTD (centile)")
