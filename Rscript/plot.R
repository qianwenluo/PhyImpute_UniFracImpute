library(scales)
library(car)
#### phyimpute ####

pv1 <- read.csv("phyImputed.csv",row.names=1)
pv1 <- log(t(pv1)+1)

pv2 <- read.csv("Unifracimputed.csv",row.names=1)
pv2 <- log(t(pv2)+1)

raw <- read.csv("otu_example.csv",row.names=1)
raw <- log(raw+1)

pdf("mean_sd_scatter_epplise.pdf",width = 6.69)
par(mfrow = c(2,2))
dataEllipse(y = colMeans(pv1[,]), x = apply(pv1[,], 2, sd),levels = 0.80,  
            main = "PhyImpute", xlim = c(0, 6), ylim = c(0,14.5), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1,grid = FALSE)
abline(v = 2.78, h=9.1)

dataEllipse(y = colMeans(pv2[,]), x = apply(pv2[,], 2, sd),levels = 0.80,
            main = "UniFracImpute", xlim = c(0, 6), ylim = c(0,14.5), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1,grid = FALSE)
abline(v = 2.81, h=9)

dataEllipse(y = colMeans(raw[,]), x = apply(raw[,], 2, sd),levels = 0.80,
            main = "NoImputation", xlim = c(0, 6), ylim = c(0,14.5), ylab = "Taxon mean", xlab = "Taxon SD", cex.lab=1,grid = FALSE)
abline(v = 4.42, h=5.85)
dev.off()
