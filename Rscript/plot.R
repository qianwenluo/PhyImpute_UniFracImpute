plot <- function(datain, main){
  input <- log(t(datain)+1)
  dataEllipse(y = colMeans(input[,]), x = apply(input[,], 2, sd),levels = 0.80,  
              main = main, xlim = c(0, max(apply(input[,], 2, sd))), ylim = c(0, max(colMeans(input[,]))), 
              ylab = "Taxon mean", xlab = "Taxon SD", 
              cex.lab=1,grid = FALSE)
  
}
#### Example ####

pv1 <- read.csv("phyImputed.csv",row.names=1)
pv2 <- read.csv("Unifracimputed.csv",row.names=1)
raw <- read.csv("otu_example.csv",row.names=1)

pdf("mean_sd_scatter_epplise.pdf",width = 6.69)
par(mfrow = c(2,2))
plot(pv1, main="phyImpute");
abline(v = 2.78, h=9.1)
plot(pv2, main = "UniFracImpute")
abline(v = 2.81, h=9)
plot(raw, main = "NoImpute")
abline(v = 4.42, h=5.85)
dev.off()
