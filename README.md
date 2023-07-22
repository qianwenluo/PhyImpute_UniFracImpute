# PhyImpute_UniFracImpute
Two imputation methods, PhyImpute and UniFracImpute, identify and impute non-biological zeros in microbial count data by borrowing information from similar samples

What we need:

•	R version 3.6.1

•	R package: phyloseq, ape, phangorn, phytools, geiger, scDoc

####### load phylogenetic tree and input data ##############

phytree <- read.tree(file="phylotree.tre")

otu.tab <- read.csv("otu_example.csv",row.names = 1,check.names = FALSE)

####### PhyImpute ##############

output1<-phyimpute(otudata=otu.tab, tree=phytree)

####### UniFracImpute ##########

output2 <- unifracimpute(otudata=otu.tab, tree=phytree)

####### Plot the Results #######
library(scales)

library(car)

pv1 <- log(t(output1)+1)

pv2 <- log(t(output2)+1)

pv3 <- log(otu.tab+1)

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
