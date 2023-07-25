# PhyImpute_UniFracImpute
Two imputation methods, PhyImpute and UniFracImpute, identify and impute non-biological zeros in microbial count data by borrowing information from similar samples

What we need:

•	R version 3.6.1

•	R package: phyloseq, ape, phangorn, phytools, geiger, scDoc, scales, car

####### load phylogenetic tree and input data ##############

phytree <- read.tree(file="phylotree.tre")

otu.tab <- read.csv("otu_example.csv",row.names = 1,check.names = FALSE)

####### PhyImpute ##############

output1<-phyimpute(otudata=otu.tab, tree=phytree)

####### UniFracImpute ##########

output2 <- unifracimpute(otudata=otu.tab, tree=phytree)

####### Plot the Results #######

plot(output1, main="PhyImpute)

abline(v = 2.78, h=9.1)
