# PhyImpute_UniFracImpute
Two imputation methods, PhyImpute and UniFracImpute, identify and impute non-biological zeros in microbial count data by borrowing information from similar samples

What we need:

•	R version 3.6.1

•	R package: phyloseq, ape, phangorn, phytools, geiger, scDoc, scales, car, metagenomeSeq

####### load phylogenetic tree and input data ##############

phytree <- read.tree(file="phylotree.tre")

otu.tab <- read.csv("otu_example.csv",row.names = 1,check.names = FALSE)

####### PhyImpute ##############

## Use PNB to estimate posterior probabilities 
output1<-phyimpute(otudata=otu.tab, tree=phytree, method = 'pnb')

## Use Zero Inflated Gamma to estimate posterior probabilities
output2<-phyimpute(otudata=otu.tab, tree=phytree,method='zig')

####### UniFracImpute ###########

output3 <- unifracimpute(otudata=otu.tab, tree=phytree, method = 'pnb')
output4 <- unifracimpute(otudata=otu.tab, tree=phytree,method='zig')

####### Plot the Results ##########

plot(output1, main="PhyImpute")

