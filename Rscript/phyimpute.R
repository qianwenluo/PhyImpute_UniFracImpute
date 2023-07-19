################### phyImpute ################

######### load packackges ####################
library(phyloseq)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(scDoc)


####### load phylogenetic tree ##############
tree <- read.tree(file="phylotree.tre")

otu.tab <- read.csv("otu_example.csv",row.names = 1,check.names = FALSE)


####### calculate the probabilities of non-biological zeros #######
input <- t(otu.tab)
offsets <- as.numeric(log(colSums(t(otu.tab))))
dmat <- prob.dropout(input = t(otu.tab), is.count = T, mcore = 3,offsets = offsets)
dprob <- 0.5

nc <- dim(input)[2]  # number of samples
ng <- dim(input)[1]  # number of otus

w.prior <- rowSums(dmat > dprob) / nc

####### calculate the sample-to-sample similarites #######

sim.o <- matrix(0, nrow=nc, ncol=nc) # output of similarity (sample x sample matrix)
rownames(sim.o) <- colnames(sim.o) <- colnames(input)

otutree <- function(ntip, edge2, otu.tab){
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i] 
    node <- edge[tip.loc, 1]            
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]   
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  return(cum)
}

for (i in 1:nc) {
  for (j in 1:nc) {
    
    if (i < j) {
      
      prob0 <- rep(1, ng)
      prob0[(dmat[,i]>dprob) | (dmat[,j]>dprob)] <- -1 
      prob0[prob0 == -1] <- w.prior[prob0 == -1]
      
      A <- input[,i] * sqrt(prob0)
      B <- input[,j] * sqrt(prob0)
      
      tab <- t(cbind(A, B))
      row.sum <- rowSums(tab)
      tab <- tab / row.sum
      n <- nrow(tab)
      
      absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
      if (length(absent) != 0) {
        tree <- drop.tip(tree, absent)
        warning("The tree has more OTU than the OTU table!")
      }
      tip.label <- tree$tip.label
      tab <- tab[, tip.label]
      
      ntip <- length(tip.label)
      nbr <- nrow(tree$edge)  
      edge <- tree$edge
      edge2 <- edge[, 2]
      br.len <- tree$edge.length
      cum <- matrix(0, nbr, n)
      cum <- otutree(ntip, edge2, tab)
      
      A_1 <- cum[, 1] * sqrt(br.len)
      B_1 <- cum[, 2] * sqrt(br.len)
      
      sim.o[i, j] <- sim.o[j, i] <- sum(A_1*B_1) / sqrt(sum(A_1^2)*sum(B_1^2))
      
    }
  }
}
out <- sim.o / max(sim.o)
sim.mat <- out

####### impute the values #######
colnames(sim.mat) <- rownames(sim.mat)
impute.mat <- impute.knn(input = t(otu.tab), dmat = dmat, sim.mat = sim.mat, 
                         k = 5, sim.cut = 1e-4)

###### save the result #######
write.csv(impute.mat$output.imp, file="phyImputed.csv")


