################### UnifracImpute ################

unifracimpute <- function(otudata, tree){
  
  ####### calculate the probabilities of non-biological zeros #######
  offsets <- as.numeric(log(colSums(t(otudata))))
  dmat <- prob.dropout(input = t(otudata), is.count = F, mcore = 3,offsets = offsets)
  dprob <- 0.5
  
  otudata <- as.matrix(otudata)
  row.sum <- rowSums(otudata)
  input <- otudata / row.sum
  
  nc <- dim(input)[1] # number of samples
  ng <- dim(input)[2] # number of otus
  
  w.prior <- rowSums(dmat > dprob) / nc
  
  ####### calculate the sample-to-sample similarity #######
  otutree <- function(ntip, edge2, otudata){
    for (i in 1:ntip) {
      tip.loc <- which(edge2 == i)
      cum[tip.loc, ] <- cum[tip.loc, ] + otudata[, i]	
      node <- edge[tip.loc, 1]						 
      node.loc <- which(edge2 == node)
      while (length(node.loc)) {
        cum[node.loc, ] <- cum[node.loc, ] + otudata[, i]		
        node <- edge[node.loc, 1]
        node.loc <- which(edge2 == node)
      }
    }
    return(cum)
  }
  
  alpha = 1
  unifracs <- array(NA, c(nc, nc))
  for (i in 1:nc){
    unifracs[i,i] <- 0
  }
  
  for (i in 2:nc) {
    for (j in 1:(i-1)) {
      prob0 <- rep(1, ng)
      prob0[(dmat[,i]>dprob) | (dmat[,j]>dprob)] <- -1 
      prob0[prob0 == -1] <- w.prior[prob0 == -1]
      
      A <- input[i,] * sqrt(prob0) # input is proportion OTU table
      B <- input[j,] * sqrt(prob0)
      tab <- t(cbind(A, B))
      row.sum <- rowSums(tab)
      tab <- tab / row.sum
      n <- nrow(tab)
      
      absent <- tree$tip.label[!(tree$tip.label %in% colnames(otudata))]
      if (length(absent) != 0) {
        tree <- drop.tip(tree, absent)
       
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
      
      cum1 <- cum[, 2]
      cum2 <- cum[, 1]
      ind <- (cum1 + cum2) != 0
      cum1 <- cum1[ind]
      cum2 <- cum2[ind]		
      br.len2 <- br.len[ind]			
      diff <- abs(cum1 - cum2) / (cum1 + cum2)		
      
      # Generalized UniFrac distance
      w <- br.len2 * (cum1 + cum2)^alpha
      unifracs[i, j] <- unifracs[j, i] <- sum(diff * w) / sum(w)
    }
  }  
  
  sim.mat <- unifracs
  sim.mat <- 1-sim.mat
  diag(sim.mat) <- 0
  
  ####### impute the values #######
  rownames(sim.mat) = colnames(sim.mat) = rownames(otudata)
  impute.mat <- impute.knn(input = t(otudata), dmat = dmat, sim.mat = sim.mat, 
                           k = 5, sim.cut = 1e-4)
  return(impute.mat$output.imp)
}
