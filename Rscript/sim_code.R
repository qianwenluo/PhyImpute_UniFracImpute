
# https://github.com/biobakery/biobakery/wiki/SparseDOSSA
library(sparseDOSSA)
library(stringi)

n.sample <- 100 #number of samples
n.microbes <- 200 #number of features
spike.perc <- 0.2
#spikeStrength <- "1.0" #effect size
n.metadata <- 1
metadata <- matrix(rbinom(n=100, size = 1, prob = 0.5), nrow = 1, ncol = 100)
# spiked-in metadata (which metadata to spike-in)

# simulate data
simulated_data <- sparseDOSSA(number_features = n.microbes, 
                              number_samples = n.sample,
                              percent_spiked = spike.perc, 
                              UserMetadata = metadata,
                              seed = 1234,
                              datasetCount = 1, 
                              noZeroInflate = FALSE)
# gather output from sparseDOSSA
simulated_result <- as.data.frame(simulated_data$OTU_count)
rownames(simulated_result) <- simulated_result$X1
simulated_result <- simulated_result[-1,-1]
colnames(simulated_result) <- paste('Sample',1:ncol(simulated_result),sep='')
data <- as.matrix(simulated_result[-c((n.metadata+1):(2*n.microbes+n.metadata)),])
data <- data.matrix(data)
class(data) <- "numeric"
truth <- c(unlist(simulated_data$truth))
truth<-truth[!stri_detect_fixed(truth,":")]
truth<-truth[(5+n.metadata):length(truth)]
truth<-as.data.frame(truth)
spikeCount <- as.character(length(1))
significant_features<-truth[seq(1, (as.numeric(spikeCount)+1)*(n.microbes*spike.perc), (as.numeric(spikeCount)+1)),]
significant_features<-as.vector(significant_features)

# Extract Features #
features <- as.data.frame(t(data[-c(1:n.metadata),]))
# Extract Metadata #
metadata<-as.data.frame(data[1,])
colnames(metadata)<-rownames(data)[1]
# Mark True Positive Features #
wh.TP = colnames(features) %in% significant_features
colnames(features)<-paste("Feature", 1:n.microbes, sep = "")
newname = paste0(colnames(features)[wh.TP], "_TP")
colnames(features)[wh.TP] <- newname;
colnames(features)[grep('TP', colnames(features))]
#### OTU_table: feature
#### metadata: metadata

otu_table <- t(features)
metadata

write.csv(otu_table, file="sim1.csv")
write.csv(metadata, file ="sim1_meta.csv")
