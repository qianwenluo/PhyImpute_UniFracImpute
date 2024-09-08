source("gamma_norm.R")

  set.seed(1234)
  T2D_coeff <- readRDS("T2D_dat2_sim_add_filter_coef.rds")
  otu_real_data_T2D <- read.csv("otu_real_data_T2D.csv", row.names = "X")
  D <- read.csv("D.csv", row.names = "X")
  meta_data_T2D <- read.csv("meta_data_T2D.csv", row.names = "X")
  y_sim = otu_real_data_T2D
  x = meta_data_T2D
  k = 5
  c1 = T2D_coeff
  c1 = c1*5
  c1[1] = T2D_coeff[1]
  c1[length(c1)-6] = 0.2
  c1[length(c1)-11] = -0.3
  
  
  filter <- lapply(X = 1:ncol(y_sim), FUN = function(col_i){
    y = y_sim[,col_i]
    n = length(y)
    nz <- sum(y <= (log10(1.01) + 1e-6))
    pz = 1 - nz/n
    test = pz - 1.96 * sqrt(pz * (1-pz)/n)
    if(nz == n || test <= 0){
      return(0)
    }else{
      return(1)
    }
  })
  
  filter_vec <- which(unlist(filter) == 1)
  y_sim = y_sim[, filter_vec]
  D = D[filter_vec,filter_vec]
  y_preserve <- y_sim
  #apply the imputation method on the simulated matrix
  m = dim(y_sim)[2]
  n = dim(y_sim)[1]
  # 53 * 193
  
  X <- as.matrix(x)
  
 
  close_taxa <- list()
  for(j in 1:m){
    
    close_dist <- D[D[,j] %in% sort(D[,j])[2:(k+1)],j]
    close_taxa_vec = which(D[,j] %in% close_dist[close_dist != max(close_dist)])
    if(length(close_taxa_vec) < k){
      close_taxa_vec <- c(close_taxa_vec, which(D[,j] == max(close_dist))[1:(k-length(close_taxa_vec))])
    }
    close_taxa[[j]] = close_taxa_vec
  }
  #generate design matrix
  p = dim(x)[2]
  if(is.null(p)){
    p = 1
  }
  #generate a row including response and a row of design matrix.
  row_length <- m * k + (n-1) * n + n*p
  idx_set <- matrix(NA, m*n, 2)
  for(i in 1:n){
    for(j in 1:m){
      idx_set[(i-1)*m+j, ] <- c(i, j)
    }
  }
  
  design_mat_gen <- matrix(0, nrow = dim(idx_set)[1], ncol = row_length)
  for(i in 1:dim(idx_set)[1]){
    design_mat_gen[i,] <- design_mat_row_gen2(y_sim, x[1:n,], idx_set[i,1], idx_set[i,2], close_taxa)
  }
  design_mat_gen <- cbind(rep(1, dim(design_mat_gen)[1]), design_mat_gen)
  imputed_value <- design_mat_gen %*% c1
  impute_mat <- y_sim
  for(i in 1:dim(idx_set)[1]){
    impute_mat[idx_set[i,1], idx_set[i,2]] = max(imputed_value[i], log10(1.01))
  }

  y_sim <- impute_mat
  # }
  ############### y_sim is the y_comp, apply mbImpute on the real data, and consider as complete data ##########
  y_sim <- y_sim - 1.5
  #write.csv(y_sim, "Karlsson_simulated.csv")
  ####### learn from learning missing rate empirically ############
  
  otable <- y_preserve
  meta_tab <- meta_data_T2D

  mean_record <- c()
  percentage_record <- c()
  D_vale_record <- c()
  beta_record <- list()
  for(j in 1:dim(otable)[2]){
    result <- gamma_norm_mix(otable[,j], data.matrix(meta_tab))
    beta_record[[j]] = result$cov_par
    mean_record <- c(mean_record, mean(otable[which(result$d < 0.5),j]))
    percentage_record <- c(percentage_record, sum(result$d > 0.5)/dim(otable)[1])
    D_vale_record <- c(D_vale_record, result$d)
  }

  
  # build a map between the mean and percentage of missing (false zeros across samples)
  missing_rate <- function(mean, emp_mean, emp_miss){
    win_len <- (range(emp_mean)[2] - range(emp_mean)[1])/3
    mean_up <- mean + win_len
    mean_lo <- mean - win_len
    sample(emp_miss[emp_mean > mean_lo & emp_mean < mean_up], 1)
  }
  # missing_rate(1, mean_record, percentage_record)
  
  col_mean <- colMeans(y_sim)
  zero_rate <- unlist(lapply(col_mean, FUN = function(x){
    print(x)
    return(missing_rate(x, mean_record, percentage_record))
  }))
  zero_mat <- matrix(NA, nrow = n, ncol = m)
  for(i in 1:m){
    zero_mat[,i] = rbinom(n,1, 1-zero_rate[i])
  }
  sim_tab_zi = y_sim * zero_mat
  sim_tab_zi[sim_tab_zi < log10(1.01)+1e-6] = log10(1.01)
  write.csv(sim_tab_zi, file="simulated.csv")
  write.csv(y_sim, file="complete.csv")
  

