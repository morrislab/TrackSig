make_data_for_vae <- function()
{
  # Count the types of cancer where there is a lot of change
  rescaled = F
  
  training_data <- c()
  training_annotation <- c()
  training_transition_data <- c()
  
  cancer_types <- as.factor(sort(unique(tumortypes[,2])))
  
  for (dir_ in list.files(DIR_RESULTS))
  {
    print(paste0("Cancer type ", dir_))
    for (tumor_id in list.files(paste0(DIR_RESULTS, dir_)))
    {
      tumor_dir <- paste0(DIR_RESULTS, dir_, "/", tumor_id, "/")
      
      if (!file.exists(tumor_dir) || 
            !file.exists(paste0(SAVED_SAMPLES_DIR, tumor_id, ".RData")) || 
            !file.exists(paste0(tumor_dir, "mixtures.mean.csv")))
      {
        next
      }
      
      load(paste0(SAVED_SAMPLES_DIR, tumor_id, ".RData"))
      
      assignments <- NA
      if (file.exists(paste0(tumor_dir, "assignments.txt")))
      {
        assignments <- read.delim(paste0(tumor_dir, "assignments.txt"), header=F, sep = " ")
      }
      
      change_points <- read.delim(paste0(tumor_dir, "changepoints.txt"), header=F, sep = " ")
      phis <- read.delim(paste0(tumor_dir, "phis.txt"), sep=" ", header=F)
      
      mixtures <- read.csv(paste0(tumor_dir, "mixtures.mean.csv"), stringsAsFactors = F)
      active_sigs <- as.numeric(gsub("S([\\d])*", "\\1", mixtures[,1]))
      
      mixtures_all_sigs <- data.frame(matrix(0, ncol=(ncol(mixtures) - 1), nrow=ncol(alex)))
      rownames(mixtures_all_sigs) <- colnames(alex)
      mixtures_all_sigs[mixtures[,1],] <- mixtures[,-1]
      
      new_item_annotation <- data.frame(tumor_id = tumor_id, type = dir_,
                                        time_point = 1:ncol(mixtures[,-1]))
      new_item_annotation <- cbind(new_item_annotation,  phi = toVerticalMatrix(unlist(phis)))
      new_item_annotation <- cbind(new_item_annotation, change_points =  toVerticalMatrix(as.numeric(1:ncol(mixtures[,-1]) %in% (change_points))))
      new_item_annotation <- cbind(new_item_annotation, assignments = toVerticalMatrix(unlist(assignments)))
      new_item_annotation <- new_item_annotation[-nrow(new_item_annotation),]
      
      new_item_data <- data.frame(tumor_id = tumor_id, t(vcf) / apply(t(vcf), 1, sum))
      new_item_data <- cbind(new_item_data, t(mixtures_all_sigs))
      
      new_item_transition_data <- cbind(new_item_data[1:(nrow(new_item_data) -1),], 
                                        new_item_data[2:nrow(new_item_data),-1],
                                        which(cancer_types == dir_))
      
      training_data <- rbind(training_data, new_item_data)
      training_transition_data <- rbind(training_transition_data, new_item_transition_data)
      training_annotation <- rbind(training_annotation, new_item_annotation)
      
      stopifnot(nrow(training_transition_data) == nrow(training_annotation))
    }
  }
  
  training_data <- training_data[,-1]
  training_transition_data <- training_transition_data[,-1]
  training_transition_data.BRCA <- training_transition_data[training_transition_data[,ncol(training_transition_data)] == which(cancer_types == "BRCA"),]
  training_transition_data.PACA <- training_transition_data[training_transition_data[,ncol(training_transition_data)] == which(cancer_types == "PACA"),]
  training_transition_data.PRAD <- training_transition_data[training_transition_data[,ncol(training_transition_data)] == which(cancer_types == "PRAD"),]
  training_transition_data.LIRI <- training_transition_data[training_transition_data[,ncol(training_transition_data)] == which(cancer_types == "LIRI"),]
  
  training_annotation.full <- training_annotation
  training_annotation$type <- sapply(training_annotation$type, function(x) which(toString(x) == cancer_types))
  training_annotation <- training_annotation[,-1]
  
  training_annotation.BRCA <- training_annotation[training_annotation$type == which(cancer_types == "BRCA"),]
  training_annotation.PACA <- training_annotation[training_annotation$type == which(cancer_types == "PACA"),]
  training_annotation.PRAD <- training_annotation[training_annotation$type == which(cancer_types == "PRAD"),]
  training_annotation.LIRI <- training_annotation[training_annotation$type == which(cancer_types == "LIRI"),]
  
  vae_dir = "/Users/yulia/Documents/=Courses=/CSC2541/autograd/vae/"
  order <- read.csv(paste0(vae_dir, "tumor_order.csv"), stringsAsFactors = F, header=F)[,1]
  #order <- sample(1:nrow(training_transition_data))
  #write.table(order, file=paste0(vae_dir, "tumor_order.csv"), col.names=F, row.names = F, sep=",")
  training_transition_data.shuffled <- training_transition_data[order,]
  training_data.shuffled <- training_data[order,]
  training_annotation.shuffled <- training_annotation[order,]
  training_annotation.full.shuffled <- training_annotation.full[order,]
  
  training_transition_data.only_mut_types.shuffled <- training_transition_data.shuffled[c(1:96, 127:222)]
  
  write.table(training_data.shuffled, file=paste0(vae_dir, "training_data.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.shuffled, file=paste0(vae_dir, "training_transition_data_w_cancer_type.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.only_mut_types.shuffled, file=paste0(vae_dir, "training_transition_data.only_mut_types.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.shuffled[,-ncol(training_transition_data)], file=paste0(vae_dir, "training_transition_data.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_annotation.shuffled, file=paste0(vae_dir, "training_annotation.csv"), row.names = F, col.names=F, sep=",")
  write.table(training_annotation.full.shuffled, file=paste0(vae_dir, "training_annotation.full.csv"), row.names = F, col.names=F, sep=",")
  
  
  write.table(training_transition_data.BRCA[,-ncol(training_transition_data)], file=paste0(vae_dir, "training_transition_data.BRCA.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.PACA[,-ncol(training_transition_data)], file=paste0(vae_dir, "training_transition_data.PACA.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.PRAD[,-ncol(training_transition_data)], file=paste0(vae_dir, "training_transition_data.PRAD.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_transition_data.LIRI[,-ncol(training_transition_data)], file=paste0(vae_dir, "training_transition_data.LIRI.csv"), col.names=F, row.names = F, sep=",")
  
  write.table(training_annotation.BRCA, file=paste0(vae_dir, "training_annotation.BRCA.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_annotation.PACA, file=paste0(vae_dir, "training_annotation.PACA.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_annotation.PRAD, file=paste0(vae_dir, "training_annotation.PRAD.csv"), col.names=F, row.names = F, sep=",")
  write.table(training_annotation.LIRI, file=paste0(vae_dir, "training_annotation.LIRI.csv"), col.names=F, row.names = F, sep=",")
  
  write(sapply(cancer_types, toString), file=paste0(vae_dir, "cancer_types.txt"), sep = " ", ncolumns=length(cancer_types))
  
  return(stats)
}



compute_mutation_probabilities <- function(dir_counts = DIR_COUNTS)
{
  MUTATION_TYPES = "pwgs/mut_types/cnvint/"
  PROBABILITIES_DIR = "mutation_probabilities/"
  
  if (!file.exists(PROBABILITIES_DIR)) {
    dir.create(PROBABILITIES_DIR)
  }
  
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 
  
  for (example in tumors)
  {
    print(paste0("Example ", example, " (", which(tumors == example), " out of ", length(tumors), ")"))
    
    pct <- proc.time()
    data_file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData")
    if (file.exists(data_file))
    {
      load(data_file)
    } else {
      print(paste0("Data file ", data_file, " does not exist"))
      next
    }
    
    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, noise_sig)
    } else {
      alex.t <- alex
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    
    mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
    changepoints <- unlist(read.table(paste0(dir_name, "changepoints.txt"), header=F))
    
    mutations <- read.delim(paste0(MUTATION_TYPES, "/", example, ".mut_types.txt"), sep="\t",
                            stringsAsFactors = F)
    colnames(mutations) <- c("chr", "pos", "phi", "ref", "alt", "context")
    exposures <-  as.data.frame(matrix(NA, ncol=ncol(alex.t), nrow=nrow(mutations)))
    colnames(exposures) <- colnames(alex.t)
    prob = data.frame( prob = matrix(NA, ncol=1, nrow=nrow(mutations)))
    
    # Compute probability of mutations based on mixtures at previous time point.
    # As mixtures are computed for sliding window, we will use mixtures that were computed before the mutation
    # first appear in a sliding window
    for (t in 1:nrow(mutations))
    {
      mixture_index <- min(max(1, round(t / (100*gap)) - 1), ncol(mixtures))
      exposures[t,] <- mixtures[,mixture_index]
      
      current_mut <- mutations[t,]
      current_mut_type <- paste(current_mut$ref, current_mut$alt, current_mut$context, sep="_")
      
      mutation_probabilities <- as.matrix(alex.t) %*% t(as.matrix(exposures[t,]))
      prob[t, ] <- mutation_probabilities[current_mut_type,]
    }
    
    mutations <- cbind(mutations, exposures, prob)
    
    write.table(mutations, file=paste0(PROBABILITIES_DIR, "/", example, ".mut_prob.txt"),
                sep = "\t", quote = F, row.names=F)
    
    print(paste("Computed example", example))
  }
}


# prediction_prob_latent_dim_10.kl_weight_decay_time_1.csv
# predictions_latent_dim_10.kl_weight_decay_time_1.csv
# pred_prob <-read.csv("/Users/yulia/Documents/=Courses=/CSC2541/autograd/vae/transition_time_point_VAE_only_mut_types/training_transition_data.only_mut_types/prediction_prob_latent_dim_50.kl_weight_decay_time_1.csv", header=F, stringsAsFactors = F)[,1]
# 
# tumor_id <- "9ddf2119-a222-4fa5-a9f3-0bec7eeea36b"
# indices <- which(training_annotation.full.shuffled$type == "BRCA" & training_annotation.full.shuffled$tumor_id == tumor_id)
# indices <- indices[indices < length(pred_prob)]
# tp_order <-  order(training_annotation.full.shuffled[indices,]$time_point)
# training_annotation.full.shuffled[indices,][tp_order,]
# y <- rep(NA, max(training_annotation.full.shuffled[indices,][tp_order,]$time_point))
# y[training_annotation.full.shuffled[indices,][tp_order,]$time_point] <- pred_prob[indices][tp_order]
# plot(y, type="l")
# 
# 98bb3025-0637-4106-8621-12df7b5d662f 
# 9ddf2119-a222-4fa5-a9f3-0bec7eeea36b
# d8c6d4b8-f279-4edc-aaa3-a1cc266aec4d 
# ddc7377d-82c3-480a-be3c-3d1da52c77d4
# fc8130df-3225-3f96-e040-11ac0d485dfe 
# fc8130df-e399-e34d-e040-11ac0c483279
# fc8130e0-096a-b991-e040-11ac0c48327d 
# fe04d042-a4cc-4a14-8197-415ea40951aa

get_simulations <- function(mixtures)
{
  simutated_counts <- sapply(1:ncol(mixtures), function(i) rmultinom(1, unlist(round(100*mixtures[1,])[i]), alex[,1]))
  for (j in 2:nrow(mixtures))
  {
    tmp <- sapply(1:ncol(mixtures), function(i) rmultinom(1, unlist(round(100*mixtures[j,])[i]), alex[,j]))
    simutated_counts <- simutated_counts + tmp
  }
  return(simutated_counts)
}

simulate_mutations_s1_s5 <- function()
{
  mutation_types <- read.delim("trinucleotide.txt", header=F)
  mutation_types <- paste(mutation_types[,1], mutation_types[,2], mutation_types[,3], sep="_")
  
  n_bootstap = 20
  n_timesteps = 20
  var = 0.03
  phis <- sort(runif(n_timesteps, min=0, max=1), decreasing=T)
  
  simulated_mixtures <- list()
  
  # Simulation 1
  sim_num = 1
  signatures <- c("S1", "S5")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  rand <- 0.8
  true_mixtures[1,] <- rand + runif(n_timesteps, min=-var, max=var)
  true_mixtures[2,] <- 1-true_mixtures[1,]
  simulated_mixtures[[sim_num]] <- true_mixtures
    
  # Simulation 2
  sim_num = 2
  signatures <- c("S1", "S5")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  rand <- 0.8
  true_mixtures[1,] <- rand + runif(n_timesteps, min=-var*2, max=var*2)
  true_mixtures[2,] <- 1-true_mixtures[1,]
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 3
  sim_num = 3
  signatures <- c("S1", "S5")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  rand <- 0.6
  true_mixtures[1,] <- rand + runif(n_timesteps, min=-var*2, max=var*2)
  true_mixtures[2,] <- 1-true_mixtures[1,]
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 4
  sim_num = 4
  signatures <- c("S1", "S5")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",1:5] <- 0.3 + runif(5, min=-var*2, max=var*2)
  true_mixtures["S1",7:12] <- 0.6 + runif(6, min=-var*2, max=var*2)
  true_mixtures["S1",14:20] <- 0.3 + runif(7, min=-var*2, max=var*2)
  true_mixtures["S1",6] <- mean(true_mixtures["S1",5], true_mixtures["S1",7])
  true_mixtures["S1",13] <- mean(true_mixtures["S1",12], true_mixtures["S1",14])
  true_mixtures[2,] <- 1-true_mixtures[1,]
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 5
  sim_num = 5
  signatures <- c("S1", "S5", "S8")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",] <- 0.1 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S5",1:10] <- 0.2 + runif(10, min=-var, max=var)
  true_mixtures["S5",14:20] <- 0.7 + runif(7, min=-var, max=var)
  true_mixtures["S5",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S5",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S5",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S8",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 6
  sim_num = 6
  signatures <- c("S1", "S5", "S8")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S5",1:10] <- 0.3 + runif(10, min=-var, max=var)
  true_mixtures["S5",14:20] <- 0.4 + runif(7, min=-var, max=var)
  true_mixtures["S5",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S5",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S5",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S8",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 7
  sim_num = 7
  signatures <- c("S1", "S5", "S8")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S5",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S1",1:10] <- 0.3 + runif(10, min=-var, max=var)
  true_mixtures["S1",14:20] <- 0.6 + runif(7, min=-var, max=var)
  true_mixtures["S1",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S1",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S1",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S8",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 8
  sim_num = 8
  signatures <- c("S1", "S5", "S3")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",] <- 0.1 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S5",1:10] <- 0.2 + runif(10, min=-var, max=var)
  true_mixtures["S5",14:20] <- 0.7 + runif(7, min=-var, max=var)
  true_mixtures["S5",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S5",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S5",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S3",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 9
  sim_num = 9
  signatures <- c("S1", "S5", "S3")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S5",1:10] <- 0.3 + runif(10, min=-var, max=var)
  true_mixtures["S5",14:20] <- 0.4 + runif(7, min=-var, max=var)
  true_mixtures["S5",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S5",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S5",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S3",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 10
  sim_num = 10
  signatures <- c("S1", "S5", "S3")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S5",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S1",1:10] <- 0.3 + runif(10, min=-var, max=var)
  true_mixtures["S1",14:20] <- 0.6 + runif(7, min=-var, max=var)
  true_mixtures["S1",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S1",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S1",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S3",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 11
  sim_num = 11
  signatures <- c("S1", "S5", "S8")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S5",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S1",] <- 0.4 + runif(n_timesteps, min=-var, max=var)
  
  true_mixtures["S8",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 12
  sim_num = 12
  signatures <- c("S1", "S5", "S3")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S5",] <- 0.2 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S1",] <- 0.4 + runif(n_timesteps, min=-var, max=var)

  true_mixtures["S3",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  # Simulation 13
  sim_num = 13
  signatures <- c("S1", "S5", "S2")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S5",] <- 0.5 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S1",] <- 0.1 + runif(n_timesteps, min=-var, max=var)
  
  true_mixtures["S2",] <- 1 - sapply(1:n_timesteps, function(i) {true_mixtures["S1",i]+ true_mixtures["S5",i]})
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  
  # Simulation 14
  sim_num = 14
  signatures <- c("S1", "S5", "S3", "S8")
  
  true_mixtures <- data.frame(matrix(NA, ncol=n_timesteps, nrow=length(signatures)))
  rownames(true_mixtures) <- signatures
  colnames(true_mixtures) <- round(phis, 3)
  
  true_mixtures["S1",] <- 0.1 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S5",] <- 0.1 + runif(n_timesteps, min=-var, max=var)
  true_mixtures["S8",1:10] <- 0.3 + runif(10, min=-var, max=var)
  true_mixtures["S8",14:20] <- 0.6 + runif(7, min=-var, max=var)
  true_mixtures["S8",12] <- mean(true_mixtures["S5",10], true_mixtures["S5",14])
  true_mixtures["S8",11] <- mean(true_mixtures["S5",12], true_mixtures["S5",14])
  true_mixtures["S8",13] <- mean(true_mixtures["S5",10], true_mixtures["S5",12])
  
  true_mixtures["S3",] <- 1 - apply(true_mixtures, 2, sum, na.rm=T)
  simulated_mixtures[[sim_num]] <- true_mixtures
  
  
  for (i in 1:length(simulated_mixtures))
  {
    suppressWarnings(dir.create(paste0(DIR_RESULTS, "simulation", i)))
    phis <- as.numeric(colnames(simulated_mixtures[[i]]))
    plot_signatures(simulated_mixtures[[i]], plot_name=paste0(DIR_RESULTS, "simulation", i, "/true_mixtures.pdf"), 
                    phis = phis)
    
    write.csv(simulated_mixtures[[i]], file=paste0(DIR_RESULTS, "simulation", i, "/true_mixtures.csv"))
    write(rownames(simulated_mixtures[[i]]), file=paste0(DIR_RESULTS, "simulation", i, "/true_mixture_sigs.txt"),
          ncolumns = nrow(simulated_mixtures[[i]]))
    
    counts <- get_simulations(simulated_mixtures[[i]])
    rownames(counts) <- mutation_types
    colnames(counts) <- phis
    counts_header <- t(rbind(paste0("Simulation_", i), round(phis, 3), counts))
    
    write.table(counts_header, file=paste0(DIR_COUNTS, "simulation", i, ".phi.txt"), 
                row.names = F, col.names=F, quote=F, sep="\t")
  }
  
  #   counts.m <- melt(counts)
  #   bootstrapped_counts <- c()
  #   for (i in 1:n_bootstap)
  #   {
  #     bootstrapped <- sample(1:nrow(counts.m), nrow(counts.m), replace=T)
  #     bootstrapped_mut <- counts.m[bootstrapped, ]
  #     bootstrapped_mut <- bootstrapped_mut[order(bootstrapped_mut$Var2, decreasing = T),]
  #     
  #     for (j in 1:n_timesteps)
  #     {
  #       #bootstrapped_mut[(100*(i-1)):(100*i),]
  #       bootstrapped_counts <- rbind()
  #     }
  #   }
  
}