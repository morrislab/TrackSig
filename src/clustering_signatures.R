cluster_by_overall_mutations <- function()
{  
  mutation_counts <- c()
  for (example in tumortypes[,1])
  {
    # Favourite example of LUSC: 1c3df485-8e75-4378-87f6-c6463a520624
    
    pct <- proc.time()
    list[vcfFile, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(example, dir_counts, tumortypes)
    
    if (is.null(vcfData))
      next
    if (nrow(vcfData) < 6) # Plots with less than 6 lines of data are meaningless so ignored
      next
    
    new_item <- c(example, acronym, apply(vcfData, 2, sum))
    mutation_counts <- rbind(mutation_counts, new_item)
  }
  
  mutation_counts.data <- matrix(as.numeric(mutation_counts[,-c(1,2)]), nrow=nrow(mutation_counts))
  rownames(mutation_counts.data) <- mutation_counts[,1]
  
  #scaling 
  mutation_counts.data.scaled <- mutation_counts.data / apply(mutation_counts.data, 1, max)
  
  mutation_counts.data.log <- log(mutation_counts.data)
  mutation_counts.data.log[mutation_counts.data.log == -Inf] <- NA
  
  mutation_counts.data.scaled.log <- log(mutation_counts.data.scaled)
  mutation_counts.data.scaled.log[mutation_counts.data.scaled.log == -Inf] <- NA
  
  pdf(file=paste0(DIR_RESULTS, "/cluster_by_log_mutation_counts.pdf"), width=20, height=20)
  col_color = toVerticalMatrix(mutation_counts[,2])
  colnames(col_color) <- "Cancer type"
  aheatmap(t(mutation_counts.data),  Rowv = NA, labCol=NULL, labRow = rownames(alex), 
           annCol=col_color)
  dev.off()

  pdf(file=paste0(DIR_RESULTS, "/cluster_by_log_mutation_counts_rescaled.pdf"), width=20, height=20)
  col_color = toVerticalMatrix(mutation_counts[,2])
  colnames(col_color) <- "Cancer type"
  aheatmap(t(mutation_counts.data.scaled),  Rowv = NA, labCol=NULL, labRow = rownames(alex), 
           annCol=col_color, annColors = list("Cancer type" = COLORS2))
  dev.off()
  
  pdf(file=paste0(DIR_RESULTS, "/cosmic_signatures_mutation_counts.pdf"), width=10, height=15)
  alex.scaled <- t(t(alex) / apply(alex, 2, max))
  aheatmap(alex.scaled,  Rowv = NA, labCol=NULL, labRow = rownames(alex), Colv=NA)
  dev.off()
}

cluster_per_tissue_by_mean_and_sd <- function()
{
  rescaled = T
  
  for (dir_ in list.files(DIR_RESULTS))
  {
    file_summary <- paste0(DIR_RESULTS, dir_, "/", "overall_summary.csv")
    if (!file.exists(file_summary))
    {
      print(paste0("No overall summary for ", dir_))
      next
    }
    
    overall_summary <- read.csv(file_summary)
    
    if (nrow(overall_summary) < 3)
      next
    
    col_name <- "mean_sigs"
    if (rescaled)
      col_name <- paste0(col_name, ".rescaled")
    overall_summary.means <- overall_summary[,grep(col_name, colnames(overall_summary), value=T)]
    
    sig_names <- gsub(paste0(col_name, ".(S[0-9]+)$"), "\\1", colnames(overall_summary.means))
    
    col_name <- "sd_sigs"
    if (rescaled)
      col_name <- paste0(col_name, ".rescaled")
    overall_summary.sd <- overall_summary[,grep(col_name, colnames(overall_summary), value=T)]
    
    mydistfun <- function(x) { dist(overall_summary.means, method = "euclidean") + dist(overall_summary.sd, method = "euclidean")}
    d <- mydistfun(NA)
    
    wss <- c()
    max.cl <- min(15, (nrow(overall_summary.means)-1))
    for (i in 1:max.cl) {
      wss[i] <- sum(kmeans(overall_summary.means, centers=i)$withinss)
    }
    n_clusters <- which(abs((wss[2:length(wss)] - wss[1:(length(wss)-1)])/wss[1:(length(wss)-1)]) < 0.1)
    if (length(n_clusters) == 0) {
      n_clusters = max.cl
    } else {
      n_clusters = min(n_clusters)
    }
    
    # determine number of clusters based on sillouette
    #     wss <- c()
    #     for (i in 1:max.cl) {
    #       wss[i] <- pam(overall_summary.means, i)$silinfo$ avg.width
    #     }
    #     n_clusters <- which.max(wss)
    
    d.clustered <- hclust(d) 
    groups <- cutree(d.clustered, k=n_clusters)
    groups <- LETTERS[groups]
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_clustered_by_means.pdf"), width=10, height=10)
    aheatmap(t(overall_summary.means),  Rowv = NA, labRow = sig_names, 
             labCol=sapply(overall_summary$tumor_id, toString), 
             annCol=list(Clusters = groups), 
             annColors=list(Clusters=brewer.pal(9,"Set1")),
             distfun = mydistfun)
    dev.off()
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_barplot_clustered_by_means.pdf"), width=7, height=7)
    exposures_in_cluster <- aggregate(overall_summary.means, by=list(groups), FUN=mean)
    barplot(t(as.matrix(exposures_in_cluster[,-1])), main="Signature exposures in clusters",
            names.arg = exposures_in_cluster[,1], col=COLORS2, ylim=c(0,1.6))
    legend("top", sig_names, fill=COLORS2, horiz = T, cex=0.6)
    dev.off()
    
  }
}


cluster_per_tissue_by_trajectories <- function()
{
  library(dtw)
  rescaled = T
  
  for (dir_ in list.files(DIR_RESULTS))
  {
    sig_indices <- which(as.logical(active_signatures[active_signatures$acronym == dir_,3:(ncol(active_signatures)-1)]))
    active_sigs_per_type <- colnames(active_signatures[,3:(ncol(active_signatures)-1)][sig_indices])
    
    if (length(active_sigs_per_type) == 0)
      next
    
    sig_timeseries <- list()
    for (sig in active_sigs_per_type)
    {
      for (file_ in  list.files(paste0(DIR_RESULTS, dir_, "/")))
      {
        if (!file.exists(paste0(DIR_RESULTS, dir_, "/", file_, "/")))
          next
        
        mixtures_file <- paste0(DIR_RESULTS, dir_, "/", file_, "/mixtures.csv")
        if (!file.exists(mixtures_file))
        {
          print(paste0("No mixtures for ", paste0(DIR_RESULTS, dir_, "/", file_, "/mixtures.csv")))
          next
        }
        
        mixtures <- read.csv(mixtures_file)
        sig_timeseries[[sig]][[file_]] <- unlist(mixtures[mixtures$X == sig,-1])
      }
    }
    
    similarities <- list()
    for (sig in active_sigs_per_type)
    {
      n_samples <- length(sig_timeseries[[sig]])
      similarities[[sig]] <- matrix(NA, ncol=n_samples, nrow=n_samples)
      for (i in 1:n_samples)
        for (j in 1:i)
        {
          similarities[[sig]][i,j] <- dtw(sig_timeseries[[sig]][[i]], sig_timeseries[[sig]][[j]], 
                                          keep.internals=T, open.end=T)$normalizedDistance
          similarities[[sig]][j,i] <- similarities[[sig]][i,j]
        }
    }
    
    W = tryCatch(SNF(similarities),  error = function(e) print(e))
    if ("error" %in% class(W))
    {
      print(paste("SNF error with", dir_))
      next
    }
    W.na <- W
    for (i in 1:nrow(W)) {
      W.na[i,i] <- NA
    }
    #       aheatmap(W.na)
    
    C = estimateNumberOfClustersGivenGraph(W, 2:15)[[1]] 
    groups = spectralClustering(W,C)
    groups <- LETTERS[groups]
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_clustered_by_SNF.pdf"), width=10, height=10)
    aheatmap(W.na[,order(groups)],  #Rowv = NA, Colv=NA,
             labCol=names(sig_timeseries[[1]])[order(groups)], labRow=NA,
             annCol=list(Clusters = groups[order(groups)]), 
             annColors=list(Clusters=brewer.pal(9,"Set1")))
    dev.off()
    
    W.na.normalized <- (W.na - min(W.na, na.rm = T)) / ( max(W.na, na.rm = T) -  min(W.na, na.rm = T))
    
    hr <- hclust(as.dist(1-W.na.normalized), method="complete")
    cl_hclust <- cutree(hr, h=max(hr$height/1.5))
    clusterCols <- rainbow(length(unique(cl_hclust)))
    myClusterSideBar <- clusterCols[cl_hclust]
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_hclust_based_on_similarities.pdf"), width=10, height=10)
    aheatmap(W.na.normalized,  #Rowv = as.dendrogram(hr), Colv=NA,
             labCol=names(sig_timeseries[[1]]), labRow=NA,
             annCol=list(Clusters=cl_hclust),
             annColors=list(Clusters=rainbow(9)))
    dev.off()
    
    
    file_summary <- paste0(DIR_RESULTS, dir_, "/", "overall_summary.csv")
    if (!file.exists(file_summary))
    {
      print(paste0("No overall summary for ", dir_))
      next
    }
    
    overall_summary <- read.csv(file_summary)
    
    if (nrow(overall_summary) < 3)
      next
    
    col_name <- "mean_sigs"
    if (rescaled)
      col_name <- paste0(col_name, ".rescaled")
    overall_summary.means <- overall_summary[,grep(col_name, colnames(overall_summary), value=T)]
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_barplot_clustered_by_SNF.pdf"), width=7, height=7)
    exposures_in_cluster <- aggregate(overall_summary.means, by=list(groups), FUN=mean)
    barplot(t(as.matrix(exposures_in_cluster[,-1])), main="Signature exposures in clusters",
            names.arg = exposures_in_cluster[,1], col=COLORS2, ylim=c(0,1.6))
    legend("top", active_sigs_per_type, fill=COLORS2, horiz = T, cex=0.6)
    dev.off()
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_barplot_hclust_based_on_similarities.pdf"), width=7, height=7)
    exposures_in_cluster <- aggregate(overall_summary.means, by=list(cl_hclust), FUN=mean)
    barplot(t(as.matrix(exposures_in_cluster[,-1])), main="Signature exposures in clusters",
            names.arg = exposures_in_cluster[,1], col=COLORS2, ylim=c(0,1.6))
    legend("top", active_sigs_per_type, fill=COLORS2, horiz = T, cex=0.6)
    dev.off()
  
    # Use clustering (affinity propagation or k-meiloids)
    library("apcluster")
    ap <- apcluster(W)
    ap_clusters <- labels(ap, type="enum")
    # ap@clusters
    # ap@exemplars
    exemplar_ids <- names(sig_timeseries[[1]])[ap@exemplars]
    
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_clustered_by_affinity_propagation.pdf"), width=10, height=10)
    aheatmap(W.na[,order(ap_clusters)],  Rowv = NA, Colv=NA,
             labCol=names(sig_timeseries[[1]]), labRow=NA,
             annCol=list(Clusters = ap_clusters[order(ap_clusters)]), 
             annColors=list(Clusters=brewer.pal(9,"Set1")))
    dev.off()
    
    # Barplots for exemplars
    pdf(file=paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_barplot_clustered_by_affinity_propagation.pdf"), width=7, height=7)
    exposures_in_cluster <- overall_summary.means[ap@exemplars, ]
    barplot(t(as.matrix(exposures_in_cluster)), main="Signature exposures in clusters",
            names.arg = 1:length(ap@exemplars), col=COLORS2, ylim=c(0,1.6))
    legend("top", active_sigs_per_type, fill=COLORS2, horiz = T, cex=0.6)
    dev.off()
    
    for (i in 1:length(exemplar_ids)) {
      exem <- exemplar_ids[i]
      exemplar_mixtures <- read_mixtures(paste0(DIR_RESULTS, dir_, "/", exem, "/mixtures.csv"))
      
      plot_name <- paste0(DIR_RESULTS, dir_, "/", dir_, "sigs_clustered_by_affinity_propagation_cluster", i, ".pdf")
      list[g,data] <- plot_signatures(exemplar_mixtures, plot_name=plot_name, save=T)
      
      cluster <- names(sig_timeseries[[1]])[ap@clusters[[i]]]
      dir_name <- paste0(DIR_RESULTS, dir_, "/", dir_, "_cluster", i)
      suppressWarnings(dir.create(dir_name))
      for (sample in cluster) {
        m <- read_mixtures(paste0(DIR_RESULTS, dir_, "/", sample, "/mixtures.csv"))
        plot_signatures(m, plot_name=paste0(dir_name, "/", sample, ".pdf"), save=T)
#         dtw_output <- tryCatch(multi_dtw(m, exemplar_mixtures), error=function(e) NULL)
#         if (is.null(dtw_output))
#           next
#         m.rescaled_to_exemplar <- rescale_to_reference(m, dtw_output)
#         
#         melted <- data.frame(rownames(m),m.rescaled_to_exemplar)
#         colnames(melted) <- colnames(data)[1:(dtw_output$jmin+1)]
#         melted <- melt(cbind(Signatures=rownames(m),melted), id.vars = "Signatures")
#         
#         g <- g + geom_line(data=melted, size=0.7, alpha=0.3,
#                            aes(x=variable, y=value, group = Signatures, color = Signatures))
      }
      #plot(g)
    }
  }
}

multi_dtw <- function(query, reference) {
  stopifnot(nrow(query) == nrow(reference))
  
  D <- data.frame(matrix(NA, ncol=ncol(query), nrow=nrow(reference)))
  for (i in 1:ncol(query))
    for (j in 1:ncol(reference))
    {
      D[j,i] <- sum(sapply(1:nrow(query), function(k) abs(query[k,i] - reference[k, j])))
    }
  return(dtw(t(D), open.end=T))
}

rescale_to_reference <- function(m, dtw_output) {
  ind <- aggregate(dtw_output$index1, by = list(dtw_output$index2), FUN=list)[,2]
  rescaled <- sapply(ind, function(i) apply(toVerticalMatrix(m[,i]), 1, mean))
  colnames(rescaled) <- 1:ncol(rescaled)
  return(rescaled)
}