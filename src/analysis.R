
# DIR_COUNTS = "counts_phis_1682samples/"
# DIR_RESULTS = "Examples_1682/"
# dir.create(DIR_RESULTS)

get_tumor_types_with_big_change <- function(sig_change_threshold = 0.2, sig_max_threshold = 0.1)
{
  # Count the types of cancer where there is a lot of change
  rescaled = F
  
  tumor_types_with_big_change <- list()
  for (dir_ in list.files(DIR_RESULTS))
  {
    file_summary <- paste0(DIR_RESULTS, dir_, "/", "overall_summary.csv")
    if (!file.exists(file_summary))
    {
      print(paste0("No overall summary for ", dir_))
      next
    }
    
    info = file.info(file_summary)
    empty = nrow(info[info$size < 5, ]) != 0
    
    if (empty)
    {
      print(paste0("No overall summary for ", dir_))
      next
    }
    
    overall_summary <- read.csv(file_summary)
    
    if (nrow(overall_summary) < 3)
      next
    
    col_name <- "max_change"
    if (rescaled)
      col_name <- paste0(col_name, ".rescaled")
    
    overall_summary.max_change <- overall_summary[,col_name]
    tumors_large_change_ids <- which(overall_summary.max_change > sig_change_threshold)
    tumors_large_change <- sapply(overall_summary$tumor_id[tumors_large_change_ids], toString)
    
    ratio_tumors_with_high_change <- length(tumors_large_change_ids) / length(overall_summary.max_change)
    if (ratio_tumors_with_high_change > 0.2 &  length(overall_summary.max_change) > 10) {
      #print(paste0(dir_, ": ratio tumors with high change: ", round(ratio_tumors_with_high_change, digits=2), ". Total samples:  ", length(overall_summary.max_change), " Mean change: ", round(mean(overall_summary.max_change),3)))
      tumor_types_with_big_change[[dir_]] <- list(tumors_large_change = tumors_large_change, ratio_tumors_with_high_change = ratio_tumors_with_high_change) 
    }

  	tumors_small_change_ids <- which(overall_summary.max_change < sig_max_threshold)
  	tumors_small_change <- sapply(overall_summary$tumor_id[tumors_small_change_ids], toString)
  	
  	ratio_tumors_with_small_change <- length(tumors_small_change_ids) / length(overall_summary.max_change)
  
  	if (ratio_tumors_with_high_change > 0.2 &  length(overall_summary.max_change) > 10) {
  	  print(paste0(dir_, ": ratio tumors with small change: ", round(ratio_tumors_with_small_change, digits=2), ". Total samples:  ", length(overall_summary.max_change), " Mean change: ", round(mean(overall_summary.max_change),3)))
  	}
  }
  
  return(tumor_types_with_big_change)
}


get_statistics <- function()
{
  # Count the types of cancer where there is a lot of change
  rescaled = F
  sig_one <- "S1"
  
  stats <- list()
  changing_sigs <- list()
  dominant_sigs <- list()
  dist_to_transition <- list()
  dist_to_transition.random <- list()
  
  for (dir_ in list.files(DIR_RESULTS))
  {
    print(paste0("Cancer type ", dir_))
    for (tumor_id in list.files(paste0(DIR_RESULTS, dir_)))
    {
      tumor_dir <- paste0(DIR_RESULTS, dir_, "/", tumor_id, "/")
      
      if (!file.exists(tumor_dir))
      {
        next
      }
      
      if (!file.exists(paste0(tumor_dir, "mixtures.err.csv")))
      {
        if (file.exists(paste0(tumor_dir, "mixtures.csv")))
        {
          print(paste(tumor_id, ": mixtures.err.csv is absent"))
        }
        next
      }
      
      stderr <- read.csv(paste0(tumor_dir, "mixtures.err.csv"))[,-1]
      stderr.max <- max(stderr)
      stderr.mean <- mean(apply(stderr, 1, mean))
      
      assignments <- transition_points <- c()
      n_subclones <- 1
      if (file.exists(paste0(tumor_dir, "assignments.txt")))
      {
        assignments <- read.delim(paste0(tumor_dir, "assignments.txt"), header=F, sep = " ")
        transition_points <- which(as.logical(diff(unlist(assignments))))
        n_subclones <- max(as.numeric(names(table(unlist(assignments)))))
      }
      
      change_points <- read.delim(paste0(tumor_dir, "changepoints.txt"), header=F, sep = " ")
      phis <- read.delim(paste0(tumor_dir, "phis.txt"), sep=" ", header=F)
      max_phi = max(phis)
      min_phi = min(phis)
      
      summary <- read.csv(paste0(tumor_dir, "summary.csv"))
      
      mixtures <- read.csv(paste0(tumor_dir, "mixtures.mean.csv"))
      rownames(mixtures) <- mixtures[,1]
      mixtures <- mixtures[,-1]
      mixtures_plus <- mixtures + stderr
      mixtures_minus <- mixtures - stderr
      
      max_indices <- apply(mixtures, 1, which.max)
      min_indices <- apply(mixtures, 1, which.min)
      distance_between_error_bars <- sapply(1:nrow(mixtures), function(x) {mixtures_plus[x,max_indices[x]]}) - 
                                      sapply(1:nrow(mixtures), function(x) {mixtures_minus[x,min_indices[x]]})
      names(distance_between_error_bars) <- rownames(mixtures)
      
      
      sigs_significant_change <- which(distance_between_error_bars > 0.1)
      n_significant_change <- length(sigs_significant_change)
      changing_sigs[[tumor_id]] <- names(sigs_significant_change)
      sigs_significant_change <- toString(names(sigs_significant_change))
      
      dominant_sigs[[tumor_id]] <- names(which(apply(mixtures,1,mean) > 0.1))
      
      var_small_sig <- mean(apply(mixtures[,-1], 2, function(x) {sum(x[x < 0.05])}))
      n_small_sig <- mean(apply(mixtures[,-1], 2, function(x) {sum(x < 0.05)}))
      n_sigs <- nrow(mixtures)   
      n_timesteps <- ncol(mixtures)  - 1
      s1_mean_exposure <- mean(unlist(mixtures[sig_one,]))
      s1_change <- max(unlist(mixtures[sig_one,])) - min(unlist(mixtures[sig_one,]))
      n_changepoints = length(change_points)
      n_transition_points <- length(transition_points)
      top_signature = names(which.max(apply(mixtures, 1, max) - apply(mixtures, 1, min)))
      max_change = max(apply(mixtures, 1, max) - apply(mixtures, 1, min))
      max_change_no_s1_s5 = max(apply(mixtures[-which(!(rownames(mixtures) %in% c(sig_one, "S5"))),], 1, max) - apply(mixtures[-which(!(rownames(mixtures) %in% c(sig_one, "S5"))),], 1, min))
      cor_s1_s5 = NA
      if (sum(c(sig_one,"S5") %in% rownames(mixtures)) == 2)
      {
        cor_s1_s5 <- cor(unlist(mixtures[c(),]), unlist(mixtures[c("S5"),]))
      }
      
      n_changepoints_unsorted = max_change_unsorted = NA
      if (file.exists(paste0(tumor_dir, "mixtures.mean.unsorted.csv")))
      {
        mixtures.unsorted <- read_mixtures(paste0(tumor_dir, "mixtures.mean.unsorted.csv"))
        max_change_unsorted = max(apply(mixtures.unsorted, 1, max) - apply(mixtures.unsorted, 1, min))
        
        n_cp.unsorted <- c()
        i <- 1
        while (file.exists(paste0(tumor_dir, "changepoints.bootstrap_", i, ".unsorted.txt")))
        {
          changepoints.unsorted <- unlist(read.table(paste0(tumor_dir, "changepoints.bootstrap_", i, ".unsorted.txt"), header=F))
          n_cp.unsorted <- c(n_cp.unsorted, length(changepoints.unsorted))
          i = i+1
        }
        n_changepoints_unsorted = mean(n_cp.unsorted)
      }
      
      
      mean_distance_to_changepoint <- -Inf
      max_distance_to_changepoint <--Inf
      transition_point_with_max_dist <- -Inf
      distance_at_transition_with_max_dist <- -Inf
      if (length(transition_points) > 0)
      {
        corr_change_vs_transition <- outer(unlist(change_points), c(transition_points), "-")
        mean_distance_to_changepoint <- mean(apply(abs(corr_change_vs_transition), 2, min))
        max_distance_to_changepoint <- max(apply(abs(corr_change_vs_transition), 2, min))
        transition_point_with_max_dist <- which.max(apply(abs(corr_change_vs_transition), 2, min))
        distance_at_transition_with_max_dist <- corr_change_vs_transition[,transition_point_with_max_dist][which.min(abs(corr_change_vs_transition[,transition_point_with_max_dist]))]  
        
        dist_to_transition[[tumor_id]] <- sapply(1:ncol(corr_change_vs_transition), function(x) {corr_change_vs_transition[apply(abs(corr_change_vs_transition), 2, which.min)[x],x]})
        
        random_change_points <- sample(1:n_timesteps, n_changepoints)
        corr_random_change_vs_transition <- outer(unlist(random_change_points), c(transition_points), "-")
        transition_point_with_max_dist.random <- which.max(apply(abs(corr_random_change_vs_transition), 2, min))
        distance_at_transition_with_max_dist.random <- corr_random_change_vs_transition[,transition_point_with_max_dist.random][which.min(abs(corr_random_change_vs_transition[,transition_point_with_max_dist.random]))]  
        
        dist_to_transition.random[[tumor_id]] <- sapply(1:ncol(corr_random_change_vs_transition), function(x) {corr_random_change_vs_transition[apply(abs(corr_random_change_vs_transition), 2, which.min)[x],x]})
        
      }
      
      new_item <- data.frame(tumor_id = tumor_id, type = dir_, 
                             stderr.max = stderr.max, stderr.mean = stderr.mean,
                             n_subclones = n_subclones, 
                             var_small_sig = var_small_sig, n_small_sig = n_small_sig,
                             n_sigs = n_sigs,
                             n_timesteps =n_timesteps,
                             n_significant_change = n_significant_change, 
                             s1_mean_exposure = s1_mean_exposure,
                             s1_change = s1_change,
                             sigs_significant_change = sigs_significant_change,
                             mean_distance_to_changepoint = mean_distance_to_changepoint,
                             max_distance_to_changepoint = max_distance_to_changepoint,
                             distance_at_transition_with_max_dist = distance_at_transition_with_max_dist,
                             distance_at_transition_with_max_dist.random = distance_at_transition_with_max_dist.random,
                             n_changepoints = n_changepoints,
                            n_transition_points = n_transition_points,
                            top_signature = top_signature,
                            max_change = max_change,
                            max_change_no_s1_s5 = max_change_no_s1_s5,
                            max_phi = max_phi,
                            min_phi = min_phi,
                            cor_s1_s5 = cor_s1_s5,
                            max_change_unsorted = max_change_unsorted,
                            n_changepoints_unsorted = n_changepoints_unsorted,
                            stringsAsFactors = F)
        
      stats <- rbind(stats, new_item)
    } 
  }
    
  write.csv(stats, file=paste0(DIR_RESULTS, "Stats.csv"), row.names = F)
  
  return(stats)
}






# Plot the histograms over time for tumors where change is more than 0.2
plot_timeseries_histograms <- function()
{
  attach(mtcars)
  tumor_types_with_big_change <- get_tumor_types_with_big_change()
  for (tumor_type in names(tumor_types_with_big_change)) {
    for (example in tumor_types_with_big_change[[tumor_type]]$tumors_large_change) {
      # Plot VAF for 2290 tumor
      # example <- "2290b078-6a5b-4c83-9dfb-b525bbf14e4e"
      list[vcfFile, vcfData, phis, acronym, dir_name] <- extract_data_for_example(example, DIR_COUNTS)
      
      data_method <- "sliding400"
      window_size=400
      shift = window_size/100
      vcf <- get_sliding_window_data(vcfData, shift=shift)
      phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
      phis_sliding_window <- phis_sliding_window / (shift+1)
      colnames(vcf) <- round(phis_sliding_window, 3)
      
      prior = NULL
      
      list[alex.t.only_known, matched_type, acronym_] <- get_signatures_for_current_sample(vcfFile, active_signatures.our_samples, alex)
      if (is.null(alex.t.only_known)) {
        alex.t.only_known <- alex
      }
      
      mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
      changepoints <- unlist(read.table(paste0(dir_name, "changepoints.txt"), header=F))
      if (length(which(diff(changepoints) == 1)) != 0) {
        changepoints <- changepoints[-(1+which(diff(changepoints) == 1))]
      }
      
      intervals <- c(1, changepoints,  ncol(vcf))
      for (i in 1:(length(intervals)-1))
      {
        # take into account that mixtures were counted on windows of 400
        make_histogram(t(vcfData[max(0, intervals[i]-3):(intervals[i+1]+3),]), paste0(DIR_RESULTS, acronym_, "/", example, "/hist_interval_", i, ".pdf"))
      }
      
      window_size=200
      shift = window_size/100
      vcf <- get_sliding_window_data(vcfData, shift, gap=2)
      phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
      phis_sliding_window <- phis_sliding_window / (shift+1)
      colnames(vcf) <- round(phis_sliding_window, 3)
      
      covered_by_changepoints <- convert_window_to_indices(nrow(vcfData), 400, 100)[changepoints, ]
      intervals_for_change_histogram <- convert_window_to_indices(nrow(vcfData), 200, 200)
      
      converted_changepoints <- find_indices_with_intersection(covered_by_changepoints, intervals_for_change_histogram)
      
      pdf(paste0(DIR_RESULTS, acronym_, "/", example, "/mutation_counts_histogram_timeline.pdf"), width=1.5* ncol(vcf), height=4)
      layout(matrix(1:ncol(vcf), 1, ncol(vcf), byrow = TRUE))
      
      for (i in 1:ncol(vcf)) {
        make_histogram_plot(toVerticalMatrix(vcf[,i]), horiz =T)
        if (i %in% converted_changepoints) {
          rect(0,-1,0.05,120, border = "red", lwd=4)
        }
      }
      dev.off()
      
      changes <- vcf[,2:ncol(vcf)] - vcf[,1:(ncol(vcf) - 1)]
      
      pdf(paste0(DIR_RESULTS, acronym_, "/", example, "/mutation_count_changes_histogram_timeline.pdf"), width=1.5* ncol(changes), height=4)
      layout(matrix(1:ncol(changes), 1, ncol(changes), byrow = TRUE))
      
      for (i in 1:ncol(changes)) {
        make_histogram_plot(toVerticalMatrix(changes[,i]), horiz =T)
        if (i %in% converted_changepoints) {
          rect(-0.05,-1,0.05,120, border = "red", lwd=4)
        }
      }
      dev.off()
      
      # make a histogram of probabilities of a mutational signatures
      #make_histogram(toVerticalMatrix(alex.t[,"S5"]), paste0(DIR_RESULTS, acronym_, "/", example, "/S5.pdf"))
    }
  }
}

get_cancer_types <- function()
{
  valid_tumors <- gsub("(.*).phi.txt", "\\1", list.files(DIR_COUNTS))
  valid_tumors <- tumortypes[tumortypes[,1] %in% valid_tumors, ]
  counts <- aggregate(valid_tumors, by=list(valid_tumors[,2]), FUN = length)
  counts <- counts[,-2]
  too_small = 25
  other <- sum(counts[counts[,2] < too_small,2])
  counts <- counts[counts[,2] >= too_small,]
  counts <- rbind(counts, c("Other", other))
  s_counts <- order(as.numeric(counts[,2]))
  st_counts <- c()
  left = 1
  middle = s_middle = length(s_counts) %/% 2 + 1
  while (left < s_middle) {
    st_counts <- c(st_counts, s_counts[left], s_counts[middle])
    left = left + 1
    middle = middle + 1
  }
  if (middle == length(s_counts)) {
    st_counts <- c(st_counts,  s_counts[length(s_counts)])
  }
  st_counts <- counts[st_counts,]
  labels <- paste0(st_counts[,1], " (", st_counts[,2], ")")
  pdf("cancer_types.pdf", width=9, height=9)
  pie(as.numeric(st_counts[,2]), labels = labels, radius = 0.7, init.angle = 315)
  dev.off()
}



get_exposure_max_changes <- function()
{
  library(vioplot)
  
  rescaled = F
  
  max_changes <- directions <- mean_exposures <- data.frame(matrix(NA, ncol=ncol(alex), nrow=0))
  colnames(max_changes) <- colnames(directions) <- colnames(mean_exposures) <- colnames(alex)
 
  for (dir_ in list.files(DIR_RESULTS))
  {
    if (!file.exists(paste0(DIR_RESULTS, dir_, "/"))) {
      next
    }
    print(paste0("Cancer type ", dir_))
    for (tumor_id in list.files(paste0(DIR_RESULTS, dir_)))
    {
      tumor_dir <- paste0(DIR_RESULTS, dir_, "/", tumor_id, "/")
      
      if (!file.exists(tumor_dir) || !file.exists(paste0(tumor_dir, "mixtures.mean.csv"))) {
        next
      }
      
      mixtures <- read.csv(paste0(tumor_dir, "mixtures.mean.csv"))
      rownames(mixtures) <- mixtures[,1]
      mixtures <- mixtures[,-1]
      
      new_item <- toHorizontalMatrix(rep(0, ncol(alex)))
      colnames(new_item) <- colnames(alex)
        
      rownames(new_item) <- tumor_id
      
      new_item_mean <- new_item_direction <- new_item
      new_item[,rownames(mixtures)] <- apply(mixtures, 1, max) - apply(mixtures, 1, min)
      
      new_item_direction[,rownames(mixtures)][apply(mixtures, 1, which.max) > apply(mixtures, 1, which.min)] <- 1
      new_item_direction[,rownames(mixtures)][apply(mixtures, 1, which.max) < apply(mixtures, 1, which.min)] <- -1
      
      new_item_mean[,rownames(mixtures)] <- apply(mixtures, 1, mean)
      
      max_changes <- rbind(max_changes, new_item)
      directions <- rbind(directions, new_item_direction)
      mean_exposures <- rbind(mean_exposures, new_item_mean)
    }
  }
  
  aggregate_by = tumortypes
  rownames(aggregate_by) <- aggregate_by[,1]
  aggregate_by <- aggregate_by[rownames(max_changes),-1]
  
  max_changes_per_tumor_type <- aggregate(max_changes, by = list(aggregate_by), FUN=mean)
  max_changes_per_tumor_type.err <- aggregate(max_changes, by = list(aggregate_by), FUN=function(x) sd(x) / sqrt(length(x)))
  directions_per_tumor_type <- aggregate(directions, by = list(aggregate_by), FUN=mean)
  mean_exposures_per_tumor_type <- aggregate(mean_exposures, by = list(aggregate_by), FUN=mean)
  
  filter_n_samples <- which(table(aggregate_by) > 10)
  max_changes_per_tumor_type <- max_changes_per_tumor_type[filter_n_samples,]
  max_changes_per_tumor_type.err <- max_changes_per_tumor_type.err[filter_n_samples,]
  directions_per_tumor_type <- directions_per_tumor_type[filter_n_samples,]
  mean_exposures_per_tumor_type <- mean_exposures_per_tumor_type[filter_n_samples,]
  mean_exposures_per_tumor_type[mean_exposures_per_tumor_type < 0.05] <- 0
  mean_exposures_per_tumor_type <- cbind(mean_exposures_per_tumor_type, 
                                         Other = 1-apply(mean_exposures_per_tumor_type[,-1], 1,sum))
  rownames(mean_exposures_per_tumor_type) <- mean_exposures_per_tumor_type[,1]
  mean_exposures_per_tumor_type <- mean_exposures_per_tumor_type[,-1]
  
  colfunc <- colorRampPalette(c("blue", "white", "red"))
  colors <- colfunc(20)
  
  pdf(paste0(DIR_RESULTS,"bar_plots_max_change.pdf"), width = 0.32 * (ncol(max_changes_per_tumor_type) - 1), height = 1.4 * nrow(max_changes_per_tumor_type))
  par(mfrow=c(nrow(max_changes_per_tumor_type),1),mar=c(2, 5, 2, 1))
  
  for (i in 1:nrow(max_changes_per_tumor_type))
  {
    barCenters <- barplot(unlist(c(max_changes_per_tumor_type[i,-1])), ylab = max_changes_per_tumor_type[i,1],
                          ylim = c(0, 0.25), col = colors[round((unlist(directions_per_tumor_type[i,-1]) + 1) * 10)])
    
    segments(barCenters, 
             toVerticalMatrix(max_changes_per_tumor_type[i,-1] - max_changes_per_tumor_type.err[i,-1]), 
             barCenters,
             toVerticalMatrix(max_changes_per_tumor_type[i,-1] + max_changes_per_tumor_type.err[i,-1]), 
             lwd = 1.5)
  }
  dev.off()

  rownames(tumortypes) <- tumortypes[,1]
  tumortypes_selected <- tumortypes[rownames(max_changes), ]
  unique_tumortypes <- names(table(tumortypes_selected[,2]))
  unique_tumortypes <- unique_tumortypes[table(tumortypes_selected[,2]) > 10]
  max_changes_directional <- max_changes * directions
  
  pdf(paste0(DIR_RESULTS,"violin_plots_max_change.pdf"), width = 0.4 * (ncol(max_changes_per_tumor_type) - 1), height = 2.0 * nrow(max_changes_per_tumor_type))
  par(mfrow=c(length(unique_tumortypes),1),mar=c(2, 5, 1.3, 1), bty='n')
  
  for (i in 1:length(unique_tumortypes))
  {
    plot(0:1,0:1,type="n",xlim=c(1,ncol(max_changes)), ylim = c(-0.5, 0.5), ann=F, xaxt="n")
    axis(side=1, labels=FALSE)
    for (j in 1:ncol(max_changes))
    {
      filter <- tumortypes_selected[,2] == unique_tumortypes[i]
      n_samples <- sum(filter)
      distr <- max_changes_directional[filter,j]
      c <- mean(directions[filter,j])
      if (sum(distr) != 0)
      {
        vioplot(distr, col = colors[round((c + 1) * 10)],add=T, at=j)
      }
      axis(side=1,at=j,labels=colnames(max_changes)[j])
    }
    title(ylab=paste0(unique_tumortypes[i], " (", n_samples, ")"))
    abline(h=0, col="gray")
  }
  dev.off()
}


plot_statistics <- function()
{
  stats <- read.csv(paste0(DIR_RESULTS, "Stats.csv"), header=T)
  
  cor(stats$n_timesteps, stats$n_changepoints)
  pdf(paste0(DIR_RESULTS,"n_timesteps_vs_n_changepoints.pdf"), width = 5, height=5)
  plot(stats$n_timesteps, stats$n_changepoints,  col=rgb(0,0,1,.2), pch=19,
       xlab = "# timesteps", ylab = "# change points")
  dev.off()
  
  cor(stats[stats$type == "BRCA",]$n_timesteps, stats[stats$type == "BRCA",]$n_changepoints)
  plot(stats[stats$type == "BRCA",]$n_timesteps, stats[stats$type == "BRCA",]$n_changepoints,  col=rgb(0,0,1,.2), pch=19)
  plot(1:11, table(stats$n_changepoints),  col=rgb(0,0,1), pch=19)

  pdf(paste0(DIR_RESULTS,"max_change.pdf"), width = 5, height=5)
  h <- hist(stats$max_change_no_s1_s5, breaks=seq(0,1,0.02), plot=F)
  plot(h$breaks[-1], cumsum(h$density) / sum(h$density), main="", 
       ylab = "Cumulative ratio of exposures", xlab = "Max signature change", col="gray", type="l",
       lwd =3)
  dev.off()
  hist(stats$max_change_no_s1_s5, breaks=seq(0,1,0.02))
  
  sum(stats$distance_at_transition_with_max_dist[stats$distance_at_transition_with_max_dist >= -1000] < 0 ) / sum(stats$distance_at_transition_with_max_dist >= -1000)
  d <- stats$distance_at_transition_with_max_dist[stats$distance_at_transition_with_max_dist >= -1000]
  sum(d <= 3 & d >= -3) / sum(stats$distance_at_transition_with_max_dist >= -1000)
  aggregate(stats$distance_at_transition_with_max_dist, by=list(stats$type), FUN=function(x) sum(x[x >= -1000] < 0 ) / sum(x >= -1000))
  aggregate(stats$distance_at_transition_with_max_dist, by=list(stats$type), FUN=function(x) sum(x[x >= -1000] == 0 ) / sum(x >= -1000))
  plot(stats$distance_at_transition_with_max_dist, stats$max_change)
  library(hexbin)
  library(RColorBrewer)
  rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  pdf(paste0(DIR_RESULTS,"Max_dist_change_transition_point.pdf"), width = 7, height=5)
  hexbinplot(max_change~distance_at_transition_with_max_dist, data=stats[stats$distance_at_transition_with_max_dist >= -1000,], colramp=rf,
             ylab="Max change", xlab="Maximum distance between predicted and real change point")
  dev.off()
  
  pdf(paste0(DIR_RESULTS,"transition_point_distance.pdf"), width = 5, height=5)
  hist(abs(unlist(dist_to_transition)), 
            xlab="Distance between predicted and real change point", xlim=c(-30,30), breaks=seq(0, 45,1), 
            main="", col="gray", freq=FALSE, plot=F)
  dev.off()
  
  pdf(paste0(DIR_RESULTS,"transition_point_distance_cum.pdf"), width = 5, height=5)
  h <- hist(abs(unlist(dist_to_transition)), 
       xlab="Distance between predicted and real change point", xlim=c(-30,30), breaks=seq(0, 45,1), 
       main="", col="gray", freq=FALSE, plot=F)
  h.random <-   hist(abs(unlist(dist_to_transition.random)), 
                     xlab="Distance between predicted and real change point", xlim=c(-30,30), breaks=seq(0, 45,1), 
                     main="", col="gray", freq=FALSE, plot=F)
  plot(h$breaks[-1], cumsum(h$density) / sum(h$density), main="", 
       ylab = "Cumulative ratio", xlab = "Distance between predicted and real change point", col="blue", type="l",
       lwd =3, xlim=c(1,30), ylim=c(0,1))
  lines(cumsum(h.random$density) / sum(h.random$density), col="red",lwd =3)
  legend("bottomright", c("Predicted", "Random"),
         lty=c(1,1), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "red"))
  dev.off()
  mean(abs(unlist(dist_to_transition)))
  
  pdf(paste0(DIR_RESULTS,"random_transition_point_distance.pdf"), width = 5, height=5)
  hist(unlist(dist_to_transition.random), 
       xlab="Distance between predicted and real change point", xlim=c(-30,30), breaks=seq(-45, 40,1), 
       main="", col="gray", freq=FALSE)
  dev.off()
  
  mean(abs(unlist(dist_to_transition.random)))
  
  sum(stats$n_changepoints > 0 & stats$n_transition_points == 0) / nrow(stats)
  sum(stats$n_changepoints == 0 & stats$n_transition_points > 0) / nrow(stats)
  plot(stats$n_changepoints, stats$n_transition_points,  col=rgb(0,0,1,.05), pch=19)
  
  
  # Plot stacked barplot of mixture composition
  # Example for presentation: "2da2b745-068d-408c-9163-3b3a78d4aaed"
  mixture_histogram <- t(alex.t) * apply(mixtures,1,mean)
  leave_sigs <- c("S1", "S5", "S8")
  other = apply(mixture_histogram[!(rownames(mixture_histogram) %in% leave_sigs),],2,sum)
  mixture_histogram_cleaned <- rbind(mixture_histogram[rownames(mixture_histogram) %in% leave_sigs,], Other = other)
  
  COLORS3 <- c("gold1", "firebrick2", "purple2", "palegreen1", "royalblue1", "maroon2")
    
  pdf(paste0(DIR_RESULTS,"example_mixture.pdf"), width=8, height=5)
  barplot(mixture_histogram_cleaned, cex.names=0.4, 
          col=COLORS3,
          #ylab = "Mutation type probability",
          ylim = ylim, xlim = xlim, las=2, yaxt="n", space=0.5)
  
  legend("topright", 
         legend = rownames(mixture_histogram_cleaned), 
         fill = COLORS3,
         cex = 0.6)
  dev.off()
  
  apply(mixtures,1,mean)[leave_sigs]
  
  make_histogram(toVerticalMatrix(apply(mixture_histogram_cleaned,2,sum) * 100), paste0(DIR_RESULTS, "example_distribution.pdf"))
  
  filter <- stats$type != "LIRI" & stats$type != "LAML" & (stats$max_phi - stats$min_phi > 0.05)
  more_than_1_subclone <- stats[filter,]$n_subclones > 1
  boxplot(stats[filter,]$max_change_no_s1_s5 ~ more_than_1_subclone, ylim=c(0,0.8), main="max_change_no_s1_s5")
  
  pdf(paste0(DIR_RESULTS,"max_change_in_subclones.pdf"), width=5, height=7)
  boxplot(stats[filter,]$max_change_no_s1_s5 ~ more_than_1_subclone, ylim=c(0,0.8),
          names=c("No subclones", "One subclone or more"))
  dev.off()
  
  pdf(paste0(DIR_RESULTS, "correlation_s1_s5.pdf"))
  hist(stats[stats$n_timesteps > 2,"cor_s1_s5"], xlim=c(-1,1), main="", xlab="Correlation between S1 and S5")
  dev.off()
  
  pdf(paste0(DIR_RESULTS, "max_changes_violin.pdf"))
  vioplot(stats$max_change[stats$max_change < 1],
          stats$max_change_no_s1_s5,
          stats$max_change_unsorted[!is.na(stats$max_change_unsorted)],
          names=c("All", "S1 and S5 change removed", "Permutated"),
          col="gold")
  dev.off()
   
  boxplot(stats$max_change_unsorted)
}


get_thresholds_on_sig_changes <-function()
{
  sig_change_threshold <- c()
  
  for (dir_ in list.files(DIR_RESULTS))
  {
    if (!file.exists(paste0(DIR_RESULTS, dir_, "/"))) {
      next
    }
    print(paste0("Cancer type ", dir_))
    for (tumor_id in list.files(paste0(DIR_RESULTS, dir_)))
    {
      tumor_dir <- paste0(DIR_RESULTS, dir_, "/", tumor_id, "/")
      
      if (!file.exists(tumor_dir) || !file.exists(paste0(tumor_dir, "mixtures.mean.csv"))) {
        next
      }
      
      if (!file.exists(paste0(tumor_dir, "mixtures.mean.unsorted.csv")) || !file.exists(paste0(tumor_dir, "mixtures.sd.unsorted.csv"))) {
        print(paste0("Mixtures on unsorted mutations missing on tumor ", tumor_id))
        next
      }
      mixtures <- read.csv(paste0(tumor_dir, "mixtures.mean.unsorted.csv"))
      
      mixtures.sd <- read.csv(paste0(tumor_dir, "mixtures.sd.unsorted.csv"))
      rownames(mixtures.sd) <- mixtures[,1]
      mixtures.sd <- mixtures.sd[,-1]
      
      threshold <- apply(mixtures.sd, 1, max)
      
      new_item <- toHorizontalMatrix(rep(0, ncol(alex)))
      colnames(new_item) <- colnames(alex)
      rownames(new_item) <- tumor_id
      
      new_item[,names(threshold)] <- threshold
      
      sig_change_threshold <- rbind(sig_change_threshold, new_item)
    }
  }
  
  write.csv(sig_change_threshold, file=paste0(DIR_RESULTS, "sig_change_threshold.csv"), row.names = F)
}
