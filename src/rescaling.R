# LUSC_files <- tumortypes[tumortypes[,2] == "LUSC",1]
# for (example in LUSC_files[1:length(LUSC_files)])
rescale_mixtures <- function(mixtures, keep_constant, beware_of_small_sigs=T)
{
  keep_constant <- toHorizontalMatrix(keep_constant)
  mixtures.rescaled <- mixtures
  mixtures_for_constant_sigs <- apply(keep_constant,1, mean)
  
  #  How important it is to rescale to this signtures
  if (beware_of_small_sigs) {
    weight = keep_constant / max(apply(mixtures,1, mean))
  } else {
    weight = rep(1, ncol(mixtures))
  }
  
  for (t in 1:ncol(mixtures))
  {
    mixtures.t <- keep_constant[,t]
    if (nrow(keep_constant) == 1) {
      rescale_coef <- (mixtures_for_constant_sigs / mixtures.t - 1) * weight[t] + 1
    } else {
      rescale_coef <- coef(lm(mixtures_for_constant_sigs ~ as.matrix(mixtures.t), weights = weight[,t]))[2]
    }
    
    mixtures.rescaled[,t] <- mixtures[,t] * rescale_coef
  }
  
  return(mixtures.rescaled)
}


determine_landscape_changing_sig <- function()
{
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(DIR_COUNTS)) 

  for (example in tumors)
  {
    print(paste0("Example ", example, " (", which(tumors == example), " out of ", length(tumors), ")"))
  
    pct <- proc.time()
    list[vcfFile, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(example, dir_counts, tumortypes)
    
    if (is.null(vcfData))
    {
      print(paste0("No data read for example ", example))
      next
    }
    if (nrow(vcfData) == 0)
      next
    if (nrow(vcfData) < 6) # Plots with less than 6 lines of data are meaningless so ignored
      next
    
    if (repeat_computation > 1) {
      repeat_dir = paste0("repeat", repeat_index, "/")
      dir.create(paste0(dir_name, repeat_dir))
      
      vcfData <- add_noise(vcfData, noise_rate = 0.05)
    }
    
    data_method <- "sliding400"
    window_size=400
    shift = window_size/100
    vcf <- get_sliding_window_data(vcfData, shift=shift)
    phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
    phis_sliding_window <- phis_sliding_window / shift
    
    phis_for_plot = phis_sliding_window
    if (sum(phis_sliding_window < 0.001) > 0.2 * length(phis_sliding_window))
    {
      phis_for_plot = NULL
    }
    
    if (!is.null(phis_for_plot))
    {
      colnames(vcf) <- round(phis_sliding_window, 3)
    } else {
      colnames(vcf) <- NULL
    }
    
    if (!is.null(assigns_phylo_nodes))
    {
      assigns_phylo_nodes_sw <- get_sliding_window_data(toVerticalMatrix(assigns_phylo_nodes), shift=shift)
      assigns_phylo_nodes_sw <- assigns_phylo_nodes_sw / shift
      assigns_phylo_nodes_sw <- round(assigns_phylo_nodes_sw)
    } else {
      assigns_phylo_nodes_sw = NULL
    }
    
    prior = NULL
    
    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(vcfFile, active_signatures.our_samples, alex)
    } else {
      alex.t <- alex
    }
    
    dir_name_ <- dir_name
    
    save(vcfData, vcf, phis, phis_sliding_window, assigns_phylo_nodes, assigns_phylo_nodes_sw, 
         file = paste0(dir_name_, "/", example, ".RData"))
    
    method_name <- "iterativeChangePoints"
    
    if (!file.exists(paste0(dir_name_, "mixtures.csv")) || !file.exists(paste0(dir_name_, "changepoints.txt")))
    {
      stop(paste0("Mixtures or changepoints are absent for sample ", example))
    } else
    {
      mixtures <- read_mixtures(paste0(dir_name_, "mixtures.csv"))
      changepoints <- unlist(read.table(paste0(dir_name_, "changepoints.txt"), header=F))
      if (ncol(mixtures) != length(phis_sliding_window))
      {
        list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t)
        
        write.csv(mixtures, file=paste0(dir_name_, "mixtures.csv"))
        write(changepoints, file=paste0(dir_name_, "changepoints.txt"), ncolumns=length(changepoints))
      }
      stopifnot(ncol(mixtures) == length(phis_sliding_window))
    }
    
    sig_changes <- apply(mixtures, 1, max) - apply(mixtures, 1, min)
    highly_changing_sigs <- rownames(mixtures)[sig_changes > 0.2]
    
    compute_proportions_to = "S1"
    
    if (length(highly_changing_sigs) > 0)
    {
      for (h_sig in highly_changing_sigs)
      {
        print(h_sig)
        mixtures_without_h_sig <- mixtures[rownames(mixtures) != h_sig,]
        
        mixtures.rescaled_to_h_sig <- rescale_mixtures(mixtures_without_h_sig, apply(mixtures_without_h_sig, 2, sum), beware_of_small_sigs=F)
        mixtures.rescaled_to_h_sig.no_h_sig <-  mixtures.rescaled_to_h_sig[rownames(mixtures.rescaled_to_h_sig) != h_sig,]
        
        proportions_before <- t(t(mixtures_without_h_sig) / unlist(mixtures_without_h_sig[compute_proportions_to,,drop=T]))
        proportions_after <- t(t(mixtures.rescaled_to_h_sig.no_h_sig) / unlist(mixtures.rescaled_to_h_sig.no_h_sig[compute_proportions_to,,drop=T]))
          
        plot(plot_signatures(mixtures_without_h_sig, plot_name=plot_name, phis = phis_sliding_window, mark_change_points=T, change_points=changepoints,
                             assigns_phylo_nodes = assigns_phylo_nodes_sw, save=F)$plot)
        
        plot(plot_signatures(mixtures.rescaled_to_h_sig.no_h_sig, plot_name=plot_name, phis = phis_sliding_window, mark_change_points=T, change_points=changepoints,
                             assigns_phylo_nodes = assigns_phylo_nodes_sw, save=F)$plot)
        
        aheatmap(proportions_before, Colv=NA, Rowv=NA)
        aheatmap(proportions_after, Colv=NA, Rowv=NA)
        
        readline()
      }
    }
    
    plot_signatures(mixtures, plot_name=plot_name, phis = phis_sliding_window, mark_change_points=T, change_points=changepoints,
                    assigns_phylo_nodes = assigns_phylo_nodes_sw, save=F)
    print(paste("Computed example", example))
  }
}

rescale_mixtures_to_make_sum_1 <- function(mixtures_)
{
  mixtures.rescaled <- sapply(1:ncol(mixtures_), function (i) mixtures_[,i] / apply(mixtures_, 2, sum)[i])
  rownames(mixtures.rescaled) <- rownames(mixtures_)
  return(mixtures.rescaled)
}