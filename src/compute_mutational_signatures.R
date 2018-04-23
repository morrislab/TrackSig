# AUTHOR: Yulia Rubanova

source("src/header.R")

group = 0
EXAMPLES_PER_GROUP <- 500

save_data_for_samples <- function(dir_counts = DIR_COUNTS,  bootstrap_counts = BOOTSTRAP_COUNTS)
{
  if (!file.exists(SAVED_SAMPLES_DIR)) {
    dir.create(SAVED_SAMPLES_DIR, recursive = T)
  }

  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 

  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)

  for (example in examples_group)
  {
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))    
    
      list[tumor_id, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(example, dir_counts, tumortypes, dir_create = F)

    if (is.null(vcfData))
    {
      print(paste0("No data read for example ", example))
      next
    }
    if (nrow(vcfData) == 0)
    {
      print(paste0("Zero rows for example " , example))
      next
    }
    if (nrow(vcfData) < 6) # Plots with less than 6 lines of data are meaningless so ignored
    {
      print(paste0("Less than 6 rows per example ", example))
      next
    }
    
    if (sliding_window) { 
      # Sliding window approach
      data_method <- "sliding400"
      window_size=400
      shift = window_size/100
      gap = 1
      vcf <- get_sliding_window_data(vcfData, shift=shift, gap = gap)
      phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
      phis_sliding_window <- phis_sliding_window / shift
      
      phis_for_plot = phis_sliding_window
    } else {
      shift <- gap <- NULL
      data_method <- "chunk100"
      vcf <- t(vcfData)
      phis_for_plot <- phis_sliding_window <- phis
    }
      
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
        assigns_phylo_nodes <-  as.factor(assigns_phylo_nodes)
        assigns_phylo_nodes_levels <- levels(assigns_phylo_nodes)
        assigns_phylo_nodes <- toVerticalMatrix(as.numeric(assigns_phylo_nodes))
        
        if (sliding_window) {
          assigns_phylo_nodes_sw <- get_sliding_window_data(assigns_phylo_nodes, shift=shift)
          assigns_phylo_nodes_sw <- assigns_phylo_nodes_sw / shift
          assigns_phylo_nodes_sw <- round(assigns_phylo_nodes_sw)
        } else {
          assigns_phylo_nodes_sw <- assigns_phylo_nodes
        } 
        
        for (l in 1:length(assigns_phylo_nodes_levels))
        {
          assigns_phylo_nodes_sw[assigns_phylo_nodes_sw == l] <- assigns_phylo_nodes_levels[l]
        }
        
        stopifnot(length(phis_for_plot) == length(assigns_phylo_nodes_sw))
      } else {
        assigns_phylo_nodes_sw = assigns_phylo_nodes
      }
      
      list[bootstrap_vcfs, bootstrap_phis] <- lapply(extract_bootstrap_data_for_example(example, bootstrap_counts), t)

      if (compute_bootstrap) {
        for (j in 1:length(bootstrap_phis))
        {
          if (sliding_window) {
            bootstrap_vcfs[[j]] <- get_sliding_window_data(bootstrap_vcfs[[j]], shift=shift, gap = gap)
            bootstrap_phis[[j]] <- get_sliding_window_data(toVerticalMatrix(bootstrap_phis[[j]]), shift=shift)
            bootstrap_phis[[j]] <- bootstrap_phis[[j]] / shift
          }
        }
      }

    save(vcfData, vcf, phis, phis_sliding_window, assigns_phylo_nodes, assigns_phylo_nodes_sw, 
         acronym, window, shift, gap, tumor_id, phis_for_plot, bootstrap_vcfs, bootstrap_phis,
         file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData"))
  }
}


compute_signatures_for_all_examples <- function(dir_counts = DIR_COUNTS)
{
  add_early_late_transition = TRUE

  age_signatures <- c("S1", "S5", "L1", "1", "5a", "5b")

    tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 

  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)
  
  mutation_types <- read.delim(trinucleotide_file, header=F)
  mutation_types <- paste(mutation_types[,1], mutation_types[,2], mutation_types[,3], sep="_")

  for (example in examples_group)
  {
    set.seed(which(examples_group == example))
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))
    
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

    if (is.null(acronym) || acronym == "") {
      print(paste("ERROR: Cancer type not found for ", example))
      next
    }
    
    if (is.null(alex.t))
    {
      print(paste0("No active signatures for sample", example, " ...."))
      next
    }

    if (is.vector(alex.t)) {
      next
    }

    if (sum(vcf[apply(alex.t,1,sum) == 0,] != 0) != 0) {
      print(paste0("Sample ", example, ": some trinucleotides have probability 0 under the model, but their count is non-zero. Sot the count vector is impossible under the model."))
      next
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    
    suppressWarnings(dir.create(dir_name, recursive = T))  
    
    if (!is.null(phis_for_plot))
    { 
      write(phis_for_plot, file=paste0(dir_name, "phis.txt"), ncolumns=length(phis_for_plot))
    }

    method_name <- "iterativeChangePoints"
    
    if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
    {
      if (changepoint_method == "PELT") {
        list[changepoints, mixtures] <- find_changepoints_pelt(vcf, alex.t)
      } else {
        list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t, n_signatures = ncol(alex.t))
      }

      write.csv(mixtures, file=paste0(dir_name, "mixtures.csv"))

      n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
      write(changepoints, file=paste0(dir_name, "changepoints.txt"), ncolumns=n_col)
    } else {
      mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
      cp_file = paste0(dir_name, "changepoints.txt")
      if (file.info(cp_file)$size == 1) {
        changepoints <- c()
      } else {
        changepoints <- unlist(read.table(cp_file, header=F))
      }
    }

    if (!is.null(assigns_phylo_nodes_sw)) {
      write(assigns_phylo_nodes_sw,  file=paste0(dir_name, "assignments.txt"), ncolumns=length(assigns_phylo_nodes_sw))
    } else  {
      n_clusters = transition_points = assigns_phylo_nodes_sw = NULL
    }

    age_signatures <- intersect(rownames(mixtures), age_signatures)

      plot_name <- paste0(dir_name, "/", acronym, "_", tumor_id, "_", sig_amount, postfix, ".pdf")

    if (PLOT_FULL_NAME)
    {
      plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix, ".pdf")
    }

    mark_cp <- !is.null(changepoints)
    plot_signatures(mixtures*100, plot_name=plot_name, phis = phis_for_plot, mark_change_points=mark_cp, change_points=changepoints, 
                    #assigns_phylo_nodes = assigns_phylo_nodes_sw, 
                    transition_points = transition_points,
                    scale=1.2)
    
    mixtures.rescaled = NULL
  
    print(paste("Computed example", example))
  }
 
}


compute_errorbars_for_all_examples <- function(bootstrap_counts = BOOTSTRAP_COUNTS)
{
  add_early_late_transition = TRUE

  tumors <- list.dirs(bootstrap_counts, recursive = F, full.names=F)
  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)

  for (example in examples_group)
  {
    set.seed(which(examples_group == example))
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))
    
    data_file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData")
    if (file.exists(data_file)) {
      load(data_file)
    } else {
      print(paste0("Data file ", data_file, " does not exist"))
      next
    }
  
    max_tp = 100
    # do not compute bootstrap for samples with number of time points > max_tp
    if (ncol(vcf) > max_tp) {
      print(paste0("Skipping ", example, ": more than ", max_tp, " timepoints"))
      next
    }

    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, noise_sig)
    } else if (sig_amount == "selectedSignatures") {
      selected = c("S1", "S3", "S5", "S8", "S9", "S16")
      alex.t <- alex[,selected]
    } else {
      alex.t <- alex
    }
    
    if (is.null(alex.t)) {
      print(paste0("No active signatures for sample", example, " ...."))
      next
    }
    
    if (acronym == "") {
      print(paste("Cancer type not found for ", tumor_id))
      next
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    suppressWarnings(dir.create(dir_name, recursive = T))  
    
    method_name <- "iterativeChangePoints"
    
    if (!file.exists(paste0(dir_name, "mixtures.csv"))) {
      next
    }   
    mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
    changepoints <- tryCatch({
      unlist(read.table(paste0(dir_name, "changepoints.txt"), header=F))
    }, error = function(e){return()})

    print("Computing bootstrapped trajectories")
    list[mixtures_bootstrap, changepoints_bootstrap] <- get_bootstrap_mixtures(bootstrap_vcfs, bootstrap_phis, alex.t, dir_name, "")
    list[mixtures.mean, mixtures.sd, mixtures.err] <- compute_mean_sd_err(mixtures_bootstrap, rownames(mixtures), dir_name)
    
    transition_points = NULL
    plot_name <- paste0(dir_name, "/", acronym, "_", example, "_", sig_amount)
    
    if (PLOT_FULL_NAME)
    {
      plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, postfix)
    }

    plot_signatures(mixtures.mean*100, plot_name=paste0(plot_name, ".mean.bootstrap_traj.pdf"), phis = phis_for_plot, 
                    mark_change_points=T,
                    change_points=changepoints,
                    transition_points = transition_points,
                    fitted_data = lapply(mixtures_bootstrap, function(x) x* 100))
                    #assigns_phylo_nodes = assigns_phylo_nodes_sw) #error_bars = mixtures.err)

    plot_signatures_real_scale(mixtures.mean*100, plot_name=paste0(plot_name, ".mean.bootstrap_traj.all.pdf"), 
                               phis = phis_for_plot, mark_change_points=T, change_points=changepoints_bootstrap, 
                               #assigns_phylo_nodes = assigns_phylo_nodes_sw, #error_bars = mixtures.err,
                               transition_points = transition_points,
                               fitted_data = lapply(mixtures_bootstrap, function(x) x* 100), remove_sigs_below = 0)

    print(paste("Computed example", example))
  }
}



save_data_for_samples()
suppressMessages(compute_signatures_for_all_examples())
if (compute_bootstrap) {compute_errorbars_for_all_examples()}
