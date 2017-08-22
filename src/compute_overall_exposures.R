# AUTHOR: Yulia Rubanova

compute_overall_exposures_for_all_examples <- function(dir_counts = DIR_COUNTS, active_signatures = active_signatures.our_samples, 
    save=F, tumors=NULL, verbose=T)
{
  if (is.null(tumors)) {
    tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 
  }

  all_exposures <- data.frame(matrix(0, nrow=0, ncol=ncol(alex)))
  colnames(all_exposures) <- colnames(alex)
  
  for (example in tumors)
  {
    if (grepl("simulation.*", example))
    {
      next
    }
    set.seed(which(tumors == example))
    if (verbose) {
      print(paste0("Example ", example, " (", which(tumors == example), " out of ", length(tumors), ")"))
    }
    pct <- proc.time()
    
    data_file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData")
    if (file.exists(data_file))
    {
      load(data_file)
    } else {
      if (verbose) {
        print(paste0("Data file ", data_file, " does not exist"))
      }
      next
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    
    #suppressWarnings(dir.create(dir_name, recursive = T))  
    if (!file.exists(dir_name)) {
      dir.create(dir_name, recursive = T)
    }
    
    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(tumor_id, active_signatures, alex, noise_sig)
    } else {
      alex.t <- alex
    }
    
    if (is.null(alex.t))
    {
      if (verbose) {
        print(paste0("No active signatures for sample", example, " ...."))
      }
      next
    }

    impossible_types <- apply(alex.t,1,sum) == 0
    mutation_counts <- apply(t(vcfData),1,sum)
    if (sum(mutation_counts[impossible_types]) > 0) {
      if (verbose) {
       print("Probability of some mutations is zero under the active set of signatures")
      }
      next
    }
    
    mixtures <- fit_mixture_of_multinomials_EM(mutation_counts, alex.t)

    mixtures <- t(data.frame(mixtures))
    colnames(mixtures) <- colnames(alex.t)
    rownames(mixtures) <- example
    
    colnames(mixtures) <- colnames(alex.t)
    
    new_item <- data.frame(matrix(0, nrow=1, ncol=ncol(alex)))
    colnames(new_item) <- colnames(alex)
    rownames(new_item) <- example
    new_item[,colnames(mixtures)] <- mixtures
    
    all_exposures <- rbind(all_exposures, new_item)
    
    if (verbose) {
      print(proc.time() - pct)
      print(paste("Computed example", example))
    }
  }
  
  if (save) {
    write.csv(all_exposures, file=paste0(DIR_RESULTS, "overall_exposures.csv"))
  }
  
  if ("SNPs" %in% colnames(all_exposures)) {
    all_exposures <- all_exposures[,-which(colnames(all_exposures) == "SNPs")]
  }
  
  return(all_exposures)
}
