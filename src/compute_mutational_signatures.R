if (file.exists("/dupa-filer/")) {
  DIR <- "/dupa-filer/yulia/mutational_signatures/"
  READ_DIR = WRITE_DIR = DIR
} else if (file.exists("/scratch/")) {
  DIR="." 
  READ_DIR = "/home/q/qmorris/yulia/"
  WRITE_DIR = "/scratch/q/qmorris/yulia/"
} else {
  DIR <- "/Users/yulia/Documents/mutational_signatures/"
  READ_DIR = WRITE_DIR = DIR
}

setwd(DIR)

library(glmnet)
library(reshape2)
library(ggplot2)
library(doParallel)
library(NMF)
library(SNFtool)
library(zoo)
nmf.options(grid.patch=TRUE)

sig_amount <- "onlyKnownSignatures"
# sig_amount <- "allSignatures"
# sig_amount <-  "selectedSignatures"

noise_sig = NULL
#noise_sig = "uniform"

compute_bootstrap = TRUE

sliding_window = TRUE
if (sliding_window)
{
  saved_dir = "sliding_window_400"
} else {
  saved_dir = "chunk100"
}

noise_sig_dir = ""
if (!is.null(noise_sig)) {
  noise_sig_dir = paste0(noise_sig, "_noise_sig/")
}

phi_source="vaf"
sig_type = "" # cosmic
# sig_type = "pcawg_"

# DIR_COUNTS = "counts_1682samples/psub/"
# DIR_RESULTS = "Examples_1682/psub/results_onlyKnownSignatures_psub/"
# SAVED_SAMPLES_DIR = "saved_data/Examples_1682/"
# mutation_assignments = "mutass_sorted_by_phi_1682/psub/"
# ANNOTATION="saved_data/annotation_data.pcawg_sigs.RData"
# BOOTSTRAP_COUNTS = "bootstrap_counts_1682samples/psub/"
# 
# ANNOTATION="saved_data/annotation_data.RData"
#ANNOTATION="saved_data/annotation_data.lydia.lydias_signatures.RData"

# DIR_COUNTS = paste0(READ_DIR, "lydia_prostate_cancer_200/counts/" ,"/")
# DIR_RESULTS = paste0(WRITE_DIR, "lydia_prostate_cancer_200/results_", sig_amount, "/")
# ANNOTATION=paste0(READ_DIR, "saved_data/annotation_data.lydia.RData")
# SAVED_SAMPLES_DIR = paste0(READ_DIR, "saved_data/Examples_lydia")
# mutation_assignments = NULL

LYDIA <- F

if (LYDIA)
{
  DIR_COUNTS = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/counts/"
  if (sliding_window) {
    DIR_RESULTS = paste0("/scratch/q/qmorris/yulia/prostate_cancer_200_fullresults/results_", sig_amount, "/", noise_sig_dir)
  } else { 
    DIR_RESULTS = paste0("/scratch/q/qmorris/yulia/prostate_cancer_200_fullresults/results_", sig_amount, "_chunk100/", noise_sig_dir)
  }
  SAVED_SAMPLES_DIR = paste0("/scratch/q/qmorris/yulia/prostate_cancer_200_fullresults/saved_data/", saved_dir, "/")
  mutation_assignments = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/prostate_cancer_200_ssms/"
  ANNOTATION="/home/q/qmorris/yulia/saved_data/annotation_data.lydia.RData"
  mutation_order = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/mut_order/"
  BOOTSTRAP_COUNTS = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/bootstrap_counts/"
} else {
  DIR_COUNTS = paste0("/home/q/qmorris/yulia/pwgs/counts/", phi_source, "/")
  # DIR_RESULTS = paste0("/scratch/q/qmorris/yulia/pwgs/samples/psub/results_", sig_amount, "/")
  if (sliding_window) {  
    DIR_RESULTS = paste0("/scratch/q/qmorris/yulia/pwgs/samples/", phi_source, "/results_", sig_amount, "_", sig_type, phi_source, "/", noise_sig_dir)
  } else {
    DIR_RESULTS = paste0("/scratch/q/qmorris/yulia/pwgs/samples/", phi_source, "/results_", sig_amount, "_", sig_type, phi_source, "_chunk100/", noise_sig_dir)
  }
  SAVED_SAMPLES_DIR = paste0("/scratch/q/qmorris/yulia/saved_data/samples/", phi_source, "/", saved_dir, "/")
  mutation_assignments = paste0("/home/q/qmorris/yulia/pwgs/mutass_sorted_by_phi/", phi_source, "/")
  if (sig_type == "pcawg_") {
  	ANNOTATION="/home/q/qmorris/yulia/saved_data/annotation_data.pcawg_sigs.RData"
  } else {
	ANNOTATION="/home/q/qmorris/yulia/saved_data/annotation_data.RData"
  }
  BOOTSTRAP_COUNTS = paste0("pwgs/bootstrap_counts/", phi_source, "/")
  mutation_order = paste0("/home/q/qmorris/yulia/pwgs/mut_order/", phi_source, "/")
}

load(ANNOTATION)

PLOT_FULL_NAME = F

if (!file.exists(DIR_RESULTS))
{
  dir.create(DIR_RESULTS, recursive=T)
}

src_files <- setdiff(grep(".*R$", list.files(paste0(READ_DIR, "src"),full.names = T), value = T), 
                     c(paste0(READ_DIR, "src/compute_mutational_signatures.R"), paste0(READ_DIR, "src/compute_mutational_signatures_copy.R")))
print(src_files)
for (file in src_files)
{
  source(file)
}
registerDoParallel(cores = detectCores()-1)


suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--group", default=0, type="integer", help="Index of example group to be evaluated")
args <- parser$parse_args()
group <- args$group

save_data_for_samples <- function(dir_counts = DIR_COUNTS,  bootstrap_counts = BOOTSTRAP_COUNTS)
{
  if (!file.exists(SAVED_SAMPLES_DIR)) {
    dir.create(SAVED_SAMPLES_DIR, recursive = T)
  }
  
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 
  EXAMPLES_PER_GROUP <- 10
  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)
  for (example in examples_group)
  {
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))    
    pct <- proc.time()
    
    list[tumor_id, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(example, dir_counts, tumortypes)
    
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
    
    list[bootstrap_vcfs, bootstrap_phis, bootstrap_vcfs_unsorted] <- lapply(extract_bootstrap_data_for_example(example, bootstrap_counts), t)
    
    for (j in 1:length(bootstrap_phis))
    {
      if (sliding_window) {
        bootstrap_vcfs[[j]] <- get_sliding_window_data(bootstrap_vcfs[[j]], shift=shift, gap = gap)
        bootstrap_vcfs_unsorted[[j]] <- get_sliding_window_data(bootstrap_vcfs_unsorted[[j]], shift=shift, gap = gap)
        bootstrap_phis[[j]] <- get_sliding_window_data(toVerticalMatrix(bootstrap_phis[[j]]), shift=shift)
        bootstrap_phis[[j]] <- bootstrap_phis[[j]] / shift
      }
    }
    
    save(vcfData, vcf, phis, phis_sliding_window, assigns_phylo_nodes, assigns_phylo_nodes_sw, 
         acronym, window, shift, gap, tumor_id, phis_for_plot, bootstrap_vcfs, bootstrap_phis, bootstrap_vcfs_unsorted,
         file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData"))
  }
}

get_overall_distribution <- function(dir_counts = DIR_COUNTS)
{
  print("Cmputing overall distribution...")
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts))
  overall_distribution <- c()
  for (example in tumors)
  {
    pct <- proc.time()

    list[tumor_id, vcfData, phis, assigns_phylo_nodes, acronym, dir_name] <- extract_data_for_example(example, dir_counts, tumortypes)

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

    overall_distribution <- rbind(overall_distribution, apply(vcfData,2,sum) / sum(apply(vcfData,2,sum)))
  }

  overall_distribution <- apply(overall_distribution,2,sum) / sum(apply(overall_distribution,2,sum))
  save(overall_distribution, file = paste0(SAVED_SAMPLES_DIR, "/overall_distribution.RData"))
}



compute_signatures_for_all_examples <- function(dir_counts = DIR_COUNTS)
{
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(dir_counts)) 
  EXAMPLES_PER_GROUP <- 10
  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)
  
  for (example in examples_group)
  {
    set.seed(which(examples_group == example))
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))
    
    pct <- proc.time()
    
    data_file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData")
    if (file.exists(data_file))
    {
      load(data_file)
    } else {
      print(paste0("Data file ", data_file, " does not exist"))
      next
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    
    #suppressWarnings(dir.create(dir_name, recursive = T))  
    dir.create(dir_name, recursive = T)

    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(tumor_id, active_signatures.our_samples, alex, noise_sig)
    } else if (sig_amount == "selectedSignatures") {
      selected = c("S1", "S3", "S5", "S8", "S9", "S16")
      alex.t <- alex[,selected]
    } else {
      alex.t <- alex
    }
    
    
    if (is.null(alex.t))
    {
      print(paste0("No active signatures for sample", example, " ...."))
      next
    }
    
    if (!is.null(phis_for_plot))
    { 
      write(phis_for_plot, file=paste0(dir_name, "phis.txt"), ncolumns=length(phis_for_plot))
    }
    if (!is.null(assigns_phylo_nodes_sw))
    {
      write(assigns_phylo_nodes_sw,  file=paste0(dir_name, "assignments.txt"), ncolumns=length(assigns_phylo_nodes_sw))
    }
    method_name <- "iterativeChangePoints"
    
    if (!file.exists(paste0(dir_name, "mixtures.csv")) || !file.exists(paste0(dir_name, "changepoints.txt")))
    {
      list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t)
      
      write.csv(mixtures, file=paste0(dir_name, "mixtures.csv"))
      write(changepoints, file=paste0(dir_name, "changepoints.txt"), ncolumns=length(changepoints))
    } else
    {
      mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
      changepoints <- unlist(read.table(paste0(dir_name, "changepoints.txt"), header=F))
    }
    
    plot_name <- paste0(dir_name, "/", acronym, "_", example, "_", sig_amount, ".pdf")
    if (PLOT_FULL_NAME)
    {
      plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, ".pdf")
    }
    plot_signatures(mixtures, plot_name=plot_name, phis = phis_for_plot, mark_change_points=T, change_points=changepoints, 
                    assigns_phylo_nodes = assigns_phylo_nodes_sw)
    
    mixtures.rescaled = NULL
    age_signatures <- c("S1", "S5", "L1", "1")
    age_signatures <- intersect(rownames(mixtures), age_signatures)
    if (sum(age_signatures %in% rownames(mixtures)) != 0)
    {
      keep_constant <- IgnoreVectorOrMatrix(mixtures[age_signatures,], function(x) {apply(x,2, sum) })
      if (!is.null(keep_constant) & (!file.exists(paste0(dir_name, "mixtures.rescaled.csv"))))
      {
        mixtures.rescaled <- rescale_mixtures(mixtures, keep_constant)
        write.csv(mixtures.rescaled, file=paste0(dir_name, "mixtures.rescaled.csv"))
      } else {
        mixtures.rescaled <- read_mixtures(paste0(dir_name, "mixtures.rescaled.csv"))
      }
      
      plot_name <- paste0(dir_name, "/", acronym, "_", example, "_", sig_amount, "_rescaled.pdf")
      if (PLOT_FULL_NAME)
      {
        plot_name <- paste0(dir_name, acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, "_rescaled.pdf")
      }
      #       plot_signatures(mixtures.rescaled, plot_name=plot_name, phis = phis_for_plot, mark_change_points=F, change_points=changepoints,
      #                       assigns_phylo_nodes = assigns_phylo_nodes_sw)
      
    }
    gather_statistics(mixtures, changepoints, tumor_id, dir_name, acronym, mixtures.rescaled)
    
    print(proc.time() - pct)
    print(paste("Computed example", example))
  }
  
  gather_summaries_per_tissue(omit_signature_information = T)
}


compute_errorbars_for_all_examples <- function(bootstrap_counts = BOOTSTRAP_COUNTS)
{
  tumors <- gsub("([^/]*)\\.phi\\.txt","\\1", list.files(bootstrap_counts)) 
  
  EXAMPLES_PER_GROUP <- 10
  examples_group <- get_examples_group(tumors, EXAMPLES_PER_GROUP, group)
  
  for (example in examples_group)
  {
    set.seed(which(examples_group == example))
    print(paste0("Example ", example, " (", which(examples_group == example), " out of ", length(examples_group), ")"))

    pct <- proc.time()
    
    data_file = paste0(SAVED_SAMPLES_DIR, "/", example, ".RData")
    if (file.exists(data_file)) {
      load(data_file)
    } else {
      print(paste0("Data file ", data_file, " does not exist"))
      next
    }
    
    dir_name <- paste0(DIR_RESULTS, acronym, "/", tumor_id, "/")
    
    suppressWarnings(dir.create(dir_name, recursive = T))  
    
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
    
    method_name <- "iterativeChangePoints"
    
    if (!file.exists(paste0(dir_name, "mixtures.csv"))) {
      next
    }   
    mixtures <- read_mixtures(paste0(dir_name, "mixtures.csv"))
    changepoints <- unlist(read.table(paste0(dir_name, "changepoints.txt"), header=F))
    
    print("Computing bootstrapped trajectories")
    list[mixtures_bootstrap, changepoints_bootstrap] <- get_bootstrap_mixtures(bootstrap_vcfs, bootstrap_phis, alex.t, dir_name, "")
    list[mixtures.mean, mixtures.sd, mixtures.err] <- compute_mean_sd_err(mixtures_bootstrap, rownames(mixtures), dir_name)
    
 if (is.null(noise_sig))
{   
    print("Computing bootstrapped trajectories on unsorted mutations")
    list[mixtures_bootstrap_unsorted, changepoints_bootstrap_unsorted] <- get_bootstrap_mixtures(bootstrap_vcfs_unsorted, bootstrap_phis, alex.t, dir_name, ".unsorted")
    list[mixtures.mean.unsorted, mixtures.sd.unsorted, mixtures.err.unsorted] <- compute_mean_sd_err(mixtures_bootstrap_unsorted, rownames(mixtures), dir_name, descr = ".unsorted")
}    
    plot_name <- paste0(dir_name, "/", acronym, "_", example, "_", sig_amount)
    if (PLOT_FULL_NAME)
    {
      plot_name <- paste0(dir_name, "/", acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name)
    }
    plot_signatures(mixtures, plot_name=paste0(plot_name, ".pdf"), phis = phis_for_plot, mark_change_points=F, change_points=changepoints, 
                       assigns_phylo_nodes = assigns_phylo_nodes_sw, error_bars = mixtures.err)
    
    plot_signatures(mixtures.mean, plot_name=paste0(plot_name, ".mean.pdf"), phis = phis_for_plot, mark_change_points=F, change_points=changepoints,     
                     assigns_phylo_nodes = assigns_phylo_nodes_sw, error_bars = mixtures.err)
    
    plot_signatures(mixtures.mean, plot_name=paste0(plot_name, ".mean.bootstrap_traj.pdf"), phis = phis_for_plot, mark_change_points=T, change_points=changepoints,     
                    assigns_phylo_nodes = assigns_phylo_nodes_sw) #error_bars = mixtures.err)
    
    p <- plot_signatures(mixtures.mean, plot_name="", phis = phis_for_plot, mark_change_points=T, change_points=changepoints, save=F)$plot   
                     #assigns_phylo_nodes = assigns_phylo_nodes_sw) #error_bars = mixtures.err)

    g <- add_change_point_distribution(p, changepoints_bootstrap, ncol(mixtures), max(mixtures) + 0.05, paste0(plot_name, ".mean.cp_distr.pdf"))

#     cp_distr_table <- data.frame(matrix(0, ncol=nrow(cp_distr), nrow=nrow(cp_distr)))
#     for (i in 1:nrow(cp_distr))
#     {
#       for (j in 1:nrow(cp_distr))
#       {
#         if (i >= j)
#         {
#           cp_distr_table[i,j] <- sum(cp_distr[i:j,1])
#         }
#       }
#     }
    plot_signatures_real_scale(mixtures.mean, plot_name=paste0(plot_name, ".mean.bootstrap_traj.all.pdf"), 
                               phis = phis_for_plot, mark_change_points=F, change_points=changepoints,     
                               #assigns_phylo_nodes = assigns_phylo_nodes_sw, #error_bars = mixtures.err,
                               fitted_data = mixtures_bootstrap, remove_sigs_below = 0.05)

if (is.null(noise_sig))
{    
    plot_signatures(mixtures.mean.unsorted, plot_name=paste0(plot_name, ".mean.unsorted.bootstrap_traj.pdf"), phis = phis_for_plot, mark_change_points=T, change_points=changepoints,     
                    assigns_phylo_nodes = assigns_phylo_nodes_sw) #error_bars = mixtures.err)
    
    plot_signatures_real_scale(mixtures.mean.unsorted, plot_name=paste0(plot_name, ".mean.unsorted.bootstrap_traj.all.pdf"), 
                               phis = phis_for_plot, mark_change_points=F, change_points=changepoints,     
                               #assigns_phylo_nodes = assigns_phylo_nodes_sw, #error_bars = mixtures.err,
                               fitted_data = mixtures_bootstrap_unsorted, remove_sigs_below = 0.05)
}    
    mixtures.rescaled <- read_mixtures(paste0(dir_name, "mixtures.rescaled.csv"))
    age_signatures <- c("S1", "S5", "L1", "1")
    age_signatures <- intersect(rownames(mixtures), age_signatures)
    if (sum(age_signatures %in% rownames(mixtures)) != 0)
    {
      keep_constant <- IgnoreVectorOrMatrix(mixtures[age_signatures,], function(x) {apply(x,2, sum) })
      
      mixtures_bootstrap.rescaled <- list()
      for (j in 1:length(mixtures_bootstrap))
      {
        if (!file.exists(paste0(dir_name,"mixtures.bootstrap_", j, ".rescaled.csv")))
        { 
          m <- rescale_mixtures(mixtures_bootstrap[[j]], keep_constant)
          mixtures_bootstrap.rescaled[[j]] <- m
          write.csv(m, file=paste0(dir_name, "mixtures.bootstrap_", j, ".rescaled.csv"))
        } else {
          m <- read_mixtures(paste0(dir_name,"mixtures.bootstrap_", j, ".rescaled.csv"))
        }
        mixtures_bootstrap.rescaled[[j]] <- m
      }
      
      list[mixtures.rescaled.mean, mixtures.rescaled.sd, mixtures.rescaled.err] <- compute_mean_sd_err(mixtures_bootstrap.rescaled, rownames(mixtures), dir_name)
      

      plot_name <- paste0(dir_name, "/", acronym, "_", example, "_", sig_amount, "_rescaled")
      if (PLOT_FULL_NAME)
      {
        plot_name <- paste0(dir_name, acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_noPrior_", method_name, "_rescaled")
      }
      #     plot_signatures(mixtures.rescaled, plot_name=paste0(plot_name, ".pdf"), phis = phis_for_plot, mark_change_points=F, change_points=changepoints,
      #                     assigns_phylo_nodes = assigns_phylo_nodes_sw, error_bars = mixtures.rescaled.err)
      # 
      #     plot_signatures(mixtures.rescaled.mean, plot_name=paste0(plot_name, ".mean.pdf"), phis = phis_for_plot, mark_change_points=F, change_points=changepoints,
      #                   assigns_phylo_nodes = assigns_phylo_nodes_sw, error_bars = mixtures.rescaled.err)
    }
    
    print(proc.time() - pct)
    print(paste("Computed example", example))
  }
}








#get_overall_distribution()

save_data_for_samples()
#suppressMessages(compute_signatures_for_all_examples())
compute_errorbars_for_all_examples()

#compute_mutation_probabilities()
