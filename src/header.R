# AUTHOR: Yulia Rubanova
# 2018

options( expressions = 5e5 )

library(reshape2)
library(ggplot2)
library(NMF)

nmf.options(grid.patch=TRUE)

# Select fitting all signatures or only signatures for the particular cancer type
#sig_amount <- "onlyKnownSignatures" # recommended
sig_amount <- "onlyKnownSignatures" # not recommended, time-consuming

# if the signatures are specified per cancer type or per sample
cancer_type_signatures = T

# if signatures trajectories need to be computed on bootstrapped signatures as well
# bootstrapping provides the uncertainty estimations on the trajectories
# warning: by default, mutations are bootstrapped 30 times and the script will run 30 time longer
compute_bootstrap = FALSE

sliding_window = FALSE
noise_sig = NULL

simulated_data = FALSE

postfix = ""

# file with cancer types of each sample
tumortype_file <- "data/tumortypes.txt"

if (simulated_data) {
  DIR_COUNTS = "./simulated_data/"
  mutation_order = NULL
  BOOTSTRAP_COUNTS = NULL
  purity_file <- NULL
} else {
  # folders with mutation counts, mutation order and bootstrapped mutations
  # don't need to be changed unless different folder were specified in make_counts.sh
  DIR_COUNTS = "data/counts/"
  mutation_order = "data/mut_order/"
  BOOTSTRAP_COUNTS = "data/bootstrap/"
  purity_file = "data/example_purity.txt"
}

# folder to write results to
DIR_RESULTS = "results_signature_trajectories/"

# file with signatures definitions
signature_file = "annotation/alexSignatures.txt"

# file with trinucleotide context
trinucleotide_file = "annotation/trinucleotide.txt"

# specifies active signatures in TCGA cancer types
active_signatures_file = "annotation/active_signatures_transposed.txt"

# specifies active signatures in each sample. Contains the active signatures for the example
# active_signatures_file = "annotation/active_in_samples.txt"

SAVED_SAMPLES_DIR = "saved_data/"

PLOT_FULL_NAME = F
mutation_assignments = ""

if (!file.exists(DIR_RESULTS))
{
  dir.create(DIR_RESULTS, recursive=T)
}

src_files <- setdiff(grep(".*R$", list.files(paste0( "src"),full.names = T), value = T), 
                     c(paste0( "src/compute_mutational_signatures.R"),
                       paste0( "src/header.R")))
for (file in src_files)
{
  source(file)
}

list[alex, tumortypes, active_signatures, active_signatures.our_samples] <- 
    load_annotation(tumortype_file, signature_file, active_signatures_file)

