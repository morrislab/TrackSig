# Trackature

A framework to infer mutational signatures over time.

## Background
Cell processes leave a unique signature of mutation types in cancer genome. Using the mutational signatures, it is possible to infer the fraction of mutations contributed by each cell process. Mutational signatures are represented as multinomial distributions over 96 mutation types. Using our framework Trackature, we can infer mutational signatures changing over time.

## Dependencies
- Python 2.7.9  
  Packages: pyvcf (pip2 install pyvcf)
  
- Perl v5.18.2  
  Packages: Bio::DB::Fasta. Please refer to http://www.cpan.org/modules/INSTALL.html on how to install packages on your machine. On Mac OS and Unix: sudo cpan Bio::DB::Fasta
  
- R 3.1.2  
  Packages: reshape2, ggplot2, NMF. R packages can be installed using install.packages("package_name") command.


## Input
- VCF file with point mutations  
  VCF file should be named as tumor_id.vcf, where tumor_id is an id of the tumor.  
  INFO column in VCF file should contain reference and alternate read counts in the following format: "t_alt_count=5;t_ref_count=20". 

Optional:  
- Sample purity  
  Please refer to the example for the format of tumor purity file. The file should be tab-delimited and contain at least two columns: "samplename" and "purity". The "samplename" column should match the name of the vcf file.

- Copy number alteration calls

## Usage 
The commands below assume starting from the 'Trackature' directory. The example.vcf is provided in the repo. Running the code as written below will compute signature trajectories for the example.

### Generating variant allele frequency estimates

We use variant allele frequency (VAF) to sort mutations by the order of their occurrence.

To generate VAF values:
```
python src/make_corrected_vaf.py --vcf data/example.vcf --output data/example_vaf.txt
```

It is recommended to correct VAF by copy number and tumor purity if those are available. You can specify file CNA calls and a file containing sample purities the following way:
```
python src/make_corrected_vaf.py --vcf data/example.vcf --cnv your_cna_call_file.txt --purity purity_file.txt --output data/example_vaf.txt
```
Please refer to the example for the format of tumor purity file. The file should be tab-delimited and contain at least two columns: "samplename" and "purity". The "samplename" column should match the name of the vcf file.

To make use of purities when correcting VAF, provide the name of the purity file to make_corrected_vaf.py with parameter "--purity".

### Making mutation counts

To make mutation counts over 96 mutation types:
```
src/make_counts.sh data/example.vcf data/example_vaf.txt
```
where the first parameter is the vcf file and the second parameter is file with VAF values generated at the previous step.

Requires hg19 reference which can be downloaded from here: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/





### Generating signature trajectories
To use cancer-type specific signatures, please provide `data/tumortypes.txt` file listing tumor IDs and their cancer types in the format (tab-delimited):
```
ID	tumortype
example	LAML
```
The names of the cancer types must match the ones in `annotation/active_signatures_transposed.txt` table which lists active signatures for TCGA cancer types.

Modify `src/header.R` to set up the path to your data
```
tumortype_file <- "data/tumortypes.txt"
```

Compute signature trajectories for all samples:
```
Rcript src/compute_mutational_signatures.R
```

Results can be found in "results_signature_trajectories" folder (by default, specified in by DIR_RESULTS in `src/header.R`) in appropriate cancer type and tumor id folders. Signature trajectories are stored in `mixtures.csv`. Rows correspond to signatures. Columns correspond to time points. The columns are named by the average cellular prevalence that corresponds to the time point.

If you wish to compute uncertainty for trajectories as well, set `compute_bootstrap parameter` in `src/header.R` to TRUE before running the script (slows down the computation).

## Other functionality

### Computing overall signature exposures across all mutations

Mean signature exposures across all mutations in the tumor can be computed the following way:
```
R
source("src/header.R")
compute_overall_exposures_for_all_examples()
```
It will compute the signature exposures (aka "mixtures") for all samples in the data/counts directory. The results can be found in the results directory under the appropriate cancer type and tumor ID ( `overall_mixtures.csv`).

Alternatively, you can call a function to compute exposures directly:
```
mixtures <- fit_mixture_of_multinomials_EM(mutation_counts, alex.t)
```
Input to fit_mixture_of_multinomials_EM:
- mutation_counts : 1x96 data frame specifying mutation counts over 96 types
- alex.t : specifies signatures to fit (96 by number of signatures). 

### Determining active signatures

**Using COSMIC per-cancer active signatures**  
Active signatures vary across cancer types. To ensure that we fit only most relevant signatures for the particular cancer type, we refer to the table of active signatures from [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures), reproduced in `annotation/active_signatures_transposed.txt` file.

To use cancer-type specific signatures, please provide `data/tumortypes.txt` file listing tumor IDs and their cancer types (see example in the repo). The names of the cancer types must match the onces in `active_signatures_transposed.txt` table. 

**Estimating active signatures from scratch**

If active signatures are unavailable, they can be estimated by computing mean exposures across all mutations (see "Computing overall signature exposures" section) and taking the signatures with highest exposures (for example, with exposure > 10%). Next, these signatures should be specified as active in `src/header.R` before running `src/compute_mutational_signatures.R`. 

We recommmend to estimate active signatures first (as described above) instead of fitting all signatures over time through Trackature, as it provides more stable results and speeds up the computation.

### Providing per-tumor signatures

To specify a separate set of active signatures for each sample:

1) in `src/header.R` set cancer_type_signatures = FALSE   
2) in `src/header.R` set active_signatures_file to the file with active signatures per sample. See example of such file in `annotation/active_in_samples.txt`.

### Providing new signatures

Signatures provided in the repo are from [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures) located in `annotation/alexSignatures.txt`. You can use your own signatures by providing path to another signature file through `signature_file` parameter in `src/header.R`.

## Notes
Please note that Trackature does not run on samples with less than 600 mutations. Less than 600 mutations will result in less than 3 time points, and there is no point to analize it as a time series. On tumors with less than 600 mutations, you can compute signature exposures without dividing mutations into time points (see "Computing overall signature exposures" section).





