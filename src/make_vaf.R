vcf_dir = "/home/q/qmorris/yulia/pwgs/consprelim/"  
vaf_dir = "/scratch/q/qmorris/yulia/pwgs/ssmvaf_full/"
phi_dir = "/home/q/qmorris/yulia/pwgs/ssmphi/consprelim.sample-5000.psub/"

get_vaf <- function(info)
{
  info <- strsplit(info, ";")[[1]]
  vaf <- as.double(gsub("VAF=(.*)", "\\1", info[grep("VAF=(.*)", info)]))
  vaf <- vaf * 2
  
  return(vaf)
}

make_vafs_for_mutations <- function()
{
  for (f in list.files(vcf_dir))
  {
    tumor_id <- gsub("(.*).annotated.snv_mnv.vcf", "\\1", f)
      
    vcf <- read.delim(paste0(vcf_dir, f), header = T, stringsAsFactors=F)
    info <- vcf[,ncol(vcf)]
    
    filter <-  vcf[,"FILTER"] == "."
    
    vaf_data <- paste(vcf[,c("X.CHROM")], vcf[,c("POS")], sep="_")
    vaf_values <- as.matrix(sapply(info, get_vaf), ncol=1)
    empty_vaf <- sapply(vaf_values, length) == 1
    
    vaf_data <- paste(vaf_data, vaf_values, sep="=")
    vaf_data <- vaf_data[filter & empty_vaf]
    
    
    write.table(vaf_data, file=paste0(vaf_dir, tumor_id, ".txt"), col.names=F, row.names = F, quote=F)
  }
} 

make_vafs_for_mutations_subset_with_phi  <- function()
{
  # Make vaf counts only for tumors and mutations where phi is available.
  for (f in list.files(vcf_dir)[2771:length(list.files(vcf_dir))])
  {
    tumor_id <- gsub("(.*).annotated.snv_mnv.vcf", "\\1", f)
    print(paste(tumor_id, " Sample ", which( list.files(vcf_dir) == f), " out of ", length( list.files(vcf_dir))))
    # Skipping tumors where we don't have phis
    if (!file.exists(paste0(phi_dir, tumor_id, ".txt"))) {
      print("Skipping...")
      next
    }
    
    vcf <- read.delim(paste0(vcf_dir, f), header = T, stringsAsFactors=F)
    info <- vcf[,ncol(vcf)]
    
    phis <- read.delim(paste0(phi_dir, tumor_id, ".txt"), header = F, stringsAsFactors=F, sep="=")
    
    filter <-  vcf[,"FILTER"] == "."
    
    vaf_data <- paste(vcf[,c("X.CHROM")], vcf[,c("POS")], sep="_")
    filter_no_phi <- vaf_data %in% phis[,1]
      
    vaf_values <- as.matrix(sapply(info, get_vaf), ncol=1)
    empty_vaf <- sapply(vaf_values, length) == 1
    
    vaf_data <- paste(vaf_data, vaf_values, sep="=")
    vaf_data <- vaf_data[filter & empty_vaf & filter_no_phi]
    
    write.table(vaf_data, file=paste0(vaf_dir, tumor_id, ".txt"), col.names=F, row.names = F, quote=F)
  }
} 