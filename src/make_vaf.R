args <- commandArgs(trailingOnly = TRUE)

vcf_dir = "data/vcf/"  
vaf_dir = "data/phi/"


if (length(args) > 0) {
	vcf_dir <- paste0(args[1], "/")
}

if (length(args) > 1) {
	vaf_dir <- paste0(args[2], "/")
}

print(paste("VCF dir:", vcf_dir))
print(paste("VAF dir:",vaf_dir))

get_vaf <- function(info)
{
  info <- strsplit(info, ";")[[1]]

  vaf <- as.double(gsub("VAF=(.*)", "\\1", info[grep("VAF=(.*)", info)]))
  vaf <- vaf * 2
  
  return(vaf)
}

make_vafs_for_mutations <- function(vcf_dir, vaf_dir)
{
  for (f in list.files(vcf_dir))
  {
    tumor_id <- gsub("(.*).vcf", "\\1", f)
      
    vcf <- read.delim(paste0(vcf_dir, f), header = T, stringsAsFactors=F)
    info <- vcf[,ncol(vcf)]

    if (sum(is.na(info) > 0)) {
      stop("VCF file does not contain VAF info")
    }

    if (sum(grepl("VAF=(.*)", info)) == 0) {
      stop("VCF file does not contain VAF info")
    }

    filter <-  vcf[,"FILTER"] == "." & grepl("VAF=(.*)", info)
    
    vaf_data <- paste(vcf[,c("X.CHROM")], vcf[,c("POS")], sep="_")
    vaf_values <- as.matrix(sapply(info, get_vaf), ncol=1)

    empty_vaf <- sapply(vaf_values, length) == 1
    
    vaf_data <- paste(vaf_data, vaf_values, sep="=")
    vaf_data <- vaf_data[filter & empty_vaf]
    
    if (!file.exists(vaf_dir)) {
    	dir.create(vaf_dir)
    }
    write.table(vaf_data, file=paste0(vaf_dir, tumor_id, ".txt"), col.names=F, row.names = F, quote=F)
  }
} 

make_vafs_for_mutations(vcf_dir, vaf_dir)