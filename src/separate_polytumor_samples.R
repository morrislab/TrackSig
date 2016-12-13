vcf_dir = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/vcfs/"
assignments_dir = "/home/q/qmorris/yulia/prostate_cancer_200_fullresults/prostate_cancer_200_ssms/"

polytumor = c("Baca_04-1243_T_wgs",
              "Baca_05-1657_T_wgs",
              "Baca_07-5021_T_wgs",
              "Baca_STID0000003042_T_wgs",
              "CPCG0040",
              "CPCG0047",
              "CPCG0124",
              "CPCG0185",
              "CPCG0190",
              "CPCG0232",
              "CPCG0259",
              "CPCG0362",
              "CPCG0369",
              "CPCG0388",
              "TCGA_7233_T1_wgs",
              "TCGA_7791_T1_wgs")
tumor1 = c("1", "3")
tumor2 = c("2")

separate_polytumor_files <- function()
{
  for (i in 1:length(polytumor))
  {
    id = polytumor[i]
    
    vcf_file = paste0(vcf_dir, "fh_", id, ".vcf")
    vcfData <- read.table(vcf_file, stringsAsFactors = F)
    rownames(vcfData) <- paste(vcfData[,1], vcfData[,2], sep = "_")
    
    assignment_file = paste0(assignments_dir, id, "_ssms.txt")
    as = read.table(assignment_file, stringsAsFactors=F, header = T)
    chr_pos <- as$gene
    chr_pos <- sapply(chr_pos, function(x) paste(strsplit(x, split=":")[[1]][2], strsplit(x, split=":")[[1]][3], sep="_"))
    rownames(as) <- chr_pos
    
    tumor1_pos <- rownames(as[which(as$node == "1" | as$node == "3"),])
    tumor2_pos <- rownames(as[which(as$node == "2"),])
    
    tumor1_vcf <- vcfData[tumor1_pos, ]
    tumor2_vcf <- vcfData[tumor2_pos, ]
    
    write.table(tumor1_vcf, file = paste0(vcf_dir, "fh_", id, ".tumor1.vcf"), row.names = F, col.names=F, quote = F, sep = "\t")
    write.table(tumor2_vcf, file = paste0(vcf_dir, "fh_", id, ".tumor2.vcf"), row.names = F, col.names=F, quote = F, sep = "\t")
  }
}
  
