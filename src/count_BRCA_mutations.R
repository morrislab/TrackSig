# Get cluster assignments of mutations at pseudo-time-points. Cluster assignments are grouped by 100. The cluster of the 100 mutations in determined by majority vote. Script uses mutation assignments, phi files and order of mutations dumped by extractFeaturesPhi.pl

count_BRCA_mutations <- function()
{
	# suppressPackageStartupMessages(library("argparse"))
	# parser <- ArgumentParser()

	# parser$add_argument("--mutass_dir")
	# parser$add_argument("--phi_dir")
	# parser$add_argument("--result_dir")

	# args <- parser$parse_args()
	mutass_dir = "pwgs/mutass/consprelim.sample-5000.cnvint/"
	phi_dir = "pwgs/ssmphi/consprelim.sample-5000.cnvint/"
	mut_order_dir="pwgs/mut_order/"
	result_dir = "pwgs/mutass_sorted_by_phi/"

	BRCA1 = data.frame(chr=17, start=41196312, end=41277500)
	BRCA2 = data.frame(chr=5, start=150522621, end=150570146)

	offset = 2000
	BRCA1_surroundings = BRCA1
	BRCA1_surroundings$start = BRCA1$start - offset
	BRCA1_surroundings$end = BRCA1$end + offset

	BRCA2_surroundings = BRCA2
        BRCA2_surroundings$start = BRCA2$start - offset
        BRCA2_surroundings$end = BRCA2$end + offset

	brca_counts <- c()
	for (mut_order_file in list.files(mut_order_dir))
	{
		counts = data.frame(matrix(0, nrow=1, ncol=4))
		colnames(counts) = c("BRCA1", "BRCA2", "BRCA1_surroundings", "BRCA2_surroundings")

		tumor_id =  gsub("([^/]*)\\.mut_order.txt","\\1", mut_order_file)

		mut <- read.delim(paste0(mut_order_dir, "/", mut_order_file), header=F)

		for (i in 1:nrow(mut))
		{
			if (is_inside(mut[i,1], mut[i,2], BRCA1))
			{
				counts$BRCA1 = counts$BRCA1 + 1
			}

			if (is_inside(mut[i,1], mut[i,2], BRCA2))
                        {
                                counts$BRCA2 = counts$BRCA2 + 1
                        }

			if (is_inside(mut[i,1], mut[i,2], BRCA1_surroundings))
                        {
                                counts$BRCA1_surroundings = counts$BRCA1_surroundings + 1
                        }

			if (is_inside(mut[i,1], mut[i,2], BRCA2_surroundings))
                        {
                                counts$BRCA2_surroundings = counts$BRCA2_surroundings + 1
                        }
		}
		
		brca_counts <- rbind(brca_counts, c(tumor_id, counts))
		print(mut_order_file)
	}

	write.table(brca_counts, file="pwgs/summary_tables/brca_counts.txt", quote=F, row.names=F)

}

is_inside <- function(chr, pos, region)
{
	result = FALSE

	if (chr == region$chr)
		if (pos >= region$start & pos <= region$end)
			result = TRUE

	return(result)
}

#count_BRCA_mutations()

