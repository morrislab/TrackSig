# Get cluster assignments of mutations at pseudo-time-points. Cluster assignments are grouped by 100. The cluster of the 100 mutations in determined by majority vote. Script uses mutation assignments, phi files and order of mutations dumped by extractFeaturesPhi.pl

get_tree_clusters_at_time_points <- function()
{
	# suppressPackageStartupMessages(library("argparse"))
	# parser <- ArgumentParser()

	# parser$add_argument("--mutass_dir")
	# parser$add_argument("--phi_dir")
	# parser$add_argument("--result_dir")

	# args <- parser$parse_args()
	mutass_dir = "pwgs/mutass/consprelim.sample-5000.psub/"
	phi_dir = "pwgs/ssmphi/consprelim.sample-5000.psub/"
	mut_order_dir="pwgs/mut_order/psub/"
	result_dir = "pwgs/mutass_sorted_by_phi/psub/"

	for (phi_file in list.files(phi_dir))
	{
		tumor_id =  gsub("([^/]*)\\.txt","\\1", phi_file)
		mutass_file = paste0(mutass_dir, "/", tumor_id, "_mutation_assignments.txt")
		mut_order_file = paste0(mut_order_dir, "/", tumor_id, ".mut_order.txt")
		
		if (!file.exists(mutass_file))
		{
			print(paste0("File ", mutass_file,  " not found"))
			next
		}
		
		if (!file.exists(mut_order_file))
                {
                        print(paste0("File ", mut_order_file,  " not found"))
                        next
                }

		mutass <- read.delim(mutass_file, header=T)
		phi <- read.delim(paste0(phi_dir, "/", phi_file), header=F)
		mut_order <- read.delim(mut_order_file, header=F)

		stopifnot(nrow(mutass) == nrow(phi))
		stopifnot(nrow(mutass) ==nrow(mut_order))

		phi <- strsplit(sapply(phi[,1], toString), "[_=]+")
	 	phi <- t(as.data.frame(matrix(unlist(phi), nrow=length(unlist(phi[1])))))
	
		
		mutass_ordered <- mutass[order(match(paste(mutass[,1], mutass[,2]),paste(mut_order[,1], mut_order[,2]))),]	
		clusters_ordered <- mutass_ordered[,3]
		
		n_hundreds <- length(clusters_ordered) %/% 100
		# + (length(clusters_ordered) %% 100 != 0)
		clusters_ordered_by_100 = c()
		for (i in 1:n_hundreds)
		{
			clusters_100 <- clusters_ordered[((i-1) * 100): min(length(clusters_ordered), (i * 100))]
			max_cluster <- names(which.max(table(clusters_100)))
			clusters_ordered_by_100 <- c(clusters_ordered_by_100, max_cluster)
		}

		result_file = paste0(result_dir, "/", tumor_id, ".tree_clusters_by100.txt")
		print(result_file)
		write(clusters_ordered_by_100, file=result_file)
	}

}

#get_tree_clusters_at_time_points()
