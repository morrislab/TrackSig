# AUTHORS: Yulia Rubanova and Nil Sahin

plotPhiSignatures <- function(map_only_known_signatures = FALSE) {
  # This function takes in the 96 trinucleotide counts of each sample and
  # plots the changes of the weights of the 30 mutational signatures during the tumor evolution.
  # A plot showing the number of times each signature is chosen with the maximum change
  # across all the samples and the tumor types are also outputted.
  
  load("saved_data/annotation_data.RData")
  
  filenames <- Sys.glob(paste0(DIR, "Counts-by-100/*.txt"))
  
  #   closeAllConnections()
  #   cl <- makeCluster(detectCores()-1)
  #   registerDoParallel(cl)
  
  for (vcfFile in filenames) {
    info <- file.info(vcfFile) # Empty files will be ignored
    if (info$size == 0) 
      next
    
    example <- gsub(".*/([^/]*)\\.phi\\.txt","\\1", vcfFile)
    list[vcfFile, vcfData, phis, acronym, dir_name] <- extract_data_for_example(example, DIR_COUNTS)
    
    if (nrow(vcfData) < 6) # Plots with less than 6 lines of data are meaningless so ignored
      next
    
    print(paste("Computing signatures for ", vcfFile))
    
    # Consecutive rows are added in a new matrix (vcf) together to enhance meaningful regression
    # vcf <- merge_data_chunks(vcfData)
    window_size=400
    shift = window_size/100
    vcf <- get_sliding_window_data(vcfData, shift=shift)
    phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
    phis_sliding_window <- phis_sliding_window / (shift+1)
    colnames(vcf) <- round(phis_sliding_window, 3)
    
    alex.t <- alex
    if (map_only_known_signatures)
    {
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(vcfFile, active_signatures.our_samples, alex)
    }
    # If mutational signatures were not found for this tumor type
    if(is.null(alex.t))
      next
    
    dd <- fit_signatures_linear_regression(vcf, alex.t, alpha = 0, lambda_type="lambda.min")
    cp <- get_changepoints_per_row(dd)
    
    plot_dir <- ifelse(map_only_known_signatures, "plots_only_known_signatures/", "plots/")
    plot_dir <- paste0(DIR, plot_dir, tumortypes[tumortypes$ID == vcfFile,]$tumor_type, "/")
    if (!file.exists(plot_dir)) {
      dir.create(plot_dir, recursive = T)
    }
    plot_name <- paste0(plot_dir, vcfFile, ".pdf")
    plot_signatures(dd, plot_name=plot_name, fitted_data = cp, mark_max_signature = F)
  }
}


plot_signatures <- function (dd, plot_name, phis = NULL, fitted_data = NULL, mark_max_signature=F, mark_change_points=F, 
                             change_points=NULL, error_bars = NULL, save=T,
                             assigns_phylo_nodes = NULL, scale = 1, ytitle = "Signature exposure (%)",
                             transition_points = NULL,
                             remove_sigs_below = 0,
                             sig_colors = NULL) {
  if (!is.null(error_bars) & sum(dim(dd) == dim(error_bars)) != 2) {
    stop("Dimentions of error bar matrix should be the same as dimentions of mixture matrix")
  }
  
  sigs_to_remove <- apply(dd, 1, mean) < remove_sigs_below
  dd <- dd[!sigs_to_remove, ]


  # Weight matrix is edited with appropriate row and column names
  signatures <- rownames(dd)
  df <- data.frame(signatures,dd)
  col_names <- c("Signatures")
  n_timepoints = ncol(dd)

  decrement <- as.integer(150/ncol(dd))
  for (n in 1:ncol(dd)) {
    col_names <- append(col_names, 150 - decrement*(n-1))
  }
  
  colnames(df) <- col_names
  
  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  # The signature with the maximum change is printed in the plot with the annotate function on ggplot (can be removed if unnecessary)
  
  df.m <- melt(df, id.vars = "Signatures")
  maxx <- apply(dd, 2, which.max)
  maxy <- apply(dd ,2,max)

  if (!is.null(assigns_phylo_nodes))
  {
    names(assigns_phylo_nodes) <- col_names[-1]

    tree_clusters = factor(assigns_phylo_nodes[sapply(df.m$variable, toString)])
    
    if ("Branch" %in% levels(tree_clusters) & 
      "Trunk" %in% levels(tree_clusters))
    {
      tree_clusters <- factor(tree_clusters, 
                  levels=c("Trunk", "Branch"), order=T) 
    }
    
    df.m$tree_clusters <- tree_clusters
  }

  alpha <- 1
  # if(!is.null(fitted_data))
  # {
  #   alpha <- 0.3
  # }
  
  size = 1
  if (!is.null(fitted_data)) {
    if (is.vector(fitted_data)) {
      size = 1.3
    }
  }

  g <- ggplot(data = df.m, aes(x = variable, y = value , group = Signatures, color = Signatures)) + 
    geom_line(alpha=alpha, size=size) + xlab("Avg number of mutant alleles per cancer cell") + ylab(ytitle) + 
    geom_point(alpha=alpha, size=size) + 
    theme_bw() + theme(text = element_text(size = 20)) +
    theme(axis.title = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 15)) 
    #scale_color_manual(values=COLORS2[-c(4,9)])
  
  if (!is.null(sig_colors)) {
    g <- g + scale_colour_manual(values=sig_colors)
  }

  if (!is.null(phis))
  {
    breaks = as.numeric(col_names[-1])
    labels=paste(round(phis,2),sep = "\n")

    # if there are too many time points, display only every other value
    if (length(col_names) > 60) {
        # Display every forth label
        labels[seq(2,length(labels),4)] <- ""
        labels[seq(3,length(labels),4)] <- ""
        labels[seq(4,length(labels),4)] <- ""
    } else if (length(col_names) > 35) {
        # Display every third label
        labels[seq(2,length(labels),3)] <- ""
        labels[seq(3,length(labels),3)] <- ""
    }
    g <- g + scale_x_discrete(breaks = breaks, labels=labels)
  }

  if (mark_max_signature)
  {
    g <- g + annotate("text", x=col_names[-1], y=maxy*1.05, label=rownames(dd)[maxx]) 
  }
  
  if (!is.null(fitted_data))
  {
    if (!is.vector(fitted_data))
    {
        fitted_d <- fitted_data[[i]][!sigs_to_remove,]

        fitted_data.m <- data.frame(Signatures = rownames(fitted_d),fitted_d)
        colnames(fitted_data.m) <- colnames(df)
        fitted_data.m <- melt(fitted_data.m, id.vars = "Signatures")
        g <- g + geom_line(data=fitted_data.m, size=0.5, aes(x=variable, y=value, group = Signatures, color = Signatures),  alpha=alpha)
    } else {
      for (i in 1:length(fitted_data))
      {
        alpha = 0.3

        fitted_d <- fitted_data[[i]][!sigs_to_remove,]

        fitted_data.m <- data.frame(Signatures = rownames(fitted_d),fitted_d)
        colnames(fitted_data.m) <- colnames(df)
        fitted_data.m <- melt(fitted_data.m, id.vars = "Signatures")
        g <- g + geom_line(data=fitted_data.m, size=0.5, aes(x=variable, y=value, group = Signatures, color = Signatures),  alpha=alpha)
      }
    }
  }
  
  if (mark_change_points)
  {
    if (is.null(change_points))
      stop("Please provide change points to mark in the plot")
    
    #g <- g + geom_vline(xintercept = change_points, size = 0.7)
    if (is.null(change_points))
      stop("Please provide change points to mark in the plot")
    
    #g <- g + geom_vline(xintercept = change_points, size = 1, show.legend = T)

    for (i in 1:length(change_points)) {
      g <- g +  annotate("rect", xmax=change_points[i]-1, 
          xmin=change_points[i], ymin=-Inf, ymax=Inf, alpha=0.3) 
    }

  }
  
  if (!is.null(error_bars)) {
    error.min <- dd - error_bars
    error.max <- dd + error_bars
    
    error.min.m <- data.frame(signatures,error.min)
    colnames(error.min.m) <- colnames(df)
    error.min.m <- melt(cbind(Signatures=df[,1],error.min.m), id.vars = "Signatures")
    
    error.max.m <- data.frame(signatures,error.max)
    colnames(error.max.m) <- colnames(df)
    error.max.m <- melt(cbind(Signatures=df[,1],error.max.m), id.vars = "Signatures")
    
    g <- g + geom_errorbar(aes(ymin=error.min.m$value, ymax=error.max.m$value), width=.2)
  }
  
  if (!is.null(assigns_phylo_nodes))
  {
    g <- g + facet_grid(. ~ tree_clusters, scales = "free_x", space="free_x")
  }

  if (!is.null(transition_points)) {
    g <- g + geom_vline(xintercept = transition_points, colour="red", size=1.5)
  }

  # add the below annotate function to print the signature with the maximum change
  # annotate("text", x=which.max(dd[maxx,]), y=max(dd[maxx,])+0.01, label=paste("S",maxx,sep=""))
  # theme(legend.position = "none")
  
  if (save) {
    # was: width = 12*scale
    suppressWarnings(ggsave(filename = plot_name, width = 12*scale, height=5))
  }
  
  return(list(plot = g, data = df))
}


get_group_colors_all_sigs <- function() {
  popular = c("1","8", "5","40","18","9")
  group.colors = c(gg_color_hue(length(popular)), sample(gg_color_hue(ncol(alex) - length(popular)),ncol(alex)- length(popular)))
  names(group.colors) <- c(popular, setdiff(colnames(alex), popular))
  return(group.colors)
}




plot_signatures_real_scale <- function (dd, plot_name, phis = NULL, fitted_data = NULL, mark_max_signature=F, mark_change_points=F, 
                             change_points=NULL, error_bars = NULL, save=T, ytitle = "Signature exposure (%)",
                             assigns_phylo_nodes = NULL,  transition_points = NULL, remove_sigs_below = 0, cut_at_range = NULL, sig_colors = NULL) {
  if (!is.null(error_bars) & sum(dim(dd) == dim(error_bars)) != 2) {
    stop("Dimentions of error bar matrix should be the same as dimentions of mixture matrix")
  }
  
  sigs_to_remove <- apply(dd, 1, mean) < remove_sigs_below
  dd <- dd[!sigs_to_remove, ]
  if (!is.null(cut_at_range)) {
    dd <- truncate_to_range(dd, cut_at_range)
  }
  
  # Weight matrix is edited with appropriate row and column names
  signatures <- rownames(dd)
  df <- data.frame(signatures,dd)
  col_names <- c("Signatures")
  
  decrement <- as.integer(150/ncol(dd))
  for (n in 1:ncol(dd)) {
    col_names <- append(col_names, 150 - decrement*(n-1))
  }
  
  #colnames(df) <- col_names
  colnames(df) <- c("Signatures", round(phis,3 ))
  
  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  # The signature with the maximum change is printed in the plot with the annotate function on ggplot (can be removed if unnecessary)
  
  df.m <- melt(df, id.vars = "Signatures")
  df.m$variable <- sapply(df.m$variable, function(x) as.numeric(toString(x)))
  maxx <- apply(dd, 2, which.max)
  maxy <- apply(dd ,2,max)
  
  if (!is.null(fitted_data))
  {
    if (is.vector(fitted_data))
    {
      all_phis <- c()
      for (i in 1:length(fitted_data)) 
      {
        all_phis <- c(all_phis, round(sapply(colnames(fitted_data[[i]]), function(x) as.numeric(toString(x))),3))
      }
      all_phis <- unique(all_phis)
    }
  }
  
  if (!is.null(assigns_phylo_nodes))
  {
    names(assigns_phylo_nodes) <- round(phis, 3) #col_names[-1]
    tree_clusters = factor(assigns_phylo_nodes[sapply(df.m$variable, toString)])
    
    order_of_clusters <- c("1", "2", "3", "4")
    
    for (i in 1:length(levels(tree_clusters)))
    {
      tree_clusters <- factor(tree_clusters)
      if (order_of_clusters[i] %in% levels(tree_clusters))
      {
        
        if (i == 1)
        {
          current_boundary <- min(as.numeric(names(tree_clusters[tree_clusters == order_of_clusters[i]])))
          phis_to_add <- all_phis[all_phis > current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t)
        }
        if (i == length(levels(tree_clusters)))
        {
          phis_to_add <- all_phis[all_phis < current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t) 
        } else {
          current_boundary_upper <- min(as.numeric(names(tree_clusters[tree_clusters == order_of_clusters[i]])))
          phis_to_add <- all_phis[all_phis > current_boundary]
          phis_to_add.t <- rep(order_of_clusters[i], length(phis_to_add))
          names(phis_to_add.t) <- phis_to_add
          tree_clusters <- c(tree_clusters,  phis_to_add.t)
          
          current_boundary <- current_boundary_upper
        }
      }
    }
    
    tree_clusters <- factor(tree_clusters)
    
    if ("Branch" %in% levels(tree_clusters) & 
          "Trunk" %in% levels(tree_clusters))
    {
      tree_clusters <- factor(tree_clusters, 
                              levels=c("Trunk", "Branch"), order=T) 
    }
    tree_clusters <- tree_clusters[order(as.numeric(names(tree_clusters)), decreasing = T)]
    names(tree_clusters) <- round(as.numeric(names(tree_clusters)),3)
    df.m$tree_clusters <- tree_clusters[sapply(round(df.m$variable,3), toString)]
  }

  alpha <- 1
  #   if(!is.null(fitted_data))
  #   {
  #     alpha <- 0.3
  #   }
  
  size = 1
  if (!is.null(fitted_data)) {
    if (is.vector(fitted_data)) {
      size = 1.5
    }
  }
  group.colors = get_group_colors_all_sigs()

  g <- ggplot(data = df.m, aes(x = variable, y = value, group = Signatures, color = Signatures)) + 
   geom_line(alpha=alpha, size=1.7) + xlab("Avg number of mutant alleles per cancer cell") + ylab(ytitle) + 
    geom_point(alpha=alpha, size=1.7) + 
    theme_bw() + theme(text = element_text(size = 20)) +
    theme(axis.title = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 15))  +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) 

  if (!is.null(sig_colors)) {
    g <- g + scale_colour_manual(values=sig_colors)
  }

  #   if (!is.null(phis))
  #   {
  #     g <- g + scale_x_discrete(breaks = as.numeric(col_names[-1]), labels=paste(round(phis,2),sep = "\n"))
  #   }
  
  g <- g + scale_x_reverse()
  
  if (mark_max_signature)
  {
    g <- g + annotate("text", x=col_names[-1], y=maxy*1.05, label=rownames(dd)[maxx]) 
  }
  
  if (!is.null(fitted_data))
  {
    if (is.vector(fitted_data))
    {
      for (i in 1:length(fitted_data))
      {
        alpha = 0.3
        
          if (!is.null(cut_at_range)) {
            fitted_data[[i]] <- truncate_to_range(fitted_data[[i]], cut_at_range)
          }
          if (ncol(fitted_data[[i]]) == 0) {
            next
          }

        fitted_data.m <- data.frame(signatures,fitted_data[[i]][!sigs_to_remove,])
        #colnames(fitted_data.m) <- colnames(df)
        colnames(fitted_data.m) <- c("Signatures", round(as.numeric(colnames(fitted_data[[i]])), 3))
        fitted_data.m <- melt(cbind(Signatures=df[,1],fitted_data.m), id.vars = "Signatures")
        fitted_data.m$variable <- sapply(fitted_data.m$variable, function(x) as.numeric(toString(x)))
        if (!is.null(assigns_phylo_nodes))
        {
          fitted_data.m$tree_clusters <- tree_clusters[sapply(round(fitted_data.m$variable,3), toString)]
        }
        g <- g + geom_line(data=fitted_data.m, size=0.7, aes(x=variable, y=value, group = Signatures, color = Signatures),  alpha=alpha)
      }
    } else {
      fitted_data.m <- data.frame(signatures,fitted_data[!sigs_to_remove,])
      colnames(fitted_data.m) <- colnames(df)
      fitted_data.m <- melt(cbind(Signatures=df[,1],fitted_data.m), id.vars = "Signatures")
      fitted_data.m$variable <- sapply(fitted_data.m$variable, function(x) as.numeric(toString(x)))
      if (!is.null(assigns_phylo_nodes))
      {
        fitted_data.m$tree_clusters <- tree_clusters[sapply(round(fitted_data.m$variable,3), toString)]
      }
      g <- g + geom_line(data=fitted_data.m, size=0.7, aes(x=variable, y=value, group = Signatures, color = Signatures)) 
    }
  }
  
  if (mark_change_points)
  {
    if (is.null(change_points))
      stop("Please provide change points to mark in the plot")
    
    #g <- g + geom_vline(xintercept = change_points, size = 1, show.legend = T)

      # if change points are in the list from various bootstrap runs, show all of them and adjust transparency
      if (class(change_points) == "list") {
          alpha = 0.5/length(change_points)
          change_points <- unlist(change_points)
        } else {
          alpha=0.3
        }
      for (i in 1:length(change_points)) {
        g <- g +  annotate("rect", xmin=round(phis,3 )[change_points[i]-1], xmax=round(phis,3 )[change_points[i]], ymin=-Inf, ymax=Inf, alpha=alpha) 
      }
  }
  
  if (!is.null(error_bars)) {
    error.min <- dd - error_bars
    error.max <- dd + error_bars
    
    error.min.m <- data.frame(signatures,error.min)
    colnames(error.min.m) <- colnames(df)
    error.min.m <- melt(cbind(Signatures=df[,1],error.min.m), id.vars = "Signatures")
    
    error.max.m <- data.frame(signatures,error.max)
    colnames(error.max.m) <- colnames(df)
    error.max.m <- melt(cbind(Signatures=df[,1],error.max.m), id.vars = "Signatures")
    
    g <- g + geom_errorbar(aes(ymin=error.min.m$value, ymax=error.max.m$value), width=.2)
  }

  if (!is.null(assigns_phylo_nodes))
  {
    g <- g + facet_grid(. ~ tree_clusters, scales = "free_x", space="fixed")
  }
  
  if (!is.null(transition_points)) {
    g <- g + geom_vline(xintercept = round(phis,3 )[transition_points], colour="red", size=1.5)
  }

  g <- g + geom_vline(xintercept = c(phis), alpha=0.3)
  
  # add the below annotate function to print the signature with the maximum change
  # annotate("text", x=which.max(dd[maxx,]), y=max(dd[maxx,])+0.01, label=paste("S",maxx,sep=""))
  # theme(legend.position = "none")
  
  if (save) {
    suppressWarnings(ggsave(filename = plot_name, width = 14, height=4))
  }
  
  return(list(plot = g, data = df))
}


run_example <- function (example) {

  list[vcfFile, vcfData, phis, acronym, dir_name] <- extract_data_for_example(example, DIR_COUNTS)
  
  if (nrow(vcfData) < 6) # Plots with less than 6 lines of data are meaningless so ignored
    return()
  
  make_histogram(t(vcfData), paste0(dir_name, "mutation_hist.pdf"))
  
  #   for (data_method in c("merged200", "sliding"))
  #   {
  data_method <- "sliding"
  if (data_method == "data_method")
  {
    # Consecutive rows are added in a new matrix (vcf) together to enhance meaningful regression
    vcf <- merge_data_chunks(vcfData)
  }
  if (data_method == "sliding")
  {
    # Make a sliding window of 400 mutations
    window_size=400
    shift = window_size/100
    vcf <- get_sliding_window_data(vcfData, shift=shift)
    phis_sliding_window <- get_sliding_window_data(toVerticalMatrix(phis), shift=shift)
    phis_sliding_window <- phis_sliding_window / (shift+1)
    colnames(vcf) <- round(phis_sliding_window, 3)
    
    data_method <- paste0(data_method, window_size)
  }
  
  for (sig_amount in c("all", "onlyKnownSignatures"))
  {
    if (sig_amount == "onlyKnownSignatures") {
      # Fit only known signatures
      list[alex.t, matched_type, acronym] <- get_signatures_for_current_sample(vcfFile, active_signatures.our_samples, alex)
    } else {
      alex.t <- alex
    }
    
    for (prior_type in c("noPrior", "overallMixturesPrior"))
    {
      if (prior_type == "noPrior")
        prior = NULL
      
      if (prior_type == "overallMixturesPrior")
      {
        # Fit linear regression and mixture model with prior of overall signature mixtures over all time points
        overall_mixtures <- fit_mixture_of_multinomials_EM(apply(vcfData, 2, sum), as.matrix(alex.t))
        prior = overall_mixtures
      }
      
      
      dd <- fit_mixture_of_multinomials_matrix(vcf, alex.t, prior = prior)
      cp <- get_changepoints_per_row(dd)
      plot_signatures(dd, plot_name=paste0(dir_name, acronym, "_", data_method, "_multMix_", sig_amount, "_", prior_type, ".pdf"), fitted_data = cp)
      
      overall_change_points <- get_changepoints_matrix_jointly(dd, lambda_type="lambda.1se")$change_points
      cp <- fit_mixture_of_multinomials_in_time_slices(vcf, overall_change_points, alex.t)
      plot_signatures(dd, plot_name=paste0(dir_name, acronym, "_", data_method, "_multMix_fittedPerTimeSlice_", sig_amount, "_", prior_type, ".pdf"), 
                      fitted_data = cp, mark_max_signature=F)
      
      for (alpha in c(0,1))
      {
        regression_type = alpha
        if (alpha == 0)
          regression_type = "ridge"
        if (alpha == 1)
          regression_type = "lasso"
        
        for (lambda_type in c("lambda.min", "lambda.1se","lambda0"))
        {
          dd <- fit_signatures_linear_regression(vcf, alex.t, alpha = alpha, lambda_type=lambda_type, prior = prior)
          cp <- get_changepoints_per_row(dd)
          plot_name <- paste0(dir_name, acronym, "_", data_method, "_", regression_type, "_", lambda_type, "_", sig_amount, "_", prior_type)
          #plot_signatures(cp, plot_name=paste0(plot_name, "_piecewise.pdf"), mark_max_signature=F)
          plot_signatures(dd, plot_name=paste0(plot_name, ".pdf"), fitted_data = cp)
          
          if (regression_type == "lasso")
          {
            dd <- fit_signatures_linear_regression(vcf, alex.t, alpha = 1, lambda_type=lambda_type, prior = prior)
            cp <- get_changepoints_per_row(dd)
            non_zero_signatures <- names(sort(apply(dd, 1 , sum), decreasing = T)[1:5])
            dd <- fit_signatures_linear_regression(vcf, alex.t[,non_zero_signatures], alpha = 0, lambda_type=lambda_type)
            cp <- get_changepoints_per_row(dd)
            plot_signatures(dd, plot_name=paste0(dir_name, acronym, "_sliding", window_size, "_lasso_lambda.min_removedZeroSignatures.pdf"), fitted_data = cp)
          }
        }
      }
    }
  }
}

make_histogram <- function(vcf, plot_name)
{ 
  pdf(plot_name, width=8, height=3)
  make_histogram_plot(vcf)
  dev.off()
}


make_histogram_plot <- function(vcf, horiz = FALSE)
{
  mutation_types <- read.delim("trinucleotide.txt", header=F)
  mutation_types.str <- as.factor(paste(mutation_types[,1], ">", mutation_types[,2]))
  
  d <- apply(vcf, 1, sum)
  d <- d / sum(abs(d))
  names(d) <- mutation_types[,3]
  
  par(mar=c(2,1,1,1))
  xlim = NULL
  ylim = NULL
  if (horiz) {
    min = 0
    if (sum(vcf < 0) > 0) {
      min = -0.05
    }
    xlim = c(min, 0.05)
  } else {
    ylim = c(0, 0.05)
  }
  b <- barplot(d, cex.axis=0.4, cex.names=0.4, 
               col=as.numeric(mutation_types.str) + 1, 
               #ylab = "Mutation type probability",
               ylim = ylim, xlim = xlim, las=2, yaxt="n", horiz=horiz)
  
  #axis(2,cex.axis=0.7)
  
  legend("topright", 
         legend = levels(mutation_types.str), 
         fill = (1:length(mutation_types.str))+1,
         cex = 0.6)
}

plot_rnn_results <- function(path)
{
  omit_first_timepoints = 3
  
  initial_data_path <- paste0(path, "/", "initial_data/")
  predictions_path <- paste0(path, "/", "predictions/")
  if (!(file.exists(initial_data_path) & file.exists(predictions_path))) {
    stop("Either initial_data or predictions dir does not exist")
  }
  
  plot_dir <- paste0(path, "/", "plots/")
  suppressWarnings(dir.create(plot_dir))
  
  for (file in list.files(initial_data_path))
  {
    if (!file.exists(paste0(predictions_path, file))) {
      stop(paste("File", file, "does not exist in predictions"))
    }
    initial_data <- read.delim(paste0(initial_data_path, file), header=F)
    predictions <- read.delim(paste0(predictions_path, file), header=F)
    rownames(initial_data) <- rownames(predictions) <- paste0("S", 1:30)
    
    non_zero_sigs <- apply(initial_data, 1, sum) != 0
    non_zero_time_points <- apply(initial_data, 2, sum) != 0
    initial_data <- initial_data[non_zero_sigs, non_zero_time_points]
    predictions <- predictions[non_zero_sigs, non_zero_time_points]
    
    predictions[,1:omit_first_timepoints] <- NA
    
    plot_name <- paste0(plot_dir, gsub("(.*).txt","\\1", file), ".pdf")
    plot_signatures(initial_data, plot_name=plot_name, fitted_data = predictions)
  }
}

add_change_point_distribution <- function(p, changepoints_bootstrap, n_timepoints, ymax, plot_name, save=T, rolsum = 1, scale = 100)
{
  cp_distr <- data.frame(rep(0, n_timepoints))
  rownames(cp_distr) <- 1:n_timepoints
  cp_distr[names(table(unlist(changepoints_bootstrap))),] <- table(unlist(changepoints_bootstrap))
  cp_distr <- cp_distr[,1]
  #pdf(paste0(plot_name, ".cp_distr.pdf"), width=10, height=4)
  #barplot(unlist(cp_distr), names.arg=round(phis_sliding_window,3), las=2)
  #dev.off()
  
  if (rolsum > 1)
  {
    cp_distr_rs <- rollapply(cp_distr, rolsum, sum)
    for (i in seq((rolsum-1),1,-1))
    {
      cp_distr_rs <-  c(sum(cp_distr[1:i]), cp_distr_rs)     
    }
    cp_distr <- cp_distr_rs
  }
  
  g_prime <- p
  for (i in 1:length(cp_distr)) { 
    g_prime <- g_prime + geom_rect(xmin = (i-0.2), xmax = (i+0.2), ymin=-0.06*scale, 
                                   ymax = (cp_distr[i]/max(cp_distr)*0.05-0.06)*scale) + guides(fill=FALSE)
  }
  g_prime <- g_prime + guides(color=guide_legend(override.aes=list(fill=NA))) + ylim(c(-0.05*scale, ymax))
  
  if (save) {
    suppressWarnings(ggsave(filename = plot_name, width = 12, height=8))
  }
  
  return(g_prime)
}

plot_cummulative <- function(vector, plot_file, axis_title = NULL)
{
  pdf(plot_file, width = 8, height=5)
  step = (max(vector) - min(vector)) / 1000
  h <- hist(vector, breaks=seq(min(vector),max(vector),step), plot=F)
  plot(h$breaks[-1], cumsum(h$density) / sum(h$density), main="", 
       ylab = "Cumulative ratio", xlab = axis_title, col="gray", type="l",
       lwd =3,  xaxt = "n")
  
  if (min(vector) < 0)
  {
    ticks <- c(sign(min(vector)) * seq(0, abs(min(vector)), 0.05), seq(0, abs(max(vector)),0.05) * sign(max(vector)))
    ticks <- round(sort(ticks), digits = 3)
  } else {
    ticks <- seq(round(min(vector),1), max(vector)+0.05,0.05)
    ticks <- round(sort(ticks), digits = 3)
  }
  
  axis(side = 1, at = ticks)
  dev.off()
}

make_cumulative_plot <- function(data, xlab, ylab, font = 1) {
  plot(1,main="", type="n", xlim=c(0,max(sapply(data,max))), ylim=c(0,1),
     ylab = ylab, xlab = xlab,  cex.axis=font, cex.lab = font, cex=font)

  colors = gg_color_hue(length(data))
  colors[seq(2,length(colors),2)] = COLORS2[(1:(length(data)/2+2))[-c(4,9)]]
  for (i in 1:length(data)) {
    h <- hist(data[[i]], breaks=seq(0,1,0.02), plot=F)
    lines(h$breaks, c(0,cumsum(h$density) / sum(h$density)), col=colors[i], type="l", lwd =2)
  }

  legend('bottomright',names(data),
         fill = colors, bty = 'n',
         border = NA, cex=font*0.7)
}

