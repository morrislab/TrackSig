# AUTHOR: Yulia Rubanova

# main function: get change points from vector of time points using change points lasso
get_changepoints_vector <- function(data, lambda_type="lambda.1se")
{
  X_dim = length(data)
  
  X <- matrix(0, ncol=X_dim, nrow=X_dim)
  for (i in 1:X_dim)
  {
    for (j in 1:i)
      X[i,j] <- 1
  }
  
  data_to_fit <- data
  
  glm <- cv.glmnet(X, data_to_fit, alpha = 1, grouped=F)
  lambda = glm[[lambda_type]]
  coefficients <- coef(glm,s=lambda)[-1,]
  intercept <-  coef(glm,s=lambda)[1,]
  #coefficients[abs(coefficients) < 2] <- 0
  return(list(change_points=which(coefficients != 0),
              fitted_values = X %*% coefficients + intercept))
}

# get change points for each row of the matrix. Each row is considered a time series
get_changepoints_per_row <- function(data_matrix)
{
  change_points <- c()
  for (i in 1:nrow(data_matrix))
  {
    if (sum(data_matrix[i,]) < 2 | length(which(data_matrix[i,] > 0)) < 3 )
      change_points <- rbind(change_points, rep(0, ncol(data_matrix)))
    else
      change_points <- rbind(change_points, t(get_changepoints_vector(data_matrix[i,])$fitted_values))
  }
  rownames(change_points) <- rownames(data_matrix)
  colnames(change_points) <- colnames(data_matrix)
  return(change_points)
}

# Get unified set of change points for all time series in the matrix using the hack:
# For each time point we compute the max difference in thise time point across signatures and fit those
get_changepoints_matrix_jointly <- function(data, lambda_type="lambda.1se")
{
  # fitting all the rows of data jointly to find common change points for all signatures
  # each row represents one time series
  
  X_dim = ncol(data)
  X <- matrix(0, ncol=X_dim, nrow=X_dim)
  for (i in 1:X_dim)
  {
    for (j in 1:i)
      X[i,j] <- 1
  }
  
  # get change values for all signatures
  change_matrix <- data[,1]
  change_matrix <- cbind(change_matrix, data[,-1] - data[,-ncol(data)])
  # get max of changes for each time points
  change_values_per_time_point <- apply(change_matrix, 2, function(x){ x[which.max(abs(x))]})
  # multiply the data by X to be able to fit it the same way as X*w, where w is the vector of change points
  data_to_fit <- X %*% change_values_per_time_point
  
  #   # Another approach: take fitted piese-wise constant values for each sig, get sum of changes of those and fit
  #   data_to_fit <- apply(get_changepoints_per_row(data), 2, sum)
  
  return(get_changepoints_vector(data_to_fit, lambda_type))
}


# Find change points iteratively. To find the next change points, consider every time point as a potential new change point
# For every time point, fit mixture of multinomials and compute log likelihood and bic criteria.
# The next change point is the one which gives lowest bic.
# To correctly evaluate the model using BIC criteria, 
# we pass n_params -- the number of ther params of the model, besides number of checkpoints
find_changepoints_over_all_signatures_one_by_one <- function(vcf, alex.t, n_signatures, parallelized=F)
{
  overall_mixtures <- fit_mixture_of_multinomials_EM(apply(vcf, 1, sum), as.matrix(alex.t))
  overall_mixtures.table <- matrix(rep(overall_mixtures, ncol(vcf)), ncol=ncol(vcf))
 
  changepoints <- c()
  mixtures <- overall_mixtures.table
  bic <- list()
  changepoints_per_iteration <- list()
  mixtures_per_iteration <- list()
  for (t in 1:ncol(vcf))
  {
    print(paste0("Searching for ", t, "-th change point"))
    pcm <- proc.time()
    next_cp <- find_next_change_point_over_all_signatures(vcf, alex.t, changepoints, mixtures, parallelized=parallelized)
    
    if (is.null(next_cp))
    {
       break
    }
    changepoints <- next_cp$changepoints
    mixtures <- next_cp$mixtures
    ll <- next_cp$ll
   
    bic[[t]] <- -2 * ll + (length(changepoints)+1)*(n_signatures + 1) * log(sum(vcf))
    changepoints_per_iteration[[t]] <- changepoints
    mixtures_per_iteration[[t]] <- mixtures
    print(paste("New changepoints", toString(changepoints)))
    print(proc.time() - pcm)
    
    if (t != 1)
      if (bic[[t]] > bic[[t-1]])
        break
  }
  
  optimal <- which.min(unlist(bic))
  optimal_changepoints <- changepoints_per_iteration[[optimal]]
  optimal_mixtures <- fit_mixture_of_multinomials_in_time_slices(vcf, optimal_changepoints, alex.t)
  
  return(list(bics = unlist(bic),
              optimal = optimal,
              changepoints = optimal_changepoints,
              mixtures = optimal_mixtures))
}

# Helper function for find_changepoints_over_all_signatures_one_by_one:
# Find the next change point by fitting every possible time split and evaluating them using BIC.
find_next_change_point_over_all_signatures <- function(vcf, alex.t, changepoints, mixtures, parallelized=F)
{
  if (parallelized)
  {
    likelihood_per_time_split <- foreach (time_point = 2:ncol(vcf)) %dopar%
    {
      if (time_point %in% changepoints)
      {
        list() 
      } else {
        find_next_change_point_over_all_signatures_job(vcf, alex.t, changepoints, mixtures, time_point)
      }
    }
  } else {
    likelihood_per_time_split <- list()
    for (time_point in 2:ncol(vcf))
    {
      if (time_point %in% changepoints) {
        result <- list()
       } else {
       result <- find_next_change_point_over_all_signatures_job(vcf, alex.t, changepoints, mixtures, time_point)
      }
      likelihood_per_time_split[[time_point]] <- result
    }
  }
 
  # Delete entries that correspond to the previous changepoints -- they are empty anyway
  if (length(changepoints) != 0)
  {
    likelihood_per_time_split <- likelihood_per_time_split[-changepoints]
  }
  
  # Get log likelihoods for each proposed time split and see which one performed best
  ll_per_time_split <- get_values_from_list(likelihood_per_time_split, "ll", concat_function="c") 
  
  if (length(ll_per_time_split) < 2)
  {
     result = NULL
  } else {
    best_change_point_index <- which.max(ll_per_time_split[2:length(ll_per_time_split)])+1

    result <- NULL
    if (length(best_change_point_index) != 0)
    {
      result <- likelihood_per_time_split[[best_change_point_index]]
    }
  }
  return(result)
}



find_next_change_point_over_all_signatures_job <- function(vcf, alex.t, changepoints, mixtures, time_point)
{
    new_changepoints <- sort(c(changepoints, time_point))
    current_change_point <- which(new_changepoints == time_point) 
    if (current_change_point > 1) {
      previos_cp <- new_changepoints[current_change_point-1]
    } else {
      previos_cp <- 1
    }
    
    if (current_change_point < length(new_changepoints)) {
      next_cp <- new_changepoints[current_change_point+1]
    } else { 
      next_cp <- ncol(vcf) + 1
    }
 
    # Compute mixtures only around new checkpoint. We don't have to recompute mixtures around time points which did not change.
    mixtures_around_new_checkpoint <- fit_mixture_of_multinomials_in_time_slices(vcf[,c((previos_cp):(next_cp-1))], c(time_point - previos_cp +1 ), alex.t)
    
    if (previos_cp == 1) {
      new_mixtures <- mixtures_around_new_checkpoint
    } else {
      new_mixtures <- mixtures[,1:(previos_cp-1)]
      new_mixtures <- cbind(new_mixtures, mixtures_around_new_checkpoint)
    }
    
    if ((next_cp+1) <= ncol(mixtures)) {
      new_mixtures <- cbind(new_mixtures, mixtures[,(next_cp):ncol(mixtures)])
    }
    
    if (ncol(new_mixtures) != ncol(mixtures))
    {
      print(ncol(new_mixtures))
      print(ncol(mixtures))
      stop("Some computations went wrong")
    }

    # Compute log likelihood for split time_point
    log_prob_per_timepoint <- c()
    for (t in 1:ncol(vcf))
    {
      log_prob_per_timepoint[t] <- log_likelihood_mixture_multinomials(vcf[,t], alex.t, new_mixtures[,t])
    }

    stopifnot(length(log_prob_per_timepoint) == ncol(vcf))
    
    return(list(ll = sum(log_prob_per_timepoint), 
         mixtures = new_mixtures,
         changepoints = new_changepoints))
}



# Find both optimal set of change points and the best set of signatures
# First compute the baseline -- change points with known set of signatures
# Then try to add new signatures -- for each signature find the set of change points and see if it yields better bic
# Then try to remove unnecessary signatures -- for each signature in the set, compute change points and see if it yields better bic
find_changepoints_and_signature_set <- function(vcf, alex.t, prior_signatures = c())
{
  if (length(prior_signatures) == 0)
    prior_signatures <- colnames(alex.t)
  
  sig_set <- prior_signatures
  alex.t.set <- alex.t[,sig_set]
  list[baseline_bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t.set, n_signatures=length(sig_set))
  baseline_min <- baseline_bics[optimal]
  print(paste("Initial baseline:", baseline_min))
  old_changepoints <- c()
  
  while(!(length(changepoints) == length(old_changepoints) &
            sum(old_changepoints == changepoints) == length(changepoints)))
  {
    quit = FALSE
    
    while (!quit)
    {
      # Try to add some signatures to sig_set
      bics_with_new_signature <- foreach (sig = setdiff(colnames(alex.t),sig_set)) %dopar%
      { 
        print(paste("Investigating", sig))
        current_sig_set <- c(sig_set, sig)
        alex.t.set <- alex.t[,current_sig_set]
  
        new_mixtures <- fit_mixture_of_multinomials_in_time_slices(vcf, changepoints, alex.t.set)
        
        # Compute log likelihood per time point if we added signature sig
        log_prob_per_timepoint <- c()
        for (t in 1:ncol(vcf))
        {
          log_prob_per_timepoint[t] <- log_likelihood_mixture_multinomials(vcf[,t], alex.t.set, new_mixtures[,t])
        }
        ll <- sum(log_prob_per_timepoint)
  
        bic <- -2 * ll + ((length(changepoints)+1) * (length(current_sig_set)-1)) * log(sum(vcf))
#         list[bic, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t.set, n_params=length(current_sig_set), parallelized=F)
#         bic <- bic[optimal]
        setNames(bic, sig)
      }
      
      bics_with_new_signature <- unlist(bics_with_new_signature)
      print(paste("Potential candidates:", toString(bics_with_new_signature)))
      new_signature <- names(bics_with_new_signature[which.min(bics_with_new_signature)])
      new_signature.min <- min(bics_with_new_signature)
      if (new_signature.min < baseline_min)
      {
        print(paste("Adding", new_signature))
        sig_set <- c(sig_set, new_signature)
        baseline_min = new_signature.min
        print(paste("New baseline:", baseline_min))
      } else {
        print("Not adding any signatures")
        quit=TRUE
      }
    }
    
    print(paste("Signatures after adding step:", toString(sig_set)))
          
    # Try to remove some signatures from sig_set
    quit = FALSE
    
    while (!quit)
    {
      # Try to add some signatures to sig_set
      bics_with_removed_signature <- foreach (sig = sig_set) %dopar%
      {
        print(paste("Investigating", sig))
        current_sig_set <- setdiff(sig_set, sig)
        alex.t.set <- alex.t[,current_sig_set]
        
        new_mixtures <- fit_mixture_of_multinomials_in_time_slices(vcf, changepoints, alex.t.set)
        
        # Compute log likelihood per time point if we remove signature sig
        log_prob_per_timepoint <- c()
        for (t in 1:ncol(vcf))
        {
          log_prob_per_timepoint[t] <- log_likelihood_mixture_multinomials(vcf[,t], alex.t.set, new_mixtures[,t])
        }
        ll <- sum(log_prob_per_timepoint)
        
        bic <- -2 * ll + (length(changepoints) + length(current_sig_set)) * log(sum(vcf))

#         list[bic, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t.set, n_params=length(current_sig_set), parallelized=F)
#         bic <- bic[optimal]
        setNames(bic, sig)
      }
      
      bics_with_removed_signature <- unlist(bics_with_removed_signature)
      print(paste("Potential signatures: ", toString(bics_with_removed_signature)))
      removed_signature <- names(bics_with_removed_signature[which.min(bics_with_removed_signature)])
      removed_signature.min <- min(bics_with_removed_signature)
      if (removed_signature.min < baseline_min)
      {
        print(paste("Removing", removed_signature))
        sig_set <- setdiff(sig_set, removed_signature)
        baseline_min = removed_signature.min
        print(paste("New baseline:", baseline_min))
      } else {
        print("Not removing any signatures")
        quit = TRUE
      }
    }
    
    print(paste("Signatures after removing step:", toString(sig_set)))
          
    old_changepoints <- changepoints
    
    alex.t.set <- alex.t[,sig_set]
    list[bics, optimal, changepoints, mixtures] <- find_changepoints_over_all_signatures_one_by_one(vcf, alex.t.set, n_signatures=length(sig_set))
    print(paste("Old baseline: ", baseline_min))
    baseline_min <- bics[optimal]
    print(paste("New baseline: ", baseline_min))
    
    print(paste("Updated changepoints: ", toString(changepoints)))
  }

  return(list(bics, optimal, changepoints, mixtures))
}
