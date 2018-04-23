# AUTHOR: Yulia Rubanova

fit_signatures_linear_regression <- function (vcf, alex.t, alpha = 0, 
                                              lambda_type=c("lambda.min","lambda.1se", "lambda0"), prior = NULL) {
  if (!is.null(prior))
    if (length(prior) != ncol(alex.t))
      stop(paste0("Length of prior should be equal to ", ncol(alex.t)))
  
  # RIDGE REGRESSION
  # The trinucleotide counts (response) are regressed against Alexandrov frequencies
  # The weights for the 30 mutational signatures are stored in a new matrix (dd)
  # GLMNET (generalized linear model) is used for this regression with the following parameters:
  # stdArray.glmnet: used for 10-fold cross-validation
  # family: "Poisson" is used because of nonnegative counts
  # lower.limits: set to 0 in order to obtain nonnegative weights
  # Ridge penalty is used with alpha = 0
  # lambda.min is used for calculating the coefficients of the regression for minimum cross-validated error
  dd <- matrix(nrow = (ncol(alex.t)+1), ncol = 0)
  for (m in 1:ncol(vcf)) {
    to_fit <- vcf[,m]
    lower.limits = 0
    
    if (!is.null(prior))
    {
      to_fit <- to_fit - as.matrix(alex.t) %*% prior
      to_fit[to_fit < 0] <- 0
      lower.limits = -prior
    }
    
    if (lambda_type == "lambda0")
    {
      glm <- glmnet(as.matrix(alex.t),to_fit, family = "poisson", lower.limits = lower.limits, alpha = alpha, lambda=0)
      lambda = 0
    } else {
      glm <- cv.glmnet(as.matrix(alex.t), to_fit, family = "poisson", lower.limits = lower.limits, alpha = alpha)
      lambda <- glm[[lambda_type]]
    }
    
    dd <- cbind(dd, as.array(coef(glm, s=lambda)))
    mse.min <- glm$cvm[glm$lambda == lambda]
  }
  
  dd <- dd[2:nrow(dd),] # First row that contains intercepts of the regression results is ignored
  #dd <- apply(dd, 2, function(x)(x/sum(x))) # Weights are normalized
  return(dd)
}



make_binary_table <- function(multinomial_vector)
{
  N = sum(multinomial_vector)
  mutation_types <- length(multinomial_vector)
  
  # Make a matrix of samples for fitting mixture of multinomials.
  # Each sample contains the one mutation. 
  # data.matrix is a binary matrix with N columns and mutation_types rows.
  data.matrix <- matrix(0, ncol=N, nrow=mutation_types)
  current_index <- 1
  for (i in 1:length(multinomial_vector))
  {
    if (multinomial_vector[i] != 0)
    {
      data.matrix[i, current_index:min(current_index+multinomial_vector[i]-1, N) ] <- 1
    }
    current_index <- current_index + multinomial_vector[i]
  }
  
  if (current_index != N + 1)
    stop("Something went wrong during binary matrix construction: current_index is off")
  
  return(data.matrix)
}

# Fit mixture of multinomials for each column of the matrix
fit_mixture_of_multinomials_matrix <- function (vcf, alex.t, prior=NULL)
{
  dd <- matrix(nrow = ncol(alex.t), ncol = 0)
  for (m in 1:ncol(vcf)) {
    mixtures <- fit_mixture_of_multinomials_EM(vcf[,m], as.matrix(alex.t))
    dd <- cbind(dd, as.array(mixtures))
  }
  
  rownames(dd) <- colnames(alex.t)
  
  return(dd)
}

# fit mixture of multinomials to the vector
fit_mixture_of_multinomials_EM <- function(multinomial_vector, composing_multinomials, prior=NULL)
{
  # Number of mutations to fit
  N = sum(multinomial_vector)
  # Number of mutation types / categories of mutinomial
  mutation_types <- length(multinomial_vector)
  # Number of multinomials/signatures to fit and to make mixture of
  M <- ncol(composing_multinomials)
  
  if (length(multinomial_vector) != nrow(composing_multinomials))
  {
	print(length(multinomial_vector))
	print(dim(composing_multinomials))
  	stop("Length of data vector is not equal to nrow of matrix to fit. Did you forget to transpose the matrix?")
  }
  data.matrix <- make_binary_table(multinomial_vector)
  
  # data_given_class[i,n] corresponds to class/signature i and sample/mutation n
  data_given_class <- matrix(0, nrow=M, ncol=N)
  for (i in 1:M)
  {
    data_given_class[i,] <- apply(composing_multinomials[,i]^data.matrix,2,prod)
  }

  # Mixtures of multinomials. Use uniform prior unless the prior is specified
  if (!is.null(prior))
  {
    if (length(prior) != M)
      stop(paste0("Length of prior should be equal to ", M))
  } else {
    prior <- rep(1/M, M)
  }
  
  pi <- prior
  pi_diff <- Inf
  iteration <- 1
  
  while (pi_diff > 0.001 & iteration < 1000)
  {
    # E-step: update posterior.     
    p_x <- apply(data_given_class * pi, 2, sum)
    
    # class_given_data[i,n] corresponds to class/signature i and sample/mutation n
    class_given_data <- t(t(data_given_class * pi) / p_x)
    
    # S-step: update mixtures
    pi_new <- 1/N * apply(class_given_data,1,sum)
    
    if (sum(pi_new > 1) != 0) {
      stop("Mixture ratio is greater than 1")
    }
    
    if (sum(pi_new < 0) != 0)
      stop("Mixture ratio is less than 0")
    
    if (sum(pi_new) > 1.5)
      stop("Sum of mixture ratios is greater than 1")
    
    pi_diff <- sum(abs(pi_new - pi))
    pi <- pi_new
    iteration <- iteration + 1
  }
  
  return(pi)
}

# fit mixture of mutinomials in each time slice specified by change_points
fit_mixture_of_multinomials_in_time_slices <- function(data, change_points, alex.t, split_data_at_change_point = T)
{
  fitted_values <- matrix(NA, ncol=ncol(data), nrow=ncol(alex.t))
  rownames(fitted_values) <- colnames(alex.t)
  
  # Get first time slice until the first check point and get sum of it
  if (length(change_points) == 0) {
    end_of_first_slice <- ncol(data)
    slice_indices <- 1:(end_of_first_slice)
  } else {
    end_of_first_slice <- change_points[1]
    slice_indices <- 1:(end_of_first_slice-1)
  }

  d <- list()

  d[[1]] <- list()
  d[[1]]$data <-  toVerticalMatrix(data[,slice_indices])
  d[[1]]$slice_indices <- slice_indices
 
  if (split_data_at_change_point)
  {
    right_side <- left_side <- floor(data[,end_of_first_slice] / 2)
    leftovers <- data[,end_of_first_slice] - right_side - left_side
    leftovers_right <- leftovers
    leftovers_right[ sort(sample(which(leftovers > 0), length(which(leftovers > 0)) / 2))] <- 0
    leftovers_left <- leftovers - leftovers_right
    
    right_side <- right_side + leftovers_right
    left_side <- left_side + leftovers_left
    
    stopifnot(sum(data[,end_of_first_slice]) == sum(right_side) + sum(left_side))
    
    d[[1]]$data <- cbind(d[[1]]$data, left_side)
    left_side <- c()
  }
  
  if (length(change_points) > 1)
  {
    for (i in 2:length(change_points))
    {
      slice_indices <- (change_points[i-1]):(change_points[i]-1)
      d[[i]] <- list()
      d[[i]]$data <- toVerticalMatrix(data[,slice_indices])
      d[[i]]$slice_indices <- slice_indices
      
      if (split_data_at_change_point)
      {
        d[[i]]$data <- toVerticalMatrix(cbind(right_side, d[[i]]$data[,-1]))
        right_side <- c()
        
        right_side <- left_side <- floor(data[,slice_indices[length(slice_indices)]] / 2)
        leftovers <- data[,slice_indices[length(slice_indices)]] - right_side - left_side
        leftovers_right <- leftovers
        leftovers_right[ sort(sample(which(leftovers > 0), length(which(leftovers > 0)) / 2))] <- 0
        leftovers_left <- leftovers - leftovers_right
        
        right_side <- right_side + leftovers_right
        left_side <- left_side + leftovers_left
        
        stopifnot(sum(data[,slice_indices[length(slice_indices)]]) == sum(right_side) + sum(left_side))
        
        d[[i]]$data <- cbind(d[[i]]$data, left_side)
        left_side <- c()
      }
    }
  }
  
  if (length(change_points) != 0)
  {
    slice_indices <- (change_points[length(change_points)]):ncol(data)
    d[[length(d) + 1]] <- list()
    d[[length(d)]]$data <- toVerticalMatrix(data[,slice_indices])
    d[[length(d)]]$slice_indices <- slice_indices
    
    if (split_data_at_change_point)
    {
      d[[length(d)]]$data <- cbind(right_side, toVerticalMatrix(d[[length(d)]]$data)[,-1])
      right_side <- c()
    }
  }
  

  for (i in 1:length(d))
  {
    current_d <- d[[i]]
    current_d.sum <- IgnoreVectorOrMatrix(current_d$data, FUN=function(x) {apply(x,1, sum)})
    fitted_for_time_slice <- fit_mixture_of_multinomials_EM(current_d.sum, alex.t)
    
    fitted_values[,current_d$slice_indices] <- matrix(rep(fitted_for_time_slice, length(current_d$slice_indices)), nrow=nrow(fitted_values))
  }
  colnames(fitted_values) <- colnames(data)
  return(fitted_values)
}


log_likelihood_mixture_multinomials <- function (multinomial_vector, composing_multinomials, mixtures) {
  mutation_probabilities_under_mixture <- compute_mutation_probabilities_under_mixture(
                                            multinomial_vector, composing_multinomials, mixtures)
  return(sum(mutation_probabilities_under_mixture))
}

compute_mutation_probabilities_under_mixture <- function(multinomial_vector, composing_multinomials, mixtures_at_timepoint) {
  mutation_binary_table <-  make_binary_table(multinomial_vector)
  
  # mutation_probabilities_under_multinomial[i,n] corresponds to class/signature i and sample/mutation n
  mutation_probabilities_under_multinomial <- matrix(0, nrow=ncol(composing_multinomials), ncol=ncol(mutation_binary_table))
  for (sig in 1:ncol(composing_multinomials)) {
    mutation_probabilities_under_multinomial[sig,] <- apply(composing_multinomials[,sig]^mutation_binary_table,2,prod)
  }
  
  mutation_probabilities_under_mixture <-  log(t(mutation_probabilities_under_multinomial) %*% as.matrix(mixtures_at_timepoint))
  stopifnot(length(mutation_probabilities_under_mixture) == sum(multinomial_vector))
  
  return(mutation_probabilities_under_mixture)
}
