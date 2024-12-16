# clear space
rm(list = ls())
# package
library("SMUT") # for efficient matrix multiplication via eigenMapMatMult()

################################################################################
# Main function using "all" U-statistics definition! i.e. (i_1 \ne i_2 \ne ... \ne i_m)
# Computes and return 4's estimators from 2-th HOIF to 5-th HOIF estimator. Higher order code is coming soon!
# Vector_1: a n-dimensional vector containing the treatment residuals (like (Aa -1) )
# Vector_2: a n-dimensional vector containing the outcome residuals (like (y - b),here no A-weighted)
# weight : a n-dimensional vector containing the treatment for Weight-ed Gram matrix ( like A or (1-A) )
# basis: a n*p matrix containing the basis transformations of the confounder.
# Split: a logistic variable, 
#        Split == 1 means using sample split to compute the Omega_hat (eHOIF),
#        Split == 0 means using whole sample to compute the Omega_hat (sHOIF).
################################################################################
compute_HOIF_general_all_U <- function(Vector_1, Vector_2, weight, basis, Split) {
  n <- length(Vector_1)
  
  if ( Split == 0) {
    basis_weight <- basis * weight
    t_basis_weight <- t(basis_weight)
    Mat <- eigenMapMatMult(t_basis_weight, basis)
    L <- chol(Mat)
    Omega_mat <- chol2inv(L) * n
    
    
    Ker <- eigenMapMatMult(eigenMapMatMult(basis, Omega_mat), t_basis_weight)
    
    U_list <- list()
    results <- calculate_u_statistics(Vector_1 = Vector_1, Vector_2 = Vector_2, A1 = Ker, A2 = Ker, A3 = Ker, A4 = Ker)
    
    for (i in 2:(5)) {
      U_list[[paste0("U_", i)]] <- (-1)^i * results[[i-1]]
    }
    
    
    BD_list <- Buildingblock_sum_HOIF(U_list, order = 5)
    
    return_list <- c(BD_list, U_list)
    return_list <- combine_results(sHOIF = return_list)
  } else if ( Split == 1) {
    
    
    # Split indices into two halves: a and b
    indices <- seq_len(n)
    set.seed(123)  # For reproducibility
    indices_a <- sample(indices, round(n / 2))
    indices_b <- setdiff(indices, indices_a)
    
    # Split basis, weight, Vector_1, Vector_2
    basis_a <- basis[indices_a, ]
    basis_b <- basis[indices_b, ]
    weight_a <- weight[indices_a]
    weight_b <- weight[indices_b]
    Vector_1_b <- Vector_1[indices_b]
    Vector_2_b <- Vector_2[indices_b]
    
    # Compute for group a
    basis_weight_a <- basis_a * weight_a
    t_basis_weight_a <- t(basis_weight_a)
    Mat_a <- eigenMapMatMult(t_basis_weight_a, basis_a)
    L_a <- chol(Mat_a)
    Omega_mat_a <- chol2inv(L_a) * (n / 2)
    
    # Compute Ker using Omega_mat_a and basis_b
    Ker_ab <- eigenMapMatMult(eigenMapMatMult(basis_b, Omega_mat_a), t(basis_weight_a))
    
    # Compute results using Ker_ab and Vector_1_b, Vector_2_b
    results_ab <- calculate_u_statistics(Vector_1 = Vector_1_b, Vector_2 = Vector_2_b, 
                                         A1 = Ker_ab, A2 = Ker_ab, A3 = Ker_ab, A4 = Ker_ab)
    
    U_list_ab <- list()
    for (i in 2:5) {
      U_list_ab[[paste0("U_", i)]] <- (-1)^i * results_ab[[i - 1]]
    }
    
    BD_list_ab <- Buildingblock_sum_HOIF(U_list_ab, order = 5)
    return_list_ab <- c(BD_list_ab, U_list_ab)
    
    # Combine results for SH-OIF (split a -> b)
    return_list_ab <- combine_results(eHOIF = return_list_ab)
    
    # Repeat the process, swapping a and b
    basis_weight_b <- basis_b * weight_b
    t_basis_weight_b <- t(basis_weight_b)
    Mat_b <- eigenMapMatMult(t_basis_weight_b, basis_b)
    L_b <- chol(Mat_b)
    Omega_mat_b <- chol2inv(L_b) * (n / 2)
    
    Ker_ba <- eigenMapMatMult(eigenMapMatMult(basis_a, Omega_mat_b), t(basis_weight_b))
    
    results_ba <- calculate_u_statistics(Vector_1 = Vector_1[indices_a], Vector_2 = Vector_2[indices_a], 
                                         A1 = Ker_ba, A2 = Ker_ba, A3 = Ker_ba, A4 = Ker_ba)
    
    U_list_ba <- list()
    for (i in 2:5) {
      U_list_ba[[paste0("U_", i)]] <- (-1)^i * results_ba[[i - 1]]
    }
    
    BD_list_ba <- Buildingblock_sum_HOIF(U_list_ba, order = 5)
    return_list_ba <- c(BD_list_ba, U_list_ba)
    
    return_list_ba <- combine_results(eHOIF = return_list_ba)
    
    # Average the two return lists
    return_list <- add_lists(return_list_ab, return_list_ba) 
    return_list <- lapply(return_list, function(x) x / 2)
    
  }
 
  return(return_list)
}


################################################################################
# Main function using "part" U-statistics definition i.e. (i_1 < i_2 < ... < i_m) 
# Computes and return any order>2 HOIF estimator,(order-1)'s estimators from 2-th HOIF to order-th HOIF estimator .
# Vector_1: a n-dimensional vector containing the treatment residuals (like (Aa -1) )
# Vector_2: a n-dimensional vector containing the outcome residuals (like (y - b) ,but not (A(y-b))! Here no A-weighted)
# weight : a n-dimensional vector containing the treatment for Weight-ed Gram matrix ( like A or (1-A) )
# basis: a n*p matrix containing the basis transformations of the confounder.
# order: number of estimators.
# Split: a logistic variable, 
#        Split == 1 means using sample split to compute the Omega_hat (eHOIF),
#        Split == 0 means using whole sample to compute the Omega_hat (sHOIF).
################################################################################
compute_HOIF_general_part_U <- function(Vector_1, Vector_2, weight, basis, order, split) {
  n <- length(Vector_1)
  
  if (split == 0) {
    # Original implementation for split == 0
    basis_weight <- basis * weight
    t_basis_weight <- t(basis_weight)
    Mat <- eigenMapMatMult(t_basis_weight, basis)
    L <- chol(Mat)
    Omega_mat <- chol2inv(L) * n
    
    Ker <- eigenMapMatMult(eigenMapMatMult(basis, Omega_mat), t_basis_weight)
    Ker_up <- Ker * upper.tri(Ker)
    t_vector_1 <- t(Vector_1)
    C <- diag(n)
    
    U_list <- list()
    
    for (i in 2:order) {
      C <- eigenMapMatMult(C, Ker_up)
      U_list[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1, C), Vector_2)) / choose(n, i)
    }
    
    BD_list <- Buildingblock_sum_HOIF(U_list, order)
    return_list <- c(BD_list, U_list)
    return_list <- combine_results(sHOIF = return_list)
    
  } else if (split == 1) {
    # Sample-split implementation for split == 1
    # Split indices into two halves: a and b
    idx_a <- sample(1:n, size = round(n / 2), replace = FALSE)
    idx_b <- setdiff(1:n, idx_a)
    
    # Split basis, weight, Vector_1, and Vector_2
    basis_a <- basis[idx_a, , drop = FALSE]
    basis_b <- basis[idx_b, , drop = FALSE]
    weight_a <- weight[idx_a]
    weight_b <- weight[idx_b]
    Vector_1_a <- Vector_1[idx_a]
    Vector_2_a <- Vector_2[idx_a]
    Vector_1_b <- Vector_1[idx_b]
    Vector_2_b <- Vector_2[idx_b]
    
    # Compute Mat and Omega_mat using subset a
    basis_weight_a <- basis_a * weight_a
    t_basis_weight_a <- t(basis_weight_a)
    Mat_a <- eigenMapMatMult(t_basis_weight_a, basis_a)
    L_a <- chol(Mat_a)
    Omega_mat_a <- chol2inv(L_a) * (n / 2)
    
    # Compute Ker using basis_b and Omega_mat_a
    t_basis_weight_b <- t(basis_b * weight_b)
    Ker_ab <- eigenMapMatMult(eigenMapMatMult(basis_b, Omega_mat_a), t_basis_weight_b)
    Ker_up_ab <- Ker_ab * upper.tri(Ker_ab)
    
    # Calculate U_list for subset b
    t_vector_1_b <- t(Vector_1_b)
    C_ab <- diag(length(idx_b))
    U_list_ab <- list()
    
    for (i in 2:order) {
      C_ab <- eigenMapMatMult(C_ab, Ker_up_ab)
      U_list_ab[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1_b, C_ab), Vector_2_b)) / choose((n / 2), i)
    }
    
    # Compute Mat and Omega_mat using subset b
    basis_weight_b <- basis_b * weight_b
    t_basis_weight_b <- t(basis_weight_b)
    Mat_b <- eigenMapMatMult(t_basis_weight_b, basis_b)
    L_b <- chol(Mat_b)
    Omega_mat_b <- chol2inv(L_b) * (n / 2)
    
    # Compute Ker using basis_a and Omega_mat_b
    t_basis_weight_a <- t(basis_a * weight_a)
    Ker_ba <- eigenMapMatMult(eigenMapMatMult(basis_a, Omega_mat_b), t_basis_weight_a)
    Ker_up_ba <- Ker_ba * upper.tri(Ker_ba)
    
    # Calculate U_list for subset a
    t_vector_1_a <- t(Vector_1_a)
    C_ba <- diag(length(idx_a))
    U_list_ba <- list()
    
    for (i in 2:order) {
      C_ba <- eigenMapMatMult(C_ba, Ker_up_ba)
      U_list_ba[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1_a, C_ba), Vector_2_a)) / choose((n / 2), i)
    }
    
    # Combine results from both splits and average
    BD_list_ab <- Buildingblock_sum_HOIF(U_list_ab, order)
    BD_list_ba <- Buildingblock_sum_HOIF(U_list_ba, order)
    
    return_list_ab <- c(BD_list_ab, U_list_ab)
    return_list_ab <- combine_results(eHOIF = return_list_ab)
    
    return_list_ba <- c(BD_list_ba, U_list_ba)
    return_list_ba <- combine_results(eHOIF = return_list_ba)
    
    # Average the two return lists
    return_list <- add_lists(return_list_ab, return_list_ba) 
    return_list <- lapply(return_list, function(x) x / 2)
  }
  
  return(return_list)
}


################################################################################
# Helper funciton for compute_HOIF_general_part_U()
# Computes and returns (order)-th HOIF estimator and IF from buildingblock U_list
################################################################################
Buildingblock_sum_HOIF <- function(U_list, order, fplugin = 0) {
  BD_list <- list()

  for (j in 2:(order)) {
    BD_value <- 0

    for (i in 0:(j - 2)) {

      combination <- choose(j - 2, i)

      U_index <- j - i

      BD_value <- BD_value + combination * U_list[[paste0("U_", U_index)]]
    }

 
    BD_list[[paste0("IF_", j)]] <- BD_value
  }


  sum_list <- list()


  for (i in 2:order) {
    sum_list[[paste0("HOIF", i)]] <- sum(unlist(BD_list[1:(i - 1)])) + fplugin
  }
  return_list <- c(sum_list, BD_list)
  return(return_list)
}
################################################################################
# Main body function in compute_HOIF_general_all_U()
# Computes and return 4 building_block estimators from 2-th HOIF to 5-th U-statistics.
# \bbU_{n, m}^{\mathsf{all}} \biggl[ f(O_{i_1}, \ldots, O_{i_m}) \biggr] = \frac{1}{\binom{n}{m} \factorial{m} } \sum_{i_1 \ne i_2 \ne \ldots \ne i_m} f(O_{i_1}, \ldots, O_{i_m})
################################################################################
calculate_u_statistics <- function(Vector_1, Vector_2, A1, A2 , A3, A4) {
  # Ensure input matrices are square and of the same size
  n <- nrow(A1)
  if (!(nrow(A2) == n && nrow(A3) == n && nrow(A4) == n)) {
    stop("All input matrices must be square and of the same size")
  }
  # Precompute no_diag matrices to avoid recomputation
  no_diag_A1 <- no_diag(A1)
  no_diag_A2 <- no_diag(A2)
  no_diag_A3 <- no_diag(A3)
  no_diag_A4 <- no_diag(A4)
  
  # Calculate Kernel matrices for different orders
  
  # Kernel matrix for 2nd order
  Ker2 <- A1
  
  # Kernel matrix for 3rd order
  Ker3 <- eigenMapMatMult(no_diag_A1, no_diag_A2)
  
  # Kernel matrix for 4th order
  Ker4 <- eigenMapMatMult(no_diag(Ker3), no_diag_A3) - 
    eigenMapMatMult(A1, diag_col_sum(t(no_diag_A2) * no_diag_A3)) + 
    hadamard(hadamard(A1, t(A2)), A3)
  
  # Kernel matrix for 5th order (more complex)
  Ker5 <- eigenMapMatMult(no_diag(Ker4), no_diag_A4) - 
    eigenMapMatMult(no_diag_A1, diag(diag(eigenMapMatMult(eigenMapMatMult(no_diag_A2, no_diag_A3), no_diag_A4)))) +
    hadamard(hadamard(no_diag_A1, t(no_diag_A2)), eigenMapMatMult(no_diag_A3, no_diag_A4)) +
    hadamard(hadamard(no_diag_A1, no_diag_A4), eigenMapMatMult(t(no_diag_A3), t(no_diag_A2))) -
    eigenMapMatMult(no_diag(eigenMapMatMult(no_diag_A1, no_diag_A2)), diag(diag(eigenMapMatMult(no_diag_A3, no_diag_A4)))) +
    eigenMapMatMult(no_diag_A1, no_diag(hadamard(hadamard(no_diag_A2, t(no_diag_A3)), no_diag_A4))) +
    hadamard(eigenMapMatMult(no_diag_A1, no_diag_A2), hadamard(t(no_diag_A3), no_diag_A4))
  
  # Calculate U-statistics
  # For each order, we use the formula from the document
  Ker2 <- no_diag(Ker2)
  Ker3 <- no_diag(Ker3)
  Ker4 <- no_diag(Ker4)
  Ker5 <- no_diag(Ker5)
  
  U2_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker2), Vector_2) / (n * (n-1))
  U3_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker3), Vector_2) / (n * (n-1) * (n-2))
  U4_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker4), Vector_2) / (n * (n-1) * (n-2) * (n-3))
  U5_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker5), Vector_2) / (n * (n-1) * (n-2) * (n-3) * (n-4))
  
  results_list <- list(
    U2 = U2_all,
    U3 = U3_all,
    U4 = U4_all,
    U5 = U5_all
  )
  # Return a list of U-statistics
  return( results_list)
}
################################################################################
# Helper funciton for compute_HOIF_general_part_U() thus for compute_HOIF_general_part_U()
# subtract, add, combine function for list
################################################################################
# Helper function to create no-diagonal matrix
no_diag <- function(mat) {
  diag(mat) <- 0
  return(mat)
}

# Helper function for Hadamard (point-wise) product
hadamard <- function(mat1, mat2) {
  return(mat1 * mat2)
}

# Helper function for Diag-Column-Sum
diag_col_sum <- function(mat) {
  diag(colSums(mat))
}

subtract_lists <- function(list1, list2) {
  if (length(list1) != length(list2)) {
    stop("Lists must be of the same length")
  }
  mapply(function(x, y) {
    x - y
  }, list1, list2, SIMPLIFY = FALSE)
}
add_lists <- function(list1, list2) {
  if (length(list1) != length(list2)) {
    stop("Lists must be of the same length")
  }
  mapply(function(x, y) {
    x + y
  }, list1, list2, SIMPLIFY = FALSE)
}
combine_results <- function(...) {

  input_lists <- list(...)
  list_names <- names(input_lists)


  if (is.null(list_names) || any(list_names == "")) {
    stop("All input lists must be named")
  }


  processed_lists <- lapply(names(input_lists), function(list_name) {
    current_list <- input_lists[[list_name]]

    if (!is.list(current_list)) {
      current_list <- list(current_list)
    }


    result <- lapply(names(current_list), function(name) {
      current_list[[name]]
    })
    names(result) <- paste0(list_name, "_", names(current_list))

    return(result)
  })


  result <- do.call(c, processed_lists)
  return(result)
}

################## Example ######################
# In this example, we estimate the bias of DML/AIPW for the potential outcomes Y(a = 1) and Y(a = 0).
# The outcome model assumes a linear relationship: Y = A + beta * X.
# The propensity score model follows a logistic regression: A = psi(alpha * X).
#
# - compute_HOIF_general_all_U(): Provides the exact formula for HOIF, but only up to the 5th order.
# - compute_HOIF_general_part_U(): Computes an approximate formula for HOIF, extendable to any order.
#       - in their return results, "_HOIF_" are the used estimator, "_IF_" or "_U_" are just middle terms.
# Note: The basis function of X is the identity in this example. For practical applications, 
#       it is recommended to apply transformations (e.g., B-splines, Fourier basis, etc.) to X.
################################################

set.seed(123)


n <- 1000
p <- 100

mu_x <- 0
alpha <- rnorm(p, mean = 0, sd = 0.1)
beta <- rnorm(p, mean = 0, sd = 0.1)

X <- matrix(rnorm(n * p, mu_x, sd = 1), nrow = n, ncol = p)

logit_prob <- X %*% alpha
prob_A <- 1 / (1 + exp(-logit_prob))
A <- as.numeric(rbinom(n, size = 1, prob = prob_A))


Y <- as.numeric(A + X %*% beta + rnorm(n, mean = 0, sd = 0.2))


propensity_model <- glm(A ~ X, family = binomial(link = "logit"))
propensity_scores <- predict(propensity_model, type = "response")
propensity_scores <- as.numeric(propensity_scores)


Y_model <- lm(Y ~ A + X)


Y_pred_1 <- predict(Y_model, newdata = data.frame(A = 1, X = X))
Y_pred_0 <- predict(Y_model, newdata = data.frame(A = 0, X = X))


summary(propensity_model)
summary(Y_model)

epsilon_A_1 <- A / propensity_scores - 1
epsilon_A_0 <- (1 - A) / (1 - propensity_scores) - 1
epsilon_Y_1 <- Y - Y_pred_1
epsilon_Y_0 <- Y - Y_pred_0

m <- 6


hoif_1_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  Split = 0
) 
hoif_1_part_shoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m,
  split = 0
)

hoif_1_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  Split = 1
)
hoif_1_part_ehoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m,
  split = 1
)


hoif_0_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  Split = 0
)
hoif_0_part_shoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m,
  split = 0
)

hoif_0_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  Split = 1
)

hoif_0_part_ehoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m,
  split = 1
)



