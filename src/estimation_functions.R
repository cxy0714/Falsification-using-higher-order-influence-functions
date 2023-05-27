
################################################################################
# Computes squared error loss
################################################################################

calc_mse <- function(true, pred){
  
  mse <- sum((true-pred)^2) / length(pred)
  
  return(mse)
}

################################################################################
# Computes cross-entropy / binary log loss
################################################################################

calc_logloss <- function(true, pred){
  
  pred <- pmax(pmin(1-1e-15, pred), 1e-15)
  
  logloss <- -sum(true*log(pred) + (1-true)*log(1-pred)) / length(pred)
  
  return(logloss)
}

################################################################################
# Computes the matrix psuedoinverse
################################################################################

compute_psuedoinverse <- function(mat){
  
  mat_svd <- svd(mat)
  
  mat_inv <- eigenMapMatMult(eigenMapMatMult(mat_svd$v, diag(mat_svd$d^(-1))), t(mat_svd$u))
  
  return(mat_inv)
}

################################################################################
# Returns a vector which contains the cross validation `k`-fold
# identifier for each row in `data`
################################################################################

get_CV_ids <- function(data, k){
  
  folds_ids <- sample(cut(1:nrow(data), breaks = k, labels = FALSE),
                      size = nrow(data),
                      replace = FALSE)
  
  return(folds_ids)
  
}

################################################################################
# Returns the AIPW estimate of the risk difference 
# trt: A vector of treatment indicators (0/1)
# outcome: A vector of outcome indicators (0/1)
# prob_trt1: A vector containing the estimated conditional treatment probabilities
# prob_outcome1_trt1: A vector containing the estimated outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt0: A vector containing the estimated outcome probabilities under no treatment (setting trt=0)
################################################################################

estimate_AIPW <- function(trt, outcome, prob_trt1, prob_outcome1_trt1, prob_outcome1_trt0){
  
  Y_1 <- 
    mean(
      ((outcome - prob_outcome1_trt1) * ( trt / prob_trt1 )) + prob_outcome1_trt1
      ) 
  
  Y_0 <-  
    mean(
      ((outcome - prob_outcome1_trt0) * ( (1-trt) / (1-prob_trt1) )) + prob_outcome1_trt0
      ) 
  
  RD <- Y_1 - Y_0
  
  return(
    list(
      Y_1 = Y_1, Y_0 = Y_0, RD = RD
    )
  )
  
}

################################################################################
# Returns an estimate of the variance of the AIPW estimate of the risk difference using a nonparametric bootstrap
# trt: A vector of treatment indicators (0/1)
# outcome: A vector of outcome indicators (0/1)
# prob_trt1: A vector containing the estimated conditional treatment probabilities
# prob_outcome1_trt1: A vector containing the estimated outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt0: A vector containing the estimated outcome probabilities under no treatment (setting trt=0)
# M: number of bootstrap replicates
################################################################################

estimate_var_AIPW_bootstrap <- function(trt, outcome, prob_trt1, prob_outcome1_trt1, prob_outcome1_trt0, M){
  
  n <- length(trt)
  
  W <- rmultinom(M, size = n, prob = rep(1/n, n))
  
  bootstraps <- 
    sapply(1:M, function(x){
      
      Y_1 <- mean(
        W[,x]*( ((outcome - prob_outcome1_trt1) * ( trt / prob_trt1 )) + prob_outcome1_trt1 )
      )
      
      Y_0 <- mean(
        W[,x]*( ((outcome - prob_outcome1_trt0) * ( (1-trt) / (1-prob_trt1) )) + prob_outcome1_trt0 )
      )
      
      RD <- Y_1 - Y_0
      
      c(Y_1 = Y_1, Y_0 = Y_0, RD = RD)
      
    })
  
  return(
    list(
      Y_1 = var(bootstraps[1,]), Y_0 = var(bootstraps[2,]), RD = var(bootstraps[3,])
    )
  )
  
}

################################################################################
# Returns an estimate of the variance of the AIPW estimate of the risk difference using the empirical influence curve-based estimator
# trt: A vector of treatment indicators (0/1)
# outcome: A vector of outcome indicators (0/1)
# prob_trt1: A vector containing the estimated conditional treatment probabilities
# prob_outcome1_trt1: A vector containing the estimated outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt0: A vector containing the estimated outcome probabilities under no treatment (setting trt=0)
################################################################################

estimate_var_AIPW <- function(trt, outcome, prob_trt1, prob_outcome1_trt1, prob_outcome1_trt0){
  
  n <- length(trt)
  
  Y_1 <- mean(
    ((outcome - prob_outcome1_trt1) * ( trt / prob_trt1 )) + prob_outcome1_trt1 
  )
  
  Y_0 <- mean(
    ((outcome - prob_outcome1_trt0) * ( (1-trt) / (1-prob_trt1) )) + prob_outcome1_trt0
  )
  
  RD <- Y_1 - Y_0
  
  var_Y_1_terms <- 
    ((outcome - prob_outcome1_trt1) * ( trt / prob_trt1 )) + prob_outcome1_trt1
  
  var_Y_0_terms <- 
    ((outcome - prob_outcome1_trt0) * ( (1-trt) / (1-prob_trt1) )) + prob_outcome1_trt0
  
  var_Y_1 <- sum((var_Y_1_terms - Y_1)^2) / (n-1)
  
  var_Y_0 <- sum((var_Y_0_terms - Y_0)^2) / (n-1)
  
  var_RD <- sum((var_Y_1_terms - var_Y_0_terms - RD)^2) / (n-1)
  
  return(
    list(
      Y_1 = var_Y_1, Y_0 = var_Y_0, RD = var_RD
    )
  )
  
}

################################################################################
# Returns the bias of a doubly robust estimate of a counterfactual mean under treatment (or no treatment)
# trt: A vector of treatment indicators (0/1)
# prob_trt1_est: A vector containing the estimated conditional treatment probabilities
# prob_trt1_true: A vector containing the known conditional treatment probabilities
# prob_outcome1_trt1_est: A vector containing the estimated outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt1_true: A vector containing the known outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt0_est: A vector containing the estimated outcome probabilities under no treatment (setting trt=0)
# prob_outcome1_trt0_true: A vector containing the known outcome probabilities under no treatment (setting trt=0)
################################################################################

compute_bias <- function(trt,
                         prob_trt1_est, prob_trt1_true, 
                         prob_outcome1_trt1_est, prob_outcome1_trt1_true, prob_outcome1_trt0_est, prob_outcome1_trt0_true){
  
  bias_Y_1 <- mean(-trt*(prob_outcome1_trt1_true - prob_outcome1_trt1_est)*((1/prob_trt1_true) - (1/prob_trt1_est)))
  
  bias_Y_0 <- mean((trt-1)*(prob_outcome1_trt0_true - prob_outcome1_trt0_est)*((1/(1-prob_trt1_true)) - (1/(1-prob_trt1_est))))
  
  return(
    list(
      Y_1 = bias_Y_1, Y_0 = bias_Y_0
    )
  )
}

################################################################################
# Returns the Cauchy-Schwarz bias of a doubly robust estimate of a counterfactual mean under treatment (or no treatment)
# trt: A vector of treatment indicators (0/1)
# prob_trt1_est: A vector containing the estimated conditional treatment probabilities
# prob_trt1_true: A vector containing the known conditional treatment probabilities
# prob_outcome1_trt1_est: A vector containing the estimated outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt1_true: A vector containing the known outcome probabilities under treatment (setting trt=1)
# prob_outcome1_trt0_est: A vector containing the estimated outcome probabilities under no treatment (setting trt=0)
# prob_outcome1_trt0_true: A vector containing the known outcome probabilities under no treatment (setting trt=0)
################################################################################

compute_CS_bias <- function(trt,
                            prob_trt1_est, prob_trt1_true, 
                            prob_outcome1_trt1_est, prob_outcome1_trt1_true, prob_outcome1_trt0_est, prob_outcome1_trt0_true){
  
  cs_bias_Y_1 <- sqrt(mean(-trt*((prob_outcome1_trt1_true - prob_outcome1_trt1_est)^2))*
                        mean(-trt*(((1/prob_trt1_true) - (1/prob_trt1_est))^2)))
  
  cs_bias_Y_0 <- sqrt(mean((trt-1)*((prob_outcome1_trt0_true - prob_outcome1_trt0_est)^2))*
                        mean((trt-1)*(((1/(1-prob_trt1_true)) - (1/(1-prob_trt1_est)))^2)))
  
  return(
    list(
      Y_1 = cs_bias_Y_1, Y_0 = cs_bias_Y_0
    )
    )
}

################################################################################
# Takes as input a dataframe and number of required splits, then returns
# a list containing randomly split dataframes
################################################################################

split_data <- function(data, num_splits){
  
  data_splits <- list()
  
  ids <- sample(cut(1:nrow(data), breaks = num_splits, labels = F),
                nrow(data),
                replace = F)
  
  
  for (split in 1:length(unique(ids))) {
    
    data_splits[[split]] <- data[which(ids == split), ]
    
  }
  
  return(data_splits)
  
}

################################################################################
# Returns the cross-validation loss for combinations of hyperparameter values for GBM
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# nfold: number of cross-validation folds
# tree_depth: a vector of values for boosted tree depth
# shrinkage_factor: a vector of shrinkage values  
# num_trees: a vector of values for number of trees (iterations) 
# num_cores: number of cores to use for parallelization
################################################################################

find_params_boosted_tree_model <- function(covariates_df, label_vector, nfold, tree_depth, 
                                           shrinkage_factor, num_trees, num_cores){
  
  params <- expand.grid(n.trees = num_trees, interaction.depth = tree_depth, shrinkage = shrinkage_factor)
  
  cv_test_logloss <- vector(length=nrow(params))
  
  test_preds <- vector(length=nrow(covariates_df))
  
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, library(gbm))
  
  clusterExport(cl, c('params', 'nfold', 'label_vector', 
                      'covariates_df', 'get_CV_ids', 'calc_logloss', 'calc_mse', 'test_preds'),
                envir=environment())
  
  gbm_out <- 
    pbsapply(1:nrow(params), function(iter){
      
      cv_ids <- get_CV_ids(covariates_df, k=nfold)
      
      for (k in 1:nfold) {
        
        train_indices <- which(cv_ids != k)
        
        test_indices <- which(cv_ids == k)
        
        model <- 
          gbm(formula = Y~., data = data.frame(Y=label_vector[train_indices], covariates_df[train_indices, ,drop=FALSE]), 
              n.trees = params$n.trees[iter], interaction.depth = params$interaction.depth[iter], shrinkage = params$shrinkage[iter],
              distribution = 'bernoulli')
        
        test_preds[which(cv_ids == k)] <- predict(model, type='response', n.trees = model$n.trees,
                                                  newdata = covariates_df[test_indices, ,drop=FALSE])
        
      }
      
      calc_logloss(label_vector, test_preds)
      
    }, cl=cl)
  
  stopCluster(cl)
  
  best_params <- as.list(c(params[resample(which(gbm_out == min(gbm_out)), 1), ]))
  best_params$attempted_params <- params
  best_params$logloss <- gbm_out
  
  return(best_params)
  
}

################################################################################
# Returns the cross-validation loss for combinations of hyperparameter values for random forest
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# nfold: number of cross-validation folds
# num_trees: a vector of values for number of trees 
# num_vars: a vector of values for number of variables randomly sampled at each tree split
# num_cores: number of cores to use for parallelization
################################################################################

find_params_random_forest_model <- function(covariates_df, label_vector, nfold, 
                                            num_trees, num_vars, num_cores){
  
  params <- expand.grid(num.trees = num_trees, mtry = num_vars)
  
  cv_test_logloss <- vector(length=nrow(params))
  
  test_preds <- vector(length=nrow(covariates_df))
  
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, library(ranger))
  
  clusterExport(cl, c('params', 'nfold', 'label_vector', 
                      'covariates_df', 'get_CV_ids', 'calc_logloss', 'test_preds'),
                envir=environment())
  
  rf_out <- 
    pbsapply(1:nrow(params), function(iter){
      
      cv_ids <- get_CV_ids(covariates_df, k=nfold)
      
      for (k in 1:nfold) {
        
        train_indices <- which(cv_ids != k)
        
        test_indices <- which(cv_ids == k)
        
        model <- 
          ranger(formula = Y~., data = data.frame(Y=as.factor(label_vector[train_indices]), covariates_df[train_indices, ,drop=FALSE]), 
                 num.trees = params$num.trees[iter], mtry = params$mtry[iter], probability = T, oob.error = F)
        
        test_preds[which(cv_ids == k)] <- predict(model, type='response',
                                                  data = covariates_df[test_indices, ,drop=FALSE])$predictions[,2]
        
      }
      
      calc_logloss(label_vector, test_preds)
      
    }, cl=cl)
  
  stopCluster(cl)
  
  best_params <- as.list(c(params[resample(which(rf_out == min(rf_out)), 1), ]))
  best_params$attempted_params <- params
  best_params$logloss <- rf_out
  
  return(best_params)
  
}

################################################################################
# Returns the cross-validation loss for different values of k for k nearest neighbors
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# k: a vector of values for k
# num_cores: number of cores to use for parallelization
################################################################################

find_params_knn <- function(covariates_df, label_vector, nfold, 
                            k, num_cores){
  
  params <- expand.grid(k = k)
  
  cv_test_logloss <- vector(length=nrow(params))
  
  test_preds <- vector(length=nrow(covariates_df))
  
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, library(caret))
  
  clusterExport(cl, c('params', 'nfold', 'label_vector', 
                      'covariates_df', 'get_CV_ids', 'calc_logloss', 'test_preds'),
                envir=environment())
  
  knn_out <- 
    pbsapply(1:nrow(params), function(iter){
      
      cv_ids <- get_CV_ids(covariates_df, k=nfold)
      
      for (k in 1:nfold) {
        
        train_indices <- which(cv_ids != k)
        
        test_indices <- which(cv_ids == k)
        
        test_preds[which(cv_ids == k)] <- attr(knn3Train(cl = label_vector[train_indices], 
                                                         train = covariates_df[train_indices, ,drop=FALSE],
                                                         test = covariates_df[test_indices, ,drop=FALSE], prob=T, k=params$k[iter]), 'prob')[,2]
        
      }
      
      calc_logloss(label_vector, test_preds)
      
    }, cl=cl)
  
  stopCluster(cl)
  
  best_params <- as.list(c(k=params[resample(which(knn_out == min(knn_out)), 1), ]))
  best_params$attempted_params <- params
  best_params$logloss <- knn_out
  
  return(best_params)
  
}

################################################################################
# Returns the cross-validation loss for different values of lambda and the model formula for LASSO
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# binary_vars: vector containing the names of the binary variables to be used in model fitting
# continuous_vars: vector containing the names of the continuous variables to be used in model fitting
# nfold: number of cross-validation folds
# continuous_var_spline_knots: number of cubic spline knots for continuous variables to be used in model fitting
# degree_of_interactions: degree of interaction terms between all variables to include in the model fitting (e.g., 0 includes no interactions, 1 includes all possible first degree interactions between model terms, etc)
################################################################################

find_params_lasso <- function(covariates_df, label_vector, binary_vars = NULL, continuous_vars = NULL, nfold,
                              continuous_var_spline_knots = 3, degree_of_interactions = 1){
  
  params <- list()
  
  continuous_vars_splined <- switch(
    is.null(continuous_vars) + 1,
    paste0('ns(', continuous_vars, ',', continuous_var_spline_knots - 1, ')'),
    NULL
  )
  
  formula_vars <- c(
    switch(is.null(binary_vars) + 1, binary_vars, NULL),
    switch(is.null(continuous_vars_splined) + 1, continuous_vars_splined, NULL)
  )
  
  if (degree_of_interactions > 1) {
    
    formula_RHS <- paste('(',paste(formula_vars, collapse = " + "), ')^', degree_of_interactions)
    
  } else {
    
    formula_RHS <- paste(formula_vars, collapse = " + ")
    
  }
  
  params$model_formula <- as.formula(paste('~', formula_RHS))
    
  X_matrix <- model.matrix(params$model_formula, covariates_df)[,-1]
  
  params$lambda <- cv.glmnet(y = as.factor(label_vector), x = X_matrix, family="binomial", alpha=1, intercept=T, nfolds=nfold)$lambda.min
  
  return(params)
  
}

################################################################################
# Returns the model formula for a GLM model
# binary_vars: vector containing the names of the binary variables to be used in model fitting
# continuous_vars: vector containing the names of the continuous variables to be used in model fitting
# continuous_var_spline_knots: number of cubic spline knots for continuous variables to be used in model fitting
# degree_of_interactions: degree of interaction terms between all variables to include in the model fitting (e.g., 0 includes no interactions, 1 includes all possible first degree interactions between model terms, etc)
################################################################################

find_params_glm <- function(binary_vars = NULL, continuous_vars = NULL, 
                            continuous_var_spline_knots = 3, degree_of_interactions = 1){
  
  params <- list()
  
  continuous_vars_splined <- switch(
    is.null(continuous_vars) + 1,
    paste0('ns(', continuous_vars, ',', continuous_var_spline_knots - 1, ')'),
    NULL
  )
  
  formula_vars <- c(
    switch(is.null(binary_vars) + 1, binary_vars, NULL),
    switch(is.null(continuous_vars_splined) + 1, continuous_vars_splined, NULL)
  )
  
  if (degree_of_interactions > 1) {
    
    formula_RHS <- paste('(',paste(formula_vars, collapse = " + "), ')^', degree_of_interactions)
    
  } else {
    
    formula_RHS <- paste(formula_vars, collapse = " + ")
    
  }
  
  params$model_formula <- as.formula(paste('label', paste('~', formula_RHS)))
  
  return(params)
  
}

################################################################################
# Returns an ensemble model combining the predictions from GBM, random forest, KNN, LASSO, and GLM
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# params_boosted_tree: a list returned by the find_params_boosted_tree_model() function
# params_random_forest: a list returned by the find_params_random_forest_model() function
# params_knn: a list returned by the find_params_knn() function
# params_lasso: a list returned by the find_params_lasso() function
# params_glm: a list returned by the find_params_glm() function
# num_spline_knots: number of cubic spline knots for continuous predictions from the individual models to be used in fitting the ensemble
# alpha: the ensemble model is fit using glmnet(). Alpha is the elastic net mixing parameter between 0 and 1. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty
# lambda: the ensemble model is fit using glmnet(). Lambda is the elastic net hyperparameter. If lambda=NULL, cross validation with default parameters will be used to choose lambda. Set lambda=0 for GLM.
################################################################################

fit_stacked_classifer_model <- function(covariates_df, label_vector,
                                        params_boosted_tree, params_random_forest, params_knn, params_lasso, params_glm,
                                        num_spline_knots, alpha, lambda=NULL){
  
  test_preds_ids <- get_CV_ids(covariates_df, k=10)
  
  input_preds <- data.frame(boosted_tree = vector(length=nrow(covariates_df)),
                            random_forest = vector(length=nrow(covariates_df)),
                            knn = vector(length=nrow(covariates_df)),
                            lasso = vector(length=nrow(covariates_df)),
                            glm = vector(length=nrow(covariates_df)))
  
  X_matrix_lasso <- model.matrix(params_lasso$model_formula, covariates_df)[,-1]
  
  for (k in 1:10) {
    
    train_indices <- which(test_preds_ids != k)
    
    test_indices <- which(test_preds_ids == k)
    
    boosted_tree_model <- 
      gbm(
        formula = Y~., 
        data = data.frame(Y=label_vector[train_indices], covariates_df[train_indices, ,drop=FALSE]),
        interaction.depth = params_boosted_tree$interaction.depth,
        shrinkage = params_boosted_tree$shrinkage,
        n.trees = params_boosted_tree$n.trees,
        distribution = "bernoulli"
      )
    
    random_forest_model <- 
      ranger(formula = Y~., data = data.frame(Y=as.factor(label_vector[train_indices]), covariates_df[train_indices, ,drop=FALSE]), 
             num.trees = params_random_forest$num.trees, mtry = params_random_forest$mtry, probability = T, oob.error = F)
    
    lasso_model <- 
      glmnet(y = as.factor(label_vector[train_indices]), x = X_matrix_lasso[train_indices, ,drop=FALSE], 
             family="binomial", alpha=1, intercept=T, lambda = params_lasso$lambda)
      
    
    glm_model <- 
      glm(formula = params_glm$model_formula, 
          data = data.frame('label' = as.factor(label_vector[train_indices]), covariates_df[train_indices, ,drop=FALSE]), 
          family = 'binomial')
    
    input_preds$boosted_tree[which(test_preds_ids == k)] <- predict(boosted_tree_model, type='response', n.trees=boosted_tree_model$n.trees,
                                                                    newdata = covariates_df[test_indices, ,drop=FALSE])
    
    input_preds$random_forest[which(test_preds_ids == k)] <- predict(random_forest_model, type='response',
                                                                     data = covariates_df[test_indices, ,drop=FALSE])$predictions[,2]
    
    input_preds$knn[which(test_preds_ids == k)] <- attr(knn3Train(cl = label_vector[train_indices], 
                                                                  train = covariates_df[train_indices, ,drop=FALSE],
                                                                  test = covariates_df[test_indices, ,drop=FALSE], prob=T, k=params_knn$k), 'prob')[,2]
    
    input_preds$lasso[which(test_preds_ids == k)] <- predict(lasso_model, type='response',
                                                             newx = X_matrix_lasso[test_indices, ,drop=FALSE]) %>% as.numeric()
    
    input_preds$glm[which(test_preds_ids == k)] <- predict(glm_model, type='response',
                                                           newdata = covariates_df[test_indices, ,drop=FALSE])
    
  }
  
  input_pred_vars_splined <- paste0('ns(', c('boosted_tree', 'random_forest', 'knn', 'lasso', 'glm'), 
                                    ',', num_spline_knots - 1, ')')
  
  formula_RHS <- paste(input_pred_vars_splined, collapse = " + ")
  
  model_formula <- as.formula(paste(" ~ ", formula_RHS))
  
  input_preds_df <- model.matrix(model_formula, input_preds)
  
  meta_model <- 
    if (is.null(lambda)){
      cv.glmnet(y=label_vector, x=input_preds_df, alpha=alpha, standardize=F, family='binomial')
    } else {
      glmnet(y=label_vector, x=input_preds_df, alpha=alpha, lambda=lambda, standardize=F, family='binomial')
    }
  
  return(list(meta_model = meta_model,
              meta_model_formula = model_formula))
}

################################################################################
# Returns predictions from a ensemble model fit using fit_stacked_classifer_model()
# covariates_df: the dataset containing 'predictor' variables for fitting the gbm model
# label_vector: the response/label variable for fitting the gbm model
# params_boosted_tree: a list returned by the find_params_boosted_tree_model() function
# params_random_forest: a list returned by the find_params_random_forest_model() function
# params_knn: a list returned by the find_params_knn() function
# params_lasso: a list returned by the find_params_lasso() function
# params_glm: a list returned by the find_params_glm() function
# meta_model: a glmnet model. Can be the meta_model object in the list returned by fit_stacked_classifer_model()
# meta_model_formula: a model formula. Can be the meta_model_formula object in the list returned by fit_stacked_classifer_model()
# predict_data: a new dataset for which predictions are desired
################################################################################

estimate_prob_stacked_classifier <- function(covariates_df, label_vector,
                                             params_boosted_tree, params_random_forest, params_knn, params_lasso, params_glm,
                                             meta_model, meta_model_formula, predict_data){
  
  X_matrix_lasso_train <- model.matrix(params_lasso$model_formula, covariates_df)[,-1]
  X_matrix_lasso_est <- model.matrix(params_lasso$model_formula, predict_data)[,-1]
  
  boosted_tree_model <- 
    gbm(
      formula = Y~., 
      data = data.frame(Y=label_vector, covariates_df),
      interaction.depth = params_boosted_tree$interaction.depth,
      shrinkage = params_boosted_tree$shrinkage,
      n.trees = params_boosted_tree$n.trees,
      distribution = "bernoulli"
    )
  
  random_forest_model <- 
    ranger(formula = Y~., data = data.frame(Y=as.factor(label_vector), covariates_df), 
           num.trees = params_random_forest$num.trees, mtry = params_random_forest$mtry, probability = T, oob.error = F)
  
  lasso_model <- 
    glmnet(y = as.factor(label_vector), x = X_matrix_lasso_train, 
           family="binomial", alpha=1, intercept=T, lambda = params_lasso$lambda)
  
  glm_model <- 
    glm(formula = params_glm$model_formula, 
        data = data.frame('label' = as.factor(label_vector), covariates_df), 
        family = 'binomial') 
  
  input_preds <- data.frame(boosted_tree = 
                              predict(boosted_tree_model, newdata=predict_data, type="response", n.trees=boosted_tree_model$n.trees),
                            random_forest = 
                              predict(random_forest_model, data=predict_data, type='response')$predictions[,2],
                            knn = attr(knn3Train(cl = label_vector, train = covariates_df, test = predict_data, 
                                                 prob=T, k=params_knn$k), 'prob')[,2],
                            lasso = predict(lasso_model, newx=X_matrix_lasso_est, type='response') %>% as.numeric(),
                            glm = predict(glm_model, newdata=predict_data, type='response'))
  
  final_preds <- 
    predict(meta_model, newx=model.matrix(meta_model_formula, input_preds), type='response', s='lambda.min')
  
  return(final_preds)
  
}

################################################################################
# Returns predictions from a gbm, random forest, knn, lasso, or glm model
# covariates_df: the dataset containing 'predictor' variables for fitting the gbm model
# label_vector: the response/label variable for fitting the gbm model
# params: a list returned by one of the following functions, find_params_boosted_tree_model(), find_params_random_forest_model(), find_params_knn(), find_params_lasso(), or find_params_glm()
# model: the model to return predictions from. One of 'boosted_tree', 'random_forest', 'knn', 'lasso', or 'glm'
# predict_data: a new dataset for which predictions are desired
################################################################################

estimate_prob_individual_model <- function(covariates_df, label_vector,
                                          params, model, predict_data){
  
  if (model == "boosted_tree") {
    boosted_tree_model <- 
      gbm(
        formula = Y~., 
        data = data.frame(Y=label_vector, covariates_df),
        interaction.depth = params$interaction.depth,
        shrinkage = params$shrinkage,
        n.trees = params$n.trees,
        distribution = 'bernoulli'
      )
    
    out <- predict(boosted_tree_model, newdata=predict_data, type="response", n.trees=boosted_tree_model$n.trees)
  }
  
  if (model == "random_forest") {
    random_forest_model <- 
      ranger(formula = Y~., data = data.frame(Y=as.factor(label_vector), covariates_df), 
             num.trees = params$num.trees, mtry = params$mtry, probability = T, oob.error = F)
    
    out <- predict(random_forest_model, data=predict_data, type='response')$predictions[,2]
  }
  
  if (model == "knn") {
    out <- attr(knn3Train(cl = label_vector, train = covariates_df, test = predict_data, 
                          prob=T, k=params$k), 'prob')[,2]
  }
  
  if (model == "lasso") {
    X_matrix_lasso <- model.matrix(params$model_formula, covariates_df)[,-1]
    
    lasso_model <- 
      glmnet(y = as.factor(label_vector), x = X_matrix_lasso, 
             family="binomial", alpha=1, intercept=T, lambda = params$lambda)
    
    new_X_matrix_lasso <- model.matrix(params$model_formula, predict_data)[,-1]
    
    out <- predict(lasso_model, newx=new_X_matrix_lasso, type='response') %>% as.numeric()
  }
  
  if (model == "glm") {
    glm_model <- glm(formula = params$model_formula,
                     data = data.frame('label' = as.factor(label_vector), covariates_df), 
                     family = 'binomial')
    
    out <- predict(glm_model, predict_data, type='response')
  }
  
  return(out)
  
}

################################################################################
# Returns a basis matrix for binary variables
# data: a dataframe containing the binary variables
# binary_vars: a vector containing the names of binary variables
# degree_of_interactions: degree of interaction terms between all binary variables (e.g., 0 includes no interactions, 1 includes all possible first degree interactions between terms, etc)
# intercept: whether an intercept should be included in the basis (TRUE/FALSE)
################################################################################

create_binary_var_basis <- function(data, binary_vars = NULL, degree_of_interactions = 2, intercept = F) {
  
  formula_vars <- binary_vars
  
  if (degree_of_interactions > 1) {
    
    formula_RHS <- paste('(',paste(formula_vars, collapse = " + "), ')^', degree_of_interactions)
    
  } else {
    
    formula_RHS <- paste(formula_vars, collapse = " + ")
    
  }
  
  model_formula <- as.formula(paste(" ~ ", formula_RHS))
  
  binary_var_basis <- model.matrix(model_formula, data)
  
  colnames(binary_var_basis) <- NULL
  
  return(binary_var_basis)
  
}

################################################################################
# Returns a basis matrix including B-spline transformations of continuous variables
# data: a dataframe containing continuous and binary variables
# binary_vars: a vector containing the names of binary variables (NULL if no binary variables are to be include in the basis)
# continuous_vars: a vector containing the names of continuous variables
# knots: a list containing, for each continuous variable in continuous_vars, the breakpoints that define the spline
# boundary_knots: a list containing, for each continuous variable in continuous_vars, the points that define the boundaries of the spline
# continuous_var_k: a vector containing, for each continuous variable in continuous_vars, the number of columns resulting from the basis transformation applied to that individual continuous variable. If less than 'degree_of_interactions' then will be set to 'degree_of_interactions'
# degree_of_interactions: degree of interaction terms between all variables, including transformations (e.g., 0 includes no interactions, 1 includes all possible first degree interactions between terms, etc)
# polynomial_degree: a vector containing, for each continuous variable in continuous_vars, the degree of the piecewise polynomial for B-spline transformation for that variable. Default is 3 for cubic splines
################################################################################

create_b_spline_basis <- function(data, binary_vars = NULL, continuous_vars = NULL, knots = NULL, boundary_knots = NULL,
                                  continuous_var_k = 5, degree_of_interactions = 2, polynomial_degree = 3) {
  
  if(!is.null(knots) & !is.null(boundary_knots)){
    
    knots <- lapply(1:length(knots), function(x) {
      paste(
        knots[[x]],
        collapse = ','
      )
    })
    
    boundary_knots <- lapply(1:length(boundary_knots), function(x) {
      paste(
        boundary_knots[[x]],
        collapse = ','
      )
    })
    
    continuous_vars_splined <- switch(
      is.null(continuous_vars) + 1,
      paste(
        sapply(1:length(continuous_vars), function(x) {
          paste0('bs(', continuous_vars[[x]], 
                 ', degree=', polynomial_degree, ', knots=c(', knots[[x]], '), Boundary.knots=c(', boundary_knots[[x]], '))')
        }), collapse = " + "),
      NULL
    )
    
  } else {
    
    continuous_vars_splined <- switch(
      is.null(continuous_vars) + 1,
      paste0('bs(', continuous_vars, ', degree=', polynomial_degree, ',', continuous_var_k, ')'),
      NULL
    )
    
  }
  
  formula_vars <- c(
    switch(is.null(binary_vars) + 1, binary_vars, NULL),
    switch(is.null(continuous_vars_splined) + 1, continuous_vars_splined, NULL)
  )
  
  if (degree_of_interactions > 1) {
    
    formula_RHS <- paste('(',paste(formula_vars, collapse = " + "), ')^', degree_of_interactions)
    
  } else {
    
    formula_RHS <- paste(formula_vars, collapse = " + ")
    
  }
  
  model_formula <- as.formula(paste(" ~ ", formula_RHS))
  
  b_spline_basis <- suppressWarnings(model.matrix(model_formula, data))
  
  colnames(b_spline_basis) <- NULL
  
  return(b_spline_basis)
  
}

################################################################################
# Returns a basis matrix including Fourier transformations of continuous variables
# data: a dataframe containing continuous and binary variables
# binary_vars: a vector containing the names of binary variables (NULL if no binary variables are to be include in the basis)
# continuous_vars: a vector containing the names of continuous variables
# nbasis: a vector containing, for each continuous variable in continuous_vars, the number of basis functions in the Fourier basis transformation applied to that individual continuous variable. 
# period: a vector containing, for each continuous variable in continuous_vars, the width of the interval over which all sine/cosine basis functions repeat themselves. If NULL, taken as the default as per the fourier() function
# degree_of_interactions: degree of interaction terms between all variables, including transformations (e.g., 0 includes no interactions, 1 includes all possible first degree interactions between terms, etc)
################################################################################

create_fourier_basis <- function(data, binary_vars = NULL, continuous_vars = NULL, nbasis, period = NULL, 
                                 degree_of_interactions = 2){
  
  if(is.null(period)){
    
    continuous_vars_fourier <- switch(
      is.null(continuous_vars) + 1,
      paste0('fourier(', continuous_vars, ', nbasis=', nbasis, ')'),
      NULL
    )
    
  } else {
    
    continuous_vars_fourier <- switch(
      is.null(continuous_vars) + 1,
      paste0('fourier(', continuous_vars, ', nbasis=', nbasis, ', period=', period, ')'),
      NULL
    )
    
  }
  
  formula_vars <- c(
    switch(is.null(binary_vars) + 1, binary_vars, NULL),
    switch(is.null(continuous_vars_fourier) + 1, continuous_vars_fourier, NULL)
  )
  
  if (degree_of_interactions > 1) {
    
    formula_RHS <- paste('(',paste(formula_vars, collapse = " + "), ')^', degree_of_interactions)
    
  } else {
    
    formula_RHS <- paste(formula_vars, collapse = " + ")
    
  }
  
  model_formula <- as.formula(paste(" ~ ", formula_RHS))
  
  fourier_basis <- model.matrix(model_formula, data)
  
  fourier_basis <- fourier_basis[,c(T,!(apply(fourier_basis, 2, var) == 0)[-1])]
  
  colnames(fourier_basis) <- NULL
  
  return(fourier_basis)
  
}

################################################################################
# Computes and returns IF22 for a counterfactual mean under treatment or no treatment
# y_resid: a vector containing the outcome residuals (see compute_resid())
# a_resid: a vector containing the treatment residuals (see compute_resid())
# basis: a matrix containing the basis transformations of the confounders
# sigma: a Gram matrix
################################################################################

compute_HOIF_22 <- function (y_resid, a_resid, basis, sigma) {
  
  n <- length(y_resid)
  
  sigma_svd <- svd(sigma)
  
  sigma_inv_sqrt <- eigenMapMatMult(eigenMapMatMult(sigma_svd$v, diag(sigma_svd$d^(-1/2))), t(sigma_svd$u))
  
  basis_sigma_inv_sqrt <- eigenMapMatMult(basis, sigma_inv_sqrt)

  y_term <- as.numeric(eigenMapMatMult(t(y_resid), basis_sigma_inv_sqrt))
  a_term <- as.numeric(eigenMapMatMult(t(a_resid), basis_sigma_inv_sqrt))
  
  stat <-  ( (sum(y_term * a_term) - 
              sum(y_resid * rowSums(basis_sigma_inv_sqrt^2) * a_resid)) / (n * (n-1)) ) * (-1)
  
  return(stat)
}

################################################################################
# Computes and returns a bootstrap estimate of the variance of IF22 for counterfactual means or risk differences
# y_resid: a vector containing the outcome residuals (see compute_resid())
# a_resid: a vector containing the treatment residuals (see compute_resid())
# y_resid2: a vector containing the outcome residuals (see compute_resid()) for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# a_resid2: a vector containing the treatment residuals (see compute_resid()) for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# basis: a matrix containing the basis transformations of the confounders
# sigma: a Gram matrix
# sigma2: a Gram matrix for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# M: number of bootstrap replicates
# num_cores: number of cores to use for parallelization
################################################################################

compute_var_IF_22 <- function (y_resid, a_resid, y_resid2 = NULL, a_resid2 = NULL, basis, sigma, sigma2 = NULL, 
                               M, num_cores){
  
  n <- length(y_resid)
  
  sigma_svd <- svd(sigma)
  
  sigma_inv_sqrt <- eigenMapMatMult(eigenMapMatMult(sigma_svd$v, diag(sigma_svd$d^(-1/2))), t(sigma_svd$u))
  
  basis_sigma_inv_sqrt <- eigenMapMatMult(basis, sigma_inv_sqrt)
  
  W <- rmultinom(M, size = n, prob = rep(1/n, n))
  
  params <- list(n=n, y_resid = y_resid, a_resid = a_resid, basis_sigma_inv_sqrt = basis_sigma_inv_sqrt,
                 estimand = estimand, M = M)
  
  if (num_cores != 1) {
    cl <- makeCluster(num_cores); cl_stop <- TRUE
  
    clusterEvalQ(cl, library(SMUT))
  
    clusterExport(cl, c('params', 'W', 'compute_var_IF_22_terms'), envir=environment())
  } else {cl <- 1; cl_stop <- FALSE}
  
  IF_22k_m <- pbsapply(1:M, function(x) {compute_var_IF_22_terms(params, W1 = W[,x], W2 = W[,x])}, cl=cl)
  
  IF_22k_m_c <- pbsapply(1:M, function(x) {compute_var_IF_22_terms(params, W1 = W[,x] - 1, W2 = W[,x] - 1)}, cl=cl)
  
  if (!is.null(y_resid2) & !is.null(a_resid2) & !is.null(sigma2)) {
    
    sigma_svd <- svd(sigma2)
    
    sigma_inv_sqrt <- eigenMapMatMult(eigenMapMatMult(sigma_svd$v, diag(sigma_svd$d^(-1/2))), t(sigma_svd$u))
    
    basis_sigma_inv_sqrt <- eigenMapMatMult(basis, sigma_inv_sqrt)
    
    params <- list(n=n, y_resid = y_resid2, a_resid = a_resid2, basis_sigma_inv_sqrt = basis_sigma_inv_sqrt,
                   estimand = estimand, M = M)
    
    if (num_cores != 1) {clusterExport(cl, c('params', 'W', 'compute_var_IF_22_terms'), envir=environment())}
    
    IF_22k_m2 <- pbsapply(1:M, function(x) {compute_var_IF_22_terms(params, W1 = W[,x], W2 = W[,x])}, cl=cl)
    
    IF_22k_m_c2 <- pbsapply(1:M, function(x) {compute_var_IF_22_terms(params, W1 = W[,x] - 1, W2 = W[,x] - 1)}, cl=cl)
    
    IF_22k_m <- IF_22k_m - IF_22k_m2
    IF_22k_m_c2 <- IF_22k_m_c - IF_22k_m_c2
  }
  
  if(cl_stop) {stopCluster(cl)}
  
  out <- ( (1 / (M-1)) * sum((IF_22k_m - mean(IF_22k_m))^2) ) - 
    ( (2 / (M-1)) * sum((IF_22k_m_c - mean(IF_22k_m_c))^2) )
  
  return(out)
  
}

################################################################################
# Helper function for compute_var_IF_22()
################################################################################

compute_var_IF_22_terms <- function (params, W1, W2){ 
  
  with(params, {
    
    y_resid <- W1 * y_resid
    a_resid <- W2 * a_resid
    
    y_term <- as.numeric(eigenMapMatMult(t(y_resid), basis_sigma_inv_sqrt))
    a_term <- as.numeric(eigenMapMatMult(t(a_resid), basis_sigma_inv_sqrt))
    
    stat <-  ( (sum(y_term * a_term) - 
                  sum(y_resid * rowSums(basis_sigma_inv_sqrt^2) * a_resid)) / (n * (n-1)) ) * (-1)
    
    return(stat)
    
  })
}

################################################################################
# Returns outcome or treatment residuals
# trt: vector containing indicators for binary treatment (0/1)
# outcome: vector containing indicators for the binary outcome (0/1)
# pred: the estimated conditional treatment or outcome probabilities
# type: the type of residuals to return ('trt' or 'outcome')
################################################################################

compute_resid <- function(trt, outcome=NULL, pred, type = 'trt') {

  if (type == 'trt') {
    
    out <- 1 - (trt/pred)
                     
  } else { if (type == 'outcome') {
    
    out <- trt * (outcome - pred)
    
  } else {
    
    stop("type must be 'trt' or 'outcome'")
    
  }
  }

  return(out)
    
}

################################################################################
# Returns a Gram matrix
# basis: a matrix containing the basis transformations of the confounders
# trt: vector containing indicators for binary treatment (0/1)
################################################################################

compute_sigma <- function(basis, trt) {
  
  out <- list()
  
  for (i in 1:length(basis)) {
    
    if (nrow(basis[[i]]) > 5000) {
      
      ids <- cut(1:nrow(basis[[i]]), breaks=20, labels=F)
      
      out[[i]] <- Reduce('+',
                         lapply(1:length(unique(ids)), function(id) eigenMapMatMult(t(trt[ids==id] * as.matrix(basis[[i]][ids==id,])), as.matrix(basis[[i]][ids==id,])) / nrow(basis[[i]][ids==id,]))
      ) / length(unique(ids))
      
    } else {
      
      out[[i]] <- eigenMapMatMult(t(trt * as.matrix(basis[[i]])), as.matrix(basis[[i]])) / nrow(basis[[i]])
      
    }
        
  }
  
  return(out)
}

################################################################################
# Computes and returns a sequence of IF22 (generally for different values of k) for a counterfactual mean under treatment or no treatment
# num_cores: number of cores to use for parallelization
# y_resid: a vector containing the outcome residuals (see compute_resid())
# a_resid: a vector containing the treatment residuals (see compute_resid())
# basis: a list containing matrices containing the basis transformations of the confounders (generally for different values of k). Each matrix in 'basis' must have dimensions compatible with the matrix of the same index in 'sigma'
# sigma: a list containing Gram matrices (generally for different values of k). 
################################################################################

compute_HOIF_sequence <- function(num_cores = length(sigma), y_resid, a_resid, basis, sigma) {
  
  if (num_cores != 1) {
    cl <- makeCluster(num_cores); cl_stop <- TRUE
    
    clusterEvalQ(cl, {
      library(SMUT) 
      library(splines)
      library(magrittr)
    })
    
    clusterExport(cl, c('compute_HOIF_22', 
                        'y_resid', 'a_resid', 'basis', 'sigma'),
                  envir=environment())
  } else {cl <- 1; cl_stop <- FALSE}
  
  out <- pbsapply(1:length(sigma), function(x) {
    
    conditional_bias_estimate <- 
      compute_HOIF_22(
        y_resid = y_resid, 
        a_resid = a_resid, 
        basis = basis[[x]], 
        sigma = sigma[[x]]
      )
    
    return(c(ncol(sigma[[x]]), conditional_bias_estimate))
    
  }, cl=cl)
  
  if(cl_stop) {stopCluster(cl)}
  
  colnames(out) <- names(sigma)
  
  return(out)
  
}

################################################################################
# Computes and returns a sequence of bootstrap estimates of the variance of IF22 (generally for different values of k) for counterfactual means or risk differences
# y_resid: a vector containing the outcome residuals (see compute_resid())
# a_resid: a vector containing the treatment residuals (see compute_resid())
# y_resid2: a vector containing the outcome residuals (see compute_resid()) for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# a_resid2: a vector containing the treatment residuals (see compute_resid()) for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# basis: a list containing matrices containing the basis transformations of the confounders (generally for different values of k). Each matrix in 'basis' must have dimensions compatible with the matrix of the same index in 'sigma'
# sigma: a list containing Gram matrices (generally for different values of k). 
# sigma2: a list containing Gram matrices for another level of treatment if the variance of IF22 for the risk difference is to be calculated. Leaving NULL will compute IF22 for the counterfactual mean under the first level of treatment
# M: number of bootstrap replicates
# num_cores: number of cores to use for parallelization
################################################################################

compute_var_HOIF_sequence <- function(y_resid, a_resid, y_resid2 = NULL, a_resid2 = NULL, basis, sigma, sigma2 = NULL, 
                                      M, num_cores) {
  
  out <- pbsapply(1:length(sigma), function(x) {
    
    if (!is.null(y_resid2) & !is.null(a_resid2) & !is.null(sigma2)){
      
      IF22_var <- 
        compute_var_IF_22(
          y_resid = y_resid, 
          y_resid2 = y_resid2,
          a_resid = a_resid, 
          a_resid2 = a_resid2,
          basis = basis[[x]], 
          sigma = sigma[[x]],
          sigma2 = sigma2[[x]],
          M = M,
          num_cores = num_cores
        )
      
    } else {
      
      IF22_var <- 
        compute_var_IF_22(
          y_resid = y_resid, 
          a_resid = a_resid, 
          basis = basis[[x]], 
          sigma = sigma[[x]],
          M = M,
          num_cores = num_cores
        )
      
    }
    
    return(c(ncol(sigma[[x]]), IF22_var))
    
  })
  
  return(out)
  
}