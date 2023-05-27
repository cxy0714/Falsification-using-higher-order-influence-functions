source('./src/dependencies.R')
source('./src/estimation_functions.R')

#################################################################################################

# Set sample size for estimation sample used to compute the conditional bias

n <- 5000
reps <- 200

#################################################################################################

load("./params/simulation_parameters_continuous_2,8,32,64,128.RData")

with(simulation_parameters_continuous, {
  
  prob_trt1_stacked_classifier <- prob_outcome1_trt1_stacked_classifier <- prob_outcome1_trt0_stacked_classifier <- 
    prob_trt1_boosted_tree <- prob_outcome1_trt1_boosted_tree <- prob_outcome1_trt0_boosted_tree <- 
    prob_trt1_random_forest <- prob_outcome1_trt1_random_forest <- prob_outcome1_trt0_random_forest <- 
    prob_trt1_knn <- prob_outcome1_trt1_knn <- prob_outcome1_trt0_knn <- 
    prob_trt1_lasso <- prob_outcome1_trt1_lasso <- prob_outcome1_trt0_lasso <- 
    prob_trt1_glm <- prob_outcome1_trt1_glm <- prob_outcome1_trt0_glm <- 
    p_X <- p_Y <- X <- Y <- vector()
  
  for(rep in 1:reps) {
  
    Z1 = runif(n, 0, 1); Z2 = runif(n, 0, 1); Z3 = runif(n, 0, 1); Z4 = runif(n, 0, 1); Z5 = runif(n, 0, 1)
    Z6 = runif(n, 0, 1); Z7 = runif(n, 0, 1); Z8 = runif(n, 0, 1); Z9 = runif(n, 0, 1); Z10 = runif(n, 0, 1)
    
    sim_data_est <- create_b_spline_basis(
      data=data.frame(Z1),
      continuous_vars=c('Z1'),
      knots=true_knots,
      boundary_knots=list(c(0,1)),
      polynomial_degree=3,
      degree_of_interactions=1
    )
    
    p_X[((n*(rep-1))+1):(n*(rep))] = sim_data_est %*% theta_X
    p_Y[((n*(rep-1))+1):(n*(rep))] = sim_data_est %*% theta_Y
    X[((n*(rep-1))+1):(n*(rep))] = rbinom(n, 1, p_X[((n*(rep-1))+1):(n*(rep))])
    Y[((n*(rep-1))+1):(n*(rep))] = rbinom(n, 1, p_Y[((n*(rep-1))+1):(n*(rep))])
    
    sim_data_est <- data.frame(Z1=Z1, Z2=Z2, Z3=Z3, Z4=Z4, Z5=Z5, Z6=Z6, Z7=Z7, Z8=Z8, Z9=Z9, Z10=Z10, 
                               X=X[((n*(rep-1))+1):(n*(rep))], 
                               Y=Y[((n*(rep-1))+1):(n*(rep))])
    
    prob_trt1_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                                label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                                params_boosted_tree = params_boosted_tree_trt, 
                                                                                                params_random_forest = params_random_forest_trt, 
                                                                                                params_knn = params_knn_trt,
                                                                                                params_lasso = params_lasso_trt,
                                                                                                params_glm = params_glm_trt,
                                                                                                meta_model = meta_model_trt$meta_model,
                                                                                                meta_model_formula = meta_model_trt$meta_model_formula,
                                                                                                predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                         params_boosted_tree = params_boosted_tree_outcome, 
                                                                                                         params_random_forest = params_random_forest_outcome, 
                                                                                                         params_knn = params_knn_outcome,
                                                                                                         params_lasso = params_lasso_outcome,
                                                                                                         params_glm = params_glm_outcome,
                                                                                                         meta_model = meta_model_outcome$meta_model,
                                                                                                         meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                         params_boosted_tree = params_boosted_tree_outcome, 
                                                                                                         params_random_forest = params_random_forest_outcome, 
                                                                                                         params_knn = params_knn_outcome,
                                                                                                         params_lasso = params_lasso_outcome,
                                                                                                         params_glm = params_glm_outcome,
                                                                                                         meta_model = meta_model_outcome$meta_model,
                                                                                                         meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                        params = params_boosted_tree_trt, 
                                                                                        model = "boosted_tree",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                 label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                 params = params_boosted_tree_outcome, 
                                                                                                 model = "boosted_tree",
                                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                 label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                 params = params_boosted_tree_outcome, 
                                                                                                 model = "boosted_tree",
                                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                         label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                         params = params_random_forest_trt, 
                                                                                         model = "random_forest",
                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                  params = params_random_forest_outcome, 
                                                                                                  model = "random_forest",
                                                                                                  predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                  params = params_random_forest_outcome, 
                                                                                                  model = "random_forest",
                                                                                                  predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                               label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                               params = params_knn_trt, 
                                                                               model = "knn",
                                                                               predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_knn_outcome, 
                                                                                        model = "knn",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_knn_outcome, 
                                                                                        model = "knn",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                 label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                 params = params_lasso_trt, 
                                                                                 model = "lasso",
                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                          params = params_lasso_outcome, 
                                                                                          model = "lasso",
                                                                                          predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                          params = params_lasso_outcome, 
                                                                                          model = "lasso",
                                                                                          predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                               label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                               params = params_glm_trt, 
                                                                               model = "glm",
                                                                               predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_glm_outcome, 
                                                                                        model = "glm",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_glm_outcome, 
                                                                                        model = "glm",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  }
  
  bias_stacked <- compute_bias(prob_trt1_est = prob_trt1_stacked_classifier, prob_trt1_true = p_X,
                               prob_outcome1_trt1_est = prob_outcome1_trt1_stacked_classifier, 
                               prob_outcome1_trt1_true = p_Y, 
                               prob_outcome1_trt0_est = prob_outcome1_trt0_stacked_classifier, 
                               prob_outcome1_trt0_true = p_Y,
                               trt = X)
  
  bias_boosted_tree <- compute_bias(prob_trt1_est = prob_trt1_boosted_tree, prob_trt1_true = p_X, 
                                    prob_outcome1_trt1_est = prob_outcome1_trt1_boosted_tree, 
                                    prob_outcome1_trt1_true = p_Y, 
                                    prob_outcome1_trt0_est = prob_outcome1_trt0_boosted_tree, 
                                    prob_outcome1_trt0_true = p_Y,
                                    trt = X)
  
  bias_random_forest <- compute_bias(prob_trt1_est = prob_trt1_random_forest, prob_trt1_true = p_X,
                                     prob_outcome1_trt1_est = prob_outcome1_trt1_random_forest, 
                                     prob_outcome1_trt1_true = p_Y, 
                                     prob_outcome1_trt0_est = prob_outcome1_trt0_random_forest, 
                                     prob_outcome1_trt0_true = p_Y,
                                     trt = X)
  
  bias_knn <- compute_bias(prob_trt1_est = prob_trt1_knn, prob_trt1_true = p_X, 
                           prob_outcome1_trt1_est = prob_outcome1_trt1_knn, 
                           prob_outcome1_trt1_true = p_Y, 
                           prob_outcome1_trt0_est = prob_outcome1_trt0_knn, 
                           prob_outcome1_trt0_true = p_Y,
                           trt = X)
  
  bias_lasso <- compute_bias(prob_trt1_est = prob_trt1_lasso, prob_trt1_true = p_X, 
                             prob_outcome1_trt1_est = prob_outcome1_trt1_lasso, 
                             prob_outcome1_trt1_true = p_Y, 
                             prob_outcome1_trt0_est = prob_outcome1_trt0_lasso, 
                             prob_outcome1_trt0_true = p_Y,
                             trt = X)
  
  bias_glm <- compute_bias(prob_trt1_est = prob_trt1_glm, prob_trt1_true = p_X,
                           prob_outcome1_trt1_est = prob_outcome1_trt1_glm, 
                           prob_outcome1_trt1_true = p_Y, 
                           prob_outcome1_trt0_est = prob_outcome1_trt0_glm, 
                           prob_outcome1_trt0_true = p_Y,
                           trt = X)
  
  bias_out <- list(bias_stacked = bias_stacked,
                   bias_boosted_tree = bias_boosted_tree,
                   bias_random_forest = bias_random_forest,
                   bias_knn = bias_knn,
                   bias_lasso = bias_lasso,
                   bias_glm = bias_glm)
  
  cs_bias_stacked <- compute_CS_bias(prob_trt1_est = prob_trt1_stacked_classifier, prob_trt1_true = p_X, 
                                     prob_outcome1_trt1_est = prob_outcome1_trt1_stacked_classifier, 
                                     prob_outcome1_trt1_true = p_Y, 
                                     prob_outcome1_trt0_est = prob_outcome1_trt0_stacked_classifier, 
                                     prob_outcome1_trt0_true = p_Y,
                                     trt = X)
  
  cs_bias_boosted_tree <- compute_CS_bias(prob_trt1_est = prob_trt1_boosted_tree, prob_trt1_true = p_X, 
                                          prob_outcome1_trt1_est = prob_outcome1_trt1_boosted_tree, 
                                          prob_outcome1_trt1_true = p_Y, 
                                          prob_outcome1_trt0_est = prob_outcome1_trt0_boosted_tree, 
                                          prob_outcome1_trt0_true = p_Y,
                                          trt = X)
  
  cs_bias_random_forest <- compute_CS_bias(prob_trt1_est = prob_trt1_random_forest, prob_trt1_true = p_X, 
                                           prob_outcome1_trt1_est = prob_outcome1_trt1_random_forest, 
                                           prob_outcome1_trt1_true = p_Y, 
                                           prob_outcome1_trt0_est = prob_outcome1_trt0_random_forest, 
                                           prob_outcome1_trt0_true = p_Y,
                                           trt = X)
  
  cs_bias_knn <- compute_CS_bias(prob_trt1_est = prob_trt1_knn, prob_trt1_true = p_X, 
                                 prob_outcome1_trt1_est = prob_outcome1_trt1_knn, 
                                 prob_outcome1_trt1_true = p_Y, 
                                 prob_outcome1_trt0_est = prob_outcome1_trt0_knn, 
                                 prob_outcome1_trt0_true = p_Y,
                                 trt = X)
  
  cs_bias_lasso <- compute_CS_bias(prob_trt1_est = prob_trt1_lasso, prob_trt1_true = p_X, 
                                   prob_outcome1_trt1_est = prob_outcome1_trt1_lasso, 
                                   prob_outcome1_trt1_true = p_Y, 
                                   prob_outcome1_trt0_est = prob_outcome1_trt0_lasso, 
                                   prob_outcome1_trt0_true = p_Y,
                                   trt = X)
  
  cs_bias_glm <- compute_CS_bias(prob_trt1_est = prob_trt1_glm, prob_trt1_true = p_X, 
                                 prob_outcome1_trt1_est = prob_outcome1_trt1_glm, 
                                 prob_outcome1_trt1_true = p_Y, 
                                 prob_outcome1_trt0_est = prob_outcome1_trt0_glm, 
                                 prob_outcome1_trt0_true = p_Y,
                                 trt = X)
  
  cs_bias_out <- list(cs_bias_stacked = cs_bias_stacked,
                      cs_bias_boosted_tree = cs_bias_boosted_tree,
                      cs_bias_random_forest = cs_bias_random_forest,
                      cs_bias_knn = cs_bias_knn,
                      cs_bias_lasso = cs_bias_lasso,
                      cs_bias_glm = cs_bias_glm)
  
  save(bias_out, 
       file=paste0("./output/bias_continuous.RData"))
  
  save(cs_bias_out, 
       file=paste0("./output/cs_bias_continuous.RData"))
  
})

#################################################################################################

load("./params/simulation_parameters_binary_1,2,3,4,5,6,7,8.RData")

with(simulation_parameters_binary, {
  
  prob_trt1_stacked_classifier <- prob_outcome1_trt1_stacked_classifier <- prob_outcome1_trt0_stacked_classifier <- 
    prob_trt1_boosted_tree <- prob_outcome1_trt1_boosted_tree <- prob_outcome1_trt0_boosted_tree <- 
    prob_trt1_random_forest <- prob_outcome1_trt1_random_forest <- prob_outcome1_trt0_random_forest <- 
    prob_trt1_knn <- prob_outcome1_trt1_knn <- prob_outcome1_trt0_knn <- 
    prob_trt1_lasso <- prob_outcome1_trt1_lasso <- prob_outcome1_trt0_lasso <- 
    prob_trt1_glm <- prob_outcome1_trt1_glm <- prob_outcome1_trt0_glm <- 
    p_X <- p_Y <- X <- Y <- vector()
  
  for(rep in 1:reps) {
  
    sim_data_est <- data.frame(
      Z1 = as.numeric(rbernoulli(n,p=0.5)), Z2 = as.numeric(rbernoulli(n,p=0.5)),
      Z3 = as.numeric(rbernoulli(n,p=0.5)), Z4 = as.numeric(rbernoulli(n,p=0.5)),
      Z5 = as.numeric(rbernoulli(n,p=0.5)), Z6 = as.numeric(rbernoulli(n,p=0.5)),
      Z7 = as.numeric(rbernoulli(n,p=0.5)), Z8 = as.numeric(rbernoulli(n,p=0.5)),
      Z9 = as.numeric(rbernoulli(n,p=0.5)), Z10 = as.numeric(rbernoulli(n,p=0.5))
    ) %>% 
      merge(., key_trt) %>%
      merge(., key_outcome) %>% 
      mutate(
        X = rbinom(n(), 1, val_trt),
        Y = rbinom(n(), 1, val_outcome)
      ) 
    
    p_X[((n*(rep-1))+1):(n*(rep))] <- sim_data_est$val_trt
    p_Y[((n*(rep-1))+1):(n*(rep))] <- sim_data_est$val_outcome
    X[((n*(rep-1))+1):(n*(rep))] <- sim_data_est$X
    Y[((n*(rep-1))+1):(n*(rep))] <- sim_data_est$Y
    
    sim_data_est <- sim_data_est %>% dplyr::select(-c(val_trt, val_outcome))
    
    prob_trt1_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                                label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                                params_boosted_tree = params_boosted_tree_trt, 
                                                                                                params_random_forest = params_random_forest_trt, 
                                                                                                params_knn = params_knn_trt,
                                                                                                params_lasso = params_lasso_trt,
                                                                                                params_glm = params_glm_trt,
                                                                                                meta_model = meta_model_trt$meta_model,
                                                                                                meta_model_formula = meta_model_trt$meta_model_formula,
                                                                                                predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                         params_boosted_tree = params_boosted_tree_outcome, 
                                                                                                         params_random_forest = params_random_forest_outcome, 
                                                                                                         params_knn = params_knn_outcome,
                                                                                                         params_lasso = params_lasso_outcome,
                                                                                                         params_glm = params_glm_outcome,
                                                                                                         meta_model = meta_model_outcome$meta_model,
                                                                                                         meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_stacked_classifier[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                         params_boosted_tree = params_boosted_tree_outcome, 
                                                                                                         params_random_forest = params_random_forest_outcome, 
                                                                                                         params_knn = params_knn_outcome,
                                                                                                         params_lasso = params_lasso_outcome,
                                                                                                         params_glm = params_glm_outcome,
                                                                                                         meta_model = meta_model_outcome$meta_model,
                                                                                                         meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                        params = params_boosted_tree_trt, 
                                                                                        model = "boosted_tree",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                 label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                 params = params_boosted_tree_outcome, 
                                                                                                 model = "boosted_tree",
                                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_boosted_tree[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                 label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                 params = params_boosted_tree_outcome, 
                                                                                                 model = "boosted_tree",
                                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                         label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                         params = params_random_forest_trt, 
                                                                                         model = "random_forest",
                                                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                  params = params_random_forest_outcome, 
                                                                                                  model = "random_forest",
                                                                                                  predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_random_forest[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                                  params = params_random_forest_outcome, 
                                                                                                  model = "random_forest",
                                                                                                  predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                               label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                               params = params_knn_trt, 
                                                                               model = "knn",
                                                                               predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_knn_outcome, 
                                                                                        model = "knn",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_knn[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_knn_outcome, 
                                                                                        model = "knn",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                                 label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                                 params = params_lasso_trt, 
                                                                                 model = "lasso",
                                                                                 predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                          params = params_lasso_outcome, 
                                                                                          model = "lasso",
                                                                                          predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_lasso[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                          params = params_lasso_outcome, 
                                                                                          model = "lasso",
                                                                                          predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
    
    prob_trt1_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                                               label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                               params = params_glm_trt, 
                                                                               model = "glm",
                                                                               predict_data = sim_data_est %>% dplyr::select(-c(Y, X)))
    
    prob_outcome1_trt1_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_glm_outcome, 
                                                                                        model = "glm",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_glm[((n*(rep-1))+1):(n*(rep))] <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                                        label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                                        params = params_glm_outcome, 
                                                                                        model = "glm",
                                                                                        predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  }
  
  bias_stacked <- compute_bias(prob_trt1_est = prob_trt1_stacked_classifier, prob_trt1_true = p_X, 
                               prob_outcome1_trt1_est = prob_outcome1_trt1_stacked_classifier, 
                               prob_outcome1_trt1_true = p_Y, 
                               prob_outcome1_trt0_est = prob_outcome1_trt0_stacked_classifier, 
                               prob_outcome1_trt0_true = p_Y,
                               trt = X)
  
  bias_boosted_tree <- compute_bias(prob_trt1_est = prob_trt1_boosted_tree, prob_trt1_true = p_X, 
                                    prob_outcome1_trt1_est = prob_outcome1_trt1_boosted_tree, 
                                    prob_outcome1_trt1_true = p_Y, 
                                    prob_outcome1_trt0_est = prob_outcome1_trt0_boosted_tree, 
                                    prob_outcome1_trt0_true = p_Y,
                                    trt = X)
  
  bias_random_forest <- compute_bias(prob_trt1_est = prob_trt1_random_forest, prob_trt1_true = p_X, 
                                     prob_outcome1_trt1_est = prob_outcome1_trt1_random_forest, 
                                     prob_outcome1_trt1_true = p_Y, 
                                     prob_outcome1_trt0_est = prob_outcome1_trt0_random_forest, 
                                     prob_outcome1_trt0_true = p_Y,
                                     trt = X)
  
  bias_knn <- compute_bias(prob_trt1_est = prob_trt1_knn, prob_trt1_true = p_X,
                           prob_outcome1_trt1_est = prob_outcome1_trt1_knn, 
                           prob_outcome1_trt1_true = p_Y, 
                           prob_outcome1_trt0_est = prob_outcome1_trt0_knn, 
                           prob_outcome1_trt0_true = p_Y,
                           trt = X)
  
  bias_lasso <- compute_bias(prob_trt1_est = prob_trt1_lasso, prob_trt1_true = p_X, 
                             prob_outcome1_trt1_est = prob_outcome1_trt1_lasso, 
                             prob_outcome1_trt1_true = p_Y, 
                             prob_outcome1_trt0_est = prob_outcome1_trt0_lasso, 
                             prob_outcome1_trt0_true = p_Y,
                             trt = X)
  
  bias_glm <- compute_bias(prob_trt1_est = prob_trt1_glm, prob_trt1_true = p_X,
                           prob_outcome1_trt1_est = prob_outcome1_trt1_glm, 
                           prob_outcome1_trt1_true = p_Y, 
                           prob_outcome1_trt0_est = prob_outcome1_trt0_glm, 
                           prob_outcome1_trt0_true = p_Y,
                           trt = X)
  
  bias_out <- list(bias_stacked = bias_stacked,
                   bias_boosted_tree = bias_boosted_tree,
                   bias_random_forest = bias_random_forest,
                   bias_knn = bias_knn,
                   bias_lasso = bias_lasso,
                   bias_glm = bias_glm)
  
  cs_bias_stacked <- compute_CS_bias(prob_trt1_est = prob_trt1_stacked_classifier, prob_trt1_true = p_X, 
                                     prob_outcome1_trt1_est = prob_outcome1_trt1_stacked_classifier, 
                                     prob_outcome1_trt1_true = p_Y, 
                                     prob_outcome1_trt0_est = prob_outcome1_trt0_stacked_classifier, 
                                     prob_outcome1_trt0_true = p_Y,
                                     trt = X)
  
  cs_bias_boosted_tree <- compute_CS_bias(prob_trt1_est = prob_trt1_boosted_tree, prob_trt1_true = p_X, 
                                          prob_outcome1_trt1_est = prob_outcome1_trt1_boosted_tree, 
                                          prob_outcome1_trt1_true = p_Y, 
                                          prob_outcome1_trt0_est = prob_outcome1_trt0_boosted_tree, 
                                          prob_outcome1_trt0_true = p_Y,
                                          trt = X)
  
  cs_bias_random_forest <- compute_CS_bias(prob_trt1_est = prob_trt1_random_forest, prob_trt1_true = p_X, 
                                           prob_outcome1_trt1_est = prob_outcome1_trt1_random_forest, 
                                           prob_outcome1_trt1_true = p_Y, 
                                           prob_outcome1_trt0_est = prob_outcome1_trt0_random_forest, 
                                           prob_outcome1_trt0_true = p_Y,
                                           trt = X)
  
  cs_bias_knn <- compute_CS_bias(prob_trt1_est = prob_trt1_knn, prob_trt1_true = p_X, 
                                 prob_outcome1_trt1_est = prob_outcome1_trt1_knn, 
                                 prob_outcome1_trt1_true = p_Y, 
                                 prob_outcome1_trt0_est = prob_outcome1_trt0_knn, 
                                 prob_outcome1_trt0_true = p_Y,
                                 trt = X)
  
  cs_bias_lasso <- compute_CS_bias(prob_trt1_est = prob_trt1_lasso, prob_trt1_true = p_X, 
                                   prob_outcome1_trt1_est = prob_outcome1_trt1_lasso, 
                                   prob_outcome1_trt1_true = p_Y, 
                                   prob_outcome1_trt0_est = prob_outcome1_trt0_lasso, 
                                   prob_outcome1_trt0_true = p_Y,
                                   trt = X)
  
  cs_bias_glm <- compute_CS_bias(prob_trt1_est = prob_trt1_glm, prob_trt1_true = p_X, 
                                 prob_outcome1_trt1_est = prob_outcome1_trt1_glm, 
                                 prob_outcome1_trt1_true = p_Y, 
                                 prob_outcome1_trt0_est = prob_outcome1_trt0_glm, 
                                 prob_outcome1_trt0_true = p_Y,
                                 trt = X)
  
  cs_bias_out <- list(cs_bias_stacked = cs_bias_stacked,
                      cs_bias_boosted_tree = cs_bias_boosted_tree,
                      cs_bias_random_forest = cs_bias_random_forest,
                      cs_bias_knn = cs_bias_knn,
                      cs_bias_lasso = cs_bias_lasso,
                      cs_bias_glm = cs_bias_glm)
  
  save(bias_out, 
       file=paste0("./output/bias_binary.RData"))
  
  save(cs_bias_out, 
       file=paste0("./output/cs_bias_binary.RData"))
  
})
