source('./src/dependencies.R')
params_save_name <- '2,8,32,128,256,2,6,10_1,1,1,1,1,2,2,2'
load(paste0("./params/simulation_parameters_CRC_GAN_", params_save_name, ".RData"))

simulation_parameters_CRC_GAN$n <- 5000
simulation_parameters_CRC_GAN$reps <- 200
simulation_parameters_CRC_GAN$save_name <- params_save_name
simulation_parameters_CRC_GAN$num_cores <- 10

df_generated <- read.csv('./df_generated_final.csv', nrows=10000000) %>% filter(!row_number() %in% simulation_parameters_CRC_GAN$oracle_tr_sample_ids)

with(simulation_parameters_CRC_GAN, {
  
  if (num_cores != 1) {
    cl <- makeCluster(num_cores); cl_stop <- TRUE
    
    clusterEvalQ(cl, {
      library(tidyverse)
      library(SMUT) 
      library(splines)
      library(gbm)
      library(ranger)
      library(caret)
      library(glmnet)
      library(pbapply)
    })
    
    clusterExport(cl, 'df_generated')
    
  } 
  else {cl <- 1; cl_stop <- FALSE}
  
  sim_out <- pbsapply (1:reps, function(rep) {
    
    source('./src/estimation_functions.R')
    
    ##########################################################
    ## Simulate estimation sample
    ##########################################################
    
    sim_data_est <- df_generated %>%
      sample_n(n, replace=F) %>% 
      dplyr::select(-c(id, death_90_cf_surg_approach_open1, death_90_cf_surg_approach_open0))
    
    basis_est <- c(
      lapply(1:length(K), function(k) {
        create_b_spline_basis(
          data = sim_data_est,
          continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
          binary_vars = names(sim_data_est)[!(names(sim_data_est) %in% 
                                                c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
          knots = knots[[k]], boundary_knots = boundary_knots[[k]],
          degree_of_interactions = degree_of_interactions[k],
          polynomial_degree = polynomial_degree[k]
        )
      })
    )
    
    sigma_est_eff1 <- c('est' = compute_sigma(basis = basis_est,
                                              trt = sim_data_est$X))
    
    sigma_est_eff0 <- c('est' = compute_sigma(basis = basis_est,
                                    trt = 1-sim_data_est$X))
    
    ##########################################################
    ## compute estimates of the nuisance parameters in estimation data
    ##########################################################
    
    prob_trt1_stacked_classifier <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, tumor_size, distance)),
                                                                     label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                                     params_boosted_tree = params_boosted_tree_trt, 
                                                                     params_random_forest = params_random_forest_trt, 
                                                                     params_knn = params_knn_trt,
                                                                     params_lasso = params_lasso_trt,
                                                                     params_glm = params_glm_trt,
                                                                     meta_model = meta_model_trt$meta_model,
                                                                     meta_model_formula = meta_model_trt$meta_model_formula,
                                                                     predict_data = sim_data_est %>% dplyr::select(c(age, days_from_dx_to_def_surg, tumor_size, distance)))
    
    prob_outcome1_trt1_stacked_classifier <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, tumor_size, distance)),
                                                                              label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome, 
                                                                              params_random_forest = params_random_forest_outcome, 
                                                                              params_knn = params_knn_outcome,
                                                                              params_lasso = params_lasso_outcome,
                                                                              params_glm = params_glm_outcome,
                                                                              meta_model = meta_model_outcome$meta_model,
                                                                              meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                              predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, tumor_size, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_stacked_classifier <- estimate_prob_stacked_classifier(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, tumor_size, distance)),
                                                                              label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome, 
                                                                              params_random_forest = params_random_forest_outcome, 
                                                                              params_knn = params_knn_outcome,
                                                                              params_lasso = params_lasso_outcome,
                                                                              params_glm = params_glm_outcome,
                                                                              meta_model = meta_model_outcome$meta_model,
                                                                              meta_model_formula = meta_model_outcome$meta_model_formula,
                                                                              predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, tumor_size, distance)) %>% mutate(X = 0))
    
    prob_trt1_boosted_tree <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)),
                                                             label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                             params = params_boosted_tree_trt, 
                                                             model = "boosted_tree",
                                                             predict_data = sim_data_est %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)))
    
    prob_outcome1_trt1_boosted_tree <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                                      label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome, 
                                                                      model = "boosted_tree",
                                                                      predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_boosted_tree <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                                      label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome, 
                                                                      model = "boosted_tree",
                                                                      predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)) %>% mutate(X = 0))
    
    prob_trt1_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)),
                                                              label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                              params = params_random_forest_trt, 
                                                              model = "random_forest",
                                                              predict_data = sim_data_est %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)))
    
    prob_outcome1_trt1_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                                       label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome, 
                                                                       model = "random_forest",
                                                                       predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                                       label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome, 
                                                                       model = "random_forest",
                                                                       predict_data = sim_data_est %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)) %>% mutate(X = 0))
    
    prob_trt1_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(days_from_dx_to_def_surg, distance)),
                                                    label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_knn_trt, 
                                                    model = "knn",
                                                    predict_data = sim_data_est %>% dplyr::select(c(days_from_dx_to_def_surg, distance)))
    
    prob_outcome1_trt1_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)),
                                                             label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome, 
                                                             model = "knn",
                                                             predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)),
                                                             label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome, 
                                                             model = "knn",
                                                             predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)) %>% mutate(X = 0))
    
    prob_trt1_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(days_from_dx_to_def_surg, tumor_size, distance)),
                                                      label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                      params = params_lasso_trt, 
                                                      model = "lasso",
                                                      predict_data = sim_data_est %>% dplyr::select(c(days_from_dx_to_def_surg, tumor_size, distance)))
    
    prob_outcome1_trt1_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, tumor_size, distance)),
                                                               label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome, 
                                                               model = "lasso",
                                                               predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, tumor_size, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, tumor_size, distance)),
                                                               label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome, 
                                                               model = "lasso",
                                                               predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, tumor_size, distance)) %>% mutate(X = 0))
    
    prob_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(days_from_dx_to_def_surg, distance)),
                                                    label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_glm_trt, 
                                                    model = "glm",
                                                    predict_data = sim_data_est %>% dplyr::select(c(days_from_dx_to_def_surg, distance)))
    
    prob_outcome1_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)),
                                                             label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome, 
                                                             model = "glm",
                                                             predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)) %>% mutate(X = 1))
    
    prob_outcome1_trt0_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)),
                                                             label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome, 
                                                             model = "glm",
                                                             predict_data = sim_data_est %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)) %>% mutate(X = 0))
    
    
    estimate_AIPW_stacked <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_stacked_classifier,
        prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier,
        prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier
      )
    
    estimate_var_AIPW_stacked <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_stacked_classifier,
        prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier,
        prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier
      )
    
    
    estimate_AIPW_boosted_tree <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_boosted_tree,
        prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree,
        prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree
      )
    
    estimate_var_AIPW_boosted_tree <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_boosted_tree,
        prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree,
        prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree
      )
    
    
    estimate_AIPW_random_forest <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_random_forest,
        prob_outcome1_trt1 = prob_outcome1_trt1_random_forest,
        prob_outcome1_trt0 = prob_outcome1_trt0_random_forest
      )
    
    estimate_var_AIPW_random_forest <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_random_forest,
        prob_outcome1_trt1 = prob_outcome1_trt1_random_forest,
        prob_outcome1_trt0 = prob_outcome1_trt0_random_forest
      )
    
    
    estimate_AIPW_knn <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_knn,
        prob_outcome1_trt1 = prob_outcome1_trt1_knn,
        prob_outcome1_trt0 = prob_outcome1_trt0_knn
      )
    
    estimate_var_AIPW_knn <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_knn,
        prob_outcome1_trt1 = prob_outcome1_trt1_knn,
        prob_outcome1_trt0 = prob_outcome1_trt0_knn
      )
    
    
    estimate_AIPW_lasso <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_lasso,
        prob_outcome1_trt1 = prob_outcome1_trt1_lasso,
        prob_outcome1_trt0 = prob_outcome1_trt0_lasso
      )
    
    estimate_var_AIPW_lasso <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_lasso,
        prob_outcome1_trt1 = prob_outcome1_trt1_lasso,
        prob_outcome1_trt0 = prob_outcome1_trt0_lasso
      )
    
    
    estimate_AIPW_glm <- 
      estimate_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_glm,
        prob_outcome1_trt1 = prob_outcome1_trt1_glm,
        prob_outcome1_trt0 = prob_outcome1_trt0_glm
      )
    
    estimate_var_AIPW_glm <- 
      estimate_var_AIPW(
        trt = sim_data_est$X,
        outcome = sim_data_est$Y,
        prob_trt1 = prob_trt1_glm,
        prob_outcome1_trt1 = prob_outcome1_trt1_glm,
        prob_outcome1_trt0 = prob_outcome1_trt0_glm
      )
    
    
    ##########################################################
    ## compute IF22 with Gram matrix computed from 
    # 1) large separate sample (oracle), 
    # 2) training sample (including non-linear shrinkage version), and
    # 3) estimation sample
    ##########################################################
    
    HOIF_stacked_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_stacked_classifier, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_stacked_classifier, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_stacked_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_stacked_classifier, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_stacked_classifier, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    HOIF_boosted_tree_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_boosted_tree, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_boosted_tree, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_boosted_tree_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_boosted_tree, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_boosted_tree, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    HOIF_random_forest_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_random_forest, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_random_forest, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_random_forest_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_random_forest, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_random_forest, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    HOIF_knn_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_knn, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_knn, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_knn_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_knn, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_knn, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    HOIF_lasso_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_lasso, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_lasso, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_lasso_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_lasso, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_lasso, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    HOIF_glm_eff1 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = sim_data_est$X, pred = prob_trt1_glm, type = 'trt'),
      y_resid = compute_resid(trt = sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
      num_cores = 1
    )
    
    HOIF_glm_eff0 <- compute_HOIF_sequence(
      a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_glm, type = 'trt'),
      y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_glm, type = 'outcome'),
      basis = rep(basis_est, 3),
      sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
      num_cores = 1
    )
    
    
    list(
      Y_1_stacked = estimate_AIPW_stacked$Y_1,
      Y_1_boosted_tree = estimate_AIPW_boosted_tree$Y_1,
      Y_1_random_forest = estimate_AIPW_random_forest$Y_1,
      Y_1_knn = estimate_AIPW_knn$Y_1,
      Y_1_lasso = estimate_AIPW_lasso$Y_1,
      Y_1_glm = estimate_AIPW_glm$Y_1,

      Y_0_stacked = estimate_AIPW_stacked$Y_0,
      Y_0_boosted_tree = estimate_AIPW_boosted_tree$Y_0,
      Y_0_random_forest = estimate_AIPW_random_forest$Y_0,
      Y_0_knn = estimate_AIPW_knn$Y_0,
      Y_0_lasso = estimate_AIPW_lasso$Y_0,
      Y_0_glm = estimate_AIPW_glm$Y_0,

      var_Y_1_stacked = estimate_var_AIPW_stacked$Y_1,
      var_Y_1_boosted_tree = estimate_var_AIPW_boosted_tree$Y_1,
      var_Y_1_random_forest = estimate_var_AIPW_random_forest$Y_1,
      var_Y_1_knn = estimate_var_AIPW_knn$Y_1,
      var_Y_1_lasso = estimate_var_AIPW_lasso$Y_1,
      var_Y_1_glm = estimate_var_AIPW_glm$Y_1,

      var_Y_0_stacked = estimate_var_AIPW_stacked$Y_0,
      var_Y_0_boosted_tree = estimate_var_AIPW_boosted_tree$Y_0,
      var_Y_0_random_forest = estimate_var_AIPW_random_forest$Y_0,
      var_Y_0_knn = estimate_var_AIPW_knn$Y_0,
      var_Y_0_lasso = estimate_var_AIPW_lasso$Y_0,
      var_Y_0_glm = estimate_var_AIPW_glm$Y_0,

      HOIF_stacked_eff1 = HOIF_stacked_eff1,
      HOIF_boosted_tree_eff1 = HOIF_boosted_tree_eff1,
      HOIF_random_forest_eff1 = HOIF_random_forest_eff1,
      HOIF_knn_eff1 = HOIF_knn_eff1,
      HOIF_lasso_eff1 = HOIF_lasso_eff1,
      HOIF_glm_eff1 = HOIF_glm_eff1,
      
      HOIF_stacked_eff0 = HOIF_stacked_eff0,
      HOIF_boosted_tree_eff0 = HOIF_boosted_tree_eff0,
      HOIF_random_forest_eff0 = HOIF_random_forest_eff0,
      HOIF_knn_eff0 = HOIF_knn_eff0,
      HOIF_lasso_eff0 = HOIF_lasso_eff0,
      HOIF_glm_eff0 = HOIF_glm_eff0
    )
  }, cl=cl)
  
  if(cl_stop) {stopCluster(cl)}
  
  save(sim_out, file=paste0("./output/CRC_GAN_sim_", 
                            simulation_parameters_CRC_GAN$n, "_", 
                            simulation_parameters_CRC_GAN$reps, "_", 
                            simulation_parameters_CRC_GAN$save_name, ".RData"))
  
})