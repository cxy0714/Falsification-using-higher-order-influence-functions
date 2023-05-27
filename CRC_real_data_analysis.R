source('./src/dependencies.R')
source('./src/estimation_functions.R')

params_save_name <- '2,8,32,128,256,512,1024,2,6,10_1,1,1,1,1,1,1,2,2,2'
load(paste0("params/parameters_CRC_real_data_", params_save_name, ".RData"))

parameters_CRC_real_data$bootstrap_M <- 200
parameters_CRC_real_data$save_name <- params_save_name
parameters_CRC_real_data$num_cores <- 1

with(parameters_CRC_real_data, {
  
  colon_ncdb <- read.csv('./colon_ncdb.csv') %>%
    mutate(X = surg_approach_open,
           Y = death_90) %>% dplyr::select(-c(surg_approach_open, death_90))
  
  set.seed(202848)
  
  colon_ncdb_split <- split_data (colon_ncdb, 2)
  
  colon_ncdb_1 <- colon_ncdb_split[[1]]
  colon_ncdb_2 <- colon_ncdb_split[[2]]
  
  ###########
  ## split 1
  ###########
  
  basis_1 <- lapply(1:length(K), function(k) {
    create_b_spline_basis(
      data = colon_ncdb_1,
      continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
      binary_vars = names(colon_ncdb_1)[!(names(colon_ncdb_1) %in% 
                                            c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  })
  
  ###########
  ## split 2
  ###########
  
  basis_2 <- lapply(1:length(K), function(k) {
    create_b_spline_basis(
      data = colon_ncdb_2,
      continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
      binary_vars = names(colon_ncdb_2)[!(names(colon_ncdb_2) %in% 
                                            c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  })
  
  ##########################################################
  ## compute estimates in estimation data
  ##########################################################
  
  ###########
  ## split 1
  ###########
  
  prob_trt1_stacked_classifier_1 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                                     label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                                     params_boosted_tree = params_boosted_tree_trt_2, 
                                                                     params_random_forest = params_random_forest_trt_2, 
                                                                     params_knn = params_knn_trt_2,
                                                                     params_lasso = params_lasso_trt_2,
                                                                     params_glm = params_glm_trt_2,
                                                                     meta_model = meta_model_trt_2$meta_model,
                                                                     meta_model_formula = meta_model_trt_2$meta_model_formula,
                                                                     predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_stacked_classifier_1 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                              label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome_2, 
                                                                              params_random_forest = params_random_forest_outcome_2, 
                                                                              params_knn = params_knn_outcome_2,
                                                                              params_lasso = params_lasso_outcome_2,
                                                                              params_glm = params_glm_outcome_2,
                                                                              meta_model = meta_model_outcome_2$meta_model,
                                                                              meta_model_formula = meta_model_outcome_2$meta_model_formula,
                                                                              predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_stacked_classifier_1 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                              label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome_2, 
                                                                              params_random_forest = params_random_forest_outcome_2, 
                                                                              params_knn = params_knn_outcome_2,
                                                                              params_lasso = params_lasso_outcome_2,
                                                                              params_glm = params_glm_outcome_2,
                                                                              meta_model = meta_model_outcome_2$meta_model,
                                                                              meta_model_formula = meta_model_outcome_2$meta_model_formula,
                                                                              predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_boosted_tree_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                             label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                             params = params_boosted_tree_trt_2, 
                                                             model = "boosted_tree",
                                                             predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_boosted_tree_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                      label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome_2, 
                                                                      model = "boosted_tree",
                                                                      predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_boosted_tree_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                      label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome_2, 
                                                                      model = "boosted_tree",
                                                                      predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_random_forest_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                              label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                              params = params_random_forest_trt_2, 
                                                              model = "random_forest",
                                                              predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_random_forest_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                       label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome_2, 
                                                                       model = "random_forest",
                                                                       predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_random_forest_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                       label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome_2, 
                                                                       model = "random_forest",
                                                                       predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_knn_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                    label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_knn_trt_2, 
                                                    model = "knn",
                                                    predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_knn_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome_2, 
                                                             model = "knn",
                                                             predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_knn_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome_2, 
                                                             model = "knn",
                                                             predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_lasso_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                      label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                      params = params_lasso_trt_2, 
                                                      model = "lasso",
                                                      predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_lasso_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                               label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome_2, 
                                                               model = "lasso",
                                                               predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_lasso_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                               label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome_2, 
                                                               model = "lasso",
                                                               predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_glm_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                    label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_glm_trt_2, 
                                                    model = "glm",
                                                    predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_glm_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome_2, 
                                                             model = "glm",
                                                             predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_glm_1 <- estimate_prob_individual_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome_2, 
                                                             model = "glm",
                                                             predict_data = colon_ncdb_1 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  
  estimate_AIPW_stacked_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_stacked_classifier_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier_1
    )
  
  estimate_var_AIPW_stacked_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_stacked_classifier_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier_1
    )
  
  
  estimate_AIPW_boosted_tree_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_boosted_tree_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree_1
    )
  
  estimate_var_AIPW_boosted_tree_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_boosted_tree_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree_1
    )
  
  
  estimate_AIPW_random_forest_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_random_forest_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_random_forest_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_random_forest_1
    )
  
  estimate_var_AIPW_random_forest_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_random_forest_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_random_forest_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_random_forest_1
    )
  
  
  estimate_AIPW_knn_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_knn_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_knn_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_knn_1
    )
  
  estimate_var_AIPW_knn_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_knn_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_knn_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_knn_1
    )
  
  
  estimate_AIPW_lasso_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_lasso_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_lasso_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_lasso_1
    )
  
  estimate_var_AIPW_lasso_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_lasso_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_lasso_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_lasso_1
    )
  
  
  estimate_AIPW_glm_1 <- 
    estimate_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_glm_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_glm_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_glm_1
    )
  
  estimate_var_AIPW_glm_1 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_1$X,
      outcome = colon_ncdb_1$Y,
      prob_trt1 = prob_trt1_glm_1,
      prob_outcome1_trt1 = prob_outcome1_trt1_glm_1,
      prob_outcome1_trt0 = prob_outcome1_trt0_glm_1
    )
  
  
  ###########
  ## split 2
  ###########
  
  prob_trt1_stacked_classifier_2 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                                     label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                                     params_boosted_tree = params_boosted_tree_trt_1, 
                                                                     params_random_forest = params_random_forest_trt_1, 
                                                                     params_knn = params_knn_trt_1,
                                                                     params_lasso = params_lasso_trt_1,
                                                                     params_glm = params_glm_trt_1,
                                                                     meta_model = meta_model_trt_1$meta_model,
                                                                     meta_model_formula = meta_model_trt_1$meta_model_formula,
                                                                     predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_stacked_classifier_2 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                              label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome_1, 
                                                                              params_random_forest = params_random_forest_outcome_1, 
                                                                              params_knn = params_knn_outcome_1,
                                                                              params_lasso = params_lasso_outcome_1,
                                                                              params_glm = params_glm_outcome_1,
                                                                              meta_model = meta_model_outcome_1$meta_model,
                                                                              meta_model_formula = meta_model_outcome_1$meta_model_formula,
                                                                              predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_stacked_classifier_2 <- estimate_prob_stacked_classifier(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                              label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                              params_boosted_tree = params_boosted_tree_outcome_1, 
                                                                              params_random_forest = params_random_forest_outcome_1, 
                                                                              params_knn = params_knn_outcome_1,
                                                                              params_lasso = params_lasso_outcome_1,
                                                                              params_glm = params_glm_outcome_1,
                                                                              meta_model = meta_model_outcome_1$meta_model,
                                                                              meta_model_formula = meta_model_outcome_1$meta_model_formula,
                                                                              predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_boosted_tree_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                             label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                             params = params_boosted_tree_trt_1, 
                                                             model = "boosted_tree",
                                                             predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_boosted_tree_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                      label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome_1, 
                                                                      model = "boosted_tree",
                                                                      predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_boosted_tree_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                      label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                      params = params_boosted_tree_outcome_1, 
                                                                      model = "boosted_tree",
                                                                      predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_random_forest_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                              label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                              params = params_random_forest_trt_1, 
                                                              model = "random_forest",
                                                              predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_random_forest_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                       label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome_1, 
                                                                       model = "random_forest",
                                                                       predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_random_forest_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                       label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                       params = params_random_forest_outcome_1, 
                                                                       model = "random_forest",
                                                                       predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_knn_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                    label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_knn_trt_1, 
                                                    model = "knn",
                                                    predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_knn_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome_1, 
                                                             model = "knn",
                                                             predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_knn_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_knn_outcome_1, 
                                                             model = "knn",
                                                             predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_lasso_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                      label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                      params = params_lasso_trt_1, 
                                                      model = "lasso",
                                                      predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_lasso_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                               label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome_1, 
                                                               model = "lasso",
                                                               predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_lasso_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                               label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                               params = params_lasso_outcome_1, 
                                                               model = "lasso",
                                                               predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  prob_trt1_glm_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                    label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                    params = params_glm_trt_1, 
                                                    model = "glm",
                                                    predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y, X)))
  
  prob_outcome1_trt1_glm_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome_1, 
                                                             model = "glm",
                                                             predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 1))
  
  prob_outcome1_trt0_glm_2 <- estimate_prob_individual_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                             label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                             params = params_glm_outcome_1, 
                                                             model = "glm",
                                                             predict_data = colon_ncdb_2 %>% dplyr::select(-c(Y)) %>% mutate(X = 0))
  
  
  estimate_AIPW_stacked_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_stacked_classifier_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier_2
    )
  
  estimate_var_AIPW_stacked_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_stacked_classifier_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_stacked_classifier_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_stacked_classifier_2
    )
  
  
  estimate_AIPW_boosted_tree_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_boosted_tree_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree_2
    )
  
  estimate_var_AIPW_boosted_tree_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_boosted_tree_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_boosted_tree_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_boosted_tree_2
    )
  
  
  estimate_AIPW_random_forest_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_random_forest_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_random_forest_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_random_forest_2
    )
  
  estimate_var_AIPW_random_forest_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_random_forest_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_random_forest_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_random_forest_2
    )
  
  
  estimate_AIPW_knn_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_knn_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_knn_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_knn_2
    )
  
  estimate_var_AIPW_knn_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_knn_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_knn_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_knn_2
    )
  
  
  estimate_AIPW_lasso_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_lasso_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_lasso_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_lasso_2
    )
  
  estimate_var_AIPW_lasso_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_lasso_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_lasso_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_lasso_2
    )
  
  
  estimate_AIPW_glm_2 <- 
    estimate_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_glm_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_glm_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_glm_2
    )
  
  estimate_var_AIPW_glm_2 <- 
    estimate_var_AIPW(
      trt = colon_ncdb_2$X,
      outcome = colon_ncdb_2$Y,
      prob_trt1 = prob_trt1_glm_2,
      prob_outcome1_trt1 = prob_outcome1_trt1_glm_2,
      prob_outcome1_trt0 = prob_outcome1_trt0_glm_2
    )
  
  
  ##########################################################
  ## compute IF22 (for risk difference) with Gram matrix computed from 
  # 1) training sample (including shrinkage versions), and
  # 3) estimation sample
  ##########################################################
  
  ###########
  ## split 1
  ###########
  
  HOIF_stacked_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_stacked_classifier_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_stacked_classifier_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_stacked_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_stacked_classifier_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_stacked_classifier_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_stacked_eff_1 <- HOIF_stacked_eff1_1-HOIF_stacked_eff0_1
  HOIF_stacked_eff_1[1,] <- HOIF_stacked_eff1_1[1,]
  
  var_HOIF_stacked_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_stacked_classifier_1, type = 'trt'),
                                                      a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_stacked_classifier_1, type = 'trt'),
                                                      y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_stacked_classifier_1, type = 'outcome'),
                                                      y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_stacked_classifier_1, type = 'outcome'),
                                                      basis = rep(basis_1, 3),
                                                      sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                      sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                      num_cores = num_cores,
                                                      M = bootstrap_M)
  
  HOIF_boosted_tree_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_boosted_tree_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_boosted_tree_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_boosted_tree_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_boosted_tree_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff_1 <- HOIF_boosted_tree_eff1_1-HOIF_boosted_tree_eff0_1
  HOIF_boosted_tree_eff_1[1,] <- HOIF_boosted_tree_eff1_1[1,]
  
  var_HOIF_boosted_tree_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_boosted_tree_1, type = 'trt'),
                                                           a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_boosted_tree_1, type = 'trt'),
                                                           y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_boosted_tree_1, type = 'outcome'),
                                                           y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_boosted_tree_1, type = 'outcome'),
                                                           basis = rep(basis_1, 3),
                                                           sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                           sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                           num_cores = num_cores,
                                                           M = bootstrap_M)
  
  HOIF_random_forest_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_random_forest_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_random_forest_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_random_forest_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_random_forest_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_random_forest_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_random_forest_eff_1 <- HOIF_random_forest_eff1_1-HOIF_random_forest_eff0_1
  HOIF_random_forest_eff_1[1,] <- HOIF_random_forest_eff1_1[1,]
  
  var_HOIF_random_forest_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_random_forest_1, type = 'trt'),
                                                            a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_random_forest_1, type = 'trt'),
                                                            y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_random_forest_1, type = 'outcome'),
                                                            y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_random_forest_1, type = 'outcome'),
                                                            basis = rep(basis_1, 3),
                                                            sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                            sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                            num_cores = num_cores,
                                                            M = bootstrap_M)
  
  HOIF_knn_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_knn_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_knn_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_knn_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_knn_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_knn_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_knn_eff_1 <- HOIF_knn_eff1_1-HOIF_knn_eff0_1
  HOIF_knn_eff_1[1,] <- HOIF_knn_eff1_1[1,]
  
  var_HOIF_knn_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_knn_1, type = 'trt'),
                                                  a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_knn_1, type = 'trt'),
                                                  y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_knn_1, type = 'outcome'),
                                                  y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_knn_1, type = 'outcome'),
                                                  basis = rep(basis_1, 3),
                                                  sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                  sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                  num_cores = num_cores,
                                                  M = bootstrap_M)
  
  HOIF_lasso_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_lasso_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_lasso_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_lasso_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_lasso_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_lasso_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_lasso_eff_1 <- HOIF_lasso_eff1_1-HOIF_lasso_eff0_1
  HOIF_lasso_eff_1[1,] <- HOIF_lasso_eff1_1[1,]
  
  var_HOIF_lasso_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_lasso_1, type = 'trt'),
                                                    a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_lasso_1, type = 'trt'),
                                                    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_lasso_1, type = 'outcome'),
                                                    y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_lasso_1, type = 'outcome'),
                                                    basis = rep(basis_1, 3),
                                                    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                    sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                    num_cores = num_cores,
                                                    M = bootstrap_M)
  
  HOIF_glm_eff1_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_glm_1, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_glm_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_glm_eff0_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_glm_1, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_glm_1, type = 'outcome'),
    basis = rep(basis_1, 3),
    sigma = c(sigma_2_eff0, sigma_1_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_glm_eff_1 <- HOIF_glm_eff1_1-HOIF_glm_eff0_1
  HOIF_glm_eff_1[1,] <- HOIF_glm_eff1_1[1,]
  
  var_HOIF_glm_eff_1 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_1$X, pred = prob_trt1_glm_1, type = 'trt'),
                                                  a_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, pred = 1-prob_trt1_glm_1, type = 'trt'),
                                                  y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_glm_1, type = 'outcome'),
                                                  y_resid2 = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_glm_1, type = 'outcome'),
                                                  basis = rep(basis_1, 3),
                                                  sigma = c(sigma_2_eff1, sigma_1_eff1[1:10]),
                                                  sigma2 = c(sigma_2_eff0, sigma_1_eff0[1:10]),
                                                  num_cores = num_cores,
                                                  M = bootstrap_M)
  
  ###########
  ## split 2
  ###########
  
  HOIF_stacked_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_stacked_classifier_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_stacked_classifier_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_stacked_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_stacked_classifier_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_stacked_classifier_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_stacked_eff_2 <- HOIF_stacked_eff1_2-HOIF_stacked_eff0_2
  HOIF_stacked_eff_2[1,] <- HOIF_stacked_eff1_2[1,]
  
  var_HOIF_stacked_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_stacked_classifier_2, type = 'trt'),
                                                              a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_stacked_classifier_2, type = 'trt'),
                                                              y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_stacked_classifier_2, type = 'outcome'),
                                                              y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_stacked_classifier_2, type = 'outcome'),
                                                              basis = c(basis_2, basis_2, basis_2),
                                                              sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                              sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                              num_cores = num_cores,
                                                              M = bootstrap_M)
  
  HOIF_boosted_tree_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_boosted_tree_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_boosted_tree_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_boosted_tree_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_boosted_tree_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff_2 <- HOIF_boosted_tree_eff1_2-HOIF_boosted_tree_eff0_2
  HOIF_boosted_tree_eff_2[1,] <- HOIF_boosted_tree_eff1_2[1,]
  
  var_HOIF_boosted_tree_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_boosted_tree_2, type = 'trt'),
                                                                   a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_boosted_tree_2, type = 'trt'),
                                                                   y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_boosted_tree_2, type = 'outcome'),
                                                                   y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_boosted_tree_2, type = 'outcome'),
                                                                   basis = c(basis_2, basis_2, basis_2),
                                                                   sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                                   sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                                   num_cores = num_cores,
                                                                   M = bootstrap_M)
  
  HOIF_random_forest_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_random_forest_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_random_forest_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_random_forest_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_random_forest_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_random_forest_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_random_forest_eff_2 <- HOIF_random_forest_eff1_2-HOIF_random_forest_eff0_2
  HOIF_random_forest_eff_2[1,] <- HOIF_random_forest_eff1_2[1,]
  
  var_HOIF_random_forest_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_random_forest_2, type = 'trt'),
                                                                    a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_random_forest_2, type = 'trt'),
                                                                    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_random_forest_2, type = 'outcome'),
                                                                    y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_random_forest_2, type = 'outcome'),
                                                                    basis = c(basis_2, basis_2, basis_2),
                                                                    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                                    sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                                    num_cores = num_cores,
                                                                    M = bootstrap_M)
  
  HOIF_knn_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_knn_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_knn_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_knn_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_knn_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_knn_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_knn_eff_2 <- HOIF_knn_eff1_2-HOIF_knn_eff0_2
  HOIF_knn_eff_2[1,] <- HOIF_knn_eff1_2[1,]
  
  var_HOIF_knn_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_knn_2, type = 'trt'),
                                                          a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_knn_2, type = 'trt'),
                                                          y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_knn_2, type = 'outcome'),
                                                          y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_knn_2, type = 'outcome'),
                                                          basis = c(basis_2, basis_2, basis_2),
                                                          sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                          sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                          num_cores = num_cores,
                                                          M = bootstrap_M)
  
  HOIF_lasso_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_lasso_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_lasso_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    tr_corr = F,
    num_cores = 1
  )
  
  HOIF_lasso_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_lasso_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_lasso_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    tr_corr = F,
    num_cores = 1
  )
  
  HOIF_lasso_eff_2 <- HOIF_lasso_eff1_2-HOIF_lasso_eff0_2
  HOIF_lasso_eff_2[1,] <- HOIF_lasso_eff1_2[1,]
  
  var_HOIF_lasso_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_lasso_2, type = 'trt'),
                                                            a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_lasso_2, type = 'trt'),
                                                            y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_lasso_2, type = 'outcome'),
                                                            y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_lasso_2, type = 'outcome'),
                                                            basis = c(basis_2, basis_2, basis_2),
                                                            sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                            sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                            num_cores = num_cores,
                                                            M = bootstrap_M)
  
  HOIF_glm_eff1_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_glm_2, type = 'trt'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_glm_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
    num_cores = 1
  )
  
  HOIF_glm_eff0_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_glm_2, type = 'trt'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_glm_2, type = 'outcome'),
    basis = c(basis_2, basis_2, basis_2),
    sigma = c(sigma_1_eff0, sigma_2_eff0[1:10]),
    num_cores = 1
  )
  
  HOIF_glm_eff_2 <- HOIF_glm_eff1_2-HOIF_glm_eff0_2
  HOIF_glm_eff_2[1,] <- HOIF_glm_eff1_2[1,]
  
  var_HOIF_glm_eff_2 <- compute_var_HOIF_sequence(a_resid = compute_resid(trt = colon_ncdb_2$X, pred = prob_trt1_glm_2, type = 'trt'),
                                                          a_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, pred = 1-prob_trt1_glm_2, type = 'trt'),
                                                          y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_glm_2, type = 'outcome'),
                                                          y_resid2 = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_glm_2, type = 'outcome'),
                                                          basis = c(basis_2, basis_2, basis_2),
                                                          sigma = c(sigma_1_eff1, sigma_2_eff1[1:10]),
                                                          sigma2 = c(sigma_1_eff0, sigma_2_eff0[1:10]),
                                                          num_cores = num_cores,
                                                          M = bootstrap_M)
  
  ##########################################################
  ## save
  ##########################################################
  
  results <- list(
    estimate_AIPW_stacked_1 = estimate_AIPW_stacked_1,
    estimate_var_AIPW_stacked_1 = estimate_var_AIPW_stacked_1,
    estimate_AIPW_boosted_tree_1 = estimate_AIPW_boosted_tree_1,
    estimate_var_AIPW_boosted_tree_1 = estimate_var_AIPW_boosted_tree_1,
    estimate_AIPW_random_forest_1 = estimate_AIPW_random_forest_1,
    estimate_var_AIPW_random_forest_1 = estimate_var_AIPW_random_forest_1,
    estimate_AIPW_knn_1 = estimate_AIPW_knn_1,
    estimate_var_AIPW_knn_1 = estimate_var_AIPW_knn_1,
    estimate_AIPW_lasso_1 = estimate_AIPW_lasso_1,
    estimate_var_AIPW_lasso_1 = estimate_var_AIPW_lasso_1,
    estimate_AIPW_glm_1 = estimate_AIPW_glm_1,
    estimate_var_AIPW_glm_1 = estimate_var_AIPW_glm_1, 
    
    estimate_AIPW_stacked_2 = estimate_AIPW_stacked_2,
    estimate_var_AIPW_stacked_2 = estimate_var_AIPW_stacked_2,
    estimate_AIPW_boosted_tree_2 = estimate_AIPW_boosted_tree_2,
    estimate_var_AIPW_boosted_tree_2 = estimate_var_AIPW_boosted_tree_2,
    estimate_AIPW_random_forest_2 = estimate_AIPW_random_forest_2,
    estimate_var_AIPW_random_forest_2 = estimate_var_AIPW_random_forest_2,
    estimate_AIPW_knn_2 = estimate_AIPW_knn_2,
    estimate_var_AIPW_knn_2 = estimate_var_AIPW_knn_2,
    estimate_AIPW_lasso_2 = estimate_AIPW_lasso_2,
    estimate_var_AIPW_lasso_2 = estimate_var_AIPW_lasso_2,
    estimate_AIPW_glm_2 = estimate_AIPW_glm_2,
    estimate_var_AIPW_glm_2 = estimate_var_AIPW_glm_2, 
    
    HOIF_stacked_eff_1 = HOIF_stacked_eff_1,
    var_HOIF_stacked_eff_1 = var_HOIF_stacked_eff_1,
    HOIF_boosted_tree_eff_1 = HOIF_boosted_tree_eff_1,
    var_HOIF_boosted_tree_eff_1 = var_HOIF_boosted_tree_eff_1,
    HOIF_random_forest_eff_1 = HOIF_random_forest_eff_1,
    var_HOIF_random_forest_eff_1 = var_HOIF_random_forest_eff_1,
    HOIF_knn_eff_1 = HOIF_knn_eff_1,
    var_HOIF_knn_eff_1 = var_HOIF_knn_eff_1,
    HOIF_lasso_eff_1 = HOIF_lasso_eff_1,
    var_HOIF_lasso_eff_1 = var_HOIF_lasso_eff_1,
    HOIF_glm_eff_1 = HOIF_glm_eff_1,
    var_HOIF_glm_eff_1 = var_HOIF_glm_eff_1,
    
    HOIF_stacked_eff_2 = HOIF_stacked_eff_2,
    var_HOIF_stacked_eff_2 = var_HOIF_stacked_eff_2,
    HOIF_boosted_tree_eff_2 = HOIF_boosted_tree_eff_2,
    var_HOIF_boosted_tree_eff_2 = var_HOIF_boosted_tree_eff_2,
    HOIF_random_forest_eff_2 = HOIF_random_forest_eff_2,
    var_HOIF_random_forest_eff_2 = var_HOIF_random_forest_eff_2,
    HOIF_knn_eff_2 = HOIF_knn_eff_2,
    var_HOIF_knn_eff_2 = var_HOIF_knn_eff_2,
    HOIF_lasso_eff_2 = HOIF_lasso_eff_2,
    var_HOIF_lasso_eff_2 = var_HOIF_lasso_eff_2,
    HOIF_glm_eff_2 = HOIF_glm_eff_2,
    var_HOIF_glm_eff_2 = var_HOIF_glm_eff_2
  )
  
  save(results, file=paste0("./output/CRC_real_data_analysis_", save_name, ".RData"))
  
})