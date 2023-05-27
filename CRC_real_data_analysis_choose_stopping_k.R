source('./src/dependencies.R')
source('./src/estimation_functions.R')

params_save_name <- '2,8,32,128,256,512,1024,2,6,10_1,1,1,1,1,1,1,2,2,2'
load(paste0("params/parameters_CRC_real_data_", params_save_name, ".RData"))

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
  ## compute estimates of the nuisance parameters and psi in estimation data
  ##########################################################
  
  ###########
  ## split 1
  ###########
  
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
  
  
  ###########
  ## split 2
  ###########

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
  
  
  ##########################################################
  ## compute IF22 for choosing the stopping value of k with Gram matrix computed from the estimation sample
  ##########################################################
  
  # a_resid is multiplied by -1 to correct for the fact that Sigma here is defined as A*t(basis)%*%basis rather than the -A*t(basis)%*%basis 
  # which is in expected in the functions used to compute IF22 for the treatment effect parameters
  
  ###########
  ## split 1
  ###########
  
  HOIF_stacked_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_stacked_classifier_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_stacked_classifier_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_stacked_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_stacked_classifier_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_stacked_classifier_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff0,
    num_cores = 1
  )

  
  HOIF_boosted_tree_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_boosted_tree_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_boosted_tree_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_boosted_tree_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_boosted_tree_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_random_forest_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_random_forest_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_random_forest_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_random_forest_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_random_forest_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_random_forest_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff0,
    num_cores = 1
  )

  
  HOIF_knn_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_knn_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_knn_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_knn_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_knn_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_knn_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff0,
    num_cores = 1
  )

  
  HOIF_lasso_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_lasso_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_lasso_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_lasso_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_lasso_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_lasso_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff0,
    num_cores = 1
  )
  
  
  HOIF_glm_eff1_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_glm_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt1_glm_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff1,
    num_cores = 1
  )
  
  HOIF_glm_eff0_outcome_1 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_1$X), outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_glm_1, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_1$X, outcome = colon_ncdb_1$Y, pred = prob_outcome1_trt0_glm_1, type = 'outcome', estimand = 'effect'),
    basis = basis_1,
    sigma = sigma_1_eff0,
    num_cores = 1
  )

  
  ###########
  ## split 2
  ###########
  
  HOIF_stacked_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_stacked_classifier_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_stacked_classifier_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_stacked_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_stacked_classifier_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_stacked_classifier_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_boosted_tree_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_boosted_tree_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_boosted_tree_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_boosted_tree_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_boosted_tree_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_boosted_tree_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_random_forest_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_random_forest_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_random_forest_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_random_forest_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_random_forest_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_random_forest_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_knn_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_knn_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_knn_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_knn_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_knn_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_knn_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_lasso_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_lasso_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_lasso_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_lasso_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_lasso_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_lasso_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  HOIF_glm_eff1_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_glm_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt1_glm_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff1,
    num_cores = 1
  )
  
  HOIF_glm_eff0_outcome_2 <- compute_HOIF_sequence(
    a_resid = compute_resid(trt = -(1-colon_ncdb_2$X), outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_glm_2, type = 'outcome', estimand = 'effect'),
    y_resid = compute_resid(trt = 1-colon_ncdb_2$X, outcome = colon_ncdb_2$Y, pred = prob_outcome1_trt0_glm_2, type = 'outcome', estimand = 'effect'),
    basis = basis_2,
    sigma = sigma_2_eff0,
    num_cores = 1
  )

  
  ##########################################################
  ## save
  ##########################################################
  
  results <- list(
    HOIF_stacked_eff1_outcome_1 = HOIF_stacked_eff1_outcome_1,
    HOIF_stacked_eff0_outcome_1 = HOIF_stacked_eff0_outcome_1,
    
    HOIF_boosted_tree_eff1_outcome_1 = HOIF_boosted_tree_eff1_outcome_1,
    HOIF_boosted_tree_eff0_outcome_1 = HOIF_boosted_tree_eff0_outcome_1,
    
    HOIF_random_forest_eff1_outcome_1 = HOIF_random_forest_eff1_outcome_1,
    HOIF_random_forest_eff0_outcome_1 = HOIF_random_forest_eff0_outcome_1,
    
    HOIF_knn_eff1_outcome_1 = HOIF_knn_eff1_outcome_1,
    HOIF_knn_eff0_outcome_1 = HOIF_knn_eff0_outcome_1,
    
    HOIF_lasso_eff1_outcome_1 = HOIF_lasso_eff1_outcome_1,
    HOIF_lasso_eff0_outcome_1 = HOIF_lasso_eff0_outcome_1,
    
    HOIF_glm_eff1_outcome_1 = HOIF_glm_eff1_outcome_1,
    HOIF_glm_eff0_outcome_1 = HOIF_glm_eff0_outcome_1,
    
    HOIF_stacked_eff1_outcome_2 = HOIF_stacked_eff1_outcome_2,
    HOIF_stacked_eff0_outcome_2 = HOIF_stacked_eff0_outcome_2,
    
    HOIF_boosted_tree_eff1_outcome_2 = HOIF_boosted_tree_eff1_outcome_2,
    HOIF_boosted_tree_eff0_outcome_2 = HOIF_boosted_tree_eff0_outcome_2,
    
    HOIF_random_forest_eff1_outcome_2 = HOIF_random_forest_eff1_outcome_2,
    HOIF_random_forest_eff0_outcome_2 = HOIF_random_forest_eff0_outcome_2,
    
    HOIF_knn_eff1_outcome_2 = HOIF_knn_eff1_outcome_2,
    HOIF_knn_eff0_outcome_2 = HOIF_knn_eff0_outcome_2,
    
    HOIF_lasso_eff1_outcome_2 = HOIF_lasso_eff1_outcome_2,
    HOIF_lasso_eff0_outcome_2 = HOIF_lasso_eff0_outcome_2,
    
    HOIF_glm_eff1_outcome_2 = HOIF_glm_eff1_outcome_2,
    HOIF_glm_eff0_outcome_2 = HOIF_glm_eff0_outcome_2
  )
  
  save(results, file=paste0("./output/CRC_real_data_analysis_choose_stopping_k_", save_name, ".RData"))
  
})