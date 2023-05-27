source('./src/dependencies.R')
source('./src/estimation_functions.R')

# parameters for sequence of IF22

K <- c(2,8,32,128,256,512,1024,2,6,10)
degree_of_interactions <- c(1,1,1,1,1,1,1,2,2,2)
polynomial_degree <- c(1,1,1,1,1,1,1,1,1,1)
params_save_name <- '2,8,32,128,256,512,1024,2,6,10_1,1,1,1,1,1,1,2,2,2'

##########################################################
## load data and create basis function
##########################################################

colon_ncdb <- read.csv('./colon_ncdb.csv') %>%
  mutate(X = surg_approach_open,
         Y = death_90) %>% dplyr::select(-c(surg_approach_open, death_90))

set.seed(202848)

colon_ncdb_split <- split_data (colon_ncdb, 2)

colon_ncdb_1 <- colon_ncdb_split[[1]]
colon_ncdb_2 <- colon_ncdb_split[[2]]

##########################################################
## Evaluate b-spline basis functions
##########################################################

knots <- lapply(1:length(K), function(k) {
  list(
    attr(bs(colon_ncdb$distance, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(colon_ncdb$age, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(colon_ncdb$tumor_size, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(colon_ncdb$days_from_dx_to_def_surg, degree=polynomial_degree[k], df=K[k]), "knots")
  )
})

boundary_knots <- lapply(1:length(K), function(k) {
  list(
    attr(bs(colon_ncdb$distance, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(colon_ncdb$age, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(colon_ncdb$tumor_size, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(colon_ncdb$days_from_dx_to_def_surg, degree=polynomial_degree[k], df=K[k]), "Boundary.knots")
  )
})

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

sigma_1_eff1 <- c(compute_sigma(basis = basis_1,
                                trt = colon_ncdb_1$X),
                  lapply(1:length(K), function(i) {
                    nlshrink_cov(colon_ncdb_1$X*basis_1[[i]], k=1)
                  }))

sigma_1_eff0 <- c(compute_sigma(basis = basis_1,
                                trt = 1-colon_ncdb_1$X),
                  lapply(1:length(K), function(i) {
                    nlshrink_cov((1-colon_ncdb_1$X)*basis_1[[i]], k=1)
                  }))

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

sigma_2_eff1 <- c(compute_sigma(basis = basis_2,
                                trt = colon_ncdb_2$X),
                  lapply(1:length(K), function(i) {
                    nlshrink_cov(colon_ncdb_2$X*basis_2[[i]], k=1)
                  }))

sigma_2_eff0 <- c(compute_sigma(basis = basis_2,
                                trt = 1-colon_ncdb_2$X),
                  lapply(1:length(K), function(i) {
                    nlshrink_cov((1-colon_ncdb_2$X)*basis_2[[i]], k=1)
                  }))

##########################################################
## estimate nuisance parameter models using training data
##########################################################

###########
## split 1
###########

params_boosted_tree_trt_1 <- find_params_boosted_tree_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                            label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                            nfold = 5,
                                                            tree_depth = c(1:5),
                                                            shrinkage_factor = seq(0.01,0.05,0.01),
                                                            num_trees = seq(250,500,50),
                                                            num_cores = 10)

params_random_forest_trt_1 <- find_params_random_forest_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                                              label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                              nfold = 5,
                                                              num_trees = seq(500,1000,50),
                                                              num_vars = seq(1,5,1),
                                                              num_cores = 10)

params_knn_trt_1 <- find_params_knn(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                    label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                    nfold = 5,
                                    k = seq(11,101,2),
                                    num_cores = 10)

params_lasso_trt_1 <- find_params_lasso(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)),
                                        label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                        continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                        binary_vars = names(colon_ncdb_1)[!(names(colon_ncdb_1) %in% 
                                                                              c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
                                        continuous_var_spline_knots = 10,
                                        degree_of_interactions = 1,
                                        nfold = 5)

params_glm_trt_1 <- find_params_glm(continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                    binary_vars = names(colon_ncdb_1)[!(names(colon_ncdb_1) %in% 
                                                                          c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
                                    continuous_var_spline_knots = 2)

meta_model_trt_1 <- fit_stacked_classifer_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y, X)), 
                                                label_vector = colon_ncdb_1 %>% dplyr::select(X) %>% {.[[1]]},
                                                params_boosted_tree = params_boosted_tree_trt_1, 
                                                params_random_forest = params_random_forest_trt_1, 
                                                params_knn = params_knn_trt_1,
                                                params_lasso = params_lasso_trt_1,
                                                params_glm = params_glm_trt_1,
                                                num_spline_knots = 4, 
                                                alpha = 0, lambda=0)

params_boosted_tree_outcome_1 <- find_params_boosted_tree_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                nfold = 5,
                                                                tree_depth = c(1:5),
                                                                shrinkage_factor = seq(0.01,0.05,0.01),
                                                                num_trees = seq(250,500,50),
                                                                num_cores = 10)

params_random_forest_outcome_1 <- find_params_random_forest_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                                                  label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                  nfold = 5,
                                                                  num_trees = seq(500,1000,50),
                                                                  num_vars = seq(1,5,1),
                                                                  num_cores = 10)

params_knn_outcome_1 <- find_params_knn(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                        label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                        nfold = 5,
                                        k = seq(11,101,2),
                                        num_cores = 10)

params_lasso_outcome_1 <- find_params_lasso(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)),
                                            label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                            continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                            binary_vars = names(colon_ncdb_1)[!(names(colon_ncdb_1) %in% 
                                                                                  c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'Y'))],
                                            continuous_var_spline_knots = 10,
                                            degree_of_interactions = 1,
                                            nfold = 5)

params_glm_outcome_1 <- find_params_glm(continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                        binary_vars = names(colon_ncdb_1)[!(names(colon_ncdb_1) %in% 
                                                                              c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'Y'))],
                                        continuous_var_spline_knots = 2)

meta_model_outcome_1 <- fit_stacked_classifer_model(covariates_df = colon_ncdb_1 %>% dplyr::select(-c(Y)), 
                                                    label_vector = colon_ncdb_1 %>% dplyr::select(Y) %>% {.[[1]]},
                                                    params_boosted_tree = params_boosted_tree_outcome_1, 
                                                    params_random_forest = params_random_forest_outcome_1, 
                                                    params_knn = params_knn_outcome_1,
                                                    params_lasso = params_lasso_outcome_1,
                                                    params_glm = params_glm_outcome_1,
                                                    num_spline_knots = 4, 
                                                    alpha = 0, lambda=0)

###########
## split 2
###########

params_boosted_tree_trt_2 <- find_params_boosted_tree_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                            label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                            nfold = 5,
                                                            tree_depth = c(1:5),
                                                            shrinkage_factor = seq(0.01,0.05,0.01),
                                                            num_trees = seq(250,500,50),
                                                            num_cores = 10)

params_random_forest_trt_2 <- find_params_random_forest_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                                              label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                              nfold = 5,
                                                              num_trees = seq(500,1000,50),
                                                              num_vars = seq(1,5,1),
                                                              num_cores = 10)

params_knn_trt_2 <- find_params_knn(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                    label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                    nfold = 5,
                                    k = seq(11,101,2),
                                    num_cores = 10)

params_lasso_trt_2 <- find_params_lasso(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)),
                                        label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                        continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                        binary_vars = names(colon_ncdb_2)[!(names(colon_ncdb_2) %in% 
                                                                              c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
                                        continuous_var_spline_knots = 10,
                                        degree_of_interactions = 1,
                                        nfold = 5)

params_glm_trt_2 <- find_params_glm(continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                    binary_vars = names(colon_ncdb_2)[!(names(colon_ncdb_2) %in% 
                                                                          c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
                                    continuous_var_spline_knots = 2)

meta_model_trt_2 <- fit_stacked_classifer_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y, X)), 
                                                label_vector = colon_ncdb_2 %>% dplyr::select(X) %>% {.[[1]]},
                                                params_boosted_tree = params_boosted_tree_trt_2, 
                                                params_random_forest = params_random_forest_trt_2, 
                                                params_knn = params_knn_trt_2,
                                                params_lasso = params_lasso_trt_2,
                                                params_glm = params_glm_trt_2,
                                                num_spline_knots = 4, 
                                                alpha = 0, lambda=0)

params_boosted_tree_outcome_2 <- find_params_boosted_tree_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                nfold = 5,
                                                                tree_depth = c(1:5),
                                                                shrinkage_factor = seq(0.01,0.05,0.01),
                                                                num_trees = seq(250,500,50),
                                                                num_cores = 10)

params_random_forest_outcome_2 <- find_params_random_forest_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                                                  label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                                  nfold = 5,
                                                                  num_trees = seq(500,1000,50),
                                                                  num_vars = seq(1,5,1),
                                                                  num_cores = 10)

params_knn_outcome_2 <- find_params_knn(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                        label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                        nfold = 5,
                                        k = seq(11,101,2),
                                        num_cores = 10)

params_lasso_outcome_2 <- find_params_lasso(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)),
                                            label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                            continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                            binary_vars = names(colon_ncdb_2)[!(names(colon_ncdb_2) %in% 
                                                                                  c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'Y'))],
                                            continuous_var_spline_knots = 10,
                                            degree_of_interactions = 1,
                                            nfold = 5)

params_glm_outcome_2 <- find_params_glm(continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
                                        binary_vars = names(colon_ncdb_2)[!(names(colon_ncdb_2) %in% 
                                                                              c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'Y'))],
                                        continuous_var_spline_knots = 2)

meta_model_outcome_2 <- fit_stacked_classifer_model(covariates_df = colon_ncdb_2 %>% dplyr::select(-c(Y)), 
                                                    label_vector = colon_ncdb_2 %>% dplyr::select(Y) %>% {.[[1]]},
                                                    params_boosted_tree = params_boosted_tree_outcome_2, 
                                                    params_random_forest = params_random_forest_outcome_2, 
                                                    params_knn = params_knn_outcome_2,
                                                    params_lasso = params_lasso_outcome_2,
                                                    params_glm = params_glm_outcome_2,
                                                    num_spline_knots = 4, 
                                                    alpha = 0, lambda=0)

##########################################################
## save parameters
##########################################################

parameters_CRC_real_data <- list(
  K = K,
  knots = knots,
  boundary_knots = boundary_knots,
  degree_of_interactions = degree_of_interactions,
  polynomial_degree = polynomial_degree,
  sigma_1_eff1 = sigma_1_eff1,
  sigma_1_eff0 = sigma_1_eff0,
  sigma_2_eff1 = sigma_2_eff1,
  sigma_2_eff0 = sigma_2_eff0,
  params_boosted_tree_trt_1 = params_boosted_tree_trt_1,
  params_random_forest_trt_1 = params_random_forest_trt_1,
  params_knn_trt_1 = params_knn_trt_1,
  params_lasso_trt_1 = params_lasso_trt_1,
  params_glm_trt_1 = params_glm_trt_1,
  meta_model_trt_1 = meta_model_trt_1,
  params_boosted_tree_outcome_1 = params_boosted_tree_outcome_1,
  params_random_forest_outcome_1 = params_random_forest_outcome_1,
  params_knn_outcome_1 = params_knn_outcome_1,
  params_lasso_outcome_1 = params_lasso_outcome_1,
  params_glm_outcome_1 = params_glm_outcome_1,
  meta_model_outcome_1 = meta_model_outcome_1,
  params_boosted_tree_trt_2 = params_boosted_tree_trt_2,
  params_random_forest_trt_2 = params_random_forest_trt_2,
  params_knn_trt_2 = params_knn_trt_2,
  params_lasso_trt_2 = params_lasso_trt_2,
  params_glm_trt_2 = params_glm_trt_2,
  meta_model_trt_2 = meta_model_trt_2,
  params_boosted_tree_outcome_2 = params_boosted_tree_outcome_2,
  params_random_forest_outcome_2 = params_random_forest_outcome_2,
  params_knn_outcome_2 = params_knn_outcome_2,
  params_lasso_outcome_2 = params_lasso_outcome_2,
  params_glm_outcome_2 = params_glm_outcome_2,
  meta_model_outcome_2 = meta_model_outcome_2
)

save(parameters_CRC_real_data, file=paste0("./params/parameters_CRC_real_data_", params_save_name, ".RData"))

rm(K, degree_of_interactions, polynomial_degree, params_save_name, colon_ncdb, colon_ncdb_split, colon_ncdb_1, colon_ncdb_2, 
   knots, boundary_knots, basis_1, basis_2,
   sigma_1_eff1, sigma_1_eff0, sigma_2_eff1, sigma_2_eff0,
   params_boosted_tree_trt_1, params_random_forest_trt_1, params_knn_trt_1, params_lasso_trt_1, params_glm_trt_1, meta_model_trt_1,
   params_boosted_tree_outcome_1, params_random_forest_outcome_1, params_knn_outcome_1, params_lasso_outcome_1, params_glm_outcome_1, meta_model_outcome_1,
   params_boosted_tree_trt_2, params_random_forest_trt_2, params_knn_trt_2, params_lasso_trt_2, params_glm_trt_2, meta_model_trt_2,
   params_boosted_tree_outcome_2, params_random_forest_outcome_2, params_knn_outcome_2, params_lasso_outcome_2, params_glm_outcome_2, meta_model_outcome_2)