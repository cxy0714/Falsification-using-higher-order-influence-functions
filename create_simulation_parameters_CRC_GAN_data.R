source('./src/dependencies.R')
source('./src/estimation_functions.R')

# parameters for sequence of IF22

n_train <- 5000
n_oracle <- 1000000
K <- c(2,8,32,128,256,2,6,10)
degree_of_interactions <- c(1,1,1,1,1,2,2,2)
polynomial_degree <- c(1,1,1,1,1,1,1,1)
params_save_name <- '2,8,32,128,256,2,6,10_1,1,1,1,1,2,2,2'

set.seed(202848)

##########################################################
## load all data, estimate ground truth, and create basis function
##########################################################

df_generated <- read.csv('./df_generated_final.csv', nrows=10000000)

oracle_tr_sample_ids <- sample(1:nrow(df_generated), n_train+n_oracle, replace=F)
oracle_sample_ids <- sample(oracle_tr_sample_ids, n_oracle, replace=F)
tr_sample_ids <- oracle_tr_sample_ids[!oracle_tr_sample_ids %in% oracle_sample_ids]

sim_data_oracle <- df_generated %>%
  filter(row_number() %in% oracle_sample_ids) %>% 
  dplyr::select(-c(id, death_90_cf_surg_approach_open1, death_90_cf_surg_approach_open0))

knots <- lapply(1:length(K), function(k) {
  list(
    attr(bs(sim_data_oracle$distance, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(sim_data_oracle$age, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(sim_data_oracle$tumor_size, degree=polynomial_degree[k], df=K[k]), "knots"),
    attr(bs(sim_data_oracle$days_from_dx_to_def_surg, degree=polynomial_degree[k], df=K[k]), "knots")
  )
})

boundary_knots <- lapply(1:length(K), function(k) {
  list(
    attr(bs(sim_data_oracle$distance, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(sim_data_oracle$age, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(sim_data_oracle$tumor_size, degree=polynomial_degree[k], df=K[k]), "Boundary.knots"),
    attr(bs(sim_data_oracle$days_from_dx_to_def_surg, degree=polynomial_degree[k], df=K[k]), "Boundary.knots")
  )
})


basis_oracle <- c(
  lapply(1:length(K), function(k) {
    create_b_spline_basis(
      data = sim_data_oracle,
      continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
      binary_vars = names(sim_data_oracle)[!(names(sim_data_oracle) %in% 
                                               c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  })
)

sigma_oracle_eff1 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = sim_data_oracle$X))

sigma_oracle_eff0 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = 1-sim_data_oracle$X))

##########################################################
## simulate training data
##########################################################

sim_data_tr <- df_generated %>%
  filter(row_number() %in% tr_sample_ids) %>% 
  dplyr::select(-c(id, death_90_cf_surg_approach_open1, death_90_cf_surg_approach_open0))

basis_tr <- c(
  lapply(1:length(K), function(k) {
    create_b_spline_basis(
      data = sim_data_tr,
      continuous_vars = c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg'),
      binary_vars = names(sim_data_tr)[!(names(sim_data_tr) %in% 
                                           c('distance', 'age', 'tumor_size', 'days_from_dx_to_def_surg', 'X', 'Y'))],
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  })
)

sigma_tr_eff1 <- c('trt' = compute_sigma(basis = basis_tr,
                                         trt = sim_data_tr$X),
                   'nlshrink' = lapply(1:length(K), function(i) {
                     nlshrink_cov(sim_data_tr$X*basis_tr[[i]], k=1)
                   }))

sigma_tr_eff0 <- c('trt' = compute_sigma(basis = basis_tr,
                                 trt = 1-sim_data_tr$X),
                   'nlshrink' = lapply(1:length(K), function(i) {
                     nlshrink_cov((1-sim_data_tr$X)*basis_tr[[i]], k=1)
                   }))

##########################################################
## estimate nuisance parameter models using training data
##########################################################

params_boosted_tree_trt <- find_params_boosted_tree_model(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)),
                                                          label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                          nfold = 5,
                                                          tree_depth = c(1:5),
                                                          shrinkage_factor = seq(0.01,0.05,0.01),
                                                          num_trees = seq(250,500,50),
                                                          num_cores = 10)

params_random_forest_trt <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, distance)),
                                                            label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                            nfold = 5,
                                                            num_trees = seq(500,1000,50),
                                                            num_vars = seq(1,3,1),
                                                            num_cores = 10)

params_knn_trt <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(c(days_from_dx_to_def_surg, distance)),
                                  label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                  nfold = 5,
                                  k = seq(11,101,2),
                                  num_cores = 10)

params_lasso_trt <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(c(days_from_dx_to_def_surg, tumor_size, distance)),
                                      label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                      continuous_vars = c('distance', 'tumor_size', 'days_from_dx_to_def_surg'),
                                      binary_vars = NULL,
                                      continuous_var_spline_knots = 10,
                                      degree_of_interactions = 1,
                                      nfold = 5)

params_glm_trt <- find_params_glm(continuous_vars = c('distance', 'days_from_dx_to_def_surg'),
                                  binary_vars = NULL,
                                  continuous_var_spline_knots = 2)

meta_model_trt <- fit_stacked_classifer_model(covariates_df = sim_data_tr %>% dplyr::select(c(age, days_from_dx_to_def_surg, tumor_size, distance)), 
                                              label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                              params_boosted_tree = params_boosted_tree_trt, 
                                              params_random_forest = params_random_forest_trt, 
                                              params_knn = params_knn_trt,
                                              params_lasso = params_lasso_trt,
                                              params_glm = params_glm_trt,
                                              num_spline_knots = 4, 
                                              alpha = 0, lambda=0)

params_boosted_tree_outcome <- find_params_boosted_tree_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                              label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                              nfold = 5,
                                                              tree_depth = c(1:5),
                                                              shrinkage_factor = seq(0.01,0.05,0.01),
                                                              num_trees = seq(250,500,50),
                                                              num_cores = 10)

params_random_forest_outcome <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, distance)),
                                                                label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                nfold = 5,
                                                                num_trees = seq(500,1000,50),
                                                                num_vars = seq(1,4,1),
                                                                num_cores = 10)

params_knn_outcome <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, distance)),
                                      label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                      nfold = 5,
                                      k = seq(11,101,2),
                                      num_cores = 10)

params_lasso_outcome <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(c(X, days_from_dx_to_def_surg, tumor_size, distance)),
                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                          continuous_vars = c('distance', 'tumor_size', 'days_from_dx_to_def_surg'),
                                          binary_vars = 'X',
                                          continuous_var_spline_knots = 10,
                                          degree_of_interactions = 1,
                                          nfold = 5)

params_glm_outcome <- find_params_glm(continuous_vars = c('distance', 'days_from_dx_to_def_surg'),
                                      binary_vars = 'X',
                                      continuous_var_spline_knots = 2)

meta_model_outcome <- fit_stacked_classifer_model(covariates_df = sim_data_tr %>% dplyr::select(c(X, age, days_from_dx_to_def_surg, tumor_size, distance)), 
                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                  params_boosted_tree = params_boosted_tree_outcome, 
                                                  params_random_forest = params_random_forest_outcome, 
                                                  params_knn = params_knn_outcome,
                                                  params_lasso = params_lasso_outcome,
                                                  params_glm = params_glm_outcome,
                                                  num_spline_knots = 4, 
                                                  alpha = 0, lambda=0)


simulation_parameters_CRC_GAN <- list(
  oracle_tr_sample_ids = oracle_tr_sample_ids,
  oracle_sample_ids = oracle_sample_ids,
  tr_sample_ids = tr_sample_ids,
  sim_data_tr = sim_data_tr,
  K = K,
  degree_of_interactions = degree_of_interactions,
  polynomial_degree = polynomial_degree,
  knots = knots,
  boundary_knots = boundary_knots,
  sigma_oracle_eff1 = sigma_oracle_eff1,
  sigma_oracle_eff0 = sigma_oracle_eff0,
  sigma_tr_eff1 = sigma_tr_eff1,
  sigma_tr_eff0 = sigma_tr_eff0,
  params_boosted_tree_trt = params_boosted_tree_trt,
  params_random_forest_trt = params_random_forest_trt,
  params_knn_trt = params_knn_trt,
  params_lasso_trt = params_lasso_trt,
  params_glm_trt = params_glm_trt,
  meta_model_trt = meta_model_trt,
  params_boosted_tree_outcome = params_boosted_tree_outcome,
  params_random_forest_outcome = params_random_forest_outcome,
  params_knn_outcome = params_knn_outcome,
  params_lasso_outcome = params_lasso_outcome,
  params_glm_outcome = params_glm_outcome,
  meta_model_outcome = meta_model_outcome
)

save(simulation_parameters_CRC_GAN, file=paste0("simulation_parameters_CRC_GAN_", params_save_name, ".RData"))

rm(n_train, n_oracle, K, polynomial_degree, degree_of_interactions, knots, boundary_knots, params_save_name,
   df_generated, sim_data_oracle, sim_data_tr, 
   oracle_tr_sample_ids, oracle_sample_ids, tr_sample_ids, 
   sigma_oracle_eff1, sigma_oracle_eff0, 
   sigma_tr_eff1, sigma_tr_eff0,
   basis_oracle, basis_tr,
   params_boosted_tree_trt, params_random_forest_trt, params_knn_trt, params_lasso_trt, params_glm_trt, meta_model_trt,
   params_boosted_tree_outcome, params_random_forest_outcome, params_knn_outcome, params_lasso_outcome, params_glm_outcome, meta_model_outcome)