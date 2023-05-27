source('./src/dependencies.R')
source('./src/estimation_functions.R')

# parameters for sequence of IF22

n_train <- 5000
n_oracle <- 1000000
degree_of_interactions <- c(1,2,3,4,5,6,7,8)
params_save_name <- '1,2,3,4,5,6,7,8'

set.seed(202848)

# parameters for true function

key_trt <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_trt=sample(0.05*1:nrow(.), replace=F))}

key_outcome <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_outcome=sample(0.05*1:nrow(.), replace=F))}

##########################################################
## simulate oracle data
##########################################################

sim_data_oracle <- data.frame(
  Z1 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z2 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z4 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z6 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z8 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z9 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z10 = as.numeric(rbernoulli(n_oracle,p=0.5))
) %>% 
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    X = rbinom(n(), 1, val_trt),
    Y = rbinom(n(), 1, val_outcome)
  ) %>% dplyr::select(-c(val_trt, val_outcome))

basis_oracle <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_oracle %>% dplyr::select(-c(X,Y)),
    binary_vars = names(sim_data_oracle %>% dplyr::select(-c(X,Y))),
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_oracle_eff1 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = sim_data_oracle$X))

sigma_oracle_eff0 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = 1-sim_data_oracle$X))

##########################################################
## simulate training data
##########################################################

sim_data_tr <- data.frame(
  Z1 = as.numeric(rbernoulli(n_train,p=0.5)), Z2 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_train,p=0.5)), Z4 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_train,p=0.5)), Z6 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_train,p=0.5)), Z8 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z9 = as.numeric(rbernoulli(n_train,p=0.5)), Z10 = as.numeric(rbernoulli(n_train,p=0.5))
) %>% 
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    X = rbinom(n(), 1, val_trt),
    Y = rbinom(n(), 1, val_outcome)
  ) %>% dplyr::select(-c(val_trt, val_outcome))

basis_tr <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_tr %>% dplyr::select(-c(X,Y)),
    binary_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_tr_eff1 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = sim_data_tr$X),
                   'nlshrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     nlshrink_cov(sim_data_tr$X*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.unequal(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.equal(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.identity(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }))

sigma_tr_eff0 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = 1-sim_data_tr$X),
                   'nlshrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     nlshrink_cov((1-sim_data_tr$X)*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.unequal(t((1-sim_data_tr$X)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.equal(t((1-sim_data_tr$X)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.identity(t((1-sim_data_tr$X)*basis_tr[[i]]), centered=T)$Sigmahat
                   }))

##########################################################
## estimate nuisance parameter models using training data
##########################################################

params_boosted_tree_trt <- find_params_boosted_tree_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                          label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                          nfold = 5,
                                                          tree_depth = c(1:5),
                                                          shrinkage_factor = seq(0.01,0.05,0.01),
                                                          num_trees = seq(250,500,50),
                                                          num_cores = 10)

params_random_forest_trt <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                                            label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                                            nfold = 5,
                                                            num_trees = seq(500,1000,50),
                                                            num_vars = seq(1,5,1),
                                                            num_cores = 10)

params_knn_trt <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                  nfold = 5,
                                  k = seq(11,101,2),
                                  num_cores = 10)

params_lasso_trt <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                      label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                      binary_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
                                      degree_of_interactions = 3,
                                      nfold = 5)

params_glm_trt <- find_params_glm(binary_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))))

meta_model_trt <- fit_stacked_classifer_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)), 
                                              label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                              params_boosted_tree = params_boosted_tree_trt, 
                                              params_random_forest = params_random_forest_trt, 
                                              params_knn = params_knn_trt,
                                              params_lasso = params_lasso_trt,
                                              params_glm = params_glm_trt,
                                              num_spline_knots = 4, 
                                              alpha = 0, lambda=0)

params_boosted_tree_outcome <- find_params_boosted_tree_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                              label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                              nfold = 5,
                                                              tree_depth = c(1:5),
                                                              shrinkage_factor = seq(0.01,0.05,0.01),
                                                              num_trees = seq(250,500,50),
                                                              num_cores = 10)

params_random_forest_outcome <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                nfold = 5,
                                                                num_trees = seq(500,1000,50),
                                                                num_vars = seq(1,5,1),
                                                                num_cores = 10)

params_knn_outcome <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                      label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                      nfold = 5,
                                      k = seq(11,101,2),
                                      num_cores = 10)

params_lasso_outcome <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                          binary_vars = names(sim_data_tr %>% dplyr::select(-c(Y))),
                                          degree_of_interactions = 3, 
                                          nfold = 5)

params_glm_outcome <- find_params_glm(binary_vars = names(sim_data_tr %>% dplyr::select(-c(Y))))

meta_model_outcome <- fit_stacked_classifer_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)), 
                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                  params_boosted_tree = params_boosted_tree_outcome, 
                                                  params_random_forest = params_random_forest_outcome, 
                                                  params_knn = params_knn_outcome,
                                                  params_lasso = params_lasso_outcome,
                                                  params_glm = params_glm_outcome,
                                                  num_spline_knots = 4, 
                                                  alpha = 0, lambda=0)

simulation_parameters_binary <- list(
  n_oracle = n_oracle,
  key_trt = key_trt,
  key_outcome = key_outcome,
  degree_of_interactions = degree_of_interactions,
  sigma_oracle_eff1 = sigma_oracle_eff1,
  sigma_oracle_eff0 = sigma_oracle_eff0,
  sim_data_tr = sim_data_tr,
  sigma_tr_eff1 = sigma_tr_eff1,
  sigma_tr_eff0 = sigma_tr_eff0,
  params_boosted_tree_trt = params_boosted_tree_trt,
  params_random_forest_trt = params_random_forest_trt,
  params_knn_trt = params_knn_trt,
  params_glm_trt = params_glm_trt,
  params_lasso_trt = params_lasso_trt,
  meta_model_trt = meta_model_trt,
  params_boosted_tree_outcome = params_boosted_tree_outcome,
  params_random_forest_outcome = params_random_forest_outcome,
  params_knn_outcome = params_knn_outcome,
  params_lasso_outcome = params_lasso_outcome,
  params_glm_outcome = params_glm_outcome,
  meta_model_outcome = meta_model_outcome
)

save(simulation_parameters_binary, file=paste0("./params/simulation_parameters_binary_", params_save_name, ".RData"))

rm(n_train, n_oracle, degree_of_interactions, key_trt, key_outcome, params_save_name,
   sim_data_oracle, sim_data_tr,
   sigma_oracle_eff1, sigma_oracle_eff0, 
   sigma_tr_eff1, sigma_tr_eff0,
   basis_oracle, basis_tr,
   params_boosted_tree_trt, params_random_forest_trt, params_knn_trt, params_lasso_trt, params_glm_trt, meta_model_trt,
   params_boosted_tree_outcome, params_random_forest_outcome, params_knn_outcome, params_lasso_outcome, params_glm_outcome, meta_model_outcome)