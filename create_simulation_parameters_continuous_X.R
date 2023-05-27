source('./src/dependencies.R')
source('./src/estimation_functions.R')

# parameters for sequence of IF22

n_train <- 5000
n_oracle <- 1000000
K_b_splines <- c(2,8,32,64,128)
K_fourier <- c(2,8,32,64,128)
degree_of_interactions <- c(1,1,1,1,1)
polynomial_degree <- c(1,1,1,1,1)
params_save_name <- '2,8,32,64,128'

set.seed(202848)

# parameters for true b-spline function

theta_X = c(0.15, 0.7, 0.3, 0.9, 0.3, 0.9, 0.6, 0.9, 0.6, 0.2, 0.4, 0.3, 0.2, -0.2, 0.1, -0.1, 0.2, 0.5, 0.1, 0.1, 0.3, 0.2, 
            0.2, -0.1, 0.2, -0.1, 0.5, -0.2, 0.5, 0.1, 0.4, 0.3, 0.9, 0.3, 0.9, 0.3, 0.9, 0.8, 0.1, 0.3, 0.2, 0.1,
            0.1, 0.8, 0.6, -0.2, 0.2, 0.9, 0.5, 0.9, 0.5, 0.1, 0.3, 0.2, -0.2, 0.1, -0.1, 0.3, 0.6, 0.2, -0.1, 0.5, 0.3, 
            0.2, -0.1, 0.2, -0.1, 0.5, 0.9, 0.6, 0.2, 0.8, 0.2, 0.5, 0.2, -0.1, 0.5, 0.9, -0.1, 0.1, 0.3, -0.1, 0.1)
theta_Y = c(0.15, 0.7, 0.3, 0.9, 0.3, 0.9, 0.6, 0.9, 0.6, 0.2, 0.4, 0.3, 0.2, -0.2, 0.1, -0.1, 0.2, 0.5, 0.1, 0.1, 0.3, 0.2, 
            0.2, -0.1, 0.2, -0.1, 0.5, -0.2, 0.5, 0.1, 0.4, 0.3, 0.9, 0.3, 0.9, 0.3, 0.9, 0.8, 0.1, 0.3, 0.2, 0.1,
            0.1, 0.8, 0.6, -0.2, 0.2, 0.9, 0.5, 0.9, 0.5, 0.1, 0.3, 0.2, -0.2, 0.1, -0.1, 0.3, 0.6, 0.2, -0.1, 0.5, 0.3, 
            0.2, -0.1, 0.2, -0.1, 0.5, 0.9, 0.6, 0.2, 0.8, 0.2, 0.5, 0.2, -0.1, 0.5, 0.9, -0.1, 0.1, 0.3, -0.1, 0.1)
true_knots = list(seq(0.0125,0.9875,0.0125))

##########################################################
## simulate oracle data
##########################################################

Z1 = runif(n_oracle, 0, 1); Z2 = runif(n_oracle, 0, 1); Z3 = runif(n_oracle, 0, 1); Z4 = runif(n_oracle, 0, 1); Z5 = runif(n_oracle, 0, 1)
Z6 = runif(n_oracle, 0, 1); Z7 = runif(n_oracle, 0, 1); Z8 = runif(n_oracle, 0, 1); Z9 = runif(n_oracle, 0, 1); Z10 = runif(n_oracle, 0, 1)

sim_data_oracle <- create_b_spline_basis(
  data=data.frame(Z1),
  continuous_vars=c('Z1'),
  knots=true_knots,
  boundary_knots=list(c(0,1)),
  polynomial_degree=3,
  degree_of_interactions=1
)

p_X = sim_data_oracle %*% theta_X
p_Y = sim_data_oracle  %*% theta_Y
X = rbinom(n_oracle, 1, p_X)
Y = rbinom(n_oracle, 1, p_Y)

sim_data_oracle <- data.frame(Z1=Z1, Z2=Z2, Z3=Z3, Z4=Z4, Z5=Z5, Z6=Z6, Z7=Z7, Z8=Z8, Z9=Z9, Z10=Z10, X=X, Y=Y)

knots <- lapply(1:length(K_b_splines), function(k) {
  lapply(1:10, function(z) {
    attr(bs(eval(parse(text=paste0('sim_data_oracle$Z',z))), degree=polynomial_degree[k], df=K_b_splines[k]), "knots")
  })
})

boundary_knots <- lapply(1:length(K_b_splines), function(k) {
  lapply(1:10, function(z) {
    attr(bs(eval(parse(text=paste0('sim_data_oracle$Z',z))), degree=polynomial_degree[k], df=K_b_splines[k]), "Boundary.knots")
  })
})

periods <- lapply(1:10, function(z) {
  max(eval(parse(text=paste0('sim_data_oracle$Z',z)))) - min(eval(parse(text=paste0('sim_data_oracle$Z',z))))
})


basis_oracle <- c(
  'oracle' = lapply(1:length(K_b_splines), function(k) {
    create_b_spline_basis(
      data = sim_data_oracle %>% dplyr::select(-c(X,Y)),
      continuous_vars = names(sim_data_oracle %>% dplyr::select(-c(X,Y))),
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  }),
  lapply(1:length(K_fourier), function(k) {
    create_fourier_basis(
      data = sim_data_oracle %>% dplyr::select(-c(X,Y)),
      continuous_vars = names(sim_data_oracle %>% dplyr::select(-c(X,Y))),
      period = periods[k],
      nbasis = K_fourier[k],
      degree_of_interactions = degree_of_interactions[k]
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

Z1 = runif(n_train, 0, 1); Z2 = runif(n_train, 0, 1); Z3 = runif(n_train, 0, 1); Z4 = runif(n_train, 0, 1); Z5 = runif(n_train, 0, 1)
Z6 = runif(n_train, 0, 1); Z7 = runif(n_train, 0, 1); Z8 = runif(n_train, 0, 1); Z9 = runif(n_train, 0, 1); Z10 = runif(n_train, 0, 1)

sim_data_tr <- create_b_spline_basis(
  data=data.frame(Z1),
  continuous_vars=c('Z1'),
  knots=true_knots,
  boundary_knots=list(c(0,1)),
  polynomial_degree=3,
  degree_of_interactions=1
)

p_X = sim_data_tr %*% theta_X
p_Y = sim_data_tr  %*% theta_Y
X = rbinom(n_train, 1, p_X)
Y = rbinom(n_train, 1, p_Y)

sim_data_tr <- data.frame(Z1=Z1, Z2=Z2, Z3=Z3, Z4=Z4, Z5=Z5, Z6=Z6, Z7=Z7, Z8=Z8, Z9=Z9, Z10=Z10, X=X, Y=Y)

basis_tr <- c(
  lapply(1:length(K_b_splines), function(k) {
    create_b_spline_basis(
      data = sim_data_tr %>% dplyr::select(-c(X,Y)),
      continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
      knots = knots[[k]], boundary_knots = boundary_knots[[k]],
      degree_of_interactions = degree_of_interactions[k],
      polynomial_degree = polynomial_degree[k]
    )
  }),
  lapply(1:length(K_fourier), function(k) {
    create_fourier_basis(
      data = sim_data_tr %>% dplyr::select(-c(X,Y)),
      continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
      period = periods[k],
      nbasis = K_fourier[k],
      degree_of_interactions = degree_of_interactions[k]
    )
  })
)

sigma_tr_eff1 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = sim_data_tr$X),
                   'nlshrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     nlshrink_cov(sim_data_tr$X*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     shrinkcovmat.unequal(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     shrinkcovmat.equal(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     shrinkcovmat.identity(t(sim_data_tr$X*basis_tr[[i]]), centered=T)$Sigmahat
                   }))
                   
sigma_tr_eff0 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = 1-sim_data_tr$X),
                   'nlshrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     nlshrink_cov((1-sim_data_tr$X)*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     shrinkcovmat.unequal(t((1-sim_data_tr$X)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
                     shrinkcovmat.equal(t((1-sim_data_tr$X)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(K_b_splines, K_fourier)), function(i) {
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

params_knn_trt <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                  label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                  nfold = 5,
                                  k = seq(11,101,2),
                                  num_cores = 10)

params_lasso_trt <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, X)),
                                      label_vector = sim_data_tr %>% dplyr::select(X) %>% {.[[1]]},
                                      continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
                                      continuous_var_spline_knots = 10,
                                      degree_of_interactions = 3,
                                      nfold = 5)

params_glm_trt <- find_params_glm(continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
                                  continuous_var_spline_knots = 2)

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
                                          continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
                                          binary_vars = c('X'),
                                          continuous_var_spline_knots = 10,
                                          degree_of_interactions = 3,
                                          nfold = 5)

params_glm_outcome <- find_params_glm(continuous_vars = names(sim_data_tr %>% dplyr::select(-c(X,Y))),
                                      binary_vars = c('X'),
                                      continuous_var_spline_knots = 2)

meta_model_outcome <- fit_stacked_classifer_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)), 
                                                  label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                  params_boosted_tree = params_boosted_tree_outcome, 
                                                  params_random_forest = params_random_forest_outcome, 
                                                  params_knn = params_knn_outcome,
                                                  params_lasso = params_lasso_outcome,
                                                  params_glm = params_glm_outcome,
                                                  num_spline_knots = 4, 
                                                  alpha = 0, lambda=0)

simulation_parameters_continuous <- list(
  n_oracle = n_oracle,
  K_b_splines = K_b_splines,
  K_fourier = K_fourier,
  theta_X = theta_X,
  theta_Y = theta_Y,
  true_knots = true_knots,
  degree_of_interactions = degree_of_interactions,
  polynomial_degree = polynomial_degree,
  periods = periods,
  knots = knots,
  boundary_knots = boundary_knots,
  sigma_oracle_eff1 = sigma_oracle_eff1,
  sigma_oracle_eff0 = sigma_oracle_eff0,
  sim_data_tr = sim_data_tr,
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

save(simulation_parameters_continuous, file=paste0("./params/simulation_parameters_continuous_", params_save_name, ".RData"))

rm(n_train, n_oracle, K_b_splines, K_fourier, polynomial_degree, degree_of_interactions, knots, boundary_knots, periods, params_save_name,
   sim_data_oracle, sim_data_tr, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, p_X, p_Y, X, Y, theta_X, theta_Y, true_knots,
   sigma_oracle_eff1, sigma_oracle_eff0, 
   sigma_tr_eff1, sigma_tr_eff0,
   basis_oracle, basis_tr,
   params_boosted_tree_trt, params_random_forest_trt, params_knn_trt, params_lasso_trt, params_glm_trt, meta_model_trt,
   params_boosted_tree_outcome, params_random_forest_outcome, params_knn_outcome, params_lasso_outcome, params_glm_outcome, meta_model_outcome)