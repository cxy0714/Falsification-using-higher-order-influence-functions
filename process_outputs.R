source('./src/dependencies.R')
options(scipen=999)

###################################################################

load("./output/continuous_sim_5000_200_2,8,32,64,128.RData");load("./output/bias_continuous.RData");load("./output/cs_bias_continuous.RData")

###################################################################

# Table 1

n <- 5000

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  EY1 = c(
    mean(unlist(sim_out['Y_1_stacked',])),
    mean(unlist(sim_out['Y_1_boosted_tree',])),
    mean(unlist(sim_out['Y_1_random_forest',])),
    mean(unlist(sim_out['Y_1_knn',])),
    mean(unlist(sim_out['Y_1_lasso',])),
    mean(unlist(sim_out['Y_1_glm',]))
  ),
  seY1 = c(
    mean(sqrt(unlist(sim_out['var_Y_1_stacked',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_boosted_tree',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_random_forest',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_knn',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_lasso',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_glm',])/n))
  ),
  BiasEY1 = c(
    bias_out$bias_stacked$Y_1,
    bias_out$bias_boosted_tree$Y_1,
    bias_out$bias_random_forest$Y_1,
    bias_out$bias_knn$Y_1,
    bias_out$bias_lasso$Y_1,
    bias_out$bias_glm$Y_1
  ),
  CSBiasa1 = c(
    cs_bias_out$cs_bias_stacked$Y_1,
    cs_bias_out$cs_bias_boosted_tree$Y_1,
    cs_bias_out$cs_bias_random_forest$Y_1,
    cs_bias_out$cs_bias_knn$Y_1,
    cs_bias_out$cs_bias_lasso$Y_1,
    cs_bias_out$cs_bias_glm$Y_1
  ),
  EY0 = c(
    mean(unlist(sim_out['Y_0_stacked',])),
    mean(unlist(sim_out['Y_0_boosted_tree',])),
    mean(unlist(sim_out['Y_0_random_forest',])),
    mean(unlist(sim_out['Y_0_knn',])),
    mean(unlist(sim_out['Y_0_lasso',])),
    mean(unlist(sim_out['Y_0_glm',]))
  ),
  seY0 = c(
    mean(sqrt(unlist(sim_out['var_Y_0_stacked',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_boosted_tree',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_random_forest',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_knn',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_lasso',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_glm',])/n))
  ),
  BiasEY0 = c(
    bias_out$bias_stacked$Y_0,
    bias_out$bias_boosted_tree$Y_0,
    bias_out$bias_random_forest$Y_0,
    bias_out$bias_knn$Y_0,
    bias_out$bias_lasso$Y_0,
    bias_out$bias_glm$Y_0
  ),
  CSBiasa0 = c(
    cs_bias_out$cs_bias_stacked$Y_0,
    cs_bias_out$cs_bias_boosted_tree$Y_0,
    cs_bias_out$cs_bias_random_forest$Y_0,
    cs_bias_out$cs_bias_knn$Y_0,
    cs_bias_out$cs_bias_lasso$Y_0,
    cs_bias_out$cs_bias_glm$Y_0
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Table 2

n_sims <- 200

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM',
                'Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  CSBiasa1 = c(
    cs_bias_out$cs_bias_stacked$Y_1,
    cs_bias_out$cs_bias_boosted_tree$Y_1,
    cs_bias_out$cs_bias_random_forest$Y_1,
    cs_bias_out$cs_bias_knn$Y_1,
    cs_bias_out$cs_bias_lasso$Y_1,
    cs_bias_out$cs_bias_glm$Y_1,
    cs_bias_out$cs_bias_stacked$Y_1,
    cs_bias_out$cs_bias_boosted_tree$Y_1,
    cs_bias_out$cs_bias_random_forest$Y_1,
    cs_bias_out$cs_bias_knn$Y_1,
    cs_bias_out$cs_bias_lasso$Y_1,
    cs_bias_out$cs_bias_glm$Y_1
  ),
  IF22oraclea1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle9']))
  ),
  IF22oraclea1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle9'])),')')
  ),
  IF22tra1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr9']))
  ),
  IF22tra1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr9'])),')')
  ),
  IF22nlshrinka1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink9']))
  ),
  IF22nlshrinka1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink9'])),')')
  ),
  IF22esta1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est9']))
  ),
  IF22esta1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est9'])),')')
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Table 3

n_sims <- 200

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM',
                'Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  CSBiasa0 = c(
    cs_bias_out$cs_bias_stacked$Y_0,
    cs_bias_out$cs_bias_boosted_tree$Y_0,
    cs_bias_out$cs_bias_random_forest$Y_0,
    cs_bias_out$cs_bias_knn$Y_0,
    cs_bias_out$cs_bias_lasso$Y_0,
    cs_bias_out$cs_bias_glm$Y_0,
    cs_bias_out$cs_bias_stacked$Y_0,
    cs_bias_out$cs_bias_boosted_tree$Y_0,
    cs_bias_out$cs_bias_random_forest$Y_0,
    cs_bias_out$cs_bias_knn$Y_0,
    cs_bias_out$cs_bias_lasso$Y_0,
    cs_bias_out$cs_bias_glm$Y_0
  ),
  IF22oraclea0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle9']))
  ),
  IF22oraclea0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle9'])),')')
  ),
  IF22tra0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr9']))
  ),
  IF22tra0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr9'])),')')
  ),
  IF22nlshrinka0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink9']))
  ),
  IF22nlshrinka0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink9'])),')')
  ),
  IF22esta0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est4'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est9'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est9']))
  ),
  IF22esta0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est4'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est9'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est9'])),')')
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Figure 2

n_sims <- 200

graph_tr_Y_1 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4) 
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')]))
             ))

graph_shrinkage_Y_1 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')]))
             ))

graph_est_Y_1 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4')]))
             ))

graph_oracle_Y_1 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')]))
             ))

graph_tr_Y_0 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4) 
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')]))
             ))

graph_shrinkage_Y_0 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4')]))
             ))

graph_est_Y_0 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4')]))
             ))

graph_oracle_Y_0 <- 
  data.frame(k=rep(c(21,81,321,641), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4')]))
             ))

col <- c("#CC6666", "#9999CC", "#66CC99")

g2_tr_Y_1_binary <- ggplot(data=graph_tr_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{tr})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.02,0.2)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_shrinkage_Y_1_binary <- ggplot(data=graph_shrinkage_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{shrinkage})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.02,0.2)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_est_Y_1_binary <- ggplot(data=graph_est_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{est})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.02,0.2)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_oracle_Y_1_binary <- ggplot(data=graph_oracle_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{oracle})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.02,0.2)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_tr_Y_0_binary <- ggplot(data=graph_tr_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{tr})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.15,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_shrinkage_Y_0_binary <- ggplot(data=graph_shrinkage_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{shrinkage})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.15,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_est_Y_0_binary <- ggplot(data=graph_est_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{est})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.15,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_oracle_Y_0_binary <- ggplot(data=graph_oracle_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(21,81,321,641)) +
  theme_bw() + ylab(expression(IF[22]^{oracle})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.12,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))


ggarrange(g2_tr_Y_1_binary, g2_tr_Y_0_binary,
          g2_shrinkage_Y_1_binary, g2_shrinkage_Y_0_binary,
          g2_est_Y_1_binary, g2_est_Y_0_binary,
          g2_oracle_Y_1_binary, g2_oracle_Y_0_binary,
          ncol=2, nrow=4, common.legend = TRUE, legend="bottom",
          font.label=list(size=20))

ggsave("./figures/cont_graph_Y_1_and_Y_0.png", dpi=600)


###################################################################

load("./output/binary_sim_5000_200_1,2,3,4,5,6,7,8.RData");load("./output/bias_binary.RData");load("./output/cs_bias_binary.RData")

###################################################################

# Table 4

n <- 5000

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  EY1 = c(
    mean(unlist(sim_out['Y_1_stacked',])),
    mean(unlist(sim_out['Y_1_boosted_tree',])),
    mean(unlist(sim_out['Y_1_random_forest',])),
    mean(unlist(sim_out['Y_1_knn',])),
    mean(unlist(sim_out['Y_1_lasso',])),
    mean(unlist(sim_out['Y_1_glm',]))
  ),
  seY1 = c(
    mean(sqrt(unlist(sim_out['var_Y_1_stacked',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_boosted_tree',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_random_forest',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_knn',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_lasso',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_1_glm',])/n))
  ),
  BiasEY1 = c(
    bias_out$bias_stacked$Y_1,
    bias_out$bias_boosted_tree$Y_1,
    bias_out$bias_random_forest$Y_1,
    bias_out$bias_knn$Y_1,
    bias_out$bias_lasso$Y_1,
    bias_out$bias_glm$Y_1
  ),
  CSBiasa1 = c(
    cs_bias_out$cs_bias_stacked$Y_1,
    cs_bias_out$cs_bias_boosted_tree$Y_1,
    cs_bias_out$cs_bias_random_forest$Y_1,
    cs_bias_out$cs_bias_knn$Y_1,
    cs_bias_out$cs_bias_lasso$Y_1,
    cs_bias_out$cs_bias_glm$Y_1
  ),
  EY0 = c(
    mean(unlist(sim_out['Y_0_stacked',])),
    mean(unlist(sim_out['Y_0_boosted_tree',])),
    mean(unlist(sim_out['Y_0_random_forest',])),
    mean(unlist(sim_out['Y_0_knn',])),
    mean(unlist(sim_out['Y_0_lasso',])),
    mean(unlist(sim_out['Y_0_glm',]))
  ),
  seY0 = c(
    mean(sqrt(unlist(sim_out['var_Y_0_stacked',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_boosted_tree',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_random_forest',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_knn',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_lasso',])/n)),
    mean(sqrt(unlist(sim_out['var_Y_0_glm',])/n))
  ),
  BiasEY0 = c(
    bias_out$bias_stacked$Y_0,
    bias_out$bias_boosted_tree$Y_0,
    bias_out$bias_random_forest$Y_0,
    bias_out$bias_knn$Y_0,
    bias_out$bias_lasso$Y_0,
    bias_out$bias_glm$Y_0
  ),
  CSBiasa0 = c(
    cs_bias_out$cs_bias_stacked$Y_0,
    cs_bias_out$cs_bias_boosted_tree$Y_0,
    cs_bias_out$cs_bias_random_forest$Y_0,
    cs_bias_out$cs_bias_knn$Y_0,
    cs_bias_out$cs_bias_lasso$Y_0,
    cs_bias_out$cs_bias_glm$Y_0
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Table 5

n_sims <- 200

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  CSBiasa1 = c(
    cs_bias_out$cs_bias_stacked$Y_1,
    cs_bias_out$cs_bias_boosted_tree$Y_1,
    cs_bias_out$cs_bias_random_forest$Y_1,
    cs_bias_out$cs_bias_knn$Y_1,
    cs_bias_out$cs_bias_lasso$Y_1,
    cs_bias_out$cs_bias_glm$Y_1
  ),
  IF22oraclea1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle5']))
  ),
  IF22oraclea1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'oracle5'])),')')
  ),
  IF22tra1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr5']))
  ),
  IF22tra1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'tr5'])),')')
  ),
  IF22nlshrinka1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink5']))
  ),
  IF22nlshrinka1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'nlshrink5'])),')')
  ),
  IF22esta1 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est5']))
  ),
  IF22esta1_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff1',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff1',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff1',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,'est5'])),')')
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Table 6 

n_sims <- 200

data.frame(
  Estimator = c('Ensemble', 'Boosted trees', 'Random forest', 'KNN', 'LASSO', 'GLM'),
  CSBiasa0 = c(
    cs_bias_out$cs_bias_stacked$Y_0,
    cs_bias_out$cs_bias_boosted_tree$Y_0,
    cs_bias_out$cs_bias_random_forest$Y_0,
    cs_bias_out$cs_bias_knn$Y_0,
    cs_bias_out$cs_bias_lasso$Y_0,
    cs_bias_out$cs_bias_glm$Y_0
  ),
  IF22oraclea0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle5']))
  ),
  IF22oraclea0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'oracle5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'oracle5'])),')')
  ),
  IF22tra0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr5']))
  ),
  IF22tra0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'tr5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'tr5'])),')')
  ),
  IF22nlshrinka0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink5']))
  ),
  IF22nlshrinka0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'nlshrink5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'nlshrink5'])),')')
  ),
  IF22esta0 = c(
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est5'])),
    mean(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est5']))
  ),
  IF22esta0_sd = c(
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_boosted_tree_eff0',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_random_forest_eff0',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_knn_eff0',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,'est5'])),')'),
    paste0('(',sd(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,'est5'])),')')
  )
) %>% knitr::kable(format='latex', booktabs = T, linesep = "")

# Figure 3

n_sims <- 200

graph_tr_Y_1 <- 
  data.frame(k=rep(c(11,56,176,386), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4) 
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')]))
             ))

graph_shrinkage_Y_1 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')]))
             ))

graph_est_Y_1 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')]))
             ))

graph_oracle_Y_1 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff1',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')]))
             ))

graph_tr_Y_0 <- 
  data.frame(k=rep(c(11,56,176,386), 3),
             Estimator=c(rep('Ensemble', 4), 
                         rep('GLM', 4),
                         rep('LASSO', 4) 
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('tr1', 'tr2', 'tr3', 'tr4')]))
             ))

graph_shrinkage_Y_0 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('nlshrink1', 'nlshrink2', 'nlshrink3', 'nlshrink4', 'nlshrink5')]))
             ))

graph_est_Y_0 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('est1', 'est2', 'est3', 'est4', 'est5')]))
             ))

graph_oracle_Y_0 <- 
  data.frame(k=rep(c(11,56,176,386,638), 3),
             Estimator=c(rep('Ensemble', 5), 
                         rep('GLM', 5),
                         rep('LASSO', 5)
             ),
             Estimate=c(rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_stacked_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_glm_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')])),
                        rowMeans(sapply(1:n_sims, function(x) sim_out['HOIF_lasso_eff0',][[x]][2,c('oracle1', 'oracle2', 'oracle3', 'oracle4', 'oracle5')]))
             ))

col <- c("#CC6666", "#9999CC", "#66CC99")

g2_tr_Y_1_binary <- ggplot(data=graph_tr_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386)) +
  theme_bw() + ylab(expression(IF[22]^{tr})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.05,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_shrinkage_Y_1_binary <- ggplot(data=graph_shrinkage_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{shrinkage})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.05,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_est_Y_1_binary <- ggplot(data=graph_est_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{est})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.05,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_oracle_Y_1_binary <- ggplot(data=graph_oracle_Y_1, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{oracle})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_1, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_1, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_1, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.05,0.03)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_tr_Y_0_binary <- ggplot(data=graph_tr_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386)) +
  theme_bw() + ylab(expression(IF[22]^{tr})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.04,0.015)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_shrinkage_Y_0_binary <- ggplot(data=graph_shrinkage_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{shrinkage})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.04,0.015)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_est_Y_0_binary <- ggplot(data=graph_est_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{est})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.04,0.015)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))

g2_oracle_Y_0_binary <- ggplot(data=graph_oracle_Y_0, aes(x=k, y=Estimate, col=Estimator)) + geom_line(size=1, position = position_dodge(width = 10)) +
  scale_x_continuous(breaks=c(11,56,176,386,638)) +
  theme_bw() + ylab(expression(IF[22]^{oracle})) + 
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  geom_hline(yintercept=bias_out$bias_stacked$Y_0, col=col[1], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_glm$Y_0, col=col[2], linetype="dashed", lwd=1.2) +
  geom_hline(yintercept=bias_out$bias_lasso$Y_0, col=col[3], linetype="dashed", lwd=1.2) +
  ylim(c(-0.04,0.015)) +
  guides(fill=guide_legend(title="Nuisance parameter estimator"),
         col = guide_legend(title="Nuisance parameter estimator"))


ggarrange(g2_tr_Y_1_binary, g2_tr_Y_0_binary,
          g2_shrinkage_Y_1_binary, g2_shrinkage_Y_0_binary,
          g2_est_Y_1_binary, g2_est_Y_0_binary,
          g2_oracle_Y_1_binary, g2_oracle_Y_0_binary,
          ncol=2, nrow=4, common.legend = TRUE, legend="bottom",
          font.label=list(size=20))

ggsave("./figures/binary_graph_Y_1_and_Y_0.png", dpi=600)
