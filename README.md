# HOIF Estimator Computation Update
---

### **Update - 2026-01-30** 📢
A fatal typo in the eHOIF implementation inside [**HOIF_function/HOIF.R**](HOIF_function/HOIF.R) has been fixed and rigorously verified via full enumeration.

Function signature remains **unchanged**.  
The only new feature is an optional random seed for sample splitting (default = 42).

For better performance, elegant API, arbitrary-order support (including >6), and faster execution, please use the new dedicated R package:  
[**HOIF**](https://github.com/cxy0714/HOIF) → `devtools::install_github("cxy0714/HOIF")`

Oh god, what kind of shit code did I write a year ago... 

To transfer to [**HOIF**](https://github.com/cxy0714/HOIF), see the below example:

```r
################## Example ######################
# In this example, we estimate the bias of DML/AIPW for the potential outcomes Y(a = 1) and Y(a = 0).
# The outcome model assumes a linear relationship: Y = A + beta * X.
# The propensity score model follows a logistic regression: A = psi(alpha * X).
#
# - compute_HOIF_general_all_U(): Provides the exact formula for HOIF, but only up to the 5th order.
# - compute_HOIF_general_part_U(): Computes an approximate formula for HOIF, extendable to any order.
#       - in their return results, "_HOIF_" are the used estimator, "_IF_" or "_U_" are just middle terms.
# Note: The basis function of X is the identity in this example. For practical applications,
#       it is recommended to apply transformations (e.g., B-splines, Fourier basis, etc.) to X.
################################################

set.seed(123)


n <- 100
p <- 10

mu_x <- 0
alpha <- rnorm(p, mean = 0, sd = 0.1)
beta <- rnorm(p, mean = 0, sd = 0.1)

X <- matrix(rnorm(n * p, mu_x, sd = 1), nrow = n, ncol = p)

logit_prob <- X %*% alpha
prob_A <- 1 / (1 + exp(-logit_prob))
A <- as.numeric(rbinom(n, size = 1, prob = prob_A))


Y <- as.numeric(A + X %*% beta + rnorm(n, mean = 0, sd = 0.2))


propensity_model <- glm(A ~ X, family = binomial(link = "logit"))
propensity_scores <- predict(propensity_model, type = "response")
propensity_scores <- as.numeric(propensity_scores)


Y_model <- lm(Y ~ A + X)


Y_pred_1 <- predict(Y_model, newdata = data.frame(A = 1, X = X))
Y_pred_0 <- predict(Y_model, newdata = data.frame(A = 0, X = X))


summary(propensity_model)
summary(Y_model)

epsilon_A_1 <- A / propensity_scores - 1
epsilon_A_0 <- (1 - A) / (1 - propensity_scores) - 1
epsilon_Y_1 <- Y - Y_pred_1
epsilon_Y_0 <- Y - Y_pred_0

m <- 6

devtools::load_all()



# old implementation for sHOIF compute_HOIF_general_all_U(), note here a minus in epsilon_A to match the library(HOIF)
hoif_1_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = -epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = 6,
  Split = 0
)
hoif_0_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = -epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = 6,
  Split = 0
)
# new  implementation for sHOIF using the library(HOIF)

library(HOIF)
library(ustats)
library(reticulate)
py_config()

# setup_ustats()

check_ustats_setup()

shoif_new <- hoif_ate(
  X = X,
  A = A,
  Y = Y,
  mu1 = Y_pred_1,
  mu0 = Y_pred_0,
  pi = propensity_scores,
  m = 6,
  transform_method = "none",
  basis_dim = p,
  sample_split = FALSE
)

# double check
shoif_old_1 <- c(
  hoif_1_all_shoif$sHOIF_HOIF2,
  hoif_1_all_shoif$sHOIF_HOIF3,
  hoif_1_all_shoif$sHOIF_HOIF4,
  hoif_1_all_shoif$sHOIF_HOIF5,
  hoif_1_all_shoif$sHOIF_HOIF6
)
shoif_new$HOIF1
shoif_old_1
all.equal(shoif_new$HOIF1, shoif_old_1)

shoif_old_0 <- c(
  hoif_0_all_shoif$sHOIF_HOIF2,
  hoif_0_all_shoif$sHOIF_HOIF3,
  hoif_0_all_shoif$sHOIF_HOIF4,
  hoif_0_all_shoif$sHOIF_HOIF5,
  hoif_0_all_shoif$sHOIF_HOIF6
)
shoif_new$HOIF0
shoif_old_0
all.equal(shoif_new$HOIF0, shoif_old_0)

# old implementation for eHOIF compute_HOIF_general_all_U(), note here a minus in epsilon_A to match the library(HOIF)

hoif_1_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = - epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = 6,
  Split= 1
)

hoif_0_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = - epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = 6,
  Split= 1
)
# new  implementation for eHOIF using the library(HOIF)

ehoif_new <- hoif_ate(
  X = X,
  A = A,
  Y = Y,
  mu1 = Y_pred_1,
  mu0 = Y_pred_0,
  pi = propensity_scores,
  m = 6,
  transform_method = "none",
  basis_dim = p,
  sample_split = TRUE,
  n_folds = 2,
  seed = 42
)

# double check
ehoif_old_1 <- c(
  hoif_1_all_ehoif$eHOIF_HOIF2,
  hoif_1_all_ehoif$eHOIF_HOIF3,
  hoif_1_all_ehoif$eHOIF_HOIF4,
  hoif_1_all_ehoif$eHOIF_HOIF5,
  hoif_1_all_ehoif$eHOIF_HOIF6
)
ehoif_new$HOIF1
ehoif_old_1
all.equal(ehoif_new$HOIF1, ehoif_old_1)

ehoif_old_0 <- c(
  hoif_0_all_ehoif$eHOIF_HOIF2,
  hoif_0_all_ehoif$eHOIF_HOIF3,
  hoif_0_all_ehoif$eHOIF_HOIF4,
  hoif_0_all_ehoif$eHOIF_HOIF5,
  hoif_0_all_ehoif$eHOIF_HOIF6
)
ehoif_new$HOIF0
ehoif_old_0
all.equal(ehoif_new$HOIF0, ehoif_old_0)
```
---

### **Update - 2024-12-19** 📢

The higher-order (now up to **the sixth** order) computation code for the HOIF (sHOIF/eHOIF) estimator has been updated and is now available in the [**HOIF_function/HOIF.R**](HOIF_function/HOIF.R). 

---

---

### **Update - 2024-12-16** 📢

The higher-order (now up to **the fifth** order) computation code for the HOIF (sHOIF/eHOIF) estimator has been updated and is now available in the [**HOIF_function/HOIF.R**](HOIF_function/HOIF.R). 

---



# Falsification using higher order influence functions for double machine learning estimators of causal effects

# Introduction
Here we provide the code to reproduce the analysis described in: 

# Organization
- `create_simulation_parameters_continuous_X.R` — simulates training and oracle samples, fits nuisance parameters models in the training sample, and computes inverse Gram matrices for the simulation with continuous X.
- `continuous_X_sim.R` — simulates estimation samples, computes estimates of counterfactual means, and computes second order bias estimates for the simulation with continuous X.
- `continuous_X_sim_choose_stopping_k.R` — uses the approach described in the supplement to choose k for simulation with continuous X.
- `create_simulation_parameters_binary_X.R` — simulates training and oracle samples, fits nuisance parameters models in the training sample, and computes inverse Gram matrices for the simulation with binary X.
- `binary_X_sim.R` — simulates estimation samples, computes estimates of counterfactual means, and computes second order bias estimates for the simulation with binary X.
- `compute_bias_all_sims.R` — computes the conditional bias, and the Cauchy-Schwarz bias for simulations with binary and continuous X.
- `create_simulation_parameters_CRC_GAN_data.R` — simulates training and oracle samples, fits nuisance parameters models in the training sample, and computes inverse Gram matrices for the simulation with artificial data based on the National Cancer Database.
- `CRC_GAN_sim.R` — simulates estimation samples, computes estimates of counterfactual means, and computes second order bias estimates for the simulation with artificial data based on the National Cancer Database.
- `create_parameters_CRC_real_data.R` — sample splits, fits nuisance parameters models, and computes inverse Gram matrices for analysis of the National Cancer Database.
- `CRC_real_data_analysis.R` — computes cross-fit estimates of counterfactual means (and the risk difference), and estimates of the second order bias.
- `CRC_real_data_analysis_choose_stopping_k.R` — uses the approach described in the supplement to choose k for analysis of the National Cancer Database.
- `process_outputs.R` — R file which reproduces the tables and figures displayed in the main text, summarizing simulation results. 
- `process_outputs_supplement.R` — R file which reproduces the tables and figures displayed in the supplement, summarizing simulation results.
- `src`  — Folder containing scripts with functions and dependencies to be called by the above files.
- `params`  — Folder should be created to which simulation parameters, which are outputs of `create_simulation_parameters...`, are saved.
- `output`  — Folder to which simulation outputs are saved.
- `figures`  — Folder containing figures, which are outputs of `process_outputs.R` and `process_outputs_supplement.R`.

# Correspondence
If you have any questions, comments, or discover an error, please contact Kerollos Wanis at knwanis@gmail.com.

## References

<sup>1</sup>  Chen, X., Zhang, R., & Liu, L. (2025). *On computing and the complexity of computing higher-order U-statistics, exactly*. arXiv:2508.12627.
