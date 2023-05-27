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
- `params`  — Folder to which simulation parameters, which are outputs of `create_simulation_parameters...`, are saved.
- `output`  — Folder to which simulation outputs are saved.
- `figures`  — Folder containing figures, which are outputs of `process_outputs.R` and `process_outputs_supplement.R`.

# Correspondence
If you have any questions, comments, or discover an error, please contact Kerollos Wanis at knwanis@gmail.com.
