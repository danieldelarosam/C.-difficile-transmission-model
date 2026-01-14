## -------------------------------
## 0. Required Packages
## -------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(dplyr,epiR)

## -------------------------------
## 1. PRCC Analysis for Reproduction Number (R)
## -------------------------------

# Calculate PRCC using model parameters and estimated reproduction number (R)
PRCC_analysis_R <- data.frame(
  delta = parameters_df$d, alfa = parameters_df$a, 
  psi_1 = parameters_df$p_1, psi_2 = parameters_df$p_2, psi_3 = parameters_df$p_3, 
  psi_4 = parameters_df$p_4, ache_1 = parameters_df$h_1, zeta = parameters_df$z,
  equis = parameters_df$x, epsilon = parameters_df$e, ve = parameters_df$v,
  gamma_1 = parameters_df$g_1, sigma_1 = parameters_df$s_1, efe_1 = parameters_df$f_1,
  R_value = parameters_df$R)

prcc_results <- epi.prcc(PRCC_analysis_R, sided.test = 2, conf.level = 0.95) %>%
  as.data.frame() %>%
  mutate(across(c(est, lower, upper), \(x) round(x, 3)))

View(prcc_results)

## -------------------------------
## 2. PRCC Analysis for Colonization Index
## -------------------------------

# Calculate PRCC using model parameters and colonization amplification index
PRCC_colonization <- data.frame(
  delta = model_simulations$delta, alfa = model_simulations$alfa, 
  psi_1 = model_simulations$psi_1, psi_2 = model_simulations$psi_2, psi_3 = model_simulations$psi_3, 
  psi_4 = model_simulations$psi_4, ache_1 = model_simulations$ache_1, zeta = model_simulations$zeta,
  equis = model_simulations$equis, epsilon = model_simulations$epsilon, ve = model_simulations$ve,
  gamma_1 = model_simulations$gamma_1, sigma_1 = model_simulations$sigma_1, efe_1 = model_simulations$efe_1,
  index = model_simulations$Index_equ)


prcc_results_2 <- epi.prcc(PRCC_colonization, sided.test = 2, conf.level = 0.95) %>%
  as.data.frame() %>%
  mutate(across(c(est, lower, upper), round, 3)) 

View(prcc_results_2)
