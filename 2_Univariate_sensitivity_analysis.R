## =========================================================
## Script 2: Univariate sensitivity analysis
## Requires Script 1
## =========================================================

## -------------------------------
## 0. Required packages
## -------------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
p_load(deSolve, dplyr, ggplot2)

## IMPORTANT: Run Script 1 first.
## This script assumes these already exist in the environment:
## - initial_conditions
## - times
## - comp_model
## - run_model(params)
## - optimize_delta(params) 
## - parameters            


## -------------------------------
## 1. Baseline 
## -------------------------------
baseline_params <- parameters
baseline_params["delta"] <- optimize_delta(baseline_params)

baseline_out  <- run_model(baseline_params)
baseline_last <- baseline_out[nrow(baseline_out), ]

cat("Baseline calibrated delta:", as.numeric(baseline_params["delta"]), "\n")
cat("Baseline diagnosed @ day 570:", as.numeric(baseline_last$dDiagnosed_CDI), "\n")


## -------------------------------
## 2. R number formula + calc_R()
## -------------------------------
R_formula <- expression(
  (d*l*(a + p_1*z)*
     (p_3*p_4^2*x + p_4^2*v*x + e*f_1*p_3*v + e*h_1*p_3*v +
        f_1*h_1*p_3*x + e*p_3*p_4*v + f_1*h_1*v*x +
        f_1*p_3*p_4*x + h_1*p_3*p_4*x + f_1*p_4*v*x +
        h_1*p_4*v*x - e*p_4^2*v*x - e*f_1*g_1*p_3*v -
        e*f_1*p_4*v*x - e*h_1*p_4*v*x - e*f_1*h_1*s_1*v*x)
  )/(p_2*p_3*(a + p_1)*(f_1 + p_4)*(h_1 + p_4)*(p_3 + v))
)

calc_R <- function(params) {
  envR <- list(
    d   = as.numeric(params["delta"]),
    l   = as.numeric(params["lambda"]),
    a   = as.numeric(params["alfa"]),
    x   = as.numeric(params["equis"]),
    e   = as.numeric(params["epsilon"]),
    v   = as.numeric(params["ve"]),
    f_1 = as.numeric(params["efe_1"]),
    h_1 = as.numeric(params["ache_1"]),
    g_1 = as.numeric(params["gamma_1"]),
    s_1 = as.numeric(params["sigma_1"]),
    p_1 = as.numeric(params["psi_1"]),
    p_2 = as.numeric(params["psi_2"]),
    p_3 = as.numeric(params["psi_3"]),
    p_4 = as.numeric(params["psi_4"]),
    z   = as.numeric(params["zeta"])
  )
  as.numeric(eval(R_formula, envir = envR))
}


## -------------------------------
## 3. Univariate grids (symbol -> parameter name)
## -------------------------------
## NOTE: names in param_map match those in `parameters` from Script 1.
param_map <- c(
  h_1 = "ache_1",
  n   = "ene",
  m   = "eme",
  z   = "zeta",
  x   = "equis",
  e   = "epsilon",
  v   = "ve",
  s_1 = "sigma_1",
  f_1 = "efe_1",
  g_1 = "gamma_1",
  a   = "alfa",
  p_1 = "psi_1",
  p_2 = "psi_2",
  p_3 = "psi_3",
  p_4 = "psi_4",
  l   = "lambda",
  p   = "pi"
)

# Wide-range sensitivity:
grid_list <- list(
  e   = c(0.04, 0.05,0.06,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,0.23),
  f_1 = c(1/2,1/1.9,1/1.8,1/1.7,1/1.6,1/1.5,1/1.4,1/1.3,1/1.2,1/1.1,1),
  g_1 = seq(0.45, 0.95, 0.05),
  h_1 = 1/(5:15),
  m   = seq(0.07, 0.17, 0.01),
  n   = seq(0.02, 0.13, 0.01),
  s_1 = seq(0.45, 0.95, 0.05),
  v   = 1/(11:1),
  x   = seq(0.45, 0.95, 0.05),
  z   = seq(0.02, 0.42, 0.04)
)

# Narrow-range sensitivity: ±20% in 4% increments
pct <- seq(-0.20, 0.20, by = 0.04)
for (sym in c("a","p_1","p_2","p_3","p_4","l","p")) {
  pname <- unname(param_map[sym])
  grid_list[[sym]] <- as.numeric(baseline_params[pname]) * (1 + pct)
}

## -------------------------------
## 4. Functions
## -------------------------------
calc_amplification_index <- function(last_row, params) {
  (as.numeric(params["psi_3"]) *
     (last_row$E + last_row$C + last_row$K_1 + last_row$K_2)) /
    (as.numeric(params["lambda"]) * as.numeric(params["ene"]))
}

run_univariate_one_param <- function(sym) {
  
  pname <- unname(param_map[sym])
  grid  <- grid_list[[sym]]
  base  <- as.numeric(baseline_params[pname])
  
  do.call(rbind, lapply(grid, function(val) {
    
    # 1) Change one parameter
    p <- baseline_params
    p[pname] <- val
    
    # 2) Recalibrate delta
    p["delta"] <- optimize_delta(p)
    
    # 3) Run model
    out  <- run_model(p)
    last <- out[nrow(out), ]
    
    # 4) Outcomes
    data.frame(
      param_symbol = sym,
      param_name   = pname,
      baseline     = base,
      param_value  = val,
      pct_change   = 100 * (val / base - 1),
      delta        = as.numeric(p["delta"]),
      R_number     = calc_R(p),
      Index_equ    = calc_amplification_index(last, p),
      Diagnosed570 = as.numeric(last$dDiagnosed_CDI),
      stringsAsFactors = FALSE
    )
  }))
}

## -------------------------------
## 5. Run analysis + summary table
## -------------------------------
symbols <- names(grid_list)
results_univariate <- do.call(rbind, lapply(symbols, run_univariate_one_param))

table_results <- results_univariate %>%
  transmute(
    Parameter = param_name,
    `Parameter value` = param_value,
    R = R_number,
    `Amplification index` = Index_equ
  ) %>%
  arrange(Parameter, `Parameter value`)
