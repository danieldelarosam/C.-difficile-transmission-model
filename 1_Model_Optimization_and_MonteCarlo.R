## -------------------------------
## 0. Load Required Packages
## -------------------------------

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

# Load additional packages
p_load(deSolve, reshape2, ggplot2, dplyr)

## -------------------------------
## 1. Define Model Structure and Parameters
## -------------------------------

### 1.1 Initial Conditions and Parameters

# Initial conditions: Number of individuals per compartment
initial_conditions <- c(
  R = 121,                    # Non-susceptible population
  S = 334,                    # Susceptible population
  E = 0,                      # Asymptomatic carriers (E type)
  C = 0,                      # Asymptomatic carriers (C type)
  K_1 = 0,                    # Asymptomatic carriers under treatment (K1, with contact precautions)
  K_2 = 0,                    # Asymptomatic carriers under treatment (K2, with contact precautions)                
  I = 0,                      # Symptomatic patients 
  A = 0,                      # Symptomatic patients under treatment with contact precautions
  dAdmission = 0,             # Cumulative number of admissions
  dDiagnosed_CDI=0)           # Cumulative number of diagnosed CDI cases

# Model parameters
parameters <- c(
  delta = 0.0007233824,    # Transmission rate of infected individuals
  equis = 0.7,             # Relative transmission for colonized patients
  alfa = 0.064,            # Antibiotic use rate
  epsilon = 0.13,          # Fraction of asymptomatic carriers developing symptoms during hospitalization
  ve = 1/4,                # Progression rate to symptomatic disease
  lambda = 76,             # Admission rate to the hospital
  ene = 0.08,              # Overall proportion of asymptomatic carriers at admission.
  pi = 0.004986150,        # Proportion of infected patients at admission
  eme = 0.12,              # Proportion of pre-symptomatic patients at admission
  zeta = 0.22,             # Proportion of susceptible patients at admission
  psi_1 = 1/6,             # Discharge rate for non-colonized population.
  psi_2 = 1/6,             # Discharge rate for susceptible population.
  psi_3 = 1/6,             # Discharge rate for asymptomatic carriers.
  psi_4 = 1/12,            # Discharge rate for symptomatic population.
  ache_1 = 1/10,           # Clearance rate due to treatment for infected patients
  ache_2 = 1/10,           # Clearance rate due to treatment for asymptomatic patients
  gamma_1 = 0.5,           # Reduction in transmission due to contact precautions for symptomatic patients
  gamma_2 = 0,             # Reduction in transmission due to contact precautions for asymptomatic patients
  efe_2 = 0,               # Proportion of screened individuals
  efe_1 = 1/1.5,           # Diagnosis rate of symptomatic carriers.
  sigma_1 = 0.7,           # Fraction cured among symptomatic carriers after treatment
  sigma_2 = 0              # Fraction cured among asymptomatic carriers after treatment
)

# Follow-up period and time steps
follow_up_duration = 570 #Study period in days.
times <- seq(from = 0, to = follow_up_duration, by = 1) # Time steps (days)


### 1.2 Model Structure: Compartmental ODEs
comp_model <- function(time, state, parameters) {  
  with(as.list(c(state, parameters)), { 
    
    # Proportions of patients from community
    kappa <- 1-eme                        # Proportion entering compartment C
    omega <- 1-(zeta+pi+ene)              # Proportion entering compartment S
    
    # Define the force of infection (beta1,beta2, beta3 and beta4) for each group capable of transmitting the infection, 
    B1 = delta                          # Symptomatic patients
    B2 = delta * (1-gamma_1)            # Symptomatic patients with contact precautions
    B3 = delta * equis                  # Asymptomatic carriers
    B4 = delta * equis * (1-gamma_2)    # Asymptomatic carriers with contact precautions
    
    # Differential equations: rate of change per compartment
    dR= lambda*omega - psi_1*R - alfa*R
    dS= alfa*R - (B1*I+B2*A+B3*E+B3*C+B4*K_1+B4*K_2)*S - psi_2*S + sigma_2*ache_2*K_1 + sigma_2*ache_2*K_2 + sigma_1*ache_1*A + lambda*zeta
  
    dE= epsilon*(B1*I+B2*A+B3*E+B3*C+B4*K_1+B4*K_2)*S + (1-efe_2)*lambda*ene*eme - ve*E - psi_3*E
    dC= (1-epsilon)*(B1*I+B2*A+B3*E+B3*C+B4*K_1+B4*K_2)*S + (1-efe_2)*lambda*ene*kappa - psi_3*C + (1-sigma_1)*ache_1*A
    
    dK_1= efe_2*lambda*ene*eme - sigma_2*ache_2*K_1 - psi_3*K_1 - ve*K_1
    dK_2= efe_2*lambda*ene*kappa - sigma_2*ache_2*K_2 - psi_3*K_2 
    
    dI= ve*E + lambda*pi - psi_4*I - efe_1*I 
    dA= efe_1*I - (1-sigma_1)*ache_1*A - sigma_1*ache_1*A - psi_4*A + ve*K_1 
    
    dAdmission= lambda
    dDiagnosed_CDI=efe_1*I + ve*K_1

    return(list(c(dR, dS, dE, dC, dK_1,dK_2, dI, dA, dAdmission, dDiagnosed_CDI))) 
  })
}

## -------------------------------
## 2. Run Baseline Model Simulation
## -------------------------------

run_model <- function(params) {
  as.data.frame(ode(
    y = initial_conditions, 
    times = times, 
    func = comp_model, 
    parms = params
  ))
}

model_output <- run_model(parameters)

View(model_output)


## -------------------------------
## 3. Optimize Transmission Rate (delta)
## -------------------------------

# Function to optimize the transmission rate (delta)
optimize_delta <- function(params, target = 639) {
  
  delta_objective <- function(delta) {
    params["delta"] <- delta
    result <- run_model(params)
    final_incidence <- tail(result$dDiagnosed_CDI, 1)
    abs(target - final_incidence)
  }
  
  optimize(f = delta_objective, interval = c(0, 1), tol = 1e-10)$minimum
}

# Optimize 'delta' for baseline parameters
parameters["delta"] <- optimize_delta(parameters)
model_output <- run_model(parameters)
View(model_output)

cat("Optimized delta:", parameters["delta"], "\n")


## -------------------------------
## 4. Monte Carlo Simulations
## -------------------------------

### 4.1 Generate Parameter Sets
uniform_values <- function(mean, n) runif(n, mean * 0.8, mean * 1.2)
fixed_values <- function(value, n) rep(value, n)

num_simulations <- 1000

parameters_list <- list( # List of random variables for simulations. 
  delta = uniform_values(0.1, num_simulations), # Any initial value; delta will be optimized within each simulation.
  equis = uniform_values(0.7, num_simulations),
  alfa = uniform_values(0.064, num_simulations),
  epsilon = uniform_values(0.13, num_simulations),
  ve = uniform_values(1/4, num_simulations),
  lambda = fixed_values(76, num_simulations),
  ene = uniform_values(0.08, num_simulations),
  pi = fixed_values(0.004986150, num_simulations),
  eme = uniform_values(0.12, num_simulations),
  zeta = uniform_values(0.22, num_simulations),
  psi_1 = uniform_values(1/6, num_simulations),
  psi_2 = uniform_values(1/6, num_simulations),
  psi_3 = uniform_values(1/6, num_simulations),
  psi_4 = uniform_values(1/12, num_simulations),
  ache_1 = uniform_values(1/10, num_simulations),
  ache_2 = uniform_values(0, num_simulations),
  gamma_1 = uniform_values(0.5, num_simulations),
  gamma_2 = uniform_values(0, num_simulations),
  efe_2 = uniform_values(0, num_simulations),
  efe_1 = uniform_values(1/1.5, num_simulations),
  sigma_1 = uniform_values(0.7, num_simulations),
  sigma_2 = uniform_values(0, num_simulations)
)

### 4.2 Run Simulations with Optimization
sim_results <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  current_parameters <- lapply(parameters_list, function(x) x[i])
  current_parameters["delta"] <- optimize_delta(current_parameters)
  optimized_output <- run_model(current_parameters)
  final_step <- optimized_output[nrow(optimized_output), ]
  
  sim_results[[i]] <- c(i, final_step["dAdmission"], final_step["dDiagnosed_CDI"],
                        final_step["E"], final_step["C"], final_step["K_1"], final_step["K_2"],
                        as.numeric(unlist(current_parameters)))
}



model_simulations <- as.data.frame(do.call(rbind, sim_results))
colnames(model_simulations)[1:7] <- c("simulation", "dAdmission", "dDiagnosed_CDI",
                                      "E_equi", "C_equi", "K1_equi", "K2_equi")
colnames(model_simulations)[8:ncol(model_simulations)] <- names(parameters_list)

model_simulations[] <- lapply(model_simulations, function(x) as.numeric(as.character(x)))

model_simulations <- model_simulations %>%
  mutate(
    Incidence = round((dDiagnosed_CDI / dAdmission) * 1000, 2),
    Index_equ = round((psi_3 * (E_equi + C_equi + K1_equi + K2_equi)) / (lambda * ene), 1)
  )

View(model_simulations)


## -------------------------------
## 5. Summarize and Visualize Outputs
## -------------------------------

summary_and_plot <- function(data, column, main_title, xlab, ylab) {
  cat("Median:", median(data[[column]], na.rm = TRUE), "\n")
  cat("Mean:", mean(data[[column]], na.rm = TRUE), "\n")
  cat("IQR:", quantile(data[[column]], probs = c(0.25, 0.75), na.rm = TRUE), "\n")
  cat("SD:", sd(data[[column]], na.rm = TRUE), "\n")
  
  freq <- table(round(data[[column]], 1))
  
  barplot(freq, main = paste("Distribution of", main_title), xlab = xlab, ylab = "Frequency",
          col = "gray", border = "black", las = 1)
  boxplot(data[[column]], main = main_title, ylab = ylab, col = "gray", border = "black")
}

summary_and_plot(model_simulations, "Index_equ", "Amplification Index", "Index", "Index")


## -------------------------------
## 6. Estimate Reproduction Number (Ri)
## -------------------------------

### 6.1 Reorganize Parameters
parameters_list_updated <- list(
  d = model_simulations$delta,
  x = model_simulations$equis,
  a = model_simulations$alfa,
  e = model_simulations$epsilon,
  v = model_simulations$ve,
  l = model_simulations$lambda,
  n = model_simulations$ene,
  pi = model_simulations$pi,
  m = model_simulations$eme,
  z = model_simulations$zeta,
  p_1 = model_simulations$psi_1,
  p_2 = model_simulations$psi_2,
  p_3 = model_simulations$psi_3,
  p_4 = model_simulations$psi_4,
  h_1 = model_simulations$ache_1,
  h_2 = model_simulations$ache_2,
  g_1 = model_simulations$gamma_1,
  g_2 = model_simulations$gamma_2,
  f_2 = model_simulations$efe_2,
  f_1 = model_simulations$efe_1,
  s_1 = model_simulations$sigma_1,
  s_2 = model_simulations$sigma_2)

parameters_df <- as.data.frame(parameters_list_updated)

### 6.2 Calculate Ri from Next-Generation Matrix Formula
### For details, see supplementary Matlab script (Script 1).
R_formula = expression((d*l*(a + p_1*z)*(p_3*p_4^2*x + p_4^2*v*x + e*f_1*p_3*v + e*h_1*p_3*v + f_1*h_1*p_3*x + e*p_3*p_4*v + f_1*h_1*v*x + f_1*p_3*p_4*x + h_1*p_3*p_4*x + f_1*p_4*v*x + h_1*p_4*v*x - e*p_4^2*v*x - e*f_1*g_1*p_3*v - e*f_1*p_4*v*x - e*h_1*p_4*v*x - e*f_1*h_1*s_1*v*x))/(p_2*p_3*(a + p_1)*(f_1 + p_4)*(h_1 + p_4)*(p_3 + v)))

# Evaluate R for each simulation and add to dataframe
R_values <- numeric(nrow(parameters_df))

for (i in seq_len(nrow(parameters_df))) {
  param_env <- list2env(as.list(parameters_df[i, ]))
  R_values[i] <- eval(R_formula, envir = param_env)
}

parameters_df$R_formula <- R_values
View(parameters_df)

summary_and_plot(parameters_df, "R_formula", "Intrinsec Reproduction Number", "R", "R")

#S1 Figure
plot_df <- data.frame(
  R = parameters_df$R_values,
  Amplification = model_simulations$Index_equ
)

ggplot(plot_df, aes(x = R, y = Amplification)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "Intrinsec Reproduction Number (Ri)",
    y = "Colonization Amplification index (Ai)"
  ) +
  theme_classic()

                              
