## -------------------------------
## 0. Required Packages
## -------------------------------

# Load necessary libraries
if (!require(pacman)) install.packages("pacman")
library(pacman)
p_load(deSolve, reshape2, ggplot2, dplyr,tidyr)


## -------------------------------
## 1. Define Model Structure and Parameters
## -------------------------------

### 1.1 Initial Conditions and Parameters
initial_conditions <- c(
  R = 121,                    # Non-susceptible population
  S = 334,                    # Susceptible population
  E = 0,                      # Asymptomatic carriers type E
  C = 0,                      # Asymptomatic carriers type C
  K_1 = 0,                    # Treated asymptomatic carriers type K1 with contact precautions
  K_2 = 0,                    # Treated asymptomatic carriers type K2 with contact precautions                  
  I = 0,                      # Symptomatic patients 
  A = 0,                      # Symptomatic patients with treatment and contact precautions
  dAdmission = 0,             # Cumulative number of admissions
  dDiagnosed_CDI=0) 

parameters <- c(
  delta = 0.0007233824,    # Transmission rate of infected individuals, to be optimized based on observed institutional cases. 
  equis = 0.7,             # Relative transmission for colonized patients.
  alfa = 0.064,            # Rate of antibiotic use.
  epsilon = 0.13,          # Fraction of asymptomatic carriers who develop symptomatic disease during hospitalization
  ve = 1/4,                # Progression rate to symptomatic disease.
  lambda = 76,             # Admission rate to the hospital or facility.
  ene = 0.08,              # Overall proportion of asymptomatic carriers at admission.
  pi = 0.004986150,        # Proportion of infected patients at admission.
  eme = 0.12,              # Proportion of pre-symptomatic patients at admission. 
  zeta = 0.22,             # Proportion of susceptible patients at admission (patient with history of antibiotic use before admission).
  psi_1 = 1/6,             # Discharge rate for non-colonized population.
  psi_2 = 1/6,             # Discharge rate for susceptible population.
  psi_3 = 1/6,             # Discharge rate for asymptomatic carriers.
  psi_4 = 1/12,            # Discharge rate for symptomatic population.
  ache_1 = 1/10,           # Bacterial clearance rate due to treatment for infected patients.
  ache_2 = 1/10,           # Bacterial clearance rate due to treatment for asymptomatic patients.
  gamma_1 = 0.5,           # Relative reduction of transmission due to contact precautions on symptomatic patients (Effectiveness of contact precautions).
  gamma_2 = 0,             # Relative reduction of transmission due to contact precautions on asymptomatic patients (Effectiveness of contact precautions).
  efe_2 = 0,               # Proportion of screened individuals
  efe_1 = 1/1.5,           # Diagnosis rate of symptomatic carriers.
  sigma_1 = 0.7,           # Fraction of cured symptomatic carriers after treatment (Effectiveness of treatment).
  sigma_2 = 0              # Fraction of cured asymptomatic carriers after treatment (Effectiveness of treatment).
)

follow_up_duration = 570 
times <- seq(from = 0, to = follow_up_duration, by = 1) 

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
    
    # Return compartments as a list
    return(list(c(dR, dS, dE, dC, dK_1,dK_2, dI, dA, dAdmission, dDiagnosed_CDI))) 
  })
}


## -------------------------------
## 2. Baseline Simulation (No Intervention)
## -------------------------------

# Baseline simulation with no interventions 
model_output <- as.data.frame(ode(
  y = initial_conditions,
  times = times,
  func = comp_model,
  parms = parameters
)) %>%
  mutate(IncidenceCDI = (dDiagnosed_CDI / dAdmission) * 1000)



## -------------------------------
## 3. Functions for Simulations
## -------------------------------
# Function for a single simulation, returns final incidence and colonization index
run_single_simulation <- function(params) {
  output <- as.data.frame(ode(
    y = initial_conditions,
    times = times,
    func = comp_model,
    parms = params
  ))
  
  last_row <- output[nrow(output), ]
  
  IncidenceCDI <- (last_row[["dDiagnosed_CDI"]] / last_row[["dAdmission"]]) * 1000
  
  index_c <- as.numeric((params["psi_3"] * (last_row["E"] + last_row["C"] + last_row["K_1"] + last_row["K_2"])) / 
                          (params["lambda"] * params["ene"]))
  
  return(list(
    IncidenceCDI = IncidenceCDI,
    index_c = index_c
  ))
}

### 3.2 Function for Intervention Scenarios
run_intervention <- function(sigma_2_val, gamma_2_val, ache_2_val, efe_2_val, intervention_label) {
  results_list <- vector("list", num_simulations)
  
  for (i in 1:num_simulations) {
    params <- parameters
    params["sigma_2"] <- ifelse(is.null(sigma_2_val), runif(1, 0.5, 1), sigma_2_val)
    params["gamma_2"] <- ifelse(is.null(gamma_2_val), runif(1, 0.5, 1), gamma_2_val)
    params["ache_2"] <- ifelse(is.null(ache_2_val), runif(1, (1/10)*0.8, (1/10)*1.2), ache_2_val)
    params["efe_2"] <- efe_2_val
    
    results <- run_single_simulation(params)
    
    results_list[[i]] <- data.frame(
      IncidenceCDI = results$IncidenceCDI,
      index_c = results$index_c,
      time = follow_up_duration,
      efe_2 = efe_2_val,
      Intervention = intervention_label
    )
  }
  
  do.call(rbind, results_list)
}


## -------------------------------
## 4. Run Simulations for Each Intervention
## -------------------------------

num_simulations <- 100
fixed_efe2_values <- seq(0, 1, by = 0.1) # Screening proportions from 0 to 1

#Simulations of interventions:

### 4.1 Isolation Only
isolation_results <- lapply(fixed_efe2_values, function(efe2) {
  run_intervention(
    sigma_2_val = 0,
    gamma_2_val = NULL,  # random
    ache_2_val = 0,
    efe_2_val = efe2,
    intervention_label = "Isolation"
  )
}) %>% bind_rows()

### 4.2 Isolation + Prophylaxis
prophylactic_results <- lapply(fixed_efe2_values, function(efe2) {
  run_intervention(
    sigma_2_val = NULL,  # random
    gamma_2_val = NULL,  # random
    ache_2_val = NULL,   # random
    efe_2_val = efe2,
    intervention_label = "Isolation + Prophylaxis"
  )
}) %>% bind_rows()

### 4.3 Prophylaxis Only
prophylaxis_only_results <- lapply(fixed_efe2_values, function(efe2) {
  run_intervention(
    sigma_2_val = NULL,  # random
    gamma_2_val = 0,
    ache_2_val = NULL,   # random
    efe_2_val = efe2,
    intervention_label = "Prophylaxis Only"
  )
}) %>% bind_rows()

# Combine results
all_results_combined <- bind_rows(isolation_results, prophylactic_results, prophylaxis_only_results)



## -------------------------------
## 5. Calculate Reductions in Incidence and Colonization
## -------------------------------

# Baseline references for reduction calculations
reference_value_incidence <- all_results_combined %>%
  filter(efe_2 == 0, Intervention == "Isolation") %>%
  summarise(mean_IncidenceCDI = mean(IncidenceCDI)) %>%
  pull(mean_IncidenceCDI)

reference_value_index_c <- all_results_combined %>%
  filter(efe_2 == 0, Intervention == "Isolation") %>%
  summarise(mean_index_c = mean(index_c)) %>%
  pull(mean_index_c)

# Apply reductions to data
all_results_combined <- all_results_combined %>%
  mutate(
    ReductionPercent_Incidence = ((reference_value_incidence - IncidenceCDI) / reference_value_incidence) * 100,
    ReductionPercent_index_c = ((reference_value_index_c - index_c) / reference_value_index_c) * 100
  )

## -------------------------------
## 6. Plots of Intervention Effects
## -------------------------------

### 6.1 Reduction in CDI Incidence
ggplot(all_results_combined %>% filter(efe_2 != 0), 
       aes(x = as.factor(efe_2), y = ReductionPercent_Incidence, fill = Intervention)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.color = "black") +
  scale_fill_manual(values = c("Isolation" = "lightgray", 
                               "Isolation + Prophylaxis" = "dimgray", 
                               "Prophylaxis Only" = "white")) +
  labs(
    title = "Reduction in CDI by Screening (efe_2) and Treatment Strategy",
    x = "Screening Proportion (efe_2)",
    y = "Reduction Percentage (%)",
    fill = "Treatment Strategy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    legend.position = "top"
  )


### 6.2 Reduction in Colonization Index
ggplot(all_results_combined %>% filter(efe_2 != 0), 
       aes(x = as.factor(efe_2), y = ReductionPercent_index_c, fill = Intervention)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.color = "black") +
  scale_fill_manual(values = c("Isolation" = "lightgray", 
                               "Isolation + Prophylaxis" = "dimgray", 
                               "Prophylaxis Only" = "white")) +
  labs(
    title = "Reduction in Colonization Index by Screening (efe_2) and Treatment Strategy",
    x = "Screening Proportion (efe_2)",
    y = "Reduction in Colonization Index (%)",
    fill = "Treatment Strategy"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black"),
    legend.position = "top"
  )

## -------------------------------
## 7. Summary Tables
## -------------------------------

table_incidence <- all_results_combined %>%
  group_by(efe_2, Intervention) %>%
  summarise(mean_Reduction_Incidence = mean(ReductionPercent_Incidence), .groups = "drop") %>%
  pivot_wider(names_from = Intervention, values_from = mean_Reduction_Incidence) %>%
  select(efe_2, Isolation, `Prophylaxis Only`, `Isolation + Prophylaxis`) %>%
  arrange(efe_2)

table_index_c <- all_results_combined %>%
  group_by(efe_2, Intervention) %>%
  summarise(mean_Reduction_index_c = mean(ReductionPercent_index_c), .groups = "drop") %>%
  pivot_wider(names_from = Intervention, values_from = mean_Reduction_index_c) %>%
  select(efe_2, Isolation, `Prophylaxis Only`, `Isolation + Prophylaxis`) %>%
  arrange(efe_2)

View(table_incidence, title = "Reduction in CDI Incidence (%)")
View(table_index_c, title = "Reduction in Colonization Index (%)")
