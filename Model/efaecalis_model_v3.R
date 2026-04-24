library(deSolve)
library(scales)

phage_tr_model <- function(parameters, init.state, times, event_dat, phage_activity = FALSE) {
  
  model_dde <- function(time, state, parameters) {
    
    #Bacterial growth parameters
    mu_v = parameters[["mu_v"]]
    mu_a = parameters[["mu_a"]]
    mu_va = parameters[["mu_va"]]
    Nmax = parameters[["Nmax"]]
    
    #Phage-related parameters
    beta = parameters[["beta"]]
    L = parameters[["L"]]
    gamma = parameters[["gamma"]]
    alpha = parameters[["alpha"]]
    tau = parameters[["tau"]]
    
    #Antibiotic effect parameters for Bv (Vancomycin)
    van_kill_max_BV = parameters[["van_kill_max_BV"]]
    amp_kill_max_BV = parameters[["amp_kill_max_BV"]]
    EC_van_BV = parameters[["EC_van_BV"]]
    EC_amp_BV = parameters[["EC_amp_BV"]]
    pow_van_BV = parameters[["pow_van_BV"]]
    pow_amp_BV = parameters[["pow_amp_BV"]]
    
    #Antibiotic effect parameters for Ba (Ampicillin)
    van_kill_max_BA = parameters[["van_kill_max_BA"]]
    amp_kill_max_BA = parameters[["amp_kill_max_BA"]]
    EC_van_BA = parameters[["EC_van_BA"]]
    EC_amp_BA = parameters[["EC_amp_BA"]]
    pow_van_BA = parameters[["pow_van_BA"]]
    pow_amp_BA = parameters[["pow_amp_BA"]]
    
    #Antibiotic effect parameters for Bva (both)
    van_kill_max_BVA = parameters[["van_kill_max_BVA"]]
    amp_kill_max_BVA = parameters[["amp_kill_max_BVA"]]
    EC_van_BVA = parameters[["EC_van_BVA"]]
    EC_amp_BVA = parameters[["EC_amp_BVA"]]
    pow_van_BVA = parameters[["pow_van_BVA"]]
    pow_amp_BVA = parameters[["pow_amp_BVA"]]
    
    #Antibiotic decay parameters
    gamma_amp = parameters[["gamma_amp"]]
    gamma_van = parameters[["gamma_van"]]
    
    #Read current state variables
    Bv = state[["Bv"]]
    Ba = state[["Ba"]]
    Bva = state[["Bva"]]
    Pl = state[["Pl"]]
    Pv = state[["Pv"]] # Phage carrying Van resistance
    Pa = state[["Pa"]] # Phage carrying Amp resistance
    amp = state[["amp"]]
    van = state[["van"]]
    
    N = Bv + Ba + Bva
    
    #Get the delayed values from tau hours ago
    if(time <= tau){
      Bv_past = 0; Ba_past = 0; Bva_past = 0
      Pl_past = 0; Pv_past = 0; Pa_past = 0
      N_past = 1; amp_past = 0; van_past = 0
    } else {
      Bv_past = lagvalue(time - tau, 1)
      Ba_past = lagvalue(time - tau, 2)
      Bva_past = lagvalue(time - tau, 3)
      Pl_past = lagvalue(time - tau, 4)
      Pv_past = lagvalue(time - tau, 5)
      Pa_past = lagvalue(time - tau, 6)
      N_past = Bv_past + Ba_past + Bva_past
      amp_past = lagvalue(time - tau, 7)
      van_past = lagvalue(time - tau, 8)
    }
    
    #Logistic growth
    link = (1 - N/Nmax)
    lambda = (1 - exp(-beta * N))
    phi_Pl = (1 - exp(-lambda * Pl/N))
    phi_Pv = (1 - exp(-lambda * Pv/N))
    phi_Pa = (1 - exp(-lambda * Pa/N))
    
    lambda_past = (1 - exp(-beta * N_past))
    phi_Pl_past = (1 - exp(-lambda_past * Pl_past/N_past))
    phi_Pv_past = (1 - exp(-lambda_past * Pv_past/N_past))
    phi_Pa_past = (1 - exp(-lambda_past * Pa_past/N_past))
    
    #Differential equations
    #Using Hill equations for antibiotic effect
    #For Vancomycin-resistant strain (Bv) ---
    van_effect_BV = van_kill_max_BV * van^pow_van_BV / (EC_van_BV^pow_van_BV + van^pow_van_BV)
    amp_effect_BV = amp_kill_max_BV * amp^pow_amp_BV / (EC_amp_BV^pow_amp_BV + amp^pow_amp_BV)
    
    #For Ampicillin-resistant strain (Ba) ---
    van_effect_BA = van_kill_max_BA * van^pow_van_BA / (EC_van_BA^pow_van_BA + van^pow_van_BA)
    amp_effect_BA = amp_kill_max_BA * amp^pow_amp_BA / (EC_amp_BA^pow_amp_BA + amp^pow_amp_BA)
    
    #For Double-resistant strain (Bva) ---
    van_effect_BVA = van_kill_max_BVA * van^pow_van_BVA / (EC_van_BVA^pow_van_BVA + van^pow_van_BVA)
    amp_effect_BVA = amp_kill_max_BVA * amp^pow_amp_BVA / (EC_amp_BVA^pow_amp_BVA + amp^pow_amp_BVA)
    
    # Burst scaling for double-resistant bacteria
    L_VA = L * max(0, (link - van_effect_BVA - amp_effect_BVA)) + 1
    L_V = L * max(0, (link - van_effect_BV - amp_effect_BV)) + 1
    L_A = L * max(0, (link - van_effect_BA - amp_effect_BA)) + 1
    
    #Dynamics for Bv, Ba, Bva 
    dBv = (mu_v * link) * (Bv - ((phi_Pl + phi_Pa - (phi_Pl*phi_Pa)) * Bv) ) -
      (phi_Pl + phi_Pa - (phi_Pl*phi_Pa)) * Bv - (van_effect_BV + amp_effect_BV)*mu_v*Bv
    
    dBa = (mu_a * link) * (Ba - ((phi_Pl + phi_Pv - (phi_Pl*phi_Pv)) * Ba) ) -
      (phi_Pl + phi_Pv - (phi_Pl*phi_Pv)) * Ba - (amp_effect_BA + van_effect_BA)*mu_a*Ba
    
    dBva = mu_va * link * (Bva - (phi_Pl*Bva) ) - phi_Pl * Bva +
      (phi_Pv - (phi_Pl*phi_Pv)) * Ba + (phi_Pa - (phi_Pl*phi_Pa)) * Bv - (amp_effect_BVA + van_effect_BVA)*mu_va*Bva
    
    dPl = phi_Pl_past * L_V * (1-alpha) * Bv_past + phi_Pl_past * 
      L_A * (1-alpha) * Ba_past + phi_Pl_past * L_VA * (1-2*alpha) * Bva_past - lambda * Pl - gamma * Pl

    dPv = phi_Pl_past * L_V * alpha * Bv_past + phi_Pl_past * 
      L_VA * alpha * Bva_past - lambda * Pv - gamma * Pv

    dPa = phi_Pl_past * L_A * alpha * Ba_past + phi_Pl_past * 
      L_VA * alpha * Bva_past - lambda * Pa - gamma * Pa 

    #Antibiotic decay equations
    damp = -amp*gamma_amp

    dvan = -van*gamma_van
    
    
    if(phage_activity){
      # Number of Bva cells currently being lysed by lytic phage
      Bva_lysis = unname(phi_Pl * Bva)
      
      # Number of new Bva cells produced by growth after accounting for lytic infection
      Bva_new = unname(mu_va * link * (Bva - (phi_Pl*Bva) ))
      
      # Return derivatives plus these extra tracked quantities
      return(list(c(dBv, dBa, dBva, dPl, dPv, dPa, damp, dvan),
                  Bva_lysis = Bva_lysis, Bva_new = Bva_new))
    } else {
      # Standard return: derivatives only
      return(list(c(dBv, dBa, dBva, dPl, dPv, dPa, damp, dvan)))
    }
  }
  
  # Solver call
  trajectory <- data.frame(dede(y = init.state, times = times, func = model_dde, 
                                parms = parameters, events = list(data = event_dat)))
  return(trajectory)
}

# Define parameters as a named vector
params <- c(
  #Bacterial Growth (mu in min^-1)
  mu_v = 0.015,     # Growth rate of Vancomycin-resistant strain
  mu_a = 0.015,     # Growth rate of Ampicillin-resistant strain
  mu_va = 0.014,    # Growth rate of double-resistant strain
  Nmax = 1e9,       # Carrying capacity (max cells/mL)
  
  #Phage Kinetics
  beta = 5e-10,     # Adsorption rate (mL/min)
  L = 100,          # Burst size (PFU/cell)
  gamma = 0.005,    # Phage decay rate (min^-1)
  alpha = 0.001,    # Probability of transduction
  tau = 0.333, # 20 min # Latent period (min)
  
  #Antibiotic Killing: Vancomycin
  van_kill_max_BV = 0.5, EC_van_BV = 2.0, pow_van_BV = 2.0, # BV strain
  van_kill_max_BA = 0.2, EC_van_BA = 5.0, pow_van_BA = 2.0, # BA strain (lower effect)
  van_kill_max_BVA = 0.1, EC_van_BVA = 10.0, pow_van_BVA = 2.0, # BVA (resistant)
  
  #Antibiotic Killing: Ampicillin
  amp_kill_max_BV = 0.2, EC_amp_BV = 5.0, pow_amp_BV = 2.0,  # BV strain
  amp_kill_max_BA = 0.5, EC_amp_BA = 2.0, pow_amp_BA = 2.0,  # BA strain
  amp_kill_max_BVA = 0.1, EC_amp_BVA = 10.0, pow_amp_BVA = 2.0, # BVA (resistant)
  
  #Decay rates (min^-1)
  gamma_amp = 0.001,
  gamma_van = 0.001
)
saveRDS(params, "../Parameters/parameters.rds")

# Define when events happen
event_dat <- data.frame(
  var = c("Pl", "van", "amp"),   # Which variable to change
  time = c(0, 8, 12),            # Time (hours)
  value = c(1e7, 5, 10),         # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

#define initial states 
init.state <- c(Bv = 1e6, Ba = 1e6, Bva = 0, 
                Pl = 0, Pv = 0, Pa = 0, 
                amp = 0, van = 0)

library(ggplot2)
# Run the model
times <- seq(0, 48, by = 0.1) # Simulate 48 hours
results <- phage_tr_model(parameters = params, init.state, times, event_dat)

# Plotting
ggplot(results, aes(x = time)) +
  geom_line(aes(y = Bv, color = "Van-Resistant (Bv)"), size = 0.8) +
  geom_line(aes(y = Ba, color = "Amp-Resistant (Ba)"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  labs(y = "CFU/mL", x = "Time (hours)", colour = "Bacteria:") +
  theme_bw()

