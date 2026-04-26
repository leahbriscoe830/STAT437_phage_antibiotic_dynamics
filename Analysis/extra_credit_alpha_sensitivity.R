library(deSolve)
library(ggplot2)
library(scales)
library(dplyr)

# Load the model
source("Model/new_model.R")

# Load parameters
abx_params <- read.csv("Parameters/abx_params.csv")
pha_params <- read.csv("Parameters/pha_params_new.csv")
bac_params <- read.csv("Parameters/bac_params.csv")

# Define model parameters
parameters <- c(
  mu_e = bac_params$mu_e[1],
  mu_t = bac_params$mu_t[1],
  mu_et = bac_params$mu_et[1],
  Nmax = bac_params$Nmax[1],
  
  beta = pha_params$beta,
  L = pha_params$L,
  tau = pha_params$tau,
  alpha = pha_params$alpha,
  P50 = pha_params$P50,
  gamma = 0,
  
  ery_kill_max_BE = abx_params$kmax[1],
  ery_kill_max_BT = abx_params$kmax[3],
  ery_kill_max_BET = abx_params$kmax[5],
  
  tet_kill_max_BE = abx_params$kmax[2],
  tet_kill_max_BT = abx_params$kmax[4],
  tet_kill_max_BET = abx_params$kmax[6],
  
  EC_ery_BE = abx_params$EC50[1],
  EC_ery_BT = abx_params$EC50[3],
  EC_ery_BET = abx_params$EC50[5],
  
  EC_tet_BE = abx_params$EC50[2],
  EC_tet_BT = abx_params$EC50[4],
  EC_tet_BET = abx_params$EC50[6],
  
  pow_ery_BE = abx_params$pow[1],
  pow_ery_BT = abx_params$pow[3],
  pow_ery_BET = abx_params$pow[5],
  
  pow_tet_BE = abx_params$pow[2],
  pow_tet_BT = abx_params$pow[4],
  pow_tet_BET = abx_params$pow[6],
  
  gamma_ery = 0,
  gamma_tet = 0
)

# Time in hours
times <- seq(0, 48, by = 0.1)

# Initial population values
yinit <- c(
  Be = 1e9,
  Bt = 1e9,
  Bet = 0,
  Pl = 0,
  Pe = 0,
  Pt = 0,
  ery = 0,
  tet = 0
)

# Add erythromycin, tetracycline, and phage at time 0
event_dat <- data.frame(
  var = c("ery", "tet", "Pl"),
  time = c(0, 0, 0),
  value = c(1, 1, 1e9),
  method = c("add", "add", "add")
)

# Alpha values for sensitivity analysis
alpha_values <- c(1e-10, 1e-8, 1e-6)

all_results <- data.frame()

for (a in alpha_values) {
  
  cat("Running alpha =", a, "\n")
  
  parameters["alpha"] <- a
  
  results <- phage_tr_model(
    parameters = parameters,
    init.state = yinit,
    times = times,
    event_dat = event_dat
  )
  
  results$alpha <- as.factor(a)
  
  all_results <- bind_rows(all_results, results)
}

# Plot double-resistant bacteria over time
alpha_plot <- ggplot(all_results, aes(x = time, y = Bet, color = alpha)) +
  geom_line(size = 1) +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  coord_cartesian(ylim = c(0.1, 1e10)) +
  theme_bw() +
  labs(
    title = "Sensitivity Analysis of Transduction Probability",
    x = "Time (hours)",
    y = expression("Double-resistant bacteria " ~ B[ET] ~ " (CFU/mL)"),
    color = expression(alpha ~ " value")
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(face = "bold")
  )

print(alpha_plot)

ggsave(
  "Figures/extra_credit_alpha_sensitivity.png",
  plot = alpha_plot,
  width = 8,
  height = 5,
  dpi = 600
)