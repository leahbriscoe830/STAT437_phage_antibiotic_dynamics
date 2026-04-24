
## SETUP ##########

library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)

source(here::here("Model", "efaecalis_model_v3.R"))

# Define parameters as a named vector
# Run the model script first to get these parameters
parameters <- readRDS("../Parameters/parameters.rds")

times <- seq(0, 24, by = 0.1)

#define initial states 
yinit <- c(Bv = 1e6, Ba = 1e6, Bva = 0, 
           Pl = 0, Pv = 0, Pa = 0, 
           amp = 0, van = 0)

all_results = data.frame()

## NOTHING ##########

# Define when events happen
event_dat <- data.frame(
  var = c("Pl", "van", "amp"),      # Which variable to change
  time = c(100, 100, 100),          # Time (hours)
  value = c(1, 1, 1e9),             # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE ONLY ##########

# Define when events happen
event_dat <- data.frame(
  var = c("Pl", "van", "amp"),      # Which variable to change
  time = c(100, 100, 0),            # Time (hours)
  value = c(1, 1, 1e9),             # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## VAN ONLY ##########

# Define when events happen
event_dat <- data.frame(
  var = c("Pl", "van", "amp"),      # Which variable to change
  time = c(0, 100, 100),            # Time (hours)
  value = c(1, 1, 1e9),             # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Vancomycin"
all_results = rbind(all_results, results)


## AMP ONLY ##########

# Define when events happen
event_dat <- data.frame(
  var = c("Pl", "van", "amp"),      # Which variable to change
  time = c(100, 0, 100),            # Time (hours)
  value = c(1, 1, 1e9),             # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Ampicillin"
all_results = rbind(all_results, results)


## VAN AND AMP ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(0, 0, 100),            # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Vancomycin + Ampicillin"
all_results = rbind(all_results, results)


## PHAGE AND VAN ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(0, 100, 0),            # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Vancomycin"
all_results = rbind(all_results, results)


## PHAGE AND AMP ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(100, 0, 0),            # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Ampicillin"
all_results = rbind(all_results, results)


## PHAGE AND VAN AND AMP ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),   # Which variable to change
  time = c(0, 0, 0),             # Time (hours)
  value = c(1, 1, 1e9),          # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Vancomycin + Ampicillin"
all_results = rbind(all_results, results)



## PHAGE ONLY, TRANSDUCTION ##########

parameters[["alpha"]] = 1.012/1e8

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),      # Which variable to change
  time = c(100, 100, 0),            # Time (hours)
  value = c(1, 1, 1e9),             # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE AND VAN ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(0, 100, 0),            # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Vancomycin"
all_results = rbind(all_results, results)


## PHAGE AND AMP ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(100, 0, 0),            # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Ampicillin"
all_results = rbind(all_results, results)


## PHAGE AND ABX, TRANSDUCTION ##########

event_dat <- data.frame(
  var = c("Pl", "van", "amp"),    # Which variable to change
  time = c(0, 0, 0),              # Time (hours)
  value = c(1, 1, 1e9),           # Concentration/Amount to add
  method = c("add", "add", "add") # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Vancomycin + Ampicillin"
all_results = rbind(all_results, results)


## FINAL PLOT ##########

all_results$abx = as.factor(all_results$abx)
all_results$abx = factor(all_results$abx, levels = levels(all_results$abx)[c(3,1,4,2)])

all_results$Pl[all_results$Pl == 0] = NA

pa = ggplot(all_results) +
  geom_line(aes(time, Bv, colour = "Bv"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Ba, colour = "Ba"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Bva, colour = "Bva"), size = 1) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ abx) +
  coord_cartesian(ylim = c(0.1, 3e11)) + # 3e11
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Bv", "Ba", "Bva", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[V]),
                                 expression(B[A]),
                                 expression(B[VA]),
                                 expression(P[L]))) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


all_results2 = all_results %>%
  filter(pha == "Phage w/ transduction") %>%
  select(-c("Pv", "Pa", "van", "amp")) %>%
  mutate(tot = Bv+Ba+Bva) %>%
  reshape2::melt(. , id.vars = c("time","pha","abx"))

all_results2$variable = factor(all_results2$variable, levels = levels(all_results2$variable)[c(5,1,2,3,4)])

bac_labs = c("Total bacteria", "Van-resistant", "Amp-resistant", "Double-resistant", "Phage")
names(bac_labs) = c("tot", "Bv", "Ba", "Bva", "Pl")

pb = ggplot(all_results2) +
  geom_line(aes(time, value, colour = abx), size = 1, alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ variable, labeller = labeller(variable = bac_labs)) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Antibiotic:") +
  scale_colour_manual(values = c("black", "royalblue3","green3", "purple3")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


plot_grid(pa,
          NULL,
          pb,
          labels = c("a)","", "b)"),
          nrow = 3, rel_heights = c(1,0.05,0.5))

ggsave(here::here("Figures", "ef_fig3.png"), width = 10, height = 11)
