# Based on Figure 3 of Moryl et al.

library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)

source(here::here("Model", "efaecalis_model_final.R"))

# Define parameters as a named vector
# Run the model script first to get these parameters
parameters <- readRDS("../Parameters/parameters.rds")

times <- seq(0, 24, by = 0.1)

## VAN-RESISTANT ##########
all_results = data.frame()

yinit <- c(Bv = 2e7, Ba = 0, Bva = 0,
           Pl = 0, Pv = 0, Pa = 0,
           amp = 0, van = 0)


## Control ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 100),          # Time (hours)
  value = c(0, 0, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Van-resistant"
results$group = "Control"
all_results = rbind(all_results, results)


## Phage ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 0.1),          # Time (hours)
  value = c(0, 0, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Van-resistant"
results$group = "Phage"
all_results = rbind(all_results, results)


## Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0.1, 100),          # Time (hours)
  value = c(0, 1, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Van-resistant"
results$group = "Ampicillin"
all_results = rbind(all_results, results)


## Phage + Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0, 0),          # Time (hours)
  value = c(0, 1, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Van-resistant"
results$group = "Phage+Ampicillin"
all_results = rbind(all_results, results)


## Plot ##

van = ggplot(all_results) +
  geom_line(aes(time, Bv, color = group), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 11e9)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "CFU/mL", x = "Time (hours)", color = "Group:", title = "Van-resistant") +
  scale_colour_manual(breaks = c("Control", "Phage", "Ampicillin", "Phage+Ampicillin"),
                      values = c("#d18b2a","#e3c530","#6db356","#63311d")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))



## Amp-resistant ##########
all_results = data.frame()

yinit <- c(Bv = 0, Ba = 2e7, Bva = 0,
           Pl = 0, Pv = 0, Pa = 0,
           amp = 0, van = 0)

## Control ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 100),          # Time (hours)
  value = c(0, 0, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Amp-resistant"
results$group = "Control"
all_results = rbind(all_results, results)


## Phage ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 0),          # Time (hours)
  value = c(0, 0, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Amp-resistant"
results$group = "Phage"
all_results = rbind(all_results, results)


## Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0.1, 100),          # Time (hours)
  value = c(0, 1, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Amp-resistant"
results$group = "Ampicillin"
all_results = rbind(all_results, results)


## Phage + Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0, 0),          # Time (hours)
  value = c(0, 1, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Amp-resistant"
results$group = "Phage+Ampicillin"
all_results = rbind(all_results, results)


## Plot ##

amp = ggplot(all_results) +
  geom_line(aes(time, Ba, color = group), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 11e9)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "CFU/mL", x = "Time (hours)", color = "Group:", title = "Amp-resistant") +
  scale_colour_manual(breaks = c("Control", "Phage", "Ampicillin", "Phage+Ampicillin"),
                      values = c("#d18b2a","#e3c530","#6db356","#63311d")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))



## Double-resistant ##########
all_results = data.frame()

yinit <- c(Bv = 0, Ba = 0, Bva = 2e7,
           Pl = 0, Pv = 0, Pa = 0,
           amp = 0, van = 0)

## Control ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 100),          # Time (hours)
  value = c(0, 0, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Double-resistant"
results$group = "Control"
all_results = rbind(all_results, results)


## Phage ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 100, 0.1),          # Time (hours)
  value = c(0, 0, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Double-resistant"
results$group = "Phage"
all_results = rbind(all_results, results)


## Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0.1, 100),          # Time (hours)
  value = c(0, 1, 0),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Double-resistant"
results$group = "Ampicillin"
all_results = rbind(all_results, results)


## Phage + Ampicillin ##
event_dat <- data.frame(
  var = c("van", "amp", "Pl"),      # Which variable to change
  time = c(100, 0, 0),          # Time (hours)
  value = c(0, 1, 2e6),               # Concentration/Amount to add
  method = c("add", "add", "add")   # "add" increases, "replace" overwrites
)

results = phage_tr_model(parameters, yinit, times, event_dat)

results$bac = "Double-resistant"
results$group = "Phage+Ampicillin"
all_results = rbind(all_results, results)


## Plot ##

doub = ggplot(all_results) +
  geom_line(aes(time, Bva, color = group), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 11e9)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "CFU/mL", x = "Time (hours)", color = "Group:", title = "Double-resistant") +
  scale_colour_manual(breaks = c("Control", "Phage", "Ampicillin", "Phage+Ampicillin"),
                      values = c("#d18b2a","#e3c530","#6db356","#63311d")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))



## FINAL PLOT ##########
plot_grid(plot_grid(
  plot_grid(van + theme(legend.position = "none"),
            amp + theme(legend.position = "none"),
            doub + theme(legend.position = "none"),
            ncol = 3,
            labels = c("a)", "b)" ,"c)")),
  get_legend(van + theme(legend.position = "bottom")),
  nrow = 2,
  rel_heights = c(1, 0.1)))

ggsave(here::here("Figures", "ef_fig4.png"), height = 9, width = 20)