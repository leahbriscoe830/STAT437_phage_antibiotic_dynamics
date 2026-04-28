# Based on Figure 6 of Leclerc et al.

library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(RColorBrewer)

source(here::here("Model", "efaecalis_model_final.R"))
set.seed(40)

palette = c(brewer.pal(n = 9, name = "BuGn")[c(3,4,6,8,9)], "#000000")

# Define parameters as a named vector
# Run the model script first to get these parameters
parameters <- readRDS("../Parameters/parameters.rds")

times <- seq(0, 24, by = 0.1)

#define initial states 
yinit <- c(Bv = 1e9, Ba = 1e9, Bva = 0, 
           Pl = 0, Pv = 0, Pa = 0, 
           amp = 0, van = 0)


event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(0, 0, 0),
                       value = c(2, 2, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bva, Bv, Ba) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pa = ggplot(all_results) +
  geom_line(aes(time, Bva, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Bv+Ba, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+4, xend = 0, y = 1e4, yend = 1e4) +
  geom_label(x = 0+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+5, xend = 0, y = 1e9, yend = 1e9) +
  geom_label(x = 0+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (CFU/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-11", "1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-11),
                                bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.01, 5e9)) + # ylim = c(0.1, 5e9)
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(10, 10, 0),
                       value = c(2, 2, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bva, Bv, Ba) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pb = ggplot(all_results) +
  geom_line(aes(time, Bva, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Bv+Ba, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+4, xend = 0, y = 1e4, yend = 1e4) +
  geom_label(x = 0+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  annotate("segment", x = 10+5, xend = 10, y = 1e9, yend = 1e9) +
  geom_label(x = 10+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (CFU/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-11", "1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-11),
                                bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.01, 5e9)) + # ylim = c(0.1, 5e9)
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(0, 0, 10),
                       value = c(2, 2, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bva, Bv, Ba) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pc = ggplot(all_results) +
  geom_line(aes(time, Bva, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Bv+Ba, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  annotate("segment", x = 10+4, xend = 10, y = 1e4, yend = 1e4) +
  geom_label(x = 10+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+5, xend = 0, y = 1e9, yend = 1e9) +
  geom_label(x = 0+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (CFU/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-11", "1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-11),
                                bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.01, 5e9)) + # ylim = c(0.1, 5e9)
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


## final plot ############

pall = plot_grid(plot_grid(plot_grid(pa + theme(legend.position = "none"),
                              pb + theme(legend.position = "none"),
                              pc + theme(legend.position = "none"),
                              ncol = 3,
                              labels = c("a)", "b)", "c)")),
                    get_legend(pa + guides(linetype = "none") + theme(legend.position = "bottom")),
                    get_legend(pa + guides(colour = "none") + theme(legend.position = "bottom",
                                                               legend.key.width = unit(2, "cm"))),
                    nrow = 3,
                    rel_heights = c(1,0.1,0.1)))#,

ggsave(here::here("Figures", "ef_fig3.png"), pall, height = 10, width = 10)
