

library(deSolve)
library(ggplot2)
library(ggtext)
library(scales)
library(cowplot)
library(RColorBrewer)
library(dplyr)

source(here::here("Model", "efaecalis_model_v3.R"))

#things to test
test_pha_con = 1e8
test_abx_con = 1
test_diffs = c(0,3,5,15)

#palette = rev(brewer.pal(n = 9, name = "RdBu"))[c(1,3,7,9)]
palette_blues = brewer.pal(n = 7, name = "Blues")
palette_greens = brewer.pal(n = 7, name = "Greens")
 

# Define parameters as a named vector
# Run the model script first to get these parameters
parameters <- readRDS("../Parameters/parameters.rds")

times <- seq(0, 48, by = 0.1)

#define initial states 
yinit <- c(Bv = 1e6, Ba = 1e6, Bva = 0, 
           Pl = 0, Pv = 0, Pa = 0, 
           amp = 0, van = 0)

all_results_pha = data.frame()


#expanded view
abx_time = test_diffs[1]
pha_time = 0
event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx1_pha1 = ggplot(results) +
  geom_line(aes(time, Bv, colour = "Bv"), size = 0.8) +
  geom_line(aes(time, Ba, colour = "Ba"), size = 0.8) +
  geom_line(aes(time, Bva, colour = "Bva"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bva, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bva, na.rm = T))+1),
           yend = max(results$Bva), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bva, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 1.5, label = "0h", size = 4) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Bv", "Ba", "Bva", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[V]),
                                 expression(B[A]),
                                 expression(B[VA]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 21, fill = "black", size = 5)


abx_time = test_diffs[2]
pha_time = 0 
event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx6_pha1 = ggplot(results) +
  geom_line(aes(time, Bv, colour = "Bv"), size = 0.8) +
  geom_line(aes(time, Ba, colour = "Ba"), size = 0.8) +
  geom_line(aes(time, Bva, colour = "Bva"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bva, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bva, na.rm = T))+1),
           yend = max(results$Bva), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bva, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 1.5, label = "3h", size = 4) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Bv", "Ba", "Bva", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[V]),
                                 expression(B[A]),
                                 expression(B[VA]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 22, fill = "black", size = 5)



abx_time = test_diffs[3]
pha_time = 0 
event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx8_pha1 = ggplot(results) +
  geom_line(aes(time, Bv, colour = "Bv"), size = 0.8) +
  geom_line(aes(time, Ba, colour = "Ba"), size = 0.8) +
  geom_line(aes(time, Bva, colour = "Bva"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bva, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bva, na.rm = T))+1),
           yend = max(results$Bva), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bva, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 2.5, label = "5h", size = 4) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Bv", "Ba", "Bva", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 8, fill = "black", size = 5)



abx_time = test_diffs[4]
pha_time = 0 
event_dat = data.frame(var = c("van", "amp", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx16_pha1 = ggplot(results) +
  geom_line(aes(time, Bv, colour = "Bv"), size = 0.8) +
  geom_line(aes(time, Ba, colour = "Ba"), size = 0.8) +
  geom_line(aes(time, Bva, colour = "Bva"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bva, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bva, na.rm = T))+1),
           yend = max(results$Bva), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bva, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 7.5, label = "15h", size = 4) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Bv", "Ba", "Bva", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[V]),
                                 expression(B[A]),
                                 expression(B[VA]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))


final_plot = plot_grid(plot_grid(
                         plot_grid(p_abx1_pha1 + theme(legend.position = "none"),
                                   p_abx6_pha1 + theme(legend.position = "none"),
                                   p_abx8_pha1 + theme(legend.position = "none"),
                                   p_abx16_pha1 + theme(legend.position = "none"),
                                   ncol = 2),
                         get_legend(p_abx1_pha1 + theme(legend.position = "right")),
                         nrow = 1,
                         rel_widths = c(1,0.1)),
                       ncol = 1,
                       rel_heights = c(1.3, 0.05, 1),
                       labels = c("", "" ,"c)"), hjust = 0)

ggsave(here::here("Figures", "ef_fig5.png"), final_plot, height = 15, width = 13)

