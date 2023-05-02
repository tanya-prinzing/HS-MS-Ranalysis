
# Oct 2, 2022 - Tanya Prinzing

# * Update ratio plots for adults only as we shouldn't compare across the full size range

# Script for plots for Horn Shark scaling paper (GSA vs MR)

library("tidyverse")
library("patchwork")
library("viridis")
library("qualpalr")


# Plot color palette
pal4mods <- qualpal(n = 4, colorspace = "pretty")
pal2mods <- qualpal(n = 2, colorspace = "pretty")


# FIGURE 1:  PLOT GSA, MMR and RMR estimates -- ###### MMR Chase, not Air

# (a) Metabolic rate PLOT

# HS_MR_Plot <- 
ggplot(HornShark_justMR_df, aes(x = MassKG_log, y = WholeEstimate_log_KG, 
                                group = MeasureType)) +
  geom_smooth(method = "lm", se = FALSE, size = 1.3, alpha = 0.8, 
              show.legend = FALSE, color = "black") + 
  geom_point(size = 3.5, alpha = 0.7, show.legend = FALSE) + 
  scale_y_continuous(name = expression('Metabolic Rate (mg'~O[2]~ h^-1~')'),
                     breaks = seq(1, 3, 1),
                     labels = c((expression(paste(10^{1}))), (expression(paste(10^{2}))),
                                (expression(paste(10^{3}))))) +
  scale_x_continuous(name = expression(paste("Body Mass (kg)")),
                     breaks = seq(-1, 0.5, 0.5),
                     labels = c((expression(paste(10^{-1}))),(expression(paste(10^{-0.5}))),
                                (expression(paste(10^{0}))), (expression(paste(10^{0.5})))),
                     limits = c(-1.5, 0.65)) +
  annotate("text", x = -1.3, y = 2.1, label = expression(paste('MMR')),
           color = "black", hjust = 0, size = 3.8) +
  annotate("text", x = -1.3, y = 1.94, 
           label = expression(paste(~italic(b) ~'= 1.073 ± 0.040')),
           color = "black", hjust = 0, size = 3.8) +
  annotate("text",x = -0.75, y = 0.56, label = expression(paste('RMR')),
           color = "black",  hjust = 0, size = 3.8) +
  annotate("text",x = -0.75, y = 0.4,  
           label = expression(paste(~italic(b) ~'= 0.966 ± 0.058')),
           color = "black",  hjust = 0, size = 3.8) +
  annotate("text", x = -1.5, y = 3.05, label = "(a)", color = "black", hjust = 0, size = 3.8) +
  theme_bw() +
  annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"),
                      short = unit(0.75, "mm"), colour = "black", sides = "l") +
  theme(axis.title.y = element_text( size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0,-2,0), "cm"),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank())
 
  
 # (b) GSA Plot
HS_GSA_Plot <- 
 ggplot(GSA_df, aes(x = MassKG_log, y = WholeEstimate_log_KG, color = "black")) +
 geom_smooth(method = "lm", se = FALSE, size = 1.5, alpha = 0.5, 
             show.legend = FALSE, color = "grey60") +
 geom_line(data = HS_GSA_BR_lines_plot_df, aes(x = X_axis_spots, y = Y_axis_spots),
           size = 1.5, color = "black") +
 geom_point(size = 3.5, alpha = 0.7, show.legend = FALSE, color = "black") +
 scale_y_continuous(name = expression('Gill Surface Area ('~cm^2~')'),
                    breaks = seq(2.5, 4, 0.5),
                    labels = c((expression(paste(10^{2.5}))), (expression(paste(10^{3}))),
                               (expression(paste(10^{3.5}))), (expression(paste(10^{4}))))) +
 scale_x_continuous(name = expression(paste("Body Mass (kg)")),
                    breaks = seq(-1, 0.5, 0.5),
                    labels = c((expression(paste(10^{-1}))),(expression(paste(10^{-0.5}))),
                               (expression(paste(10^{0}))), (expression(paste(10^{0.5})))),
                    limits = c(-1.5, 0.65)) +
  # annotate("text", x = -0.1, y = 3, label = expression(paste('GSA')),
  #          color = "dark grey",  hjust = 0, size = 3.8) +
  annotate("text", x = -1.3, y = 2.95,  
           label = expression(paste(''~italic(b) ~'= 0.877 ± 0.067')),
           color = "grey60",  hjust = 0, size = 3.9) +
 annotate("text", x = -0.1, y = 3., label = expression(paste('> 0.203 kg')),
          color = "black",  hjust = 0, size = 3.8) +
 annotate("text", x = -0.1, y = 2.89,  
            label = expression(paste(''~italic(b) ~'= 1.012 ± 0.113')),
            color = "black",  hjust = 0, size = 3.8) +
 annotate("text", x = -1.0, y = 2.31, label = expression(paste('< 0.203 kg')),
          color = "black",  hjust = 0, size = 3.8) +
 annotate("text", x = -1.0, y = 2.2, 
            label = expression(paste(''~italic(b) ~'= 0.564 ± 0.261')),
            color = "black",  hjust = 0, size = 3.8) +
 annotate("text", x = -1.5, y = 4.05, label = "(b)",
          color = "black", hjust = 0, size = 3.8) +
 theme_bw() +
 theme(axis.title = element_text( size = 12),
       axis.text.x = element_text(size = 12),
       axis.text.y = element_text(size = 12),
       plot.margin = unit(c(-0.5,0,0,0), "cm"),
       panel.grid = element_blank()) +
 annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"),
                     short = unit(0.75, "mm"), colour = "black")

#Stack MR and GSA plots on top of each other
HS_MR_GSA_double_plot <- HS_MR_Plot  +  HS_GSA_Plot +  
  plot_layout(ncol = 1)
# Save plot
#ggsave(HS_MR_GSA_double_plot, filename = "HS_MR_GSA_double_plot_Apr30.png", bg = "transparent", 
#       height = 9.5, width = 5.7) 
      




######################## --------------------- ###############################
######################## --------------------- ###############################
######################## --------------------- ###############################


# FIGURE 2 - GSA:MMR ratio ~ Body Mass linear regression 

GSA_MR_LOG_ratio_plot <- 
ggplot(HS_calculate_ratios_adult_df, aes(x = MassKG_log, y = GSA_MMR_ratio_log10)) +
 geom_smooth(method = "lm", color = pal2mods$hex[1], se = FALSE) +
 geom_smooth(data = HS_calculate_ratios_adult_df, aes(x = MassKG_log, y = GSA_RMR_ratio_log10),
             method = "lm", color = pal2mods$hex[2], se = FALSE) +
  geom_point(show.legend = FALSE, size = 3.5, color = pal2mods$hex[1], alpha = 0.8) + 
geom_point(data = HS_calculate_ratios_adult_df, aes(x = MassKG_log, y = GSA_RMR_ratio_log10),
           color = pal2mods$hex[2], size = 3.5, alpha = 0.8) +
scale_y_continuous(name = expression('Ratio of GSA'~(cm^2)~ 'to Metabolic Rate (mg' ~O[2]~ h^-1~')'),
                   breaks = seq(0.5, 2, 0.5),
                   labels = c((expression(paste(10^{0.5}))), (expression(paste(10^{1}))),
                              (expression(paste(10^{1.5}))), (expression(paste(10^{2}))))) +
scale_x_continuous(name = expression(paste("Body Mass (kg)")),
                   breaks = seq(-1.5, 0.5, 0.5),
                   labels = c((expression(paste(10^{-1.5}))), (expression(paste(10^{-1}))),
                              (expression(paste(10^{-0.5}))),
                              (expression(paste(10^{0}))), (expression(paste(10^{0.5}))))) +
  annotate("text", x = -0.62, y = 1.95, label = "GSA/RMR",
           color = pal2mods$hex[2], hjust = 0, size = 3.8) +
  annotate("text", x = -0.62, y = 1.87, label = expression(paste(''~italic(b) ~'= 0.04 ± 0.188')),
           color = pal2mods$hex[2], hjust = 0, size = 3.8) +
  annotate("text", x = -0.62, y = 1.28, label = "GSA/MMR",
           color = pal2mods$hex[1], hjust = 0, size = 3.8) +
  annotate("text", x = -0.62, y = 1.2, label =  expression(paste(''~italic(b) ~'= -0.102 ± 0.122')),
           color = pal2mods$hex[1], hjust = 0, size = 3.8) +
  theme_bw() +
  theme(axis.title.y = element_text( size = 12),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank()) +
annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"),
                    short = unit(0.75, "mm"), colour = "black")
# Save plot
#ggsave(GSA_MR_LOG_ratio_plot, filename = "GSA_MR_LOG_ratio_plot.png", bg = "transparent",
#       height = 5.1, width = 5.7) 

 

######################## --------------------- ###############################
######################## --------------------- ###############################
######################## --------------------- ###############################


# FIGURE 3 - GILL COMPONENTS STACKED PLOT OF THREE

# ---- FILAMENT LENGTH --- #
 HS_GSA_FilLength_Plot <- 
  ggplot(GSA_raw_df, aes(x = log10(MassKG), y = log10(TotalFilamentLength2Sides))) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5, alpha = 0.7, show.legend = FALSE, color = "black") +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(name = expression('Total Fillament Length (cm)'),
                     breaks = seq(3.6, 4.2, 0.2),
                     labels = c((expression(paste(10^{3.6}))), (expression(paste(10^{3.8}))),
                                (expression(paste(10^{4}))), (expression(paste(10^{4.2}))))) +
  annotate("text", x = -0.75, y = 4.15, color = "black", hjust = 0, size = 3.8,
           label = expression(paste(''~italic(b) ~'= 0.40 ± 0.017'))) +
  annotate("text", x = -1.48, y = 4.29, label = "(a)",
           color = "black", hjust = 0, size = 3.8) +
  theme_bw() +
  theme(axis.title.y = element_text( size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0,0,-2,0), "cm"),
        panel.grid = element_blank()) +
  annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"), sides = "l",
                      short = unit(0.75, "mm"), colour = "black")

# --- LAMELLAR FREQUENCY --- #
 HS_GSA_LamFreq_Plot <- 
ggplot(GSA_raw_df, aes(x = log10(MassKG), y = log10(LamellarFrequency))) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5, alpha = 0.5, 
              show.legend = FALSE, color = "black") +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(name = expression('Lamellar Frequency ('~mm^-1~')'),
                     breaks = seq(0.9, 1.1, 0.1),
                     labels = c((expression(paste(10^{0.9}))), 
                                (expression(paste(10^{1.0}))), (expression(paste(10^{1.1}))))) +
   annotate("text", x = -0.52, y = 1.08, color = "black", hjust = 0, size = 3.8,
            label = expression(paste(''~italic(b) ~'= -0.103 ± 0.011'))) +
  annotate("text", x = -1.48, y = 0.905, 
           label = "(b)",
           color = "black", hjust = 0, size = 3.8) +
  theme_bw() +
  theme(axis.title.y = element_text( size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(-0.5,0,-2,0), "cm"),
        panel.grid = element_blank()) +
  annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"), sides = "l",
                      short = unit(0.75, "mm"), colour = "black")


# ---- LAMELLAR SURFACE AREA ---- #
HS_GSA_LamArea_Plot <- 
ggplot(GSA_raw_df, aes(x = log10(MassKG), y = log10(BilateralLamellarSA))) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5, alpha = 0.5, 
              show.legend = FALSE, color = "grey60") +
  geom_line(data = HS_LamSA_BR_lines_plot_df, aes(x = X_axis_spots, y = Y_axis_spots),
            size = 1.5, color = "black") +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(name = expression('Lamellar Surface Area ('~mm^2~')'),
                     breaks = seq(-0.4, 0.4, 0.4),
                     labels = c((expression(paste(10^{-0.4}))),  (expression(paste(10^{0}))),
                                (expression(paste(10^{0.4}))))) +
  scale_x_continuous(name = expression(paste("Body Mass (kg)")),
                     breaks = seq(-1, 0.5, 0.5),
                     labels = c((expression(paste(10^{-1}))), (expression(paste(10^{-0.5}))), (expression(paste(10^{0}))),
                                (expression(paste(10^{0.5}))))) +
  annotate("text", x = -1.3, y = -0.2,
           label = ''~italic(b) ~'= 0.577 ± 0.073', 
           color = "grey60",  hjust = 0, size = 3.9) +
  annotate("text", x = -0.05, y = -0.18, label = '> 0.203 kg',
           color = "black",  hjust = 0, size = 3.8) +
  annotate("text", x = -0.05, y = -0.26,  
           label = expression(''~italic(b) ~'= 0.742 ± 0.111'),
           color = "black",  hjust = 0, size = 3.8) +
  annotate("text", x = -1.0, y = -0.64, label = '< 0.203 kg',
           color = "black",  hjust = 0, size = 3.8) +
  annotate("text", x = -1.0, y = -0.72, 
           label = expression(''~italic(b) ~'= 0.196 ± 0.257'),
           color = "black",  hjust = 0, size = 3.8) +
annotate("text", x = -1.48, y = 0.52, 
            label = "(c)",
            color = "black", hjust = 0, size = 3.8) +
  theme_bw() +
  theme(axis.title = element_text( size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(-0.5,0,0,0), "cm"),
        panel.grid = element_blank()) +
  annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"), 
                      short = unit(0.75, "mm"), colour = "black") 

# Join all 4 GSA components plots together into one long plot for the MS
GSA_components_strip_plot <- HS_GSA_FilLength_Plot + 
  HS_GSA_LamFreq_Plot + HS_GSA_LamArea_Plot +  plot_layout(ncol = 1)

# Save plot
#ggsave(GSA_components_strip_plot, filename = "GSA_components_plots.png", bg = "transparent", 
#       height = 13, width = 5.7) 


######################## --------------------- ###############################
######################## --------------------- ###############################
######################## --------------------- ###############################

# Figure S-1 - MMR methods comparison 

HS_MMR_Plot <- 
ggplot(HS_MMR_chaseVSair_df, aes(x = MassKG_log, y = WholeEstimate_log_KG, 
                                 group = MeasureType, color = MeasureType)) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5, alpha = 0.3, show.legend = FALSE) + 
  geom_point(size = 3.5, alpha = 0.7, show.legend = FALSE) +  
  scale_y_continuous(name = expression('Metabolic Rate (mg'~O[2]~ h^-1~')'),
                     breaks = seq(1, 3, 1),
                     labels = c((expression(paste(10^{1}))), (expression(paste(10^{2}))),
                                (expression(paste(10^{3}))))) +
  scale_x_continuous(name = expression(paste("Body Mass (kg)")),
                     breaks = seq(-1, 0.5, 0.5),
                     labels = c((expression(paste(10^{-1}))),(expression(paste(10^{-0.5}))),
                                (expression(paste(10^{0}))), (expression(paste(10^{0.5})))),
                     limits = c(-1.5, 0.65)) +
  annotate("text", x = -1.3, y = 2.1, label = expression(paste('MMR Chase')),
           color = "black", hjust = 0, size = 3.8) +
  annotate("text", x = -1.3, y = 1.95, 
           label = expression(paste(~italic(b) ~'= 1.073 ± 0.040')),
           color = "black", hjust = 0, size = 3.8) +
  
  annotate("text", x = -0.8, y = 1.2, label = expression(paste('MMR Chase + Air')),
           color = "black", hjust = 0, size = 3.8) +
  annotate("text", x = -0.8, y = 1.05, 
           label = expression(paste(~italic(b) ~'= 1.119 ± 0.045')),
           color = "black", hjust = 0, size = 3.8) +
  theme_bw() +
  theme(axis.title = element_text( size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.05, vjust = -8, color = "grey30", size = 12),
        panel.grid = element_blank()) +
  annotation_logticks(long = unit(2.5, "mm"), mid = unit(1.5, "mm"), 
                      short = unit(0.75, "mm"), colour = "black") 
# Save plot
#ggsave(HS_MMR_Plot, filename = "HS_MMR_Plot.png", height = 5.2, width = 5.7) 




