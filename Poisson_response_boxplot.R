# This script was used to investigate the pattern of collinearity shift across different levels
# of training collinearity, model complexity, and model algorithm
# Range shift = 0
# Color palette needs to be interpolated because the number of categories is greater than 9/11.
# Source link of color palette interpolation: https://www.r-bloggers.com/2013/09/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/

# Load packages
library(RColorBrewer)
library(gapminder)
library(ggplot2)
library(cowplot)

# remotes::install_github('rpkgs/gg.layers')
# for element_markdown
library(gg.layers) 

# Color palettee
library(paletteer)

library(khroma)

# Load data
load("data_analysis_poisson_GLM_no_RF.RData")
data_analysis_poisson <- data_analysis_poisson

# Color palette interpolation
colourCount <- length(unique(data_analysis_poisson$Collinearity_shift))
getPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))
getPalette <- colorRampPalette(brewer.pal(11, "PuOr"))

# getPalette <- colorRampPalette(brewer.pal(11, "BrBG"))
# getPalette <- colorRampPalette(brewer.pal(11, "RdBu"))
#getPalette <- colorRampPalette(brewer.pal(9, "Set1"))


library(paletteer)
a  <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", n = 19, direction = 1)


library(khroma)
smooth_rainbow <- colour("smooth rainbow")
plot(smooth_rainbow(256, range = c(.25, 1)))
plot(smooth_rainbow(19, range = c(.25, 1)))

# End at brown
b <- smooth_rainbow(19, range = c(.25, 1))
plot(b)

# End at red
c <- smooth_rainbow(19, range = c(0.25, 0.7))
plot(c)

# Plotting collinearity shift when range shift = 0
library(scales)

# Change the decimal places before creating boxplots
# Label formatter for free scales in ggplot2
# Source: https://stackoverflow.com/questions/37669064/facet-wrap-title-wrapping-decimal-places-on-free-y-axis-ggplot2

formatter <- function(...){
  function(x) format(round(x, 2), ...)
}

# Function to produce summary statistics (mean and +/- sd)
data_summary_sd <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x) #quantile(x, 0.95)
  ymax <- m + sd(x) #quantile(x, 0.05)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}  

# Function to produce summary statistics for lower and upper 5%
data_summary <- function(x) {
  m <- median(x)
  ymin <- quantile(x, 0.05)
  ymax <- quantile(x, 0.95)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}

# Function to produce summary statistics for 1st and 3rd quantile
data_summary_1 <- function(x) {
  m <- median(x)
  ymin <- quantile(x, 0.25)
  ymax <- quantile(x, 0.75)
  data.frame(y = m, ymin = ymin, ymax = ymax)
}


# One quick example for RMSE vs Collinearity shift when training collinearity
data_analysis_poisson_high_interaction_GLM <- data_analysis_poisson[data_analysis_poisson$Predictor_coli == "High"&data_analysis_poisson$Range_shift %in% c(0)&
                                                                    data_analysis_poisson$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                                    (data_analysis_poisson$Algorithm == "GLM")&
                                                                    data_analysis_poisson$Model == "Interaction", ]

ggplot(data = data_analysis_poisson_high_interaction_GLM, aes(Range_shift, RMSE, fill = Collinearity_shift)) + 
  geom_boxplot(colour = "black", position=position_dodge(1), coef = 10) +
  scale_y_continuous(name = "RMSE", trans = 'log10', labels = formatter(nsmall = 1)) +
  scale_fill_manual(name = "Correlation between\nX1 and X2", values = b) +
  theme(legend.position="bottom")

ggplot(data = data_analysis_poisson_high_interaction_GLM, aes(Collinearity_shift, RMSE, fill = Collinearity_shift)) + 
  stat_summary(fun.data = data_summary, color = b, geom = "errorbar", width = 0.5, linewidth = 1) +
  geom_boxplot(colour = "black", outlier.shape = NA) +
  scale_y_continuous(name = "RMSE", trans = 'log10', labels = formatter(nsmall = 1)) +
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2))) +
  scale_fill_manual(name = "Correlation between\nX1 and X2", values = b) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
# Jitter plots might be better than box plot this time
ggplot(data = data_analysis_poisson_high_interaction_GLM, aes(Collinearity_shift, RMSE, color = Collinearity_shift)) +
  stat_summary(fun.data = data_summary, color = "blue", size = 0.5, shape = 20, stroke = 1.5, linewidth = 1) +
  geom_jitter(width = 0.35, size = 0.001, alpha = 0.3) +
  scale_y_continuous(name = "RRMSE", trans = 'log10', labels = formatter(nsmall = 1)) +
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2)), name = "Correlation between\nX1 and X2") +
  scale_color_manual(name = "Correlation between\nX1 and X2", values = b) +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  stat_summary(fun.data = data_summary, color = "blue", size = 0.5, shape = 20, stroke = 1.5, linewidth = 1)

# Only use stat_summary with (linerange, errorbar, and pointrange geom) or (errorbar and crossbar geom)
ggplot(data = data_analysis_poisson_high_interaction_GLM, aes(Collinearity_shift, RMSE, color = Collinearity_shift)) +
  #stat_summary(fun.data = data_summary, geom = "linerange", size = 1) +
  stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.5, color = b, linewidth = 1) +
  stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.5, fill = b, color = "black") + 
  scale_y_continuous(name = "RRMSE", trans = 'log', labels = formatter(nsmall = 1), breaks = c(1.0, 3.0, 5.0)) +
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2)), name = "Correlation between X<sub>1</sub> and X<sub>2</sub>") +
  scale_color_manual(name = "Correlation between\nX1 and X2", values = b) +
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_markdown(size = 9),
        axis.title.y = element_markdown(size = 9))


# Apply the stat_summary with errorbar and cross bar geom to all collinearity shift

# Create the dataframe for different training collinearity
vline <- data.frame(Predictor_coli = factor(rep(levels(data_analysis_poisson$Predictor_coli), 1), levels = c(levels(data_analysis_poisson$Predictor_coli))),
                    vline = factor(rep(c("0.9", "0.6", "0.3"), 1), levels = c("0.9", "0.6", "0.3")))


Model_labs <- c("Product complexity", "Quadratic complexity", "Linear complexity")
names(Model_labs) <- c("Product", "Quadratic", "Linear")
  
Predictor_coli_labs <- c("High training collinearity<br />(*r* = 0.9)", 
                         "Mid training collinearity<br />(*r* = 0.6)", 
                         "Low training collinearity<br />(*r* = 0.3)")
names(Predictor_coli_labs) <- c("High", "Mid", "Low")

GLM_col <- ggplot(data = data_analysis_poisson[data_analysis_poisson$Range_shift %in% c(0)&
                                                data_analysis_poisson$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                (data_analysis_poisson$Algorithm == "GLM"), ],
                  aes(Collinearity_shift, RMSE, color = Collinearity_shift, fill = Collinearity_shift)) +
  stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.5, linewidth = .5) +
  stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.5, color = "black") + 
  scale_y_continuous(name = "RMSE<br>", labels = formatter(nsmall = 1)) +
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2)), name = "Correlation between X<sub>1</sub> and X<sub>2</sub>") +
  scale_color_manual(name = "Correlation between <br /> X<sub>1</sub> and X<sub>2</sub>", values = rep("black", 19)) +
  scale_fill_manual(name = "Correlation between <br /> X<sub>1</sub> and X<sub>2</sub>", values = rep("darkgrey", 19)) +
  facet_grid(Model~Predictor_coli, scales = "free_y", 
             labeller = labeller(Model = Model_labs, Predictor_coli = Predictor_coli_labs)) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_markdown(size = 9),
        axis.title.y = element_markdown(size = 9)) +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text = element_markdown(face = "bold"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_markdown(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  geom_vline(data = vline, mapping = aes(xintercept = vline), 
             linetype = "dashed", linewidth = 0.5, color = "black")

RF_col <- ggplot(data = data_analysis_poisson[data_analysis_poisson$Range_shift %in% c(0)&
                                               data_analysis_poisson$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                               (data_analysis_poisson$Algorithm == "RF"), ],
                 aes(Collinearity_shift, RMSE, color = Collinearity_shift, fill = Collinearity_shift)) +
  stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.5, linewidth = 1) +
  stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.5, color = "black") + 
  scale_y_continuous(name = "Fold change of RMSE<br>", labels = formatter(nsmall = 0)) +
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2)), name = "Correlation between X<sub>1</sub> and X<sub>2</sub>") +
  scale_color_manual(name = "Correlation between <br /> X<sub>1</sub> and X<sub>2</sub>", values = c) +
  scale_fill_manual(name = "Correlation between <br /> X<sub>1</sub> and X<sub>2</sub>", values = c) +
  facet_grid(Model~Predictor_coli, scales = "free_y", 
             labeller = labeller(Model = Model_labs, Predictor_coli = Predictor_coli_labs)) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_markdown(size = 9),
        axis.title.y = element_markdown(size = 9)) +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text = element_markdown(face = "bold"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_markdown(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  geom_vline(data = vline, mapping = aes(xintercept = vline), 
             linetype = "dotted", linewidth = 1, color = "darkgrey")


Poisson_GLM_RF_col <- plot_grid(GLM_col, RF_col,
                               labels = c("(a)", "(b)"),
                               label_size = 10,
                               rows = 2,
                               rel_heights = c(1, 1))

# GLM only
# Save the plot for collinearity shift
ggsave("Poisson_RMSE_GLM_col.png", width = 16, height = 14, units = "cm", dpi = 200)
#unlink("Poisson_RMSE_GLM_col.png")


# Apply the stat_summary with errorbar and cross bar geom to all range shift
# Color palette interpolation
# colourCount_range <- length(unique(data_analysis_poisson$Range_shift))
# getPalette_range <- colorRampPalette(brewer.pal(7, "RdYlGn"))
d  <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", n = 10, direction = 1)
e  <- paletteer_d("RColorBrewer::PRGn", n = 11, direction = -1)

Model_labs <- c("Product complexity", "Quadratic complexity", "Linear complexity")
names(Model_labs) <- c("Product", "Quadratic", "Linear")

Predictor_coli_labs <- c("High training collinearity<br />(*r* = 0.9)", 
                         "Mid training collinearity<br />(*r* = 0.6)", 
                         "Low training collinearity<br />(*r* = 0.3)")
names(Predictor_coli_labs) <- c("High", "Mid", "Low")


formatter <- function(...){
  function(x) format(round(x, 4), ...)
}

range_GLM <- data_analysis_poisson[data_analysis_poisson$Range_shift %in% c(0, 0.2, 0.4, 0.6, 2, 4, 6)&
                                    data_analysis_poisson$Collinearity_shift %in% c(0.9, 0.6, 0.3)&
                                    (data_analysis_poisson$Algorithm == "GLM"), ]

GLM_nol <- ggplot(data = data_analysis_poisson[data_analysis_poisson$Range_shift %in% c(0, 0.2, 0.4, 0.6, 2, 4, 6)&
                                                data_analysis_poisson$Collinearity_shift %in% c(0.9, 0.6, 0.3)&
                                                (data_analysis_poisson$Algorithm == "GLM"), ],
                  aes(Range_shift, RMSE, fill = Range_shift)) + 
  stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.5, linewidth = 0.5) +
  stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.5, color = "black") + 
  scale_y_continuous(name = "RMSE<br>", labels = formatter(nsmall = 0)) +
  scale_x_discrete(breaks = c(0, .2, .4, .6, 2, 4, 6), name = "<br>Factor of increased predictor novelty") +
  scale_color_manual(name = "Factor of increased<br />predictor novelty", values = rep("black", 7)) +
  scale_fill_manual(name = "Factor of increased<br />predictor novelty", values = rep("darkgrey", 7)) +
  facet_grid(Model~Predictor_coli, scales = "free_y", 
             labeller = labeller(Model = Model_labs, Predictor_coli = Predictor_coli_labs)) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_markdown(size = 9),
        axis.title.y = element_markdown(size = 9)) +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text = element_markdown(face = "bold"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_markdown(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


RF_nol <- ggplot(data = data_analysis_poisson[data_analysis_poisson$Range_shift %in% c(0, 0.2, 0.4, 0.6, 2, 4, 6)&
                                               data_analysis_poisson$Collinearity_shift %in% c(0.9, 0.6, 0.3)&
                                               (data_analysis_poisson$Algorithm == "RF"), ],
                 aes(Range_shift, RMSE, fill = Range_shift)) + 
  stat_summary(fun.data = data_summary, geom = "errorbar", width = 0.5, linewidth = 1) +
  stat_summary(fun.data = data_summary_1, geom = "crossbar", width = 0.5, color = "black") + 
  scale_y_continuous(name = "RMSE<br>", labels = formatter(nsmall = 0)) +
  scale_x_discrete(breaks = c(0, .2, .4, .6, 2, 4, 6), name = "<br>Factor of increased predictor novelty") +
  scale_color_manual(name = "Factor of increased<br />predictor novelty", values = d[c(1:2, 4, 6, 7:9)]) +
  scale_fill_manual(name = "Factor of increased<br />predictor novelty", values = d[c(1:2, 4, 6, 7:9)]) +
  facet_grid(Model~Predictor_coli, scales = "free_y", 
             labeller = labeller(Model = Model_labs, Predictor_coli = Predictor_coli_labs)) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_markdown(size = 9),
        axis.title.y = element_markdown(size = 9)) +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "grey95", color = "white"),
        strip.text = element_markdown(face = "bold"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_markdown(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))


Poisson_GLM_RF_mol <-  plot_grid(GLM_nol, RF_nol,
                                labels = c("(a)", "(b)"),
                                label_size = 10,
                                rows = 2,
                                rel_heights = c(1, 1))

# GLM only
# Save the plot for predictor novelty
ggsave("Poisson_RMSE_GLM_nol.png", width = 16, height = 14, units = "cm", dpi = 200)
#unlink("Poisson_RMSE_GLM_nol.png")



