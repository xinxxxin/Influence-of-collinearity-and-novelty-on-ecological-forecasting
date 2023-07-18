# This script is used to create heat maps for visualizing GLMactions

# Load packages
library(ggplot2)
library(metR)

# remotes::install_github('rpkgs/gg.layers')
# for element_markdown
library(gg.layers)

library(cowplot)

# Load data
load("C:/FSU/Collinearity_extra/Revision/Normal_response/data_analysis_normal_GLM_no_RF_raw_RMSE.RData")

# Color palettee
library(paletteer)
d  <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", n = 46, direction = 1)

# Calculate the median
data_heatmap_GLM <- data_analysis_normal[data_analysis_normal$Algorithm == "GLM", ]

data_heatmap_GLM <- aggregate(cbind(RMSE) ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift,
                                        data = data_heatmap_GLM,
                                        FUN = median)


# Create the starting points of the predictor collinearity
data_training_point <- data.frame(Predictor_coli = unique(data_heatmap_GLM$Predictor_coli),
                                  x = c(0.9, 0.6, 0.3),
                                  y = c(0, 0, 0),
                                  RMSE = c(1.0, 1.0, 1.0))

Model_labs <- c("Product complexity", "Quadratic complexity", "Linear complexity")
names(Model_labs) <- c("Product", "Quadratic", "Linear")

Predictor_coli_labs <- c("High training collinearity<br />(*r* = 0.9)", 
                         "Mid training collinearity<br />(*r* = 0.6)", 
                         "Low training collinearity<br />(*r* = 0.3)")
names(Predictor_coli_labs) <- c("High", "Mid", "Low")

heatmap_normal_GLM <- ggplot(data_heatmap_GLM[data_heatmap_GLM$Range_shift %in% c(0, 0.2, 0.4, 0.6), ], mapping = aes(x = as.numeric(as.character(Collinearity_shift)), y = as.numeric(as.character(Range_shift)), z = RMSE)) +
  # geom_tile(aes(fill = RMSE_ratio)) +
  geom_contour_filled(aes(z = RMSE), binwidth = 0.05) +
  scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = d, name = "RMSE") +
  # metR::geom_text_contour(aes(z = RMSE), stroke = 0.15) +
  metR::geom_contour2(aes(z = RMSE, label = ..level..), binwidth = 0.25, skip = 0, label_size = 3,
                      label.placer = label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()
                      )) +
  scale_x_continuous(name = "Correlation between X<sub>1</sub> and X<sub>2</sub>", breaks = c(-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9),
                     labels = c("-0.9", "-0.6", "-0.3", "0.0", "0.3", "0.6", "0.9")) +
  scale_y_continuous(name = "Factor of increased predictor novelty<br />") +
  #facet_grid(Algorithm~Predictor_coli) +
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
        strip.text.y = element_markdown(face = "bold"),
        panel.border = element_blank(),
        legend.position = "none",
        legend.title = element_markdown(size = 8, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.4, 'cm')) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + 
  geom_point(data = data_training_point, aes(x = x, y = y), color = "blue", shape = 3, stroke = 1.5)


# Save the heatmap for GLM
ggsave("Normal_RMSE_GLM_heatmap.png", width = 14, height = 14, units = "cm", dpi = 200)
#unlink("Normal_RMSE_heatmap.png")


# Calculate the median
data_heatmap_RF <- data_analysis_normal[data_analysis_normal$Algorithm == "RF", ]

data_heatmap_RF <- aggregate(cbind(RMSE) ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift,
                             data = data_heatmap_RF,
                             FUN = median)

# Color palettee
e  <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", n = 446, direction = 1)

heatmap_normal_RF <- ggplot(data_heatmap_RF[data_heatmap_RF$Range_shift %in% c(0, 0.2, 0.4, 0.6), ], mapping = aes(x = as.numeric(as.character(Collinearity_shift)), y = as.numeric(as.character(Range_shift)), z = RMSE)) +
  # geom_tile(aes(fill = RMSE_ratio)) +
  geom_contour_filled(aes(z = RMSE), binwidth = 0.05) +
  scale_colour_manual(aesthetics = 'fill', drop = FALSE, values = e, name = "RMSE") +
  # metR::geom_text_contour(aes(z = RMSE), stroke = 0.15) +
  metR::geom_contour2(aes(z = RMSE, label = ..level..), binwidth = 2.5, skip = 0, label_size = 3,
                      label.placer = label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()
                      )) +
  scale_x_continuous(name = "Correlation between X<sub>1</sub> and X<sub>2</sub>", breaks = c(-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9),
                     labels = c("-0.9", "-0.6", "-0.3", "0.0", "0.3", "0.6", "0.9")) +
  scale_y_continuous(name = "Factor of increased predictor novelty<br />") +
  #facet_grid(Algorithm~Predictor_coli) +
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
  geom_point(data = data_training_point, aes(x = x, y = y), color = "blue", shape = 3, stroke = 1.5)

# Save the heatmap for RF
ggsave("Normal_RMSE_RF_heatmap.png", width = 14, height = 14, units = "cm", dpi = 200)
#unlink("Normal_RMSE_heatmap.png")

heatmap_normal_GLM_RF <- plot_grid(heatmap_normal_GLM, heatmap_normal_RF,
                                   labels = c("(a)", "(b)"),
                                   label_size = 10,
                                   rows = 2,
                                   rel_heights = c(1, 1))

# Save the heatmap for GLM+RF
ggsave("Normal_RMSE_GLM_RF_heatmap.png", width = 16, height = 28, units = "cm", dpi = 200)
#unlink("Normal_RMSE_heatmap.png")
