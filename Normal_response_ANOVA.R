# This script was used to perform ANOVA for log-transformed RMSE change that was equal to RMSE(test)/RMSE(control)
# Close scientific notation for a while
options(scipen=999)
# Turn it back on
# options(scipen=0) 

# Load package
library(car)
library(ggpubr)
library(nortest)
library(heplots)


# Load data
load("data_analysis_normal_GLM_no_RF_raw_RMSE.RData")

# Log-transformed RMSE
data_anova <- data_analysis_normal
data_anova$RMSE <- log10(data_anova$RMSE)

# ANOVA that violates the assumption but return Type III sum of squares
# [1] "2023-05-07 17:30:16 EDT"
# [1] "2023-05-07 18:11:20 EDT"
# Over 11GB object for 1251000 observations

# Need 80GB RAM to be safe if allowed
Sys.time()
anova_all <-  lm(RMSE ~ Model + Predictor_coli + Range_shift + Collinearity_shift +
                   Model:Predictor_coli + Model:Range_shift + Predictor_coli:Range_shift + Model:Collinearity_shift +
                   Predictor_coli:Collinearity_shift + Range_shift:Collinearity_shift,
                 contrasts=list(Model = contr.sum, Predictor_coli = contr.sum, 
                                Range_shift = contr.sum, Collinearity_shift = contr.sum),
                 data = data_anova)
Sys.time()


# Check out the assumption of homogeneity of residuals
hist(resid(anova_all), xlab = "Residuals", main = "")

# Run ANOVA using type III sum of squares with an interaction to get partial R squares
options(contrasts = c("contr.sum","contr.poly")) # This has to be run before fitting lm or glm
anova_all_type3 <- Anova(anova_all, type = "III")
anova_all_type3

options(contrasts = c("contr.sum","contr.poly"))
anova_all_glm_type3 <- Anova(anova_all_glm, type = "III")
anova_all_glm_type3

Sys.time()
library(heplots)
anova_all_type3_check <- etasq(anova_all, type = 3, anova = TRUE, partial = TRUE)
anova_all_type3_check

anova_all_type2_check <- etasq(anova_all, type = 2, anova = TRUE, partial = TRUE)
anova_all_type2_check
Sys.time()



# Perform ANOVA for lower range shifts
# Run ANOVA using type III sum of squares with an interaction to get partial R squares
options(contrasts = c("contr.sum","contr.poly")) # This has to be run before fitting lm or glm

Sys.time()
anova_low_range_shift <- lm(RMSE ~ Model + Predictor_coli + Range_shift + Collinearity_shift +
                              Model:Predictor_coli + Model:Range_shift + Predictor_coli:Range_shift + Model:Collinearity_shift +
                              Predictor_coli:Collinearity_shift + Range_shift:Collinearity_shift,
                            contrasts=list(Model = contr.sum, Predictor_coli = contr.sum, 
                                           Range_shift = contr.sum, Collinearity_shift = contr.sum),
                            data = data_anova[data_anova$Range_shift %in% c(0, 0.2, 0.4, 0.6), ]) 
Sys.time()

anova_low_range_shift_type3 <- Anova(anova_low_range_shift, type = 3)
anova_low_range_shift_type3_check <- etasq(anova_low_range_shift, type = 3, anova = TRUE, partial = TRUE)

library(effectsize)
options(es.use_symbols = TRUE)
eta_squared(car::Anova(anova_low_range_shift, type = 3))

library(sjstats)
options(es.use_symbols = TRUE)
anova_stats(car::Anova(anova_low_range_shift, type = 3))

Sys.time()
anova_low_range_shift_1 <- lm(RMSE ~ Model + Predictor_coli + Range_shift + Collinearity_shift +
                                Model:Predictor_coli + Model:Range_shift + Predictor_coli:Range_shift + Model:Collinearity_shift +
                                Predictor_coli:Collinearity_shift + Range_shift:Collinearity_shift,
                              contrasts=list(Model = contr.sum, Predictor_coli = contr.sum, 
                                             Range_shift = contr.sum, Collinearity_shift = contr.sum),
                              data = data_anova[data_anova$Range_shift %in% c(0, 0.2, 0.4), ]) 

Sys.time()

anova_low_range_shift_1_type3 <- Anova(anova_low_range_shift_1, type = 3)
anova_low_range_shift_1_type3_check <- etasq(anova_low_range_shift_1, type = 3, anova = TRUE, partial = TRUE)


Sys.time()
anova_low_range_shift_2 <- lm(RMSE ~ Model + Predictor_coli + Range_shift + Collinearity_shift +
                                Model:Predictor_coli + Model:Range_shift + Predictor_coli:Range_shift + Model:Collinearity_shift +
                                Predictor_coli:Collinearity_shift + Range_shift:Collinearity_shift,
                              contrasts=list(Model = contr.sum, Predictor_coli = contr.sum, 
                                             Range_shift = contr.sum, Collinearity_shift = contr.sum),
                              data = data_anova[data_anova$Range_shift %in% c(0, 0.2), ]) 

Sys.time()

anova_low_range_shift_2_type3 <- Anova(anova_low_range_shift_2, type = 3)
anova_low_range_shift_2_type3_check <- etasq(anova_low_range_shift_2, type = 3, anova = TRUE, partial = TRUE)






# High training collinearity and low collinearity shift
Sys.time()
anova_low_range_shift_3 <- lm(RMSE ~ Range_shift + Collinearity_shift + Range_shift:Collinearity_shift,
                              contrasts=list(Range_shift = contr.sum, Collinearity_shift = contr.sum),
                              data = data_anova[data_anova$Range_shift %in% c(0, 0.2)&
                                                  data_anova$Predictor_coli %in% c("High")&
                                                  data_anova$Model %in% c("Product"), ]) 

Sys.time()

anova_low_range_shift_3_type3 <- Anova(anova_low_range_shift_3, type = 3)
anova_low_range_shift_3_type3_check <- etasq(anova_low_range_shift_3, type = 3, anova = TRUE, partial = TRUE)


# Export data
write.csv(anova_all_type3_check, "anova_all_type3_check.csv")
write.csv(anova_low_range_shift_type3_check, "anova_low_range_shift_type3_check.csv")
write.csv(anova_low_range_shift_1_type3_check, "anova_low_range_shift_1_type3_check.csv")
write.csv(anova_low_range_shift_2_type3_check, "anova_low_range_shift_2_type3_check.csv")

# Save the object
save(anova_all, file = "data_analysis_normal_GLM_no_RF_anova_2_way.RData")




# Try glm + permutation test
Sys.time()
anova_all_glm <-  glm(RMSE ~ Model + Predictor_coli + Range_shift + Collinearity_shift + 
                        Model:Predictor_coli + Model:Range_shift + Predictor_coli:Range_shift + Model:Collinearity_shift +
                        Predictor_coli:Collinearity_shift + Range_shift:Collinearity_shift,
                      family = gaussian,
                      contrasts=list(Model = contr.sum, Predictor_coli = contr.sum, 
                                     Range_shift = contr.sum, Collinearity_shift = contr.sum),
                      data = data_anova)
Sys.time()

# Permutation test on all coefficients
library(prettyglm)
Sys.time()
pretty_coefficients(model_object = anova_all,
                    type_iii = 'LR', 
                    significance_level = 0.05, 
                    vimethod = 'permute', 
                    target = 'RMSE', 
                    metric = 'auc',
                    return_data = FALSE,
                    pred_wrapper = predict.glm, 
                    reference_class = 0)
Sys.time()



data(Soils)  # from car package
soils.mod <- lm(cbind(pH,N,Dens,P,Ca,Mg,K,Na,Conduc) ~ Block + Contour*Depth, data=Soils)
#Anova(soils.mod)
etasq(Anova(soils.mod))
etasq(soils.mod) # same
etasq(Anova(soils.mod), anova=TRUE)

# Example 1: one-way ANOVA

outcome <- c( 1.4,2.1,3.0,2.1,3.2,4.7,3.5,4.5,5.4 )  # data
treatment1 <- factor( c( 1,1,1,2,2,2,3,3,3 ))        # grouping variable
anova1 <- aov( outcome ~ treatment1 )                # run the ANOVA
summary( anova1 )                                    # print the ANOVA table
etaSquared( anova1 )                                 # effect size

# Example 2: two-way ANOVA

treatment2 <- factor( c( 1,2,3,1,2,3,1,2,3 ))      # second grouping variable
anova2 <- aov( outcome ~ treatment1 + treatment2 + treatment1*treatment2 ) # run the ANOVA
summary( anova2 )                                  # print the ANOVA table
etaSquared( anova2 )         


# Perform ANOVA for lower range shifts
Sys.time()
anova_low_range_shift <- lm(RMSE ~ Model*Predictor_coli*Range_shift*Collinearity_shift, 
                            data = data_anova[data_anova$Range_shift %in% c(0, 0.2, 0.4, 0.6), ]) 
  
Sys.time()

anova_low_range_shift_type3 <- Anova(anova_low_range_shift, type = 3)
anova_low_range_shift_type3_check <- etasq(anova_low_range_shift)

# Perform ANOVA for 0 and 0.2 range shifts
Sys.time()
anova_0_0_2_range_shift <- lm(RMSE ~ Model*Predictor_coli*Range_shift*Collinearity_shift,
                              data = data_anova[data_anova$Range_shift %in% c(0, 0.2), ]) 

Sys.time()

anova_0_0_2_range_shift_type3 <- Anova(anova_0_0_2_range_shift, type = 3)
anova_0_0_2_range_shift_type3_check <- etasq(anova_0_0_2_range_shift)


#####################################################################################################################################
#####################################################################################################################################
# Run ANOVA in each level of model complexity and gradually added one level of novelty in each run
# because novelty accounted for the most variation and the interaction between model complexity and novelty came to the 2nd
# model complexity was the 3rd
# Therefore, within each model complexity, other the variation explained by factors of interest may be further investigated

# GLM models
# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_GLM_A_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- as.data.frame(etasq(anova_log10_RMSE_factors_GLM_A_product))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_GLM_B_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_B_product))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_GLM_C_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_C_product))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_GLM_D_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_D_product))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_GLM_E_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_E_product))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_GLM_F_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                               Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                               Range_shift*Collinearity_shift,
                                             data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                                 data_anova$Model == "Product"&
                                                                 data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_F_product))

# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_GLM_A_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_A_quadratic))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_GLM_B_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_B_quadratic))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_GLM_C_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_C_quadratic))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_GLM_D_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_D_quadratic))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_GLM_E_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_E_quadratic))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_GLM_F_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                 Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                 Range_shift*Collinearity_shift,
                                               data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                                   data_anova$Model == "Quadratic"&
                                                                   data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_F_quadratic))

# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_GLM_A_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_A_Linear))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_GLM_B_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_B_Linear))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_GLM_C_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_C_Linear))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_GLM_D_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_D_Linear))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_GLM_E_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_E_Linear))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_GLM_F_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                                data_anova$Model == "Linear"&
                                                                data_anova$Algorithm == "GLM", ])
GLM_partial_R <- rbind(GLM_partial_R, etasq(anova_log10_RMSE_factors_GLM_F_Linear))

# Removed duplicated row names and NAs
GLM_partial_R$Factors <- rep(c("Training collinearity", "Predictor novelty", "Collinearity shift", 
                               "Training collinearity x Predictor novelty",
                               "Training collinearity x Collinearity shift",
                               "Predictor novelty x Collinearity shift",
                               "Residuals"), 18)
GLM_partial_R <- GLM_partial_R[!is.na(GLM_partial_R$`Partial eta^2`), ]
row.names(GLM_partial_R) <- 1:108
GLM_partial_R$ANOVA <- c(rep(c(rep("A", 6),
                               rep("B", 6),
                               rep("C", 6),
                               rep("D", 6),
                               rep("E", 6),
                               rep("F", 6)), 3))

GLM_partial_R$Model_complexity <- c(rep("Product relationship", 36), 
                                    rep("Quadratic relationship", 36), 
                                    rep("Linear relationship", 36))
names(GLM_partial_R)[1] <- "value" # For ggplot syntax afterwards only

# RF models
# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_RF_A_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- as.data.frame(etasq(anova_log10_RMSE_factors_RF_A_product))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_RF_B_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_B_product))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_RF_C_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_C_product))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_RF_D_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_D_product))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_RF_E_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_E_product))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_RF_F_product <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                              Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                              Range_shift*Collinearity_shift,
                                            data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                                data_anova$Model == "Product"&
                                                                data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_F_product))

# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_RF_A_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_A_quadratic))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_RF_B_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_B_quadratic))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_RF_C_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_C_quadratic))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_RF_D_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_D_quadratic))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_RF_E_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_E_quadratic))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_RF_F_quadratic <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                                Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                                Range_shift*Collinearity_shift,
                                              data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                                  data_anova$Model == "Quadratic"&
                                                                  data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_F_quadratic))

# Warning: subset Range_shift using %in% to make sure all the cases were selected
# A: Range_shift == c(0, 0.2)
anova_log10_RMSE_factors_RF_A_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_A_Linear))

# B: Range_shift == c(0, 0.2, 0.4)
anova_log10_RMSE_factors_RF_B_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2, .4)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_B_Linear))

# C: Range_shift == c(0, 0.2, 0.4, 0.6)
anova_log10_RMSE_factors_RF_C_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_C_Linear))

# D: Range_shift == c(0, 0.2, 0.4, 0.6, 2)
anova_log10_RMSE_factors_RF_D_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_D_Linear))

# E: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4)
anova_log10_RMSE_factors_RF_E_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_E_Linear))

# F: Range_shift == c(0, 0.2, 0.4, 0.6, 2, 4, 6)
anova_log10_RMSE_factors_RF_F_Linear <- lm(RMSE ~ Predictor_coli + Range_shift + Collinearity_shift +
                                             Predictor_coli*Range_shift + Predictor_coli*Collinearity_shift +
                                             Range_shift*Collinearity_shift,
                                           data = data_anova[data_anova$Range_shift %in% c(0, .2, .4, .6, 2, 4, 6)&
                                                               data_anova$Model == "Linear"&
                                                               data_anova$Algorithm == "RF", ])
RF_partial_R <- rbind(RF_partial_R, etasq(anova_log10_RMSE_factors_RF_F_Linear))

# Removed duplicated row names and NAs
RF_partial_R$Factors <- rep(c("Training collinearity", "Predictor novelty", "Collinearity shift", 
                              "Training collinearity x Predictor novelty",
                              "Training collinearity x Collinearity shift",
                              "Predictor novelty x Collinearity shift",
                              "Residuals"), 18)
RF_partial_R <- RF_partial_R[!is.na(RF_partial_R$`Partial eta^2`), ]
row.names(RF_partial_R) <- 1:108
RF_partial_R$ANOVA <- c(rep(c(rep("A", 6),
                              rep("B", 6),
                              rep("C", 6),
                              rep("D", 6),
                              rep("E", 6),
                              rep("F", 6)), 3))

RF_partial_R$Model_complexity <- c(rep("Product relationship", 36), 
                                   rep("Quadratic relationship", 36), 
                                   rep("Linear relationship", 36))
names(RF_partial_R)[1] <- "value" # For ggplot syntax afterwards only

