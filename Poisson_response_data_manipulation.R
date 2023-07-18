# Load the package just in case
library(reshape2)

# Load the transpose function
source("Poisson_response_transpose_function.R")


# Load the results of experiments for normal response
load("Poisson_response_theta_GLM_no_noise_no_RF.RData")

#######################################################################################################################################
#######################################################################################################################################
# Part1: Manipulate raw RMSE of control and the results of both GLM and RF

# High training collinearity + Interaction
data_high_inter_rrmse <- data_transpose(data = high_inter_rrmse/high_inter_rrmse$test_same, algorithm = "GLM", model = "Product", 
                                        training_collinearity = "High", training_correlation = 0.9)
data_high_inter_rrmse_rf <- data_transpose(data = high_inter_rrmse_rf/high_inter_rrmse_rf$test_same, algorithm = "RF", model = "Product", 
                                           training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Quadratic
data_high_quad_rrmse <- data_transpose(data = high_quad_rrmse/high_quad_rrmse$test_same, algorithm = "GLM", model = "Quadratic", 
                                       training_collinearity = "High", training_correlation = 0.9)
data_high_quad_rrmse_rf <- data_transpose(data = high_quad_rrmse_rf/high_quad_rrmse_rf$test_same, algorithm = "RF", model = "Quadratic", 
                                          training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Linear
data_high_linear_rrmse <- data_transpose(data = high_linear_rrmse/high_linear_rrmse$test_same, algorithm = "GLM", model = "Linear", 
                                         training_collinearity = "High", training_correlation = 0.9)
data_high_linear_rrmse_rf <- data_transpose(data = high_linear_rrmse_rf/high_linear_rrmse_rf$test_same, algorithm = "RF", model = "Linear", 
                                            training_collinearity = "High", training_correlation = 0.9)

# Mid training collinearity + Interaction
data_mid_inter_rrmse <- data_transpose(data = mid_inter_rrmse/mid_inter_rrmse$test_same, algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Mid", training_correlation = 0.6)
data_mid_inter_rrmse_rf <- data_transpose(data = mid_inter_rrmse_rf/mid_inter_rrmse_rf$test_same, algorithm = "RF", model = "Product", 
                                          training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Quadratic
data_mid_quad_rrmse <- data_transpose(data = mid_quad_rrmse/mid_quad_rrmse$test_same, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Mid", training_correlation = 0.6)
data_mid_quad_rrmse_rf <- data_transpose(data = mid_quad_rrmse_rf/mid_quad_rrmse_rf$test_same, algorithm = "RF", model = "Quadratic", 
                                         training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Linear
data_mid_linear_rrmse <- data_transpose(data = mid_linear_rrmse/mid_linear_rrmse$test_same, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Mid", training_correlation = 0.6)
data_mid_linear_rrmse_rf <- data_transpose(data = mid_linear_rrmse_rf/mid_linear_rrmse_rf$test_same, algorithm = "RF", model = "Linear", 
                                           training_collinearity = "Mid", training_correlation = 0.6)

# Low training collinearity + Interaction
data_low_inter_rrmse <- data_transpose(data = low_inter_rrmse/low_inter_rrmse$test_same, algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Low", training_correlation = 0.3)
data_low_inter_rrmse_rf <- data_transpose(data = low_inter_rrmse_rf/low_inter_rrmse_rf$test_same, algorithm = "RF", model = "Product", 
                                          training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Quadratic
data_low_quad_rrmse <- data_transpose(data = low_quad_rrmse/low_quad_rrmse$test_same, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Low", training_correlation = 0.3)
data_low_quad_rrmse_rf <- data_transpose(data = low_quad_rrmse_rf/low_quad_rrmse_rf$test_same, algorithm = "RF", model = "Quadratic", 
                                         training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Linear
data_low_linear_rrmse <- data_transpose(data = low_linear_rrmse/low_linear_rrmse$test_same, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Low", training_correlation = 0.3)
data_low_linear_rrmse_rf <- data_transpose(data = low_linear_rrmse_rf/low_linear_rrmse_rf$test_same, algorithm = "RF", model = "Linear", 
                                           training_collinearity = "Low", training_correlation = 0.3)

# Combine all data
data_analysis_poisson <- rbind(data_high_inter_rrmse, data_high_inter_rrmse_rf, data_high_quad_rrmse,
                               data_high_quad_rrmse_rf, data_high_linear_rrmse, data_high_linear_rrmse_rf,
                               data_mid_inter_rrmse, data_mid_inter_rrmse_rf, data_mid_quad_rrmse, 
                               data_mid_quad_rrmse_rf, data_mid_linear_rrmse, data_mid_linear_rrmse_rf,
                               data_low_inter_rrmse, data_low_inter_rrmse_rf, data_low_quad_rrmse, 
                               data_low_quad_rrmse_rf, data_low_linear_rrmse, data_low_linear_rrmse_rf)
                               
rm(data_high_inter_rrmse, data_high_inter_rrmse_rf, data_high_quad_rrmse, 
   data_high_quad_rrmse_rf, data_high_linear_rrmse, data_high_linear_rrmse_rf,
   data_mid_inter_rrmse, data_mid_inter_rrmse_rf, data_mid_quad_rrmse, 
   data_mid_quad_rrmse_rf, data_mid_linear_rrmse, data_mid_linear_rrmse_rf,
   data_low_inter_rrmse, data_low_inter_rrmse_rf, data_low_quad_rrmse, 
   data_low_quad_rrmse_rf, data_low_linear_rrmse, data_low_linear_rrmse_rf)

# Export data as a RData file
save(data_analysis_poisson, file = "data_analysis_poisson_GLM_RF.RData")
#######################################################################################################################################









#######################################################################################################################################
#######################################################################################################################################
# Part2: Calculate relative RMSE of GLMs for making boxplots
# High training collinearity + Interaction
data_high_inter_rrmse <- data_transpose(data = high_inter_rrmse#/high_inter_rrmse$test_same
                                        ,algorithm = "GLM", model = "Product", 
                                        training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Quadratic
data_high_quad_rrmse <- data_transpose(data = high_quad_rrmse#/high_quad_rrmse$test_same
                                       , algorithm = "GLM", model = "Quadratic", 
                                       training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Linear
data_high_linear_rrmse <- data_transpose(data = high_linear_rrmse#/high_linear_rrmse$test_same
                                         , algorithm = "GLM", model = "Linear", 
                                         training_collinearity = "High", training_correlation = 0.9)

# Mid training collinearity + Interaction
data_mid_inter_rrmse <- data_transpose(data = mid_inter_rrmse#/mid_inter_rrmse$test_same
                                       , algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Quadratic
data_mid_quad_rrmse <- data_transpose(data = mid_quad_rrmse#/mid_quad_rrmse$test_same
                                      , algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Linear
data_mid_linear_rrmse <- data_transpose(data = mid_linear_rrmse#/mid_linear_rrmse$test_same
                                        , algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Mid", training_correlation = 0.6)

# Low training collinearity + Interaction
data_low_inter_rrmse <- data_transpose(data = low_inter_rrmse#/low_inter_rrmse$test_same
                                       , algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Quadratic
data_low_quad_rrmse <- data_transpose(data = low_quad_rrmse#/low_quad_rrmse$test_same
                                      , algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Linear
data_low_linear_rrmse <- data_transpose(data = low_linear_rrmse#/low_linear_rrmse$test_same
                                        , algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Low", training_correlation = 0.3)

# Combine all data
data_analysis_poisson_GLM <- rbind(data_high_inter_rrmse, data_high_quad_rrmse, data_high_linear_rrmse,
                                   data_mid_inter_rrmse, data_mid_quad_rrmse, data_mid_linear_rrmse, 
                                   data_low_inter_rrmse, data_low_quad_rrmse, data_low_linear_rrmse)

rm(data_high_inter_rrmse, data_high_quad_rrmse, data_high_linear_rrmse, 
   data_mid_inter_rrmse,  data_mid_quad_rrmse, data_mid_linear_rrmse, 
   data_low_inter_rrmse,  data_low_quad_rrmse, data_low_linear_rrmse)

# Export data as a RData file
save(data_analysis_poisson_GLM, file = "data_analysis_poisson_GLM.RData")
#######################################################################################################################################




######################################################################################################################################
#######################################################################################################################################
# Part3: Manipulate raw RMSE of control and the results of GLM only

# High training collinearity + Interaction
data_high_inter_rrmse <- data_transpose(data = high_inter_rrmse, algorithm = "GLM", model = "Product", 
                                        training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Quadratic
data_high_quad_rrmse <- data_transpose(data = high_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                       training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Linear
data_high_linear_rrmse <- data_transpose(data = high_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                         training_collinearity = "High", training_correlation = 0.9)
# Mid training collinearity + Interaction
data_mid_inter_rrmse <- data_transpose(data = mid_inter_rrmse, algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Quadratic
data_mid_quad_rrmse <- data_transpose(data = mid_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Linear
data_mid_linear_rrmse <- data_transpose(data = mid_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Mid", training_correlation = 0.6)
# Low training collinearity + Interaction
data_low_inter_rrmse <- data_transpose(data = low_inter_rrmse, algorithm = "GLM", model = "Product", 
                                       training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Quadratic
data_low_quad_rrmse <- data_transpose(data = low_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Linear
data_low_linear_rrmse <- data_transpose(data = low_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Low", training_correlation = 0.3)

# Combine all data
data_analysis_poisson <- rbind(data_high_inter_rrmse, data_high_quad_rrmse, data_high_linear_rrmse, 
                               data_mid_inter_rrmse, data_mid_quad_rrmse, data_mid_linear_rrmse,
                               data_low_inter_rrmse, data_low_quad_rrmse, data_low_linear_rrmse)

rm(data_high_inter_rrmse, data_high_quad_rrmse, data_high_linear_rrmse, 
   data_mid_inter_rrmse, data_mid_quad_rrmse, data_mid_linear_rrmse,
   data_low_inter_rrmse, data_low_quad_rrmse, data_low_linear_rrmse)

# Export data as a RData file
save(data_analysis_poisson, file = "data_analysis_poisson_GLM_no_RF.RData")
#######################################################################################################################################




