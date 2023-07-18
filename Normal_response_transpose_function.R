library(reshape2)

data_transpose <- function(data, algorithm = "GLM", 
                           model = "Product", 
                           training_collinearity = "High",
                           training_correlation = 0.9){
  # Transpose all collinearity shift
  nrow_analysis <- nrow(data)
  a <- data[1:nrow_analysis, 3:21]
  
  a$Algorithm <- as.character(algorithm)
  a$Model <- as.character(model)
  a$Predictor_coli <- as.character(training_collinearity)
  # The range shift column represent the ratio of increased range in the test data
  # The range shift for all collinearity shift scenarios is equal to 0
  a$Range_shift <- 0
  a$Collinearity_shift <- 0
  
  a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))
  
  a$Algorithm <- as.factor(a$Algorithm)
  a$Model <- as.factor(a$Model)
  a$Predictor_coli <- as.factor(a$Predictor_coli)
  a$Range_shift <- as.factor(a$Range_shift)
  a$Collinearity_shift <- as.factor(a$Collinearity_shift)
  
  # Rename
  names(a)[7] <- "RMSE"
  
  # Assign values to each level of correlation in the test data
  a$Collinearity_shift <- c(rep(0.9, nrow_analysis),
                            rep(0.8, nrow_analysis),
                            rep(0.7, nrow_analysis),
                            rep(0.6, nrow_analysis),
                            rep(0.5, nrow_analysis),
                            rep(0.4, nrow_analysis),
                            rep(0.3, nrow_analysis),
                            rep(0.2, nrow_analysis),
                            rep(0.1, nrow_analysis),
                            rep(0, nrow_analysis),
                            rep(-0.1, nrow_analysis),
                            rep(-0.2, nrow_analysis),
                            rep(-0.3, nrow_analysis),
                            rep(-0.4, nrow_analysis),
                            rep(-0.5, nrow_analysis),
                            rep(-0.6, nrow_analysis),
                            rep(-0.7, nrow_analysis),
                            rep(-0.8, nrow_analysis),
                            rep(-0.9, nrow_analysis))
  a$Collinearity_shift <- as.factor(a$Collinearity_shift)
  
  # Remove the redundant column created by melt function
  a <- a[, -6]
  
  # Transpose all range shift
  b <- data[1:nrow_analysis, 22:27]
  b$Algorithm <- as.character(algorithm)
  b$Model <- as.character(model)
  b$Predictor_coli <- as.character(training_collinearity)
  b$Range_shift <- 0
  # The collinearity shift column represent the correlation in the test data
  # The correlation for range shift only scenarios is equal to the training correlation
  b$Collinearity_shift <- training_correlation # This is the training collinearity
  
  b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))
  
  b$Algorithm <- as.factor(b$Algorithm)
  b$Model <- as.factor(b$Model)
  b$Predictor_coli <- as.factor(b$Predictor_coli)
  b$Range_shift <- as.factor(b$Range_shift)
  b$Collinearity_shift <- as.factor(b$Collinearity_shift)
  
  # Rename
  names(b)[7] <- "RMSE"
  
  # Assign values to each level of increased range in the test data
  b$Range_shift <- c(rep(0.2, nrow_analysis),
                     rep(0.4, nrow_analysis),
                     rep(0.6, nrow_analysis),
                     rep(2, nrow_analysis),
                     rep(4, nrow_analysis),
                     rep(6, nrow_analysis))
  b$Range_shift <- as.factor(b$Range_shift)
  
  # Remove the redundant column created by melt function
  b <- b[, -6]
  
  # Transpose all interactions
  c <- data[1:nrow_analysis, 28:141]
  
  c$Algorithm <- as.character(algorithm)
  c$Model <- as.character(model)
  c$Predictor_coli <- as.character(training_collinearity)
  c$Range_shift <- 0
  c$Collinearity_shift <- 0
  
  c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
  c$Algorithm <- as.factor(c$Algorithm)
  c$Model <- as.factor(c$Model)
  c$Predictor_coli <- as.factor(c$Predictor_coli)
  c$Range_shift <- as.factor(c$Range_shift)
  c$Collinearity_shift <- as.factor(c$Collinearity_shift)
  
  # Rename
  names(c)[7] <- c("RMSE")
  
  # Assign values to each level of increased range in the test data
  c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                     rep(0.4, nrow_analysis*19),
                     rep(0.6, nrow_analysis*19),
                     rep(2, nrow_analysis*19),
                     rep(4, nrow_analysis*19),
                     rep(6, nrow_analysis*19))
  
  # Assign values to each level of correlation in the test data
  c$Collinearity_shift <- rep(c(rep(0.9, nrow_analysis),
                                rep(0.8, nrow_analysis),
                                rep(0.7, nrow_analysis),
                                rep(0.6, nrow_analysis),
                                rep(0.5, nrow_analysis),
                                rep(0.4, nrow_analysis),
                                rep(0.3, nrow_analysis),
                                rep(0.2, nrow_analysis),
                                rep(0.1, nrow_analysis),
                                rep(0, nrow_analysis),
                                rep(-0.1, nrow_analysis),
                                rep(-0.2, nrow_analysis),
                                rep(-0.3, nrow_analysis),
                                rep(-0.4, nrow_analysis),
                                rep(-0.5, nrow_analysis),
                                rep(-0.6, nrow_analysis),
                                rep(-0.7, nrow_analysis),
                                rep(-0.8, nrow_analysis),
                                rep(-0.9, nrow_analysis)), 6)
  c$Range_shift <- as.factor(c$Range_shift)
  c$Collinearity_shift <- as.factor(c$Collinearity_shift)
  
  # Remove the redundant column created by melt function
  c <- c[, -6]
  
  rbind(a, b, c)
}
