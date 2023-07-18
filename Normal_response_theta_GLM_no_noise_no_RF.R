# This script was used to simulate the response variable that followed normal distributions

# Data generator
library(faux)
library(mvtnorm)
library(data.table)
library(ggplot2)
library(data.table)

library(foreach)
library(doParallel)


# Parallel stochastic 
library(doRNG)

# Data generation function
######################################################################################################################################################################################################################
# Date: 03/17/2023
# Change: create a new column called theta to save the parameter of distribution
######################################################################################################################################################################################################################

collineariser_extra <- function(N = 1000, 
                                varnames = c("X1", "X2", "X3"),
                                mu.training = c(0, 0, 0),
                                sd.training = c(1, 1, 1),
                                rho.training = c(0.5, 0.3, 0.2),
                                empirical = FALSE,
                                rho.test = seq(0.9, -0.9, by = -0.1),
                                range.vec = c(.2, .4, .6, 2, 4, 6),
                                method = "normal",
                                seed = 1001,
                                FUN = function(X) 2*X[,1] - 1.5*X[,2] - 2*X[,3], ...){
  
  library(mvtnorm)
  library(faux)
  
  # helper functions: --------------------------------------------------------
  Xgen <- function(N = N, mu = mu.training, sd = sd.training, r = rho.training,
                   varnames = varnames, empirical = empirical, ...){
    
    X_multinorm <- rnorm_multi(n = N, 
                               mu = mu.training,
                               sd = sd.training,
                               r = rho.training, 
                               varnames = varnames,
                               empirical = empirical)
    X_multinorm
  }
  
  # Create a list object for each test scenario
  # Date: 07/21/2022
  # Change: Added r = 0 between r =  0.1 and -0.1 to make a full gradient of collinearity shift
  # 1-2: test and same
  # 3-11: a gradient of collinearity with 9 levels (positive collinearity)
  # 12: r = 0
  # 12+1(13)-20+1(21): a gradient of collinearity with 9 levels (negative collinearity)
  # 21+1(22)-26+1(27): a gradient of range shift with 6 levels
  
  # 27+1(28)-44+2(46): a gradient of interaction between collinearity and range shift, 2*range x 18 collinearity levels
  # 45+2(47)-62+3(65): a gradient of interaction between collinearity and range shift, 3*range x 18 collinearity levels
  # 63+3(66)-80+4(84): a gradient of interaction between collinearity and range shift, 4*range x 18 collinearity levels
  # 81+4(85)-98+5(103): a gradient of interaction between collinearity and range shift, 5*range x 18 collinearity levels
  # 99+5(104)-116+6(122): a gradient of interaction between collinearity and range shift, 6*range x 18 collinearity 
  # 117+6(123)-134+7(141): a gradient of interaction between collinearity and range shift, 7*range x 18 collinearity levels
  
  # Create a list of training and all test data sets
  Data <- vector("list", 141)
  
  # Name each scenario
  # Name training and test same data sets
  names(Data)[1:2] <- c("train", "test_same")
  
  # Name collinearity shift data
  for (i in 1:19) {names(Data)[i+2] = paste("r", rho.test[i], sep = "=")}
  
  # for (i in 1:9) {names(Data)[i+2] = paste("r", rho.test[i], sep = "=")}
  # for (i in 1:9) {names(Data)[i+11] = paste("r", rev(rho.test)[i], sep = "=-")}
  
  # Name range shift data sets
  names(Data)[22:27] <- c("1Range", "2Range", "3Range", "4Range", "5Range", "6Range")
  
  # Name range and collinearity shift data
  for (j in 1:6) {
    for (i in 1:19) {
      names(Data)[i+27+19*(j-1)] = paste(names(Data)[22:27][j], names(Data)[i+2], sep = "*")
    }
  }
  
  # generate train/test.same data:
  if (!is.null(seed)) set.seed(seed)
  X <- Xgen(N = 2*N, mu = mu.training, sd = sd.training, r = rho.training,
            varnames = varnames, empirical = empirical)
  
  colnames(X) <- paste("X", 1:3, sep="")
  Data[[1]] <- X[1:N,]
  # Create quadratic and interaction terms to calculate multicollinearity
  # Create (X1)^2
  Data[[1]][, 4] <- (Data[[1]][, 1])^2
  # Create (X2)^2
  Data[[1]][, 5] <- (Data[[1]][, 2])^2
  # Create (X3)^2
  Data[[1]][, 6] <- (Data[[1]][, 3])^2
  # Create X1*X2
  Data[[1]][, 7] <- (Data[[1]][, 1])*(Data[[1]][, 2])
  colnames(Data[[1]])[4:6] <- paste("X", 1:3,"^2", sep = "")
  colnames(Data[[1]])[7] <- paste("X1","X2", sep = "*")
  
  Data[[2]] <-  X[(N+1):(2*N),]
  # Create quadratic and interaction terms to calculate multicollinearity
  # Create (X1)^2
  Data[[2]][, 4] <- (Data[[2]][, 1])^2
  # Create (X2)^2
  Data[[2]][, 5] <- (Data[[2]][, 2])^2
  # Create (X3)^2
  Data[[2]][, 6] <- (Data[[2]][, 3])^2
  # Create X1*X2
  Data[[2]][, 7] <- (Data[[2]][, 1])*(Data[[2]][, 2])
  colnames(Data[[2]])[4:6] <- paste("X", 1:3,"^2", sep = "")
  colnames(Data[[2]])[7] <- paste("X1","X2", sep = "*")
  
  # generate test data with collinearity shift only
  # X1 and X2 are positively correlated starting from 0.9 to 0.1
  
  # Date: 07/21/2022
  # Change: X1 and X2 are correlated starting from 0.9 to -0.9
  for (i in 1:19) {
    Data[[i+2]] <- Data[[1]]
    Data[[i+2]][, 2] <- rnorm_pre(Data[[1]][, 1], mu = 0, sd = 1, r = rho.test[i], empirical = F)
    colnames(Data[[i+2]])[1:3] <- paste("X", 1:3, sep="")
    
    # Create quadratic and interaction terms to calculate multicollinearity
    # Create (X1)^2
    Data[[i+2]][, 4] <- (Data[[i+2]][, 1])^2
    # Create (X2)^2
    Data[[i+2]][, 5] <- (Data[[i+2]][, 2])^2
    # Create (X3)^2
    Data[[i+2]][, 6] <- (Data[[i+2]][, 3])^2
    # Create X1*X2
    Data[[i+2]][, 7] <- (Data[[i+2]][, 1])*(Data[[i+2]][, 2])
    colnames(Data[[i+2]])[4:6] <- paste("X", 1:3,"^2", sep = "")
    colnames(Data[[i+2]])[7] <- paste("X1","X2", sep = "*")
  }
  
  # # X1 and X2 are negatively correlated starting from -0.1 to -0.9
  # for (i in 1:9) {
  #   Data[[i+11]] <- Data[[1]]
  #   # This formula change the sign of correlation but does not change the range
  #   Data[[i+11]][, 2] <- rnorm_pre(Data[[1]][, 1], mu = 0, sd = 1, r = -rev(rho.test)[i], empirical = F)
  #   colnames(Data[[i+11]]) <- paste("X", 1:3, sep="")
  #   
  #   # Create quadratic and interaction terms to calculate multicollinearity
  #   # Create (X1)^2
  #   Data[[i+11]][, 4] <- (Data[[i+11]][, 1])^2
  #   # Create (X2)^2
  #   Data[[i+11]][, 5] <- (Data[[i+11]][, 2])^2
  #   # Create (X3)^2
  #   Data[[i+11]][, 6] <- (Data[[i+11]][, 3])^2
  #   # Create X1*X2
  #   Data[[i+11]][, 7] <- (Data[[i+11]][, 1])*(Data[[i+11]][, 2])
  #   colnames(Data[[i+11]])[4:6] <- paste("X", 1:3,"^2", sep = "")
  #   colnames(Data[[i+11]])[7] <- paste("X1","X2", sep = "*")
  #   
  # }
  
  # generate test data with range shift only
  # 01/24/2022 Fixed range shift for X1,X2,X3. Originally all range shift were based on max(X1) - min(X1)
  for (i in 1:6) {
    Data[[i+21]] <- Data[[1]]
    # Data[[i+21]][, 1] <- range.vec[i]*(max(Data[[1]][, 1])-min(Data[[1]][, 1])) + Data[[i+21]][, 1]
    # 08/08/2022
    # Only change the novelty of X2
    Data[[i+21]][, 2] <- range.vec[i]*(max(Data[[1]][, 2])-min(Data[[1]][, 2])) + Data[[i+21]][, 2]
    # Data[[i+21]][, 2] <- range.vec[i]*(max(Data[[1]][, 2])-min(Data[[1]][, 2])) + Data[[i+21]][, 2]
    # Data[[i+21]][, 3] <- range.vec[i]*(max(Data[[1]][, 3])-min(Data[[1]][, 3])) + Data[[i+21]][, 3]
    colnames(Data[[i+21]]) <- paste("X", 1:3, sep="")
    
    # Create quadratic and interaction terms to calculate multicollinearity
    # Create (X1)^2
    Data[[i+21]][, 4] <- (Data[[i+21]][, 1])^2
    # Create (X2)^2
    Data[[i+21]][, 5] <- (Data[[i+21]][, 2])^2
    # Create (X3)^2
    Data[[i+21]][, 6] <- (Data[[i+21]][, 3])^2
    # Create X1*X2
    Data[[i+21]][, 7] <- (Data[[i+21]][, 1])*(Data[[i+21]][, 2])
    colnames(Data[[i+21]])[4:6] <- paste("X", 1:3,"^2", sep = "")
    colnames(Data[[i+21]])[7] <- paste("X1","X2", sep = "*")
    
  }
  
  # generate test data with collinearity and range shift
  # 01/24/2022 Fixed interaction between range and col shift for X1,X2,X3. 
  # Originally all range shift were based on collinearity shifted Data[i](max(X1) - min(X1)) [i = 3-10]
  for (j in 1:6) {
    for (i in 1:19) {
      Data[[i+27+19*(j-1)]] <- Data[[i+2]]
      # 08/08/2022
      # Only change the novelty of X2
      # Data[[i+27+19*(j-1)]][, 1] <- range.vec[j]*(max(Data[[1]][, 1])-min(Data[[1]][, 1])) + Data[[i+27+19*(j-1)]][, 1]
      Data[[i+27+19*(j-1)]][, 2] <- range.vec[j]*(max(Data[[1]][, 2])-min(Data[[1]][, 2])) + Data[[i+27+19*(j-1)]][, 2]
      # Data[[i+27+19*(j-1)]][, 3] <- range.vec[j]*(max(Data[[1]][, 3])-min(Data[[1]][, 3])) + Data[[i+27+19*(j-1)]][, 3]
      
      colnames(Data[[i+27+19*(j-1)]]) <- paste("X", 1:3, sep="")
      
      # Create quadratic and interaction terms to calculate multicollinearity
      # Create (X1)^2
      Data[[i+27+19*(j-1)]][, 4] <- (Data[[i+27+19*(j-1)]][, 1])^2
      # Create (X2)^2
      Data[[i+27+19*(j-1)]][, 5] <- (Data[[i+27+19*(j-1)]][, 2])^2
      # Create (X3)^2
      Data[[i+27+19*(j-1)]][, 6] <- (Data[[i+27+19*(j-1)]][, 3])^2
      # Create X1*X2
      Data[[i+27+19*(j-1)]][, 7] <- (Data[[i+27+19*(j-1)]][, 1])*(Data[[i+27+19*(j-1)]][, 2])
      colnames(Data[[i+27+19*(j-1)]])[4:6] <- paste("X", 1:3,"^2", sep = "")
      colnames(Data[[i+27+19*(j-1)]])[7] <- paste("X1","X2", sep = "*")
    }
  }
  
  # make y for each Data set:
  
  
  # X <- Xmaker(N=2*N, ...)
  #curve(2*x-1.5x*x, from=-2, to=2)
  f <- FUN
  if (method == "normal"){
    # Calculate Y and theta for Normal response
    Thetas <- lapply(Data, function(X) f(X))
    
    #  Relative SD
    Ys <- lapply(Data, function(X) rnorm(n = nrow(X), mean = f(X), sd = sqrt(var(Thetas[[1]])*0.3)))
    # Absolute SD
    #Ys <- lapply(Data.scaled, function(X) rnorm(n = nrow(X), mean = f(X), sd = 10))
    
  } else if (method == "poisson") {
    # Calculate Y and theta for Poisson response
    Ys <- lapply(Data, function(X) rpois(n = nrow(X), lambda = exp(f(X))))
    
    # Theta was calculated on link level
    # Calculate log(lambda)
    Thetas <- lapply(Data, function(X) f(X))
    
  } else if (method == "binomial") {
    # Calculate Y and theta for Binomial response
    Ys <- lapply(Data, function(X) rbinom(n = nrow(X), size = 1, prob = plogis(f(X))))
    
    # Calculate logit(p)
    Thetas <- lapply(Data, function(X) f(X))
    
  }
  
  Data.final <- Data
  attr(Data.final, "function") <- attr(f, "source")
  
  # Some options of noise tested
  # noise <- 1/rgamma(length(Ys[[1]]), 10, 0.01)
  # noise <- rnorm(length(Ys[[1]]), mean = 0, sd = 0.15*mean(Ys[[1]]))
  # noise <- runif(length(Ys[[1]]), -10, 10)
  
  # Added the noise as 0.05 to 0.35 of response variable
  # noise <- round(runif(length(Ys[[1]]), min = abs(0.05*Ys[[1]]), max = abs(0.35*Ys[[1]])))
  noise <- 0
  
  # Only add the noise in the training data would not work because the absolute RMSE between training and test
  # would be this noise if there is no any other effects
  Data.final[[1]] <- cbind(y = Ys[[1]] + noise, theta = Thetas[[1]], Data[[1]])
  
  for (i in 1:(length(Data)-1))
    # No noise in the test data
    Data.final[[i+1]] <- cbind(y = Ys[[i+1]], theta = Thetas[[i+1]], Data[[i+1]])
  
  Data.final
  
}



# Simulation function
# RMSE was calculated for true theta not true Y
# For each simulator function, mean.mu, mean.lamda, and mean.detection/occurrence should be changed based on the response
# For each simulator function, the link function in GLM should be changed based on the response

linear_model_1 <- function(vars = c("X1", "X2", "X3"),
                           # f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2]) - 3*X[,3], 
                           mu = c(0, 0, 0),
                           sd = c(1, 1, 1),
                           rho = c(0.5, 0.3, 0.2),
                           response = "y",
                           predictors = c("X1", "X2", "I(X1^2)", "X1*X2", "X3"),
                           predictors_index = c("X1", "X2", "X1^2", "X1*X2", "X3")){
  
  s <- sample(1:10000, 1)
  
  # Make all the coefficients as random
  set.seed(s)
  beta0 <- runif(n = 1, min = -10, max = 10)
  beta1 <- runif(n = 1, min = -10, max = 10)
  beta2 <- runif(n = 1, min = -10, max = 10)
  beta3 <- runif(n = 1, min = -10, max = 10)
  beta4 <- runif(n = 1, min = -10, max = 10)
  beta5 <- runif(n = 1, min = -10, max = 10)
  
  # Define functional relationship based on the number of predictors
  if (length(predictors) == 5) {
    # Linear, quadratic, and product
    f <- function(X) beta1*X[,1] + beta2*X[,1]^2 + beta3*X[,2] + beta4*(X[,1])*(X[,2]) + beta5*X[,3] + beta0
  } else if (length(predictors) == 4) {
    # Linear and quadratic
    f <- function(X) beta1*X[,1] + beta2*X[,1]^2 + beta3*X[,2] + beta5*X[,3] + beta0
  } else {
    # Linear
    f <- function(X) beta1*X[,1] + beta3*X[,2] + beta5*X[,3] + beta0
  }
  
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        method = "normal",
                                        # varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = seq(0.9, -0.9, by = -0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        N = 1000)
  length(simulated_data$train$y / simulated_data$train$theta > 0)
  
  summary(simulated_data$train$y)
  summary(simulated_data$`6Range`$y)
  summary(simulated_data$train$theta)
  sd(simulated_data$train$y)
  sd(simulated_data$train$theta)
  
  summary(simulated_data$`r=0.9`$y)
  summary(simulated_data$`r=-0.9`$y)
  sd(simulated_data$`r=0.9`$y)
  sd(simulated_data$`r=-0.9`$y)
  
  sum(simulated_data$train$y > 0)
  summary(f(simulated_data$train[, 3:7]))
  summary(f(simulated_data$`6Range`[, 3:7]))
  
  # Fit a linear model for the simulated data
  formula <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  fitted <- glm(data = as.data.frame(simulated_data$train), formula, family = "gaussian")
  
  # # Fit a RF model for the simulated data
  # Based on Random Forest Prediction Intervals. Haozhe Zhang, Joshua Zimmerman, Dan Nettleton, and Dan Nordman.
  # The American Statistician, 74(4):392-406, 2020, specifying the relationship between Y and X may increase the
  # model extrapolation of Random Forests
  
  # 08/09/2022
  # Although adding squared terms and interaction term does not considerably help model fitting,
  # keeping the same model structure as GLM makes some sense.
  # formula_rf <- as.formula(
  #   paste(response,
  #         paste(predictors, collapse = "+"),
  #         sep = "~")
  # )
  
  # # Keep the same formula used in developing GLMs for RF
  # Use alternative interface to data because the illegal columns issue
  # fitted_rf <- ranger::ranger(y = simulated_data$train[, "y"],
  #                             x = simulated_data$train[, c(predictors_index)],
  #                             data = simulated_data$train[, c("y", predictors_index)])

  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                   #
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                   #
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                   #
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  # RMSE
  # Date: 03/09/202
  # Note: RMSE is not scale-dependent and is related to the scale of original measurement. RSE might be also considered since it is not scale-dependent
  # Change: Change RMSE to RRMSE
  # Note: sum((predict(fitted, as.data.frame(X), type = "response")-as.data.frame(X)$y)^2)/nrow(X) <- mean((predict(fitted, as.data.frame(X), type = "response")-as.data.frame(X)$y)^2)
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                       #
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                       #
  #                                                                     Make sure type = response for poisson and logistic regression!!!!!!!!!!!!                                                                       #
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  ###################################################################################################################################################################################################################
  
  rmse <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted, as.data.frame(X), type = "link")-as.data.frame(X)$theta)^2))))
  
  # RMSE1 = sqrt((sum((Ypred - Yi)/Yi)^2)/n)
  rrmse_1 <- lapply(simulated_data, function(X) (sqrt(mean(((predict(fitted, as.data.frame(X), type = "link")-as.data.frame(X)$theta)/as.data.frame(X)$theta)^2))))
  
  # RMSE2 = sqrt((sum(Ypred - Yi)^2)/n)/sd(Yi)
  rrmse_2 <- lapply(simulated_data, function(X) ((sqrt(mean((predict(fitted, as.data.frame(X), type = "link")-as.data.frame(X)$theta)^2)))/(sd(as.data.frame(X)$theta))))
  
  # /sum((predict(fitted, as.data.frame(X), type = "response"))^2)
  # Note: ranger package returns predictions in predict(ranger.object, newdata)$predictions instead of predict(RandomForest.object, newdata, type = "response)
  # rmse_rf <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted_rf, as.data.frame(X))$predictions-as.data.frame(X)$theta)^2))))
  
  # RMSE1 = sqrt((sum((Ypred - Yi)/Yi)^2)/n)
  # rrmse_rf_1 <- lapply(simulated_data, function(X) (sqrt(mean(((predict(fitted_rf, as.data.frame(X))$predictions-as.data.frame(X)$theta)/as.data.frame(X)$theta)^2))))
  
  # RMSE2 = sqrt((sum(Ypred - Yi)^2)/n)/sd(Yi)
  # rrmse_rf_2 <- lapply(simulated_data, function(X) ((sqrt(mean((predict(fitted_rf, as.data.frame(X))$predictions-as.data.frame(X)$theta)^2)))/(sd(as.data.frame(X)$theta))))
  
  # # Variance
  # variance <- lapply(simulated_data, function(X) (mean((mean(predict(fitted, as.data.frame(X))) - predict(fitted, as.data.frame(X)))^2)))
  
  # Calculate the condition number for all predictors
  cn_all <-  lapply(simulated_data, function(X)(kappa(cor(X[, predictors_index]), exact = TRUE)))
  
  # Retrieve R-square and standard errors of coefficients for X2 that induced collinearity shift
  # For lm object
  # model_r_square <- summary(fitted)$r.square
  # For glm object
  model_r_square <- with(summary(fitted), 1 - deviance/null.deviance)
  std_error_X2 <- summary(fitted)$coefficients[2, 2]
  
  #rf_r_square <- fitted_rf$rsq[length(fitted_rf$rsq)]
  
  # Calculate Mahalanobis distance
  center_train <- colMeans(simulated_data$train[, predictors_index])
  cov_train <- cov(simulated_data$train[, predictors_index])
  
  # calculate Mahalanobis distance using the covariance in the training data
  mal_dist_mean <- lapply(simulated_data, function(X)(mean(mahalanobis(x = X[, predictors_index], 
                                                                       center = center_train, 
                                                                       cov = cov_train))))
  
  mal_dist_median <- lapply(simulated_data, function(X)(median(mahalanobis(x = X[, predictors_index],
                                                                           center = center_train, 
                                                                           cov = cov_train))))
  
  # result <- list(RMSE = NULL,
  #                CN_all = NULL,
  #                R_Square = NULL,
  #                Std_error = NULL,
  #                M_distance_mean = NULL,
  #                M_distance_median = NULL)
  # result$RRMSE <- rrmse
  # result$CN_all <- cn_all
  # result$R_Square <- model_r_square
  # result$Std_error <- std_error_X2
  # result$M_distance_mean <- mal_dist_mean
  # result$M_distance_median <- mal_dist_median
  # return(result)
  
  return(list(#"RRMSE" = rrmse, #"RMSE_RF" = rmse_rf,
    "RMSE" = rmse,
    "RRMSE1" = rrmse_1,
    "RRMSE2" = rrmse_2,
    #"RMSERF" = rmse_rf,
    #"RRMSE1_RF" = rrmse_rf_1,
    #"RRMSE2_RF" = rrmse_rf_2,
    "CN_all" = cn_all,
    "R_square" = model_r_square, 
    #"R_square_rf" = rf_r_square,
    "Std_error" = std_error_X2,
    "M_distance_mean" = mal_dist_mean,
    "M_distance_median" = mal_dist_median))
}

# Run simulation for high training collinearity using foreach
parallel::detectCores()
n.cores <- 14

# Create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Reproducibility for parallel processing
registerDoRNG(1234)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

start_time <- Sys.time()
# set.seed(1001) no use set.seed only works for sequential processing
high_train_coli_interaction <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.9, 0.8, 0.7),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X1*X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X1*X2", "X3"))
}

end_time <- Sys.time()

# Stop the cluster
parallel::stopCluster(cl = my.cluster)

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
high_inter_rrmse <- setDF(rbindlist(high_train_coli_interaction[names(high_train_coli_interaction) == "RMSE"]))
high_inter_rrmse

high_inter_m_dist <- setDF(rbindlist(high_train_coli_interaction[names(high_train_coli_interaction) == "M_distance_mean"]))
high_inter_m_dist

# Calculate mean RMSE over all simulations
high_inter_rrmse_mean <- colMeans(high_inter_rrmse)
high_inter_rrmse_median <- apply(high_inter_rrmse, 2, median)

# Calculate mean M distance over all simulations
high_inter_m_dist_mean <- colMeans(high_inter_m_dist)
high_inter_m_dist_median <- apply(high_inter_m_dist, 2, median)

# Visualize novelty = 0
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_mean)[3:21], 
     xlab = "Correlation between X1 and X2", ylab = "RRMSE", main = "Novelty = 0")
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_median)[3:21], 
     xlab = "Correlation between X1 and X2", ylab = "RRMSE", main = "Novelty = 0")

# Check R-squares
summary(unlist(high_train_coli_inter[names(high_train_coli_inter) =="R_square"]))


# Versus M distance novelty = 0
plot((high_inter_m_dist_mean)[3:21], (high_inter_rrmse_mean)[3:21], 
     xlab = "Mahalanobis distance", ylab = "rrmse", main = "Novelty = 0")
plot((high_inter_m_dist_median)[3:21], (high_inter_rrmse_median)[3:21], 
     xlab = "Mahalanobis distance", ylab = "rrmse", main = "Novelty = 0")

# Visualize novelty = 0.2
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_mean)[28:46],#/(high_inter_rrmse_mean)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 0.2")
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_median)[28:46],#/(high_inter_rrmse_median)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 0.2")

# Versus M distance novelty = 0.2
plot((high_inter_m_dist_mean)[28:46], (high_inter_rrmse_mean)[28:46],#/(high_inter_rrmse_mean)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 0.2")

# Visualize novelty = 4
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_mean)[104:122],#/(high_inter_rrmse_mean)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 4")
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rrmse_median)[104:122]/(high_inter_rrmse_median)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 4")

# Versus M distance novelty = 4
plot((high_inter_m_dist_mean)[104:122], (high_inter_rrmse_mean)[104:122]/(high_inter_rrmse_mean)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 4")
plot((high_inter_m_dist_median)[104:122], (high_inter_rrmse_median)[104:122]/(high_inter_rrmse_median)[2], 
     xlab = "Correlation between X1 and X2", ylab = "rrmse", main = "Novelty = 4")


# Visualize novelty only
plot(c(0, 0.2, 0.4, 0.6, 2, 4, 6), high_inter_rrmse_mean[c(2, 22:27)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "rrmse", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

plot(c(0, 0.2, 0.4, 0.6, 2, 4), high_inter_rrmse_mean[c(2, 22:26)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "rrmse", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

plot(c(0, 0.2, 0.4, 0.6, 2), high_inter_rrmse_mean[c(2, 22:25)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "rrmse", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

plot(c(0, 0.2, 0.4, 0.6, 2), high_inter_rrmse_median[c(2, 22:25)]/(high_inter_rrmse_median)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "RMSE", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

plot(c(0.2, 0.4, 0.6), high_inter_rrmse_mean[c(22:24)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "RRMSE", main = "Collinearity shift = 0")

plot(c(2, 4), high_inter_rrmse_mean[c(25:26)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "rrmse", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

# Versus M distance
plot(high_inter_m_dist_m)

# Versus M distance
plot(high_inter_m_dist_mean[c(3, 22:27)], high_inter_rrmse_mean[c(3, 22:27)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "RMSE", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)

plot(high_inter_m_dist_mean[c(3, 22:25)], high_inter_rrmse_mean[c(3, 22:25)]/(high_inter_rrmse_mean)[2], col = "red", 
     xlab = "Increased range of X2", ylab = "RMSE", main = "Collinearity shift = 0")
abline(h = 2, v = 0.4)


###############################################################################################################################################
# Date:06/23/2023
# Change: 1. Run simulation for the gradients of training collinearity and model complexity
#         2. Combine all responses variables to one data generator function
#         3. Fit Random Forests using ranger package in simulator function
###############################################################################################################################################

# Run simulation for high training collinearity using foreach
parallel::detectCores()
n.cores <- 14

# Create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Reproducibility for parallel processing
registerDoRNG(1234)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

start_time <- Sys.time()
# set.seed(1001) set.seed(1001) no use set.seed only works for sequential processing
# high training collinearity + Interaction
high_train_coli_inter <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.9, 0.8, 0.7),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X1*X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X1*X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
high_inter_rrmse <- setDF(rbindlist(high_train_coli_inter[names(high_train_coli_inter) == "RMSE"]))
high_inter_rrmse

high_inter_m_dist <- setDF(rbindlist(high_train_coli_inter[names(high_train_coli_inter) == "M_distance_mean"]))
high_inter_m_dist

# high training collinearity + Quadratic
high_train_coli_quad <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.9, 0.8, 0.7),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
high_quad_rrmse <- setDF(rbindlist(high_train_coli_quad[names(high_train_coli_quad) == "RMSE"]))
high_quad_rrmse

high_quad_m_dist <- setDF(rbindlist(high_train_coli_quad[names(high_train_coli_quad) == "M_distance_mean"]))
high_quad_m_dist

# high training collinearity + linear
high_train_coli_linear <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.9, 0.8, 0.7),
                 response = "y",
                 predictors = c("X1", "X2", "X3"),
                 predictors_index = c("X1", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
high_linear_rrmse <- setDF(rbindlist(high_train_coli_linear[names(high_train_coli_linear) == "RMSE"]))
high_linear_rrmse

high_linear_m_dist <- setDF(rbindlist(high_train_coli_linear[names(high_train_coli_linear) == "M_distance_mean"]))
high_linear_m_dist

# mid training collinearity + Interaction
mid_train_coli_inter <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.6, 0.5, 0.4),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X1*X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X1*X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
mid_inter_rrmse <- setDF(rbindlist(mid_train_coli_inter[names(mid_train_coli_inter) == "RMSE"]))
mid_inter_rrmse

mid_inter_m_dist <- setDF(rbindlist(mid_train_coli_inter[names(mid_train_coli_inter) == "M_distance_mean"]))
mid_inter_m_dist

# mid training collinearity + Quadratic
mid_train_coli_quad <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.6, 0.5, 0.4),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
mid_quad_rrmse <- setDF(rbindlist(mid_train_coli_quad[names(mid_train_coli_quad) == "RMSE"]))
mid_quad_rrmse

mid_quad_m_dist <- setDF(rbindlist(mid_train_coli_quad[names(mid_train_coli_quad) == "M_distance_mean"]))
mid_quad_m_dist

# mid training collinearity + linear
mid_train_coli_linear <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.6, 0.5, 0.4),
                 response = "y",
                 predictors = c("X1", "X2", "X3"),
                 predictors_index = c("X1", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
mid_linear_rrmse <- setDF(rbindlist(mid_train_coli_linear[names(mid_train_coli_linear) == "RMSE"]))
mid_linear_rrmse

mid_linear_m_dist <- setDF(rbindlist(mid_train_coli_linear[names(mid_train_coli_linear) == "M_distance_mean"]))
mid_linear_m_dist

# low training collinearity + Interaction
low_train_coli_inter <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.3, 0.2, 0.1),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X1*X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X1*X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
low_inter_rrmse <- setDF(rbindlist(low_train_coli_inter[names(low_train_coli_inter) == "RMSE"]))
low_inter_rrmse

low_inter_m_dist <- setDF(rbindlist(low_train_coli_inter[names(low_train_coli_inter) == "M_distance_mean"]))
low_inter_m_dist

# low training collinearity + Quadratic
low_train_coli_quad <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.3, 0.2, 0.1),
                 response = "y",
                 predictors = c("X1", "I(X1^2)", "X2", "X3"),
                 predictors_index = c("X1", "X1^2", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
low_quad_rrmse <- setDF(rbindlist(low_train_coli_quad[names(low_train_coli_quad) == "RMSE"]))
low_quad_rrmse

low_quad_m_dist <- setDF(rbindlist(low_train_coli_quad[names(low_train_coli_quad) == "M_distance_mean"]))
low_quad_m_dist

# low training collinearity + linear
low_train_coli_linear <- foreach(icount(1000), .combine = 'c') %dopar% {
  linear_model_1(vars = c("X1", "X2", "X3"),
                 mu = c(0, 0, 0),
                 sd = c(1, 1, 1),
                 rho = c(0.3, 0.2, 0.1),
                 response = "y",
                 predictors = c("X1", "X2", "X3"),
                 predictors_index = c("X1", "X2", "X3"))
}

# Based on the structure of what foreach returns, subset the lists using list names and rbind all the lists
low_linear_rrmse <- setDF(rbindlist(low_train_coli_linear[names(low_train_coli_linear) == "RMSE"]))
low_linear_rrmse

low_linear_m_dist <- setDF(rbindlist(low_train_coli_linear[names(low_train_coli_linear) == "M_distance_mean"]))
low_linear_m_dist

end_time <- Sys.time()

# Stop the cluster
parallel::stopCluster(cl = my.cluster)

# Save the output
save.image("C:/FSU/Collinearity_extra/Revision/Normal_response/Normal_response_theta_GLM_no_noise_no_RF.RData")

# Check some results real quick
# Visualize novelty = 0
plot(seq(0.9, -0.9, by = -0.1), (apply(high_inter_rrmse, 2, median))[3:21]/(apply(high_inter_rrmse, 2, median))[2], 
     xlab = "Correlation between X1 and X2", ylab = "RMSE", main = "Novelty = 0")

####################################################################################################################################
# Date: 04/17/2023
# Change: Calculate the minimum of the normal response to make sure the mu's are all positive after adding the beta0

# Product functional relationship
normal_minimum_product  <- function(){
  s <- sample(1:10000, 1)
  
  # Make all the coefficients as random
  set.seed(s)
  beta0 <- runif(n = 1, min = 10, max = 50)
  beta1 <- runif(n = 1, min = 0.001, max = 0.1)
  beta2 <- runif(n = 1, min = 0.001, max = 0.1)
  beta3 <- runif(n = 1, min = 0.001, max = 0.1)
  beta4 <- runif(n = 1, min = 0.0001, max = 0.001)
  beta5 <- runif(n = 1, min = 0.001, max = 0.1)
  
  f <- function(X) beta1*X[,1] + beta2*X[,1]^2 - beta3*X[,2] + beta4*(X[,1])*(X[,2]) + beta5*X[,3]
  
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        method = "normal",
                                        # varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = seq(0.9, -0.9, by = -0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5,
                                        N = 1000)
  min_train <- min(simulated_data$train$y)
  min_novelty <- min(simulated_data$`6Range`$y)
  min <- min(min_train, min_novelty)
  min
}

# Quadratic functional relationship
normal_minimum_quad <- function(){
  s <- sample(1:10000, 1)
  
  # Make all the coefficients as random
  set.seed(s)
  beta0 <- runif(n = 1, min = 10, max = 50)
  beta1 <- runif(n = 1, min = 0.001, max = 0.1)
  beta2 <- runif(n = 1, min = 0.001, max = 0.1)
  beta3 <- runif(n = 1, min = 0.001, max = 0.1)
  beta4 <- runif(n = 1, min = 0.0001, max = 0.001)
  beta5 <- runif(n = 1, min = 0.001, max = 0.1)
  
  f <- function(X) beta1*X[,1] + beta2*X[,1]^2 - beta3*X[,2] + beta5*X[,3]
  
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        method = "normal",
                                        # varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = seq(0.9, -0.9, by = -0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5,
                                        N = 1000)
  min_train <- min(simulated_data$train$y)
  min_novelty <- min(simulated_data$`6Range`$y)
  min <- min(min_train, min_novelty)
  min
}

# Linear functional relationship
normal_minimum_linear <- function(){
  s <- sample(1:10000, 1)
  
  # Make all the coefficients as random
  set.seed(s)
  beta0 <- runif(n = 1, min = 10, max = 50)
  beta1 <- runif(n = 1, min = 0.001, max = 0.1)
  beta2 <- runif(n = 1, min = 0.001, max = 0.1)
  beta3 <- runif(n = 1, min = 0.001, max = 0.1)
  beta4 <- runif(n = 1, min = 0.0001, max = 0.001)
  beta5 <- runif(n = 1, min = 0.001, max = 0.1)
  
  f <- function(X) beta1*X[,1] - beta3*X[,2] + beta5*X[,3]
  
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        method = "normal",
                                        # varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = seq(0.9, -0.9, by = -0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5,
                                        N = 1000)
  min_train <- min(simulated_data$train$y)
  min_novelty <- min(simulated_data$`6Range`$y)
  min <- min(min_train, min_novelty)
  min
}


parallel::detectCores()
n.cores <- 14

# Create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Reproducibility for parallel processing
registerDoRNG(1)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

normal_minimum_simulation_product <- foreach(icount(100), .combine = 'c') %dopar% {
  normal_minimum_product()
}

normal_minimum_simulation_quad <- foreach(icount(100), .combine = 'c') %dopar% {
  normal_minimum_quad()
}

normal_minimum_simulation_linear <- foreach(icount(100), .combine = 'c') %dopar% {
  normal_minimum_linear()
}

# Stop the cluster
parallel::stopCluster(cl = my.cluster)

min(normal_minimum_simulation_product)
hist(normal_minimum_simulation_product, main = "Product functional relationship")

min(normal_minimum_simulation_quad)
hist(normal_minimum_simulation_quad, main = "Quadratic functional relationship")

min(normal_minimum_simulation_linear)
hist(normal_minimum_simulation_linear, main = "Linear functional relationship")
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# Date: 04/23/2023
# Change: Find the range of y in training data

normal_minimum_product_sd  <- function(){
  s <- sample(1:10000, 1)
  
  # Make all the coefficients as random
  set.seed(s)
  beta0 <- runif(n = 1, min = -10, max = 10)
  beta1 <- runif(n = 1, min = -10, max = 10)
  beta2 <- runif(n = 1, min = -10, max = 10)
  beta3 <- runif(n = 1, min = -10, max = 10)
  beta4 <- runif(n = 1, min = -10, max = 10)
  beta5 <- runif(n = 1, min = -10, max = 10)
  
  f <- function(X) beta1*X[,1] + beta2*X[,1]^2 + beta3*X[,2] + beta4*(X[,1])*(X[,2]) + beta5*X[,3] + beta0
  
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        method = "normal",
                                        # varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = seq(0.9, -0.9, by = -0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        N = 1000)
  # min_train <- min(simulated_data$train$y)
  # min_novelty <- min(simulated_data$`6Range`$y)
  # min <- min(min_train, min_novelty)
  # min
  hist(simulated_data$train$y)
  mean(simulated_data$train$y)
}

#################################################################################################################################
#################################################################################################################################
# Data manipulation for the downstream analysis
# Treat collinearity shift as a factor first
# The number of rows for the data is 2, 376, 000
#################################################################################################################################
#################################################################################################################################
# Small test
# Transpose all collinearity shift
library(reshape2)

data_transpose <- function(data, algorithm = "GLM", 
                           model = "Interaction", 
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

# High training collinearity + Interaction
data_high_inter_rrmse <- data_transpose(data = high_inter_rrmse, algorithm = "GLM", model = "Interaction", 
                                        training_collinearity = "High", training_correlation = 0.9)
data_high_inter_rrmse_rf <- data_transpose(data = high_inter_rrmse_rf, algorithm = "RF", model = "Interaction", 
                                           training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Quadratic
data_high_quad_rrmse <- data_transpose(data = high_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                       training_collinearity = "High", training_correlation = 0.9)
data_high_quad_rrmse_rf <- data_transpose(data = high_quad_rrmse_rf, algorithm = "RF", model = "Quadratic", 
                                          training_collinearity = "High", training_correlation = 0.9)
# High training collinearity + Linear
data_high_linear_rrmse <- data_transpose(data = high_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                         training_collinearity = "High", training_correlation = 0.9)
data_high_linear_rrmse_rf <- data_transpose(data = high_linear_rrmse_rf, algorithm = "RF", model = "Linear", 
                                            training_collinearity = "High", training_correlation = 0.9)

# Mid training collinearity + Interaction
data_mid_inter_rrmse <- data_transpose(data = mid_inter_rrmse, algorithm = "GLM", model = "Interaction", 
                                       training_collinearity = "Mid", training_correlation = 0.6)
data_mid_inter_rrmse_rf <- data_transpose(data = mid_inter_rrmse_rf, algorithm = "RF", model = "Interaction", 
                                          training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Quadratic
data_mid_quad_rrmse <- data_transpose(data = mid_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Mid", training_correlation = 0.6)
data_mid_quad_rrmse_rf <- data_transpose(data = mid_quad_rrmse_rf, algorithm = "RF", model = "Quadratic", 
                                         training_collinearity = "Mid", training_correlation = 0.6)
# Mid training collinearity + Linear
data_mid_linear_rrmse <- data_transpose(data = mid_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Mid", training_correlation = 0.6)
data_mid_linear_rrmse_rf <- data_transpose(data = mid_linear_rrmse_rf, algorithm = "RF", model = "Linear", 
                                           training_collinearity = "Mid", training_correlation = 0.6)

# Low training collinearity + Interaction
data_low_inter_rrmse <- data_transpose(data = low_inter_rrmse, algorithm = "GLM", model = "Interaction", 
                                       training_collinearity = "Low", training_correlation = 0.3)
data_low_inter_rrmse_rf <- data_transpose(data = low_inter_rrmse_rf, algorithm = "RF", model = "Interaction", 
                                          training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Quadratic
data_low_quad_rrmse <- data_transpose(data = low_quad_rrmse, algorithm = "GLM", model = "Quadratic", 
                                      training_collinearity = "Low", training_correlation = 0.3)
data_low_quad_rrmse_rf <- data_transpose(data = low_quad_rrmse_rf, algorithm = "RF", model = "Quadratic", 
                                         training_collinearity = "Low", training_correlation = 0.3)
# Low training collinearity + Linear
data_low_linear_rrmse <- data_transpose(data = low_linear_rrmse, algorithm = "GLM", model = "Linear", 
                                        training_collinearity = "Low", training_correlation = 0.3)
data_low_linear_rrmse_rf <- data_transpose(data = low_linear_rrmse_rf, algorithm = "RF", model = "Linear", 
                                           training_collinearity = "Low", training_correlation = 0.3)

# Combine all data
data_analysis_normal <- rbind(data_high_inter_rrmse, data_high_inter_rrmse_rf, data_high_quad_rrmse, 
                              data_high_quad_rrmse_rf, data_high_linear_rrmse, data_high_linear_rrmse_rf,
                              data_mid_inter_rrmse, data_mid_inter_rrmse_rf, data_mid_quad_rrmse, 
                              data_mid_quad_rrmse_rf, data_mid_linear_rrmse, data_mid_linear_rrmse_rf,
                              data_low_inter_rrmse, data_low_inter_rrmse_rf, data_low_quad_rrmse, 
                              data_low_quad_rrmse_rf, data_low_linear_rrmse, data_low_linear_rrmse_rf)




