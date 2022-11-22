###################################################################################################################################
# Collinearity shift to -0.9 to 0.9                                                                                               #
# Range shift might be limited to less than 2 times of the range of training data                                                 #
# Change the sign of correlation in the test data by using Range(X) - X if using binomial + normal noise with decay for X's       #
# Predictor variables are generated from multinomial distributions with fixed correlation, mu, and sigma                          #
###################################################################################################################################

# 02/22/2022
library(faux)
library(mvtnorm)
library(data.table)
library(randomForest)
library(ggplot2)
library(data.table)

collineariser_extra <- function(N = 1000, 
                                varnames = c("X1", "X2", "X3"),
                                mu.training = c(0, 0, 0),
                                sd.training = c(1, 1, 1),
                                rho.training = c(0.5, 0.3, 0.2),
                                empirical = FALSE,
                                rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                                range.vec = c(.2, .4, .6, 2, 4, 6),
                                method = "normal",
                                seed = 1001,
                                y.noise = 0.5, scale = FALSE,
                                FUN = function(X) 25 + 2*X[,1] - 1.5*X[,2] - 2*X[,3], ...){
  
  # Function to generate collinear data with the following properties:
  # * 4 clusters with k variables each
  # * 1 uncorrelated variable
  # * y as a function of Xs
  #
  # k       number of variables for each of the four clusters
  # N       number of data points
  # decay   how fast shall the strength of collinearity decay in the most
  #         highly correlated cluster? good value range is 0.01 to 0.25:
  #         0.01 is extremely high, with no correlation below 0.64!
  #         0.02 min r = 0.5 (in least collinear cluster), max r = 0.99 (in most collinear cluster)
  #         0.1 is good collinearity: max around 0.7 and min (low cluster) of absent
  #         0.25 is practically very low correlation!
  # plotit  returns an image of the correlation matrix ("image") or a
  #         histogram ("hist") or nothing ("no")
  # method  "normal" produces normally distributed y-values, any other binary
  # seed    set random seed; or NULL
  # scale   shall the data be scaled to mean=0 and sd=1? Defaults to FALSE.
  
  # returns a list with 6 data sets:
  # 1. train, 2. test.same, 3. test.more, 4. test.less, 5. test.nonlin and 6.
  # test.none with appropriately adapted collinearity structure.
  # The first variable, y, is the response and it is NOT scaled!
  # Variables X1:X5, X6:X10, ... form a cluster each; X21 is uncorrelated
  # with all variables.
  # The functional relationship between y and X is depicted in the attributes
  # of the list.
  
  
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
  
  # scale all X-data sets:
  if (scale) {Data.scaled <- lapply(Data, scale)} else {Data.scaled = Data}
  
  # make y for each Data set:
  
  
  # X <- Xmaker(N=2*N, ...)
  #curve(2*x-1.5x*x, from=-2, to=2)
  f <- FUN
  if (method=="normal"){
    Ys <- lapply(Data.scaled, function(X) rnorm(n = nrow(X), mean = f(X), sd = y.noise))
  } else {
    Ys <- lapply(Data.scaled,
                 function(X) rbinom(nrow(X), size=1, plogis(rnorm(n=nrow(X), mean=f(X)-25, s = y.noise)) ))
  }
  
  Data.final <- Data.scaled
  attr(Data.final, "function") <- attr(f, "source")
  
  noise <- rnorm(length(Ys[[1]]), mean = 0, sd = 0.15*mean(Ys[[1]]))
  
  # Only add the noise in the training data would not work because the absolute RMSE between training and test
  # would be this noise if there is no any other effects
  Data.final[[1]] <- cbind(y = Ys[[1]] + noise, Data.scaled[[1]])
  
  for (i in 1:(length(Data)-1))
    # No noise in the test data
    Data.final[[i+1]] <- cbind(y = Ys[[i+1]], Data.scaled[[i+1]])
  
  # Only add the noise in the training data would not work because the absolute RMSE between training and test
  # would be this noise if there is no any other effects
  
  #for (i in 1:length(Data))
  # Should add the same noise to test data
  #Data.final[[i]] <- cbind(y = Ys[[i]] + noise, Data.scaled[[i]])
  
  Data.final
  
}


# Check the scatter plot of X1 and X2 under each scenario
# Testing the modified function with a user defined functional relationship using X1, X2, and X3
set.seed(1001)
data_try <- collineariser_extra(N = 1000,
                                varnames = c("X1", "X2", "X3"),
                                mu.training = c(0, 0, 0),
                                sd.training = c(1, 1, 1),
                                rho.training = c(0.3, 0.2, 0.1),
                                empirical = FALSE,
                                rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0,
                                             -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9),
                                range.vec = c(.2, .4, .6, 2, 4, 6),
                                method="normal",
                                seed = 1001,
                                y.noise = 0.5, scale = FALSE,
                                FUN=function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 1.5*(X[,1])*(X[,2])- 3*X[,3])

# Plot how X1,X2 change in different scenarios
train_try <- data.frame(data_try$train)
train_try$Group <- "A"

r_0.9_try <- data.frame(data_try$`r=0.9`)
r_0.9_try$Group <- "B"
r_0.9_try <- rbind(train_try, r_0.9_try)

r_0.8_try <- data.frame(data_try$`r=0.8`)
r_0.8_try$Group <- "B"
r_0.8_try <- rbind(train_try, r_0.8_try)

r_0.7_try <- data.frame(data_try$`r=0.7`)
r_0.7_try$Group <- "B"
r_0.7_try <- rbind(train_try, r_0.7_try)

r_0.6_try <- data.frame(data_try$`r=0.6`)
r_0.6_try$Group <- "B"
r_0.6_try <- rbind(train_try, r_0.6_try)

r_0.5_try <- data.frame(data_try$`r=0.5`)
r_0.5_try$Group <- "B"
r_0.5_try <- rbind(train_try, r_0.5_try)

r_0.4_try <- data.frame(data_try$`r=0.4`)
r_0.4_try$Group <- "B"
r_0.4_try <- rbind(train_try, r_0.4_try)

r_0.3_try <- data.frame(data_try$`r=0.3`)
r_0.3_try$Group <- "B"
r_0.3_try <- rbind(train_try, r_0.3_try)

r_0.2_try <- data.frame(data_try$`r=0.2`)
r_0.2_try$Group <- "B"
r_0.2_try <- rbind(train_try, r_0.2_try)

r_0.1_try <- data.frame(data_try$`r=0.1`)
r_0.1_try$Group <- "B"
r_0.1_try <- rbind(train_try, r_0.1_try)

r_0_try <- data.frame(data_try$`r=0`)
r_0_try$Group <- "B"
r_0_try <- rbind(train_try, r_0_try)

r_neg_0.9_try <- data.frame(data_try$`r=-0.9`)
r_neg_0.9_try$Group <- "B"
r_neg_0.9_try <- rbind(train_try, r_neg_0.9_try)

r_neg_0.8_try <- data.frame(data_try$`r=-0.8`)
r_neg_0.8_try$Group <- "B"
r_neg_0.8_try <- rbind(train_try, r_neg_0.8_try)

r_neg_0.7_try <- data.frame(data_try$`r=-0.7`)
r_neg_0.7_try$Group <- "B"
r_neg_0.7_try <- rbind(train_try, r_neg_0.7_try)

r_neg_0.6_try <- data.frame(data_try$`r=-0.6`)
r_neg_0.6_try$Group <- "B"
r_neg_0.6_try <- rbind(train_try, r_neg_0.6_try)

r_neg_0.5_try <- data.frame(data_try$`r=-0.5`)
r_neg_0.5_try$Group <- "B"
r_neg_0.5_try <- rbind(train_try, r_neg_0.5_try)

r_neg_0.4_try <- data.frame(data_try$`r=-0.4`)
r_neg_0.4_try$Group <- "B"
r_neg_0.4_try <- rbind(train_try, r_neg_0.4_try)

r_neg_0.3_try <- data.frame(data_try$`r=-0.3`)
r_neg_0.3_try$Group <- "B"
r_neg_0.3_try <- rbind(train_try, r_neg_0.3_try)

r_neg_0.2_try <- data.frame(data_try$`r=-0.2`)
r_neg_0.2_try$Group <- "B"
r_neg_0.2_try <- rbind(train_try, r_neg_0.2_try)

r_neg_0.1_try <- data.frame(data_try$`r=-0.1`)
r_neg_0.1_try$Group <- "B"
r_neg_0.1_try <- rbind(train_try, r_neg_0.1_try)

plot_col_try <- rbind(r_0.9_try, r_0.8_try, r_0.7_try, r_0.6_try, r_0.5_try,
                      r_0.4_try, r_0.3_try, r_0.2_try, r_0.1_try, r_0_try,
                      r_neg_0.9_try, r_neg_0.8_try, r_neg_0.7_try, r_neg_0.6_try, r_neg_0.5_try,
                      r_neg_0.4_try, r_neg_0.3_try, r_neg_0.2_try, r_neg_0.1_try)

plot_col_try$Shift <- "Col shift"
plot_col_try$label <- c(
  rep("r = 0.9\n", nrow(r_0.9_try)),
  rep("r = 0.8\n", nrow(r_0.8_try)),
  rep("r = 0.7\n", nrow(r_0.7_try)),
  rep("r = 0.6\n", nrow(r_0.6_try)),
  rep("r = 0.5\n", nrow(r_0.5_try)),
  rep("r = 0.4\n", nrow(r_0.4_try)),
  rep("r = 0.3\nin training data", nrow(r_0.3_try)),
  rep("r = 0.2\n", nrow(r_0.2_try)),
  rep("r = 0.1\n", nrow(r_0.1_try)),
  rep("r = 0\n", nrow(r_0.1_try)),
  rep("r = -0.9\n", nrow(r_neg_0.9_try)),
  rep("r = -0.8\n", nrow(r_neg_0.8_try)),
  rep("r = -0.7\n", nrow(r_neg_0.7_try)),
  rep("r = -0.6\n", nrow(r_neg_0.6_try)),
  rep("r = -0.5\n", nrow(r_neg_0.5_try)),
  rep("r = -0.4\n", nrow(r_neg_0.4_try)),
  rep("r = -0.3\n", nrow(r_neg_0.3_try)),
  rep("r = -0.2\n", nrow(r_neg_0.2_try)),
  rep("r = -0.1\n", nrow(r_neg_0.1_try))
)

plot_col_try$Shift <- c(rep("Positive\ncollinearity", 2000*10), rep("Negative\ncollinearity", 2000*9))

Range4_r_0.9_try <- data.frame(data_try$`4Range*r=0.9`)
Range4_r_0.9_try$Group <- "B"
Range4_r_0.9_try <- rbind(train_try, Range4_r_0.9_try)

Range4_r_0.8_try <- data.frame(data_try$`4Range*r=0.8`)
Range4_r_0.8_try$Group <- "B"
Range4_r_0.8_try <- rbind(train_try, Range4_r_0.8_try)

Range4_r_0.7_try <- data.frame(data_try$`4Range*r=0.7`)
Range4_r_0.7_try$Group <- "B"
Range4_r_0.7_try <- rbind(train_try, Range4_r_0.7_try)

Range4_r_0.6_try <- data.frame(data_try$`4Range*r=0.6`)
Range4_r_0.6_try$Group <- "B"
Range4_r_0.6_try <- rbind(train_try, Range4_r_0.6_try)

Range4_r_0.5_try <- data.frame(data_try$`4Range*r=0.5`)
Range4_r_0.5_try$Group <- "B"
Range4_r_0.5_try <- rbind(train_try, Range4_r_0.5_try)

Range4_r_0.4_try <- data.frame(data_try$`4Range*r=0.4`)
Range4_r_0.4_try$Group <- "B"
Range4_r_0.4_try <- rbind(train_try, Range4_r_0.4_try)

Range4_r_0.3_try <- data.frame(data_try$`4Range*r=0.3`)
Range4_r_0.3_try$Group <- "B"
Range4_r_0.3_try <- rbind(train_try, Range4_r_0.3_try)

Range4_r_0.2_try <- data.frame(data_try$`4Range*r=0.2`)
Range4_r_0.2_try$Group <- "B"
Range4_r_0.2_try <- rbind(train_try, Range4_r_0.2_try)

Range4_r_0.1_try <- data.frame(data_try$`4Range*r=0.1`)
Range4_r_0.1_try$Group <- "B"
Range4_r_0.1_try <- rbind(train_try, Range4_r_0.1_try)

Range4_r_0_try <- data.frame(data_try$`4Range*r=0`)
Range4_r_0_try$Group <- "B"
Range4_r_0_try <- rbind(train_try, Range4_r_0_try)

Range4_r_neg_0.9_try <- data.frame(data_try$`4Range*r=-0.9`)
Range4_r_neg_0.9_try$Group <- "B"
Range4_r_neg_0.9_try <- rbind(train_try, Range4_r_neg_0.9_try)

Range4_r_neg_0.8_try <- data.frame(data_try$`4Range*r=-0.8`)
Range4_r_neg_0.8_try$Group <- "B"
Range4_r_neg_0.8_try <- rbind(train_try, Range4_r_neg_0.8_try)

Range4_r_neg_0.7_try <- data.frame(data_try$`4Range*r=-0.7`)
Range4_r_neg_0.7_try$Group <- "B"
Range4_r_neg_0.7_try <- rbind(train_try, Range4_r_neg_0.7_try)

Range4_r_neg_0.6_try <- data.frame(data_try$`4Range*r=-0.6`)
Range4_r_neg_0.6_try$Group <- "B"
Range4_r_neg_0.6_try <- rbind(train_try, Range4_r_neg_0.6_try)

Range4_r_neg_0.5_try <- data.frame(data_try$`4Range*r=-0.5`)
Range4_r_neg_0.5_try$Group <- "B"
Range4_r_neg_0.5_try <- rbind(train_try, Range4_r_neg_0.5_try)

Range4_r_neg_0.4_try <- data.frame(data_try$`4Range*r=-0.4`)
Range4_r_neg_0.4_try$Group <- "B"
Range4_r_neg_0.4_try <- rbind(train_try, Range4_r_neg_0.4_try)

Range4_r_neg_0.3_try <- data.frame(data_try$`4Range*r=-0.3`)
Range4_r_neg_0.3_try$Group <- "B"
Range4_r_neg_0.3_try <- rbind(train_try, Range4_r_neg_0.3_try)

Range4_r_neg_0.2_try <- data.frame(data_try$`4Range*r=-0.2`)
Range4_r_neg_0.2_try$Group <- "B"
Range4_r_neg_0.2_try <- rbind(train_try, Range4_r_neg_0.2_try)

Range4_r_neg_0.1_try <- data.frame(data_try$`4Range*r=-0.1`)
Range4_r_neg_0.1_try$Group <- "B"
Range4_r_neg_0.1_try <- rbind(train_try, Range4_r_neg_0.1_try)

# Visualize all the scatter plots for Range 2
Range2_r_0.9_try <- data.frame(data_try$`2Range*r=0.9`)
Range2_r_0.9_try$Group <- "B"
Range2_r_0.9_try <- rbind(train_try, Range2_r_0.9_try)

Range2_r_0.8_try <- data.frame(data_try$`2Range*r=0.8`)
Range2_r_0.8_try$Group <- "B"
Range2_r_0.8_try <- rbind(train_try, Range2_r_0.8_try)

Range2_r_0.7_try <- data.frame(data_try$`2Range*r=0.7`)
Range2_r_0.7_try$Group <- "B"
Range2_r_0.7_try <- rbind(train_try, Range2_r_0.7_try)

Range2_r_0.6_try <- data.frame(data_try$`2Range*r=0.6`)
Range2_r_0.6_try$Group <- "B"
Range2_r_0.6_try <- rbind(train_try, Range2_r_0.6_try)

Range2_r_0.5_try <- data.frame(data_try$`2Range*r=0.5`)
Range2_r_0.5_try$Group <- "B"
Range2_r_0.5_try <- rbind(train_try, Range2_r_0.5_try)

Range2_r_0.4_try <- data.frame(data_try$`2Range*r=0.4`)
Range2_r_0.4_try$Group <- "B"
Range2_r_0.4_try <- rbind(train_try, Range2_r_0.4_try)

Range2_r_0.3_try <- data.frame(data_try$`2Range*r=0.3`)
Range2_r_0.3_try$Group <- "B"
Range2_r_0.3_try <- rbind(train_try, Range2_r_0.3_try)

Range2_r_0.2_try <- data.frame(data_try$`2Range*r=0.2`)
Range2_r_0.2_try$Group <- "B"
Range2_r_0.2_try <- rbind(train_try, Range2_r_0.2_try)

Range2_r_0.1_try <- data.frame(data_try$`2Range*r=0.1`)
Range2_r_0.1_try$Group <- "B"
Range2_r_0.1_try <- rbind(train_try, Range2_r_0.1_try)

Range2_r_0_try <- data.frame(data_try$`2Range*r=0.0`)
Range2_r_0_try$Group <- "B"
Range2_r_0_try <- rbind(train_try, Range2_r_0_try)

Range2_r_neg_0.9_try <- data.frame(data_try$`2Range*r=-0.9`)
Range2_r_neg_0.9_try$Group <- "B"
Range2_r_neg_0.9_try <- rbind(train_try, Range2_r_neg_0.9_try)

Range2_r_neg_0.8_try <- data.frame(data_try$`2Range*r=-0.8`)
Range2_r_neg_0.8_try$Group <- "B"
Range2_r_neg_0.8_try <- rbind(train_try, Range2_r_neg_0.8_try)

Range2_r_neg_0.7_try <- data.frame(data_try$`2Range*r=-0.7`)
Range2_r_neg_0.7_try$Group <- "B"
Range2_r_neg_0.7_try <- rbind(train_try, Range2_r_neg_0.7_try)

Range2_r_neg_0.6_try <- data.frame(data_try$`2Range*r=-0.6`)
Range2_r_neg_0.6_try$Group <- "B"
Range2_r_neg_0.6_try <- rbind(train_try, Range2_r_neg_0.6_try)

Range2_r_neg_0.5_try <- data.frame(data_try$`2Range*r=-0.5`)
Range2_r_neg_0.5_try$Group <- "B"
Range2_r_neg_0.5_try <- rbind(train_try, Range2_r_neg_0.5_try)

Range2_r_neg_0.4_try <- data.frame(data_try$`2Range*r=-0.4`)
Range2_r_neg_0.4_try$Group <- "B"
Range2_r_neg_0.4_try <- rbind(train_try, Range2_r_neg_0.4_try)

Range2_r_neg_0.3_try <- data.frame(data_try$`2Range*r=-0.3`)
Range2_r_neg_0.3_try$Group <- "B"
Range2_r_neg_0.3_try <- rbind(train_try, Range2_r_neg_0.3_try)

Range2_r_neg_0.2_try <- data.frame(data_try$`2Range*r=-0.2`)
Range2_r_neg_0.2_try$Group <- "B"
Range2_r_neg_0.2_try <- rbind(train_try, Range2_r_neg_0.2_try)

Range2_r_neg_0.1_try <- data.frame(data_try$`2Range*r=-0.1`)
Range2_r_neg_0.1_try$Group <- "B"
Range2_r_neg_0.1_try <- rbind(train_try, Range2_r_neg_0.1_try)

plot_range_try <- rbind(Range4_r_0.9_try, Range4_r_0.8_try, Range4_r_0.7_try, Range4_r_0.6_try, Range4_r_0.5_try,
                        Range4_r_0.4_try, Range4_r_0.3_try, Range4_r_0.2_try, Range4_r_0.1_try, Range4_r_0_try,
                        Range4_r_neg_0.9_try, Range4_r_neg_0.8_try, Range4_r_neg_0.7_try, Range4_r_neg_0.6_try, Range4_r_neg_0.5_try,
                        Range4_r_neg_0.4_try, Range4_r_neg_0.3_try, Range4_r_neg_0.2_try, Range4_r_neg_0.1_try)

plot_range_try$Shift <- "4xRange + Col shift"

# check Range 2
plot_range_try <- rbind(Range2_r_0.9_try, Range2_r_0.8_try, Range2_r_0.7_try, Range2_r_0.6_try, Range2_r_0.5_try,
                        Range2_r_0.4_try, Range2_r_0.3_try, Range2_r_0.2_try, Range2_r_0.1_try, Range2_r_0_try,
                        Range2_r_neg_0.9_try, Range2_r_neg_0.8_try, Range2_r_neg_0.7_try, Range2_r_neg_0.6_try, Range2_r_neg_0.5_try,
                        Range2_r_neg_0.4_try, Range2_r_neg_0.3_try, Range2_r_neg_0.2_try, Range2_r_neg_0.1_try)

plot_range_try$Shift <- "2xRange + Col shift"

plot_range_try$label <- c(
  rep("r = 0.9", nrow(Range4_r_0.9_try)),
  rep("r = 0.8", nrow(Range4_r_0.8_try)),
  rep("r = 0.7", nrow(Range4_r_0.7_try)),
  rep("r = 0.6", nrow(Range4_r_0.6_try)),
  rep("r = 0.5", nrow(Range4_r_0.5_try)),
  rep("r = 0.4", nrow(Range4_r_0.4_try)),
  rep("r = 0.3", nrow(Range4_r_0.3_try)),
  rep("r = 0.2", nrow(Range4_r_0.2_try)),
  rep("r = 0.1", nrow(Range4_r_0.1_try)),
  rep("r = 0", nrow(Range4_r_0_try)),
  rep("r = -0.9", nrow(Range4_r_neg_0.9_try)),
  rep("r = -0.8", nrow(Range4_r_neg_0.8_try)),
  rep("r = -0.7", nrow(Range4_r_neg_0.7_try)),
  rep("r = -0.6", nrow(Range4_r_neg_0.6_try)),
  rep("r = -0.5", nrow(Range4_r_neg_0.5_try)),
  rep("r = -0.4", nrow(Range4_r_neg_0.4_try)),
  rep("r = -0.3", nrow(Range4_r_neg_0.3_try)),
  rep("r = -0.2", nrow(Range4_r_neg_0.2_try)),
  rep("r = -0.1", nrow(Range4_r_neg_0.1_try))
)

# Check Range 2
plot_range_try$label <- c(
  rep("2xRange\n*r=0.9", nrow(Range2_r_0.9_try)),
  rep("2xRange\n*r=0.8", nrow(Range2_r_0.8_try)),
  rep("2xRange\n*r=0.7", nrow(Range2_r_0.7_try)),
  rep("2xRange\n*r=0.6", nrow(Range2_r_0.6_try)),
  rep("2xRange\n*r=0.5", nrow(Range2_r_0.5_try)),
  rep("2xRange\n*r=0.4", nrow(Range2_r_0.4_try)),
  rep("2xRange\n*r=0.3", nrow(Range2_r_0.3_try)),
  rep("2xRange\n*r=0.2", nrow(Range2_r_0.2_try)),
  rep("2xRange\n*r=0.1", nrow(Range2_r_0.1_try)),
  rep("2xRange\n*r=0", nrow(Range2_r_0_try)),
  rep("2xRange\n*r=-0.9", nrow(Range2_r_neg_0.9_try)),
  rep("2xRange\n*r=-0.8", nrow(Range2_r_neg_0.8_try)),
  rep("2xRange\n*r=-0.7", nrow(Range2_r_neg_0.7_try)),
  rep("2xRange\n*r=-0.6", nrow(Range2_r_neg_0.6_try)),
  rep("2xRange\n*r=-0.5", nrow(Range2_r_neg_0.5_try)),
  rep("2xRange\n*r=-0.4", nrow(Range2_r_neg_0.4_try)),
  rep("2xRange\n*r=-0.3", nrow(Range2_r_neg_0.3_try)),
  rep("2xRange\n*r=-0.2", nrow(Range2_r_neg_0.2_try)),
  rep("2xRange\n*r=-0.1", nrow(Range2_r_neg_0.1_try))
)

plot_range_try$Shift <- c(rep("Positive\ncollinearity", 2000*10), rep("Negative\ncollinearity", 2000*9))

library(ggplot2)
library(cowplot)
# Less_collinearity_range_shift
plot_col_try_plot <- plot_col_try[plot_col_try$Shift == "Positive\ncollinearity"&plot_col_try$Group == "B",]

pos_col_shift <- ggplot(plot_col_try_plot[plot_col_try_plot$label == "r = 0.3\nin training data"|
                                            plot_col_try_plot$label == "r = 0.6\n"|
                                            plot_col_try_plot$label == "r = 0.9\n", ], aes(X1, X2)) + 
  geom_point(alpha = .2, color = "blue", shape = 19) + 
  # scale_color_manual(labels = c("test data"), values = c("blue")) +
  facet_grid(.~label) +
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(axis.title.y = element_text(size = 10, vjust = 3),
        axis.title.x = element_text(size = 9, vjust = 3, margin =margin(t = 10, r = 0, b = 0, l = 0)),
        #panel.border = element_blank(),
        #panel.grid.minor = element_blank(), # This could be deleted if the information is too much
        legend.position = "bottom",
        legend.margin=margin(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5))


zero_col_shift <- ggplot(plot_col_try_plot[plot_col_try_plot$label == "r = 0\n", ], aes(X1, X2)) + 
  geom_point(alpha = .2, color = "blue", shape = 19) + 
  # scale_color_manual(labels = c("test data"), values = c("blue")) +
  facet_grid(.~label) +
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(axis.title.y = element_text(size = 10, vjust = 3),
        axis.title.x = element_text(size = 9, vjust = 3, margin =margin(t = 10, r = 0, b = 0, l = 0)),
        #panel.border = element_blank(),
        #panel.grid.minor = element_blank(), # This could be deleted if the information is too much
        legend.position = "bottom",
        legend.margin=margin(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5))


# More_collinearity_range_shift
plot_col_try_plot_neg <- plot_col_try[plot_col_try$Shift == "Negative\ncollinearity"&plot_col_try$Group == "B",]

neg_col_shift <- ggplot(plot_col_try_plot_neg[plot_col_try_plot_neg$label == "r = -0.3\n"|
                                                plot_col_try_plot_neg$label == "r = -0.6\n"|
                                                plot_col_try_plot_neg$label == "r = -0.9\n", ], aes(X1, X2)) + 
  geom_point(alpha = .2, color = "red", shape = 19) + 
  # scale_color_manual(labels = c("training data", "test data"), values = c("none", "blue")) +
  facet_grid(.~label) +
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(axis.title.y = element_text(size = 10, vjust = 3),
        axis.title.x = element_text(size = 9, vjust = 3, margin =margin(t = 10, r = 0, b = 0, l = 0)),
        #panel.border = element_blank(),
        #panel.grid.minor = element_blank(), # This could be deleted if the information is too much
        legend.position = "bottom",
        legend.margin=margin(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5))

# Visualize X1,X2 along the gradients of collinearity shift
plot_grid(plot_grid(NULL, pos_col_shift, NULL, neg_col_shift, nrow = 4, rel_heights = c(0.14, 1, 0.14, 1)),
          plot_grid(NULL, zero_col_shift, NULL, NULL, nrow = 4, rel_heights = c(0.14, 1, 0.14, 1)),
          ncol = 2,
          rel_widths = c(1, .4),
          label_size = 10,
          label_x = 0, label_y = 1
)


# Visualize X1,X2 along the gradients of collinearity shift
plot_grid(NULL, 
          ggplot(plot_col_try[plot_col_try$Shift == "Positive\ncollinearity", ], aes(X1, X2.2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
            labs(y = expression("X2^2")) +
            facet_grid(Shift~label) +  
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "none"),
          NULL,
          ggplot(plot_col_try[plot_col_try$Shift == "Negative\ncollinearity", ], aes(X1, X2.2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
            labs(y = expression("X2^2")) +
            facet_grid(Shift~label) +
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "bottom"),
          nrow = 4,
          labels = c("(a)" ,"", "(b)", ""),
          label_size = 10,
          rel_heights = c(0.15, 1, 0.15, 1.3),
          label_x = 0, label_y = 1
)

# Visualize X1,X1X2 along the gradients of collinearity shift
plot_grid(NULL, 
          ggplot(plot_col_try[plot_col_try$Shift == "Positive\ncollinearity", ], aes(X1, X1.X2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
            labs(y = expression("X1*X2")) +
            facet_grid(Shift~label) +  
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "none"),
          NULL,
          ggplot(plot_col_try[plot_col_try$Shift == "Negative\ncollinearity", ], aes(X1, X1.X2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
            labs(y = expression("X1*X2")) +
            facet_grid(Shift~label) +
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "bottom"),
          nrow = 4,
          labels = c("(a)" ,"", "(b)", ""),
          label_size = 10,
          rel_heights = c(0.15, 1, 0.15, 1.3),
          label_x = 0, label_y = 1
)

# Visualize X1,X2 along the gradients of collinearity and range shift 
ggplot(plot_range_try[plot_range_try$Shift == "Positive\ncollinearity", ], aes(X1, X3, color = Group)) + 
  geom_point(alpha = .3) + 
  scale_color_manual(labels = c("training data", "test data"), values = c("grey", "red")) +
  labs(y = expression("X2")) +
  facet_grid(Shift~label) +  
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none")

# More_collinearity_range_col_shift 
ggplot(plot_range_try[plot_range_try$Shift == "Negative\ncollinearity", ], aes(X1, X2, color = Group)) + 
  geom_point(alpha = .3) + 
  scale_color_manual(labels = c("training data", "test data"), values = c("grey", "red")) +
  labs(y = expression("X2")) +
  facet_grid(Shift~label) +
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "bottom")

plot_grid(NULL, 
          ggplot(plot_range_try[plot_range_try$Shift == "Positive\ncollinearity", ], aes(X1, X2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "red")) +
            labs(y = expression("X1*X2")) +
            facet_grid(Shift~label) +  
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "none"),
          NULL,
          ggplot(plot_range_try[plot_range_try$Shift == "Negative\ncollinearity", ], aes(X1, X2, color = Group)) + 
            geom_point(alpha = .3) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "red")) +
            labs(y = expression("X1*X2")) +
            facet_grid(Shift~label) +
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8),
                  legend.position = "bottom"),
          nrow = 4,
          labels = c("(a)" ,"", "(b)", ""),
          label_size = 10,
          rel_heights = c(0.15, 1, 0.15, 1.3),
          label_x = 0, label_y = 1
)




# Create the figure in method section
# Load package
library(ggplot2)
library(cowplot)

# Subset the data
plot_col_try_method <- plot_col_try[plot_col_try$label == "r = -0.9"|plot_col_try$label == "r = -0.6"|
                                      plot_col_try$label == "r = -0.3"|plot_col_try$label == "r = 0.3"|
                                      plot_col_try$label == "r = 0.6"|plot_col_try$label == "r = 0.9", ]

plot_col_try_method$cor[plot_col_try_method$label == "r = -0.9"] <- "A"
plot_col_try_method$cor[plot_col_try_method$label == "r = -0.6"] <- "B"
plot_col_try_method$cor[plot_col_try_method$label == "r = -0.3"] <- "C"
plot_col_try_method$cor[plot_col_try_method$label == "r = 0.3"] <- "D"
plot_col_try_method$cor[plot_col_try_method$label == "r = 0.6"] <- "E"
plot_col_try_method$cor[plot_col_try_method$label == "r = 0.9"] <- "F"

plot_col_try_method$cor <- factor(plot_col_try_method$cor, levels = c("A", "B", "C", "D", "E", "F"), 
                                  labels = c("r = -0.9", "r = -0.6", "r = -0.3", "r = 0.3", "r = 0.6", "r = 0.9"))

plot_col_try_method$range <- "A"

plot_range_try_method <- plot_range_try[plot_range_try$label == "r = -0.9"|plot_range_try$label == "r = -0.6"|
                                          plot_range_try$label == "r = -0.3"|plot_range_try$label == "r = 0.3"|
                                          plot_range_try$label == "r = 0.6"|plot_range_try$label == "r = 0.9", ]

plot_range_try_method$cor[plot_range_try_method$label == "r = -0.9"] <- "A"
plot_range_try_method$cor[plot_range_try_method$label == "r = -0.6"] <- "B"
plot_range_try_method$cor[plot_range_try_method$label == "r = -0.3"] <- "C"
plot_range_try_method$cor[plot_range_try_method$label == "r = 0.3"] <- "D"
plot_range_try_method$cor[plot_range_try_method$label == "r = 0.6"] <- "E"
plot_range_try_method$cor[plot_range_try_method$label == "r = 0.9"] <- "F"

plot_range_try_method$cor <- factor(plot_range_try_method$cor, levels = c("A", "B", "C", "D", "E", "F"), 
                                    labels = c("r = -0.9", "r = -0.6", "r = -0.3", "r = 0.3", "r = 0.6", "r = 0.9"))

plot_range_try_method$range <- "B"


plot_extra_method <- rbind(plot_col_try_method, plot_range_try_method)
plot_extra_method$range <- factor(plot_extra_method$range, levels = c("A", "B"), 
                                  labels = c("Same\nrange", "1.6 times\nrange"))

ggplot(plot_extra_method, aes(X1, X2, color = Group)) + 
  geom_point(alpha = .3, size = .8) + 
  scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
  facet_grid(range~cor) +
  theme_light() +
  theme(strip.background = element_rect(color = "black", fill = "light grey", size = .5, linetype = "solid"),
        strip.text.x = element_text(color = "black", size = 8),
        strip.text.y = element_text(color = "black", size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "light grey"),
        panel.border = element_rect(colour = "black", size = 1),
        legend.position = "none")

# Visualize X1,X2 along the gradients of collinearity shift
plot_grid(ggplot(plot_extra_method, aes(X1, X2, color = Group)) + 
            geom_point(alpha = .3, size = .8) + 
            scale_color_manual(labels = c("training data", "test data"), values = c("grey", "blue")) +
            facet_grid(range~cor) +
            theme_light() +
            theme(strip.background = element_rect(color = "black", fill = "light grey", size = .5, linetype = "solid"),
                  strip.text.x = element_text(color = "black", size = 8),
                  strip.text.y = element_text(color = "black", size = 8),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(color = "light grey"),
                  panel.border = element_rect(colour = "black", size = 1),
                  legend.position = "none"),
          nrow = 1,
          labels = c("(a)"),
          label_size = 10,
          rel_heights = c(0.15, 1),
          label_x = 0, label_y = 1
)



# set.seed(1001) 
# data_try <- collineariser_4(N = 1000,
#                             varnames = c("X1", "X2", "X3"),
#                             mu.training = c(0, 0, 0),
#                             sd.training = c(1, 1, 1),
#                             rho.training = c(0.9, 0.3, 0.2),
#                             empirical = FALSE,
#                             rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
#                             range.vec = c(2, 3, 4),
#                             method="normal",
#                             seed = round(max(rho.training)*1000),y.noise = 0.5, scale = FALSE,
#                             FUN=function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 1.5*(X[,1])*(X[,2])- 3*X[,3])


# Simulator based on collineariser_extra and keep track of RMSE, coefficients, and other attributes
linear_model_1 <- function(vars = c("X1", "X2", "X3"),
                           f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] - 1.5*(X[,1])*(X[,2]) - 3*X[,3], 
                           mu = c(0, 0, 0),
                           sd = c(1, 1, 1),
                           rho = c(0.5, 0.3, 0.2),
                           response = "y",
                           predictors = c("X1", "X2", "I(X1^2)", "X1*X2", "X3")){
  
  s <- sample(1:10000, 1)
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0,
                                                     -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5)
  
  # Fit a linear model for the simulated data
  formula <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  fitted <- lm(data = as.data.frame(simulated_data$train), formula)
  
  # Calculate metrics
  
  # RMSE
  rmse <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted, as.data.frame(X))-as.data.frame(X)$y)^2))))
  
  # Calculate the condition number for all predictors
  cn_all <-  lapply(simulated_data, function(X)(kappa(cor(X[,2:8]), exact = TRUE)))
  
  # Save the training data and test data with test_same and col = -0.9
  data_kept <- simulated_data[c("train", "test_same", "r=-0.9")]
  
  # Extract estimates of betas and R-squares
  model_r_square <- summary(fitted)$r.square
  
  coef_intercept <- fitted$coefficients[1]
  coef_X1 <- fitted$coefficients[2]
  coef_X2 <- fitted$coefficients[3]
  coef_X1_square <- fitted$coefficients[4]
  coef_X3 <- fitted$coefficients[5]
  coef_X1X2 <- fitted$coefficients[6]
  
  return(list("RMSE" = rmse, "CN_all" = cn_all, "data_kept" = data_kept, 
              "R_square" = model_r_square,
              "Coef_intercept" = coef_intercept,
              "Coef_X1" = coef_X1,
              "Coef_X2" = coef_X2,
              "Coef_X1_square" = coef_X1_square,
              "Coef_X3" = coef_X3,
              "Coef_X1X2" = coef_X1X2))
}

# Simulator based on collineariser_4 with randomForests
linear_model_2 <- function(vars = c("X1", "X2", "X3"),
                           f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2]) - 3*X[,3], 
                           mu = c(0, 0, 0),
                           sd = c(1, 1, 1),
                           rho = c(0.5, 0.3, 0.2),
                           response = "y",
                           predictors = c("X1", "X2", "I(X1^2)", "X1*X2", "X3")){
  
  s <- sample(1:10000, 1)
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5)
  
  # Fit a linear model for the simulated data
  formula <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  fitted <- lm(data = as.data.frame(simulated_data$train), formula)
  
  # # Fit a RF model for the simulated data
  # Based on Random Forest Prediction Intervals. Haozhe Zhang, Joshua Zimmerman, Dan Nettleton, and Dan Nordman.
  # The American Statistician, 74(4):392-406, 2020, specifying the relationship between Y and X may increase the
  # model extrapolation of Random Forests
  # formula_rf <- as.formula(
  #   paste(response,
  #         paste(c("X1", "X2", "X3"), collapse = "+"),
  #         sep = "~")
  # )
  
  # Keep the same formula used in developing GLMs for RF
  fitted_rf <- randomForest(formula, data = as.data.frame(simulated_data$train))
  
  # Calculate metrics
  
  # RMSE
  rmse <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted, as.data.frame(X))-as.data.frame(X)$y)^2))))
  rmse_rf <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted_rf, as.data.frame(X))-as.data.frame(X)$y)^2))))
  
  # Variance
  variance <- lapply(simulated_data, function(X) (mean((mean(predict(fitted, as.data.frame(X))) - predict(fitted, as.data.frame(X)))^2)))
  
  # Calculate the correlation between X1 and X2
  cor_pair <- lapply(simulated_data, function(X)(cor(X[,c(2,4)])))
  
  # Calculate the upper triangle of all predictors
  cor_uppertri <- lapply(simulated_data, function(X)(sum(cor(X[,2:4])[upper.tri(cor(X[,2:4]), diag = FALSE)])))
  
  # Calculate the condition number for X1,X2,X3
  cn_essential <-  lapply(simulated_data, function(X)(kappa(cor(X[,2:4]), exact = TRUE)))
  
  # Calculate the condition number for other terms
  cn_non_essential <-  lapply(simulated_data, function(X)(kappa(cor(X[,5:8]), exact = TRUE)))
  
  # Calculate the condition number for all predictors
  cn_all <-  lapply(simulated_data, function(X)(kappa(cor(X[,2:8]), exact = TRUE)))
  
  # Retrieve R-square and standard errors of coefficients
  model_r_square <- summary(fitted)$r.square
  std_error_X2 <- summary(fitted)$coefficients[2, 2]
  
  rf_r_square <- fitted_rf$rsq[length(fitted_rf$rsq)]
  
  #rmse_variance_list <- list("Rmse" = rmse, "Variance" = variance)
  
  return(list("RMSE" = rmse, "RMSE_RF" = rmse_rf, "Variance" = variance, "Cor_pair" = cor_pair, "Cor_uppertri" = cor_uppertri,
              "CN_essential" = cn_essential, "CN_non_essential" = cn_non_essential, "CN_all" = cn_all,
              "R_square" = model_r_square, "Std_error" = std_error_X2, "R_square_rf" = rf_r_square))
}

# Run simulation for all training collinearity
#  Sys.time()
#  "2022-03-16 18:03:04 EDT"

# Sys.time()
#  "2022-03-16 22:43:09 EDT"


# Date: 07-15-2022
# Change: Add the M distance to the returned list


# RF took much time so making all RF related commands as comments first
linear_model_2 <- function(vars = c("X1", "X2", "X3"),
                           f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2]) - 3*X[,3], 
                           mu = c(0, 0, 0),
                           sd = c(1, 1, 1),
                           rho = c(0.5, 0.3, 0.2),
                           response = "y",
                           predictors = c("X1", "X2", "I(X2^2)", "X1*X2", "X3"),
                           predictors_index = c("X1", "X2", "X2^2", "X1*X2", "X3")){
  
  s <- sample(1:10000, 1)
  # Create simulated data set
  simulated_data <- collineariser_extra(FUN = f,
                                        varnames = vars,
                                        sd.training = sd,
                                        mu.training = mu,
                                        rho.training = rho,
                                        rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0,
                                                     -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9),
                                        range.vec = c(.2, .4, .6, 2, 4, 6),
                                        seed = s,
                                        empirical = FALSE,
                                        y.noise = 0.5)
  
  # Fit a linear model for the simulated data
  formula <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  fitted <- lm(data = as.data.frame(simulated_data$train), formula)
  
  # # Fit a RF model for the simulated data
  # Based on Random Forest Prediction Intervals. Haozhe Zhang, Joshua Zimmerman, Dan Nettleton, and Dan Nordman.
  # The American Statistician, 74(4):392-406, 2020, specifying the relationship between Y and X may increase the
  # model extrapolation of Random Forests
  
  # 08/09/2022
  # Although adding squared terms and interaction term does not considerably help model fitting,
  # keeping the same model structure as GLM makes some sense.
  formula_rf <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  # Keep the same formula used in developing GLMs for RF
  fitted_rf <- randomForest(formula, data = as.data.frame(simulated_data$train))
  
  # Calculate metrics
  
  # RMSE
  rmse <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted, as.data.frame(X))-as.data.frame(X)$y)^2))))
  rmse_rf <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted_rf, as.data.frame(X))-as.data.frame(X)$y)^2))))
  
  # # Variance
  # variance <- lapply(simulated_data, function(X) (mean((mean(predict(fitted, as.data.frame(X))) - predict(fitted, as.data.frame(X)))^2)))
  
  # Calculate the condition number for all predictors
  cn_all <-  lapply(simulated_data, function(X)(kappa(cor(X[, predictors_index]), exact = TRUE)))
  
  # Retrieve R-square and standard errors of coefficients for X2 that induced collinearity shift
  model_r_square <- summary(fitted)$r.square
  std_error_X2 <- summary(fitted)$coefficients[2, 2]
  
  rf_r_square <- fitted_rf$rsq[length(fitted_rf$rsq)]
  
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
  
  return(list("RMSE" = rmse, "RMSE_RF" = rmse_rf, "CN_all" = cn_all,
              "R_square" = model_r_square, "R_square_rf" = rf_r_square, "Std_error" = std_error_X2,
              "M_distance_mean" = mal_dist_mean,
              "M_distance_median" = mal_dist_median))
}



# Run simulation for high training collinearity
start_time <- Sys.time()
set.seed(1001)
high_train_coli_interaction <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                              f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2])- 3*X[,3],
                                                              mu = c(0, 0, 0),
                                                              sd = c(1, 1, 1),
                                                              rho = c(0.9, 0.8, 0.7),
                                                              response = "y",
                                                              predictors = c("X1", "X2", "I(X1^2)","X1*X2", "X3"),
                                                              predictors_index = c("X1", "X2", "X1^2","X1*X2", "X3")))

high_inter_rmse <- setDF(rbindlist(high_train_coli_interaction["RMSE", ]))
high_inter_rmse_rf <- setDF(rbindlist(high_train_coli_interaction["RMSE_RF", ]))
# high_inter_cor_pair <- setDF(rbindlist(high_train_coli_interaction["Cor_pair", ]))
# high_inter_cor_uppertri <- setDF(rbindlist(high_train_coli_interaction["Cor_uppertri", ]))
# high_inter_CN_essential <- setDF(rbindlist(high_train_coli_interaction["CN_essential", ]))
# high_inter_CN_non_essential <- setDF(rbindlist(high_train_coli_interaction["CN_non_essential", ]))
high_inter_CN_all <- setDF(rbindlist(high_train_coli_interaction["CN_all", ]))
high_inter_r_square <- unlist(high_train_coli_interaction["R_square", ])
high_inter_r_square_rf <- unlist(high_train_coli_interaction["R_square_rf", ])

high_inter_std_error_X2 <- unlist(high_train_coli_interaction["Std_error", ])

high_inter_M_distance_mean <- setDF(rbindlist(high_train_coli_interaction["M_distance_mean", ]))
high_inter_M_distance_median <- setDF(rbindlist(high_train_coli_interaction["M_distance_median", ]))

high_inter_M_distance_mean_mean <- colMeans(high_inter_M_distance_mean)
high_inter_M_distance_median_mean <- colMeans(high_inter_M_distance_median)


high_inter_rmse_mean <- colMeans(high_inter_rmse)
high_inter_rmse_rf_mean <- colMeans(high_inter_rmse_rf)
# high_inter_CN_essential_mean <- colMeans(high_inter_CN_essential)
# high_inter_CN_non_essential_mean <- colMeans(high_inter_CN_non_essential)
high_inter_CN_all_mean <- colMeans(high_inter_CN_all)

plot(seq(0.9, -0.9, by = -0.1), (high_inter_rmse_mean/0.5833850)[3:21])
plot(seq(0.9, -0.9, by = -0.1), (high_inter_rmse_rf_mean/1.549172)[3:21])



plot(c(0, 0.2, 0.4, 0.6, 2, 4, 6), high_inter_rmse_mean[c(3, 28, 47, 66, 85, 104, 123)], col = "red")
points(c(0, 0.2, 0.4, 0.6, 2, 4, 6), high_inter_rmse_mean[c(21, 46, 65, 84, 103, 122, 141)], col = "blue")

plot(c(0, 0.2, 0.4, 0.6), high_inter_rmse_mean[c(3, 28, 47, 66)], col = "red", ylim = c(0, 2))
points(c(0, 0.2, 0.4, 0.6), high_inter_rmse_mean[c(21, 46, 65, 84)], col = "blue")

plot(c(2, 4, 6), high_inter_rmse_mean[c(85, 104, 123)], col = "red")
points(c(2, 4, 6), high_inter_rmse_mean[c(103, 122, 141)], col = "blue")



set.seed(1002)
high_train_coli_quad <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                       f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 3*X[,3],
                                                       mu = c(0, 0, 0),
                                                       sd = c(1, 1, 1),
                                                       rho = c(0.9, 0.8, 0.7),
                                                       response = "y",
                                                       predictors = c("X1", "X2", "I(X1^2)", "X3"),
                                                       predictors_index = c("X1", "X2", "X1^2", "X3")))

high_quad_rmse <- setDF(rbindlist(high_train_coli_quad["RMSE", ]))
high_quad_rmse_rf <- setDF(rbindlist(high_train_coli_quad["RMSE_RF", ]))
# high_quad_cor_pair <- setDF(rbindlist(high_train_coli_quad["Cor_pair", ]))
# high_quad_cor_uppertri <- setDF(rbindlist(high_train_coli_quad["Cor_uppertri", ]))
# high_quad_CN_essential <- setDF(rbindlist(high_train_coli_quad["CN_essential", ]))
# high_quad_CN_non_essential <- setDF(rbindlist(high_train_coli_quad["CN_non_essential", ]))
high_quad_CN_all <- setDF(rbindlist(high_train_coli_quad["CN_all", ]))
high_quad_r_square <- unlist(high_train_coli_quad["R_square", ])
high_quad_r_square_rf <- unlist(high_train_coli_quad["R_square_rf", ])

high_quad_std_error_X2 <- unlist(high_train_coli_quad["Std_error", ])

high_quad_M_distance_mean <- setDF(rbindlist(high_train_coli_quad["M_distance_mean", ]))
high_quad_M_distance_median <- setDF(rbindlist(high_train_coli_quad["M_distance_median", ]))

high_quad_M_distance_mean_mean <- colMeans(high_quad_M_distance_mean)
high_quad_M_distance_median_mean <- colMeans(high_quad_M_distance_median)

high_quad_rmse_mean <- colMeans(high_quad_rmse)
high_quad_rmse_rf_mean <- colMeans(high_quad_rmse_rf)
# high_quad_CN_essential_mean <- colMeans(high_quad_CN_essential)
# high_quad_CN_non_essential_mean <- colMeans(high_quad_CN_non_essential)
high_quad_CN_all_mean <- colMeans(high_quad_CN_all)

plot(seq(0.9, -0.9, by = -0.1), high_quad_rmse_mean[3:21])
plot(seq(0.9, -0.9, by = -0.1), high_quad_rmse_mean[28:46])

plot(seq(0.9, -0.9, by = -0.1), high_quad_rmse_rf_mean[3:21])

set.seed(1003)
high_train_coli_linear <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                         f = function(X) 25 + 2*X[,1] + 1.5*X[,2]- 3*X[,3],
                                                         mu = c(0, 0, 0),
                                                         sd = c(1, 1, 1),
                                                         rho = c(0.9, 0.8, 0.7),
                                                         response = "y",
                                                         predictors = c("X1", "X2", "X3"),
                                                         predictors_index = c("X1", "X2", "X3")))

high_linear_rmse <- setDF(rbindlist(high_train_coli_linear["RMSE", ]))
high_linear_rmse_rf <- setDF(rbindlist(high_train_coli_linear["RMSE_RF", ]))
# high_linear_cor_pair <- setDF(rbindlist(high_train_coli_linear["Cor_pair", ]))
# high_linear_cor_uppertri <- setDF(rbindlist(high_train_coli_linear["Cor_uppertri", ]))
# high_linear_CN_essential <- setDF(rbindlist(high_train_coli_linear["CN_essential", ]))
# high_linear_CN_non_essential <- setDF(rbindlist(high_train_coli_linear["CN_non_essential", ]))
high_linear_CN_all <- setDF(rbindlist(high_train_coli_linear["CN_all", ]))
high_linear_r_square <- unlist(high_train_coli_linear["R_square", ])
high_linear_r_square_rf <- unlist(high_train_coli_linear["R_square_rf", ])

high_linear_std_error_X2 <- unlist(high_train_coli_linear["Std_error", ])

high_linear_M_distance_mean <- setDF(rbindlist(high_train_coli_linear["M_distance_mean", ]))
high_linear_M_distance_median <- setDF(rbindlist(high_train_coli_linear["M_distance_median", ]))

high_linear_M_distance_mean_mean <- colMeans(high_linear_M_distance_mean)
high_linear_M_distance_median_mean <- colMeans(high_linear_M_distance_median)

high_linear_rmse_mean <- colMeans(high_linear_rmse)
high_linear_rmse_rf_mean <- colMeans(high_linear_rmse_rf)
# high_linear_CN_essential_mean <- colMeans(high_linear_CN_essential)
# high_linear_CN_non_essential_mean <- colMeans(high_linear_CN_non_essential)
high_linear_CN_all_mean <- colMeans(high_linear_CN_all)


# Run simulation for mid training collinearity
set.seed(1004)
mid_train_coli_interaction <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                             f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2])- 3*X[,3],
                                                             mu = c(0, 0, 0),
                                                             sd = c(1, 1, 1),
                                                             rho = c(0.6, 0.5, 0.4),
                                                             response = "y",
                                                             predictors = c("X1", "X2", "I(X1^2)","X1*X2", "X3"),
                                                             predictors_index = c("X1", "X2", "X1^2","X1*X2", "X3")))

mid_inter_rmse <- setDF(rbindlist(mid_train_coli_interaction["RMSE", ]))
mid_inter_rmse_rf <- setDF(rbindlist(mid_train_coli_interaction["RMSE_RF", ]))
# mid_inter_cor_pair <- setDF(rbindlist(mid_train_coli_interaction["Cor_pair", ]))
# mid_inter_cor_uppertri <- setDF(rbindlist(mid_train_coli_interaction["Cor_uppertri", ]))
# mid_inter_CN_essential <- setDF(rbindlist(mid_train_coli_interaction["CN_essential", ]))
# mid_inter_CN_non_essential <- setDF(rbindlist(mid_train_coli_interaction["CN_non_essential", ]))
mid_inter_CN_all <- setDF(rbindlist(mid_train_coli_interaction["CN_all", ]))
mid_inter_r_square <- unlist(mid_train_coli_interaction["R_square", ])
mid_inter_r_square_rf <- unlist(mid_train_coli_interaction["R_square_rf", ])

mid_inter_std_error_X2 <- unlist(mid_train_coli_interaction["Std_error", ])

mid_inter_M_distance_mean <- setDF(rbindlist(mid_train_coli_interaction["M_distance_mean", ]))
mid_inter_M_distance_median <- setDF(rbindlist(mid_train_coli_interaction["M_distance_median", ]))

mid_inter_M_distance_mean_mean <- colMeans(mid_inter_M_distance_mean)
mid_inter_M_distance_median_mean <- colMeans(mid_inter_M_distance_median)

mid_inter_rmse_mean <- colMeans(mid_inter_rmse)
mid_inter_rmse_rf_mean <- colMeans(mid_inter_rmse_rf)
# mid_inter_CN_essential_mean <- colMeans(mid_inter_CN_essential)
# mid_inter_CN_non_essential_mean <- colMeans(mid_inter_CN_non_essential)
mid_inter_CN_all_mean <- colMeans(mid_inter_CN_all)

plot(seq(0.9, -0.9, by = -0.1), mid_inter_rmse_mean[3:21])
plot(seq(0.9, -0.9, by = -0.1), mid_inter_rmse_rf_mean[3:21])

set.seed(1005)
mid_train_coli_quad <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                      f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 3*X[,3],
                                                      mu = c(0, 0, 0),
                                                      sd = c(1, 1, 1),
                                                      rho = c(0.6, 0.5, 0.4),
                                                      response = "y",
                                                      predictors = c("X1", "X2", "I(X1^2)", "X3"),
                                                      predictors_index = c("X1", "X2", "X1^2", "X3")))

mid_quad_rmse <- setDF(rbindlist(mid_train_coli_quad["RMSE", ]))
mid_quad_rmse_rf <- setDF(rbindlist(mid_train_coli_quad["RMSE_RF", ]))
# mid_quad_cor_pair <- setDF(rbindlist(mid_train_coli_quad["Cor_pair", ]))
# mid_quad_cor_uppertri <- setDF(rbindlist(mid_train_coli_quad["Cor_uppertri", ]))
# mid_quad_CN_essential <- setDF(rbindlist(mid_train_coli_quad["CN_essential", ]))
# mid_quad_CN_non_essential <- setDF(rbindlist(mid_train_coli_quad["CN_non_essential", ]))
mid_quad_CN_all <- setDF(rbindlist(mid_train_coli_quad["CN_all", ]))
mid_quad_r_square <- unlist(mid_train_coli_quad["R_square", ])
mid_quad_r_square_rf <- unlist(mid_train_coli_quad["R_square_rf", ])

mid_quad_std_error_X2 <- unlist(mid_train_coli_quad["Std_error", ])

mid_quad_M_distance_mean <- setDF(rbindlist(mid_train_coli_quad["M_distance_mean", ]))
mid_quad_M_distance_median <- setDF(rbindlist(mid_train_coli_quad["M_distance_median", ]))

mid_quad_M_distance_mean_mean <- colMeans(mid_quad_M_distance_mean)
mid_quad_M_distance_median_mean <- colMeans(mid_quad_M_distance_median)

mid_quad_rmse_mean <- colMeans(mid_quad_rmse)
mid_quad_rmse_rf_mean <- colMeans(mid_quad_rmse_rf)
# mid_quad_CN_essential_mean <- colMeans(mid_quad_CN_essential)
# mid_quad_CN_non_essential_mean <- colMeans(mid_quad_CN_non_essential)
mid_quad_CN_all_mean <- colMeans(mid_quad_CN_all)

set.seed(1006)
mid_train_coli_linear <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                        f = function(X) 25 + 2*X[,1] + 1.5*X[,2]- 3*X[,3],
                                                        mu = c(0, 0, 0),
                                                        sd = c(1, 1, 1),
                                                        rho = c(0.6, 0.5, 0.4),
                                                        response = "y",
                                                        predictors = c("X1", "X2", "X3"),
                                                        predictors_index = c("X1", "X2", "X3")))

mid_linear_rmse <- setDF(rbindlist(mid_train_coli_linear["RMSE", ]))
mid_linear_rmse_rf <- setDF(rbindlist(mid_train_coli_linear["RMSE_RF", ]))
# mid_linear_cor_pair <- setDF(rbindlist(mid_train_coli_linear["Cor_pair", ]))
# mid_linear_cor_uppertri <- setDF(rbindlist(mid_train_coli_linear["Cor_uppertri", ]))
# mid_linear_CN_essential <- setDF(rbindlist(mid_train_coli_linear["CN_essential", ]))
# mid_linear_CN_non_essential <- setDF(rbindlist(mid_train_coli_linear["CN_non_essential", ]))
mid_linear_CN_all <- setDF(rbindlist(mid_train_coli_linear["CN_all", ]))
mid_linear_r_square <- unlist(mid_train_coli_linear["R_square", ])
mid_linear_r_square_rf <- unlist(mid_train_coli_linear["R_square_rf", ])

mid_linear_std_error_X2 <- unlist(mid_train_coli_linear["Std_error", ])

mid_linear_M_distance_mean <- setDF(rbindlist(mid_train_coli_linear["M_distance_mean", ]))
mid_linear_M_distance_median <- setDF(rbindlist(mid_train_coli_linear["M_distance_median", ]))

mid_linear_M_distance_mean_mean <- colMeans(mid_linear_M_distance_mean)
mid_linear_M_distance_median_mean <- colMeans(mid_linear_M_distance_median)

mid_linear_rmse_mean <- colMeans(mid_linear_rmse)
mid_linear_rmse_rf_mean <- colMeans(mid_linear_rmse_rf)
# mid_linear_CN_essential_mean <- colMeans(mid_linear_CN_essential)
# mid_linear_CN_non_essential_mean <- colMeans(mid_linear_CN_non_essential)
mid_linear_CN_all_mean <- colMeans(mid_linear_CN_all)


# Run simulation for low training collinearity
set.seed(1007)
low_train_coli_interaction <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                             f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] + 1.5*(X[,1])*(X[,2])- 3*X[,3],
                                                             mu = c(0, 0, 0),
                                                             sd = c(1, 1, 1),
                                                             rho = c(0.3, 0.2, 0.1),
                                                             response = "y",
                                                             predictors = c("X1", "X2", "I(X1^2)","X1*X2", "X3"),
                                                             predictors_index = c("X1", "X2", "X1^2","X1*X2", "X3")))

low_inter_rmse <- setDF(rbindlist(low_train_coli_interaction["RMSE", ]))
low_inter_rmse_rf <- setDF(rbindlist(low_train_coli_interaction["RMSE_RF", ]))
# low_inter_cor_pair <- setDF(rbindlist(low_train_coli_interaction["Cor_pair", ]))
# low_inter_cor_uppertri <- setDF(rbindlist(low_train_coli_interaction["Cor_uppertri", ]))
# low_inter_CN_essential <- setDF(rbindlist(low_train_coli_interaction["CN_essential", ]))
# low_inter_CN_non_essential <- setDF(rbindlist(low_train_coli_interaction["CN_non_essential", ]))
low_inter_CN_all <- setDF(rbindlist(low_train_coli_interaction["CN_all", ]))
low_inter_r_square <- unlist(low_train_coli_interaction["R_square", ])
low_inter_r_square_rf <- unlist(low_train_coli_interaction["R_square_rf", ])

low_inter_std_error_X2 <- unlist(low_train_coli_interaction["Std_error", ])

low_inter_M_distance_mean <- setDF(rbindlist(low_train_coli_interaction["M_distance_mean", ]))
low_inter_M_distance_median <- setDF(rbindlist(low_train_coli_interaction["M_distance_median", ]))

low_inter_M_distance_mean_mean <- colMeans(low_inter_M_distance_mean)
low_inter_M_distance_median_mean <- colMeans(low_inter_M_distance_median)

low_inter_rmse_mean <- colMeans(low_inter_rmse)
low_inter_rmse_rf_mean <- colMeans(low_inter_rmse_rf)
# low_inter_CN_essential_mean <- colMeans(low_inter_CN_essential)
# low_inter_CN_non_essential_mean <- colMeans(low_inter_CN_non_essential)
low_inter_CN_all_mean <- colMeans(low_inter_CN_all)


plot(seq(0.9, -0.9, by = -0.1), low_inter_rmse_mean[3:21])
plot(seq(0.9, -0.9, by = -0.1), low_inter_rmse_rf_mean[3:21])

set.seed(1008)
low_train_coli_quad <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                      f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 3*X[,3],
                                                      mu = c(0, 0, 0),
                                                      sd = c(1, 1, 1),
                                                      rho = c(0.3, 0.2, 0.1),
                                                      response = "y",
                                                      predictors = c("X1", "X2", "I(X1^2)", "X3"),
                                                      predictors_index = c("X1", "X2", "X1^2", "X3")))

low_quad_rmse <- setDF(rbindlist(low_train_coli_quad["RMSE", ]))
low_quad_rmse_rf <- setDF(rbindlist(low_train_coli_quad["RMSE_RF", ]))
# low_quad_cor_pair <- setDF(rbindlist(low_train_coli_quad["Cor_pair", ]))
# low_quad_cor_uppertri <- setDF(rbindlist(low_train_coli_quad["Cor_uppertri", ]))
# low_quad_CN_essential <- setDF(rbindlist(low_train_coli_quad["CN_essential", ]))
# low_quad_CN_non_essential <- setDF(rbindlist(low_train_coli_quad["CN_non_essential", ]))
low_quad_CN_all <- setDF(rbindlist(low_train_coli_quad["CN_all", ]))
low_quad_r_square <- unlist(low_train_coli_quad["R_square", ])
low_quad_r_square_rf <- unlist(low_train_coli_quad["R_square_rf", ])

low_quad_std_error_X2 <- unlist(low_train_coli_quad["Std_error", ])

low_quad_M_distance_mean <- setDF(rbindlist(low_train_coli_quad["M_distance_mean", ]))
low_quad_M_distance_median <- setDF(rbindlist(low_train_coli_quad["M_distance_median", ]))

low_quad_M_distance_mean_mean <- colMeans(low_quad_M_distance_mean)
low_quad_M_distance_median_mean <- colMeans(low_quad_M_distance_median)

low_quad_rmse_mean <- colMeans(low_quad_rmse)
low_quad_rmse_rf_mean <- colMeans(low_quad_rmse_rf)
# low_quad_CN_essential_mean <- colMeans(low_quad_CN_essential)
# low_quad_CN_non_essential_mean <- colMeans(low_quad_CN_non_essential)
low_quad_CN_all_mean <- colMeans(low_quad_CN_all)

set.seed(1009)
low_train_coli_linear <- replicate(1000, linear_model_2(vars = c("X1", "X2", "X3"),
                                                        f = function(X) 25 + 2*X[,1] + 1.5*X[,2]- 3*X[,3],
                                                        mu = c(0, 0, 0),
                                                        sd = c(1, 1, 1),
                                                        rho = c(0.3, 0.2, 0.1),
                                                        response = "y",
                                                        predictors = c("X1", "X2", "X3"),
                                                        predictors_index = c("X1", "X2", "X3")))

low_linear_rmse <- setDF(rbindlist(low_train_coli_linear["RMSE", ]))
low_linear_rmse_rf <- setDF(rbindlist(low_train_coli_linear["RMSE_RF", ]))
# low_linear_cor_pair <- setDF(rbindlist(low_train_coli_linear["Cor_pair", ]))
# low_linear_cor_uppertri <- setDF(rbindlist(low_train_coli_linear["Cor_uppertri", ]))
# low_linear_CN_essential <- setDF(rbindlist(low_train_coli_linear["CN_essential", ]))
# low_linear_CN_non_essential <- setDF(rbindlist(low_train_coli_linear["CN_non_essential", ]))
low_linear_CN_all <- setDF(rbindlist(low_train_coli_linear["CN_all", ]))
low_linear_r_square <- unlist(low_train_coli_linear["R_square", ])
low_linear_r_square_rf <- unlist(low_train_coli_linear["R_square_rf", ])

low_linear_std_error_X2 <- unlist(low_train_coli_linear["Std_error", ])

low_linear_M_distance_mean <- setDF(rbindlist(low_train_coli_linear["M_distance_mean", ]))
low_linear_M_distance_median <- setDF(rbindlist(low_train_coli_linear["M_distance_median", ]))

low_linear_M_distance_mean_mean <- colMeans(low_linear_M_distance_mean)
low_linear_M_distance_median_mean <- colMeans(low_linear_M_distance_median)


low_linear_rmse_mean <- colMeans(low_linear_rmse)
low_linear_rmse_rf_mean <- colMeans(low_linear_rmse_rf)
# low_linear_CN_essential_mean <- colMeans(low_linear_CN_essential)
# low_linear_CN_non_essential_mean <- colMeans(low_linear_CN_non_essential)
low_linear_CN_all_mean <- colMeans(low_linear_CN_all)

end_time <- Sys.time()

save.image("C:/FSU/Collinearity_extra/colli_mod_w_M_distance.RData")


# Export RMSE for each simulation
write.table(high_linear_rmse_mean, "high_linear_rmse_mean.csv", sep = ",", row.names = T)
write.table(high_quad_rmse_mean, "high_quad_rmse_mean.csv", sep = ",", row.names = T)
write.table(high_inter_rmse_mean, "high_inter_rmse_mean.csv", sep = ",", row.names = T)

write.table(mid_linear_rmse_mean, "mid_linear_rmse_mean.csv", sep = ",", row.names = T)
write.table(mid_quad_rmse_mean, "mid_quad_rmse_mean.csv", sep = ",", row.names = T)
write.table(mid_inter_rmse_mean, "mid_inter_rmse_mean.csv", sep = ",", row.names = T)

write.table(low_linear_rmse_mean, "low_linear_rmse_mean.csv", sep = ",", row.names = T)
write.table(low_quad_rmse_mean, "low_quad_rmse_mean.csv", sep = ",", row.names = T)
write.table(low_inter_rmse_mean, "low_inter_rmse_mean.csv", sep = ",", row.names = T)


#######################################################################################################################################################################
#######################################################################################################################################################################
# ANOVA
# Treat collinearity shift as a factor first
# The number of rows for the data is 2, 376, 000
#######################################################################################################################################################################
#######################################################################################################################################################################

# Small test
# Transpose all collinearity shift
library(reshape2)
a <- high_inter_rmse[1:3, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
a$Collinearity_shift <- c(rep(0.9, 3),
                          rep(0.8, 3),
                          rep(0.7, 3),
                          rep(0.6, 3),
                          rep(0.5, 3),
                          rep(0.4, 3),
                          rep(0.3, 3),
                          rep(0.2, 3),
                          rep(0.1, 3),
                          rep(0, 3),
                          rep(-0.1, 3),
                          rep(-0.2, 3),
                          rep(-0.3, 3),
                          rep(-0.4, 3),
                          rep(-0.5, 3),
                          rep(-0.6, 3),
                          rep(-0.7, 3),
                          rep(-0.8, 3),
                          rep(-0.9, 3))
a$Collinearity_shift <- as.factor(a$Collinearity_shift)
a <- a[, -6]

# Transpose all range shift
b <- high_inter_rmse[1:3, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, 3),
                   rep(0.4, 3),
                   rep(0.6, 3),
                   rep(2, 3),
                   rep(4, 3),
                   rep(6, 3))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]



# Transpose all interactions
c <- high_inter_rmse[1:3, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, 3*19),
                   rep(0.4, 3*19),
                   rep(0.6, 3*19),
                   rep(2, 3*19),
                   rep(4, 3*19),
                   rep(6, 3*19))

c$Collinearity_shift <- rep(c(rep(0.9, 3),
                              rep(0.8, 3),
                              rep(0.7, 3),
                              rep(0.6, 3),
                              rep(0.5, 3),
                              rep(0.4, 3),
                              rep(0.3, 3),
                              rep(0.2, 3),
                              rep(0.1, 3),
                              rep(0, 3),
                              rep(-0.1, 3),
                              rep(-0.2, 3),
                              rep(-0.3, 3),
                              rep(-0.4, 3),
                              rep(-0.5, 3),
                              rep(-0.6, 3),
                              rep(-0.7, 3),
                              rep(-0.8, 3),
                              rep(-0.9, 3)), 6)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)
c <- c[, -6]

d_analysis <- rbind(a, b, c)

# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_inter_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_inter_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_inter_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_name_analysis <- rbind(a, b, c)


# Manipulate high_inter_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_inter_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_inter_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_inter_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(a, b, c)

# Manipulate high_inter_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_inter_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_inter_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_inter_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate high_quad_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_quad_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_quad_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_quad_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate high_quad_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_quad_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
                          rep(0.-0.1, nrow_analysis),
                          rep(0.-0.2, nrow_analysis),
                          rep(0.-0.3, nrow_analysis),
                          rep(0.-0.4, nrow_analysis),
                          rep(0.-0.5, nrow_analysis),
                          rep(0.-0.6, nrow_analysis),
                          rep(0.-0.7, nrow_analysis),
                          rep(0.-0.8, nrow_analysis),
                          rep(0.-0.9, nrow_analysis))
a$Collinearity_shift <- as.factor(a$Collinearity_shift)
a <- a[, -6]

# Transpose all range shift
b <- high_quad_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_quad_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
                              rep(0.-0.1, nrow_analysis),
                              rep(0.-0.2, nrow_analysis),
                              rep(0.-0.3, nrow_analysis),
                              rep(0.-0.4, nrow_analysis),
                              rep(0.-0.5, nrow_analysis),
                              rep(0.-0.6, nrow_analysis),
                              rep(0.-0.7, nrow_analysis),
                              rep(0.-0.8, nrow_analysis),
                              rep(0.-0.9, nrow_analysis)), 6)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate high_linear_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_linear_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_linear_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_linear_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate high_linear_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_linear_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- high_linear_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_linear_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate mid_inter_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_inter_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_inter_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_inter_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate mid_inter_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_inter_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_inter_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_inter_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate mid_quad_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_quad_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_quad_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_quad_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate mid_quad_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_quad_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_quad_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_quad_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate mid_linear_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_linear_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_linear_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_linear_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate mid_linear_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_linear_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_linear_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_linear_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate low_inter_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_inter_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
                          rep(0.-0.1, nrow_analysis),
                          rep(0.-0.2, nrow_analysis),
                          rep(0.-0.3, nrow_analysis),
                          rep(0.-0.4, nrow_analysis),
                          rep(0.-0.5, nrow_analysis),
                          rep(0.-0.6, nrow_analysis),
                          rep(0.-0.7, nrow_analysis),
                          rep(0.-0.8, nrow_analysis),
                          rep(0.-0.9, nrow_analysis))
a$Collinearity_shift <- as.factor(a$Collinearity_shift)
a <- a[, -6]

# Transpose all range shift
b <- low_inter_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_inter_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate low_inter_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_inter_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- low_inter_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_inter_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate low_quad_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_quad_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- low_quad_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_quad_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)


# Manipulate low_quad_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_quad_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- low_quad_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_quad_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate low_linear_rmse
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_linear_rmse[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- low_linear_rmse[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_linear_rmse[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)

# Manipulate low_linear_rmse_rf
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_linear_rmse_rf[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "RMSE"
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
a <- a[, -6]

# Transpose all range shift
b <- low_linear_rmse_rf[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "RMSE"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_linear_rmse_rf[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("RMSE")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis <- rbind(data_analysis, a, b, c)



# Some sample structures of anova

# Loop each range gradient
for (i in c(0, .2, .4, .6, 2, 4 ,6)){
  model.i <- lm(RMSE ~ collinearity_shift, data = data[data$range_shift ==i, ])
  anova(model.i)
}

# Double loop range and model gradients
for (j in c("Interaction", "Quadratic", "Linear")){
  for (i in c(0, .2, .4, .6, 2, 4 ,6)){
    model.j.i <- lm(RMSE ~ collinearity_shift + range_shift, data = data[data$Model = j&data$range_shift ==i, ])
    anova(model.j.i)
  } 
}





# Check all the factors
str(data_analysis)

anova_model <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, data = data_analysis)

anova(anova_model)

anova_model_coli <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Collinearity_shift, data = data_analysis[data_analysis$Range_shift == 0, ])
anova(anova_model_coli)


anova_model_coli_lowrange <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, data = data_analysis[data_analysis$Range_shift %in% c(0, 0.2, 0.4, 0.6), ])
anova(anova_model_coli_lowrange)


levels(data_analysis$Range_shift)


# Anova for the subsets of the data
# For range shift (.2, .4, .6) and (0, 2, 4, 6)
# w/o interaction
anova_model_rangeshift_1 <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift, 
                               data = data_analysis[data_analysis$Range_shift %in% c(0.2, 0.4, 0.6), ])
anova(anova_model_rangeshift_1)

anova_model_rangeshift_2 <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift, 
                               data = data_analysis[!(data_analysis$Range_shift %in% c(0.2, 0.4, 0.6)), ]) 
anova(anova_model_rangeshift_2)

# w/ interaction
anova_model_rangeshift_3 <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                               data = data_analysis[data_analysis$Range_shift %in% c(0.2, 0.4, 0.6), ])
anova(anova_model_rangeshift_3)

anova_model_rangeshift_4 <- lm(RMSE ~ Model + Algorithm + Predictor_coli + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                               data = data_analysis[!(data_analysis$Range_shift %in% c(0.2, 0.4, 0.6)), ]) 
anova(anova_model_rangeshift_4)



# Run three GLMs for different training collinearity
# Only consider interaction model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_model_training_1 <- lm(RMSE ~ Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                                    data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                    data_analysis$Algorithm == "GLM"&
                                                    data_analysis$Model == "Interaction", ])
anova(anova_model_training_1)


ggplot(data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == .2&
                              data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                              data_analysis$Algorithm == "GLM"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()

library(PMCMRplus)
Sys.time()
set.seed(1234)
jonckheereTest(RMSE ~ Collinearity_shift, 
               data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                      data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                      data_analysis$Algorithm == "GLM"&
                                      data_analysis$Model == "Interaction", ], alternative = "less")
Sys.time()
library(DescTools)
a <- data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                     data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                     data_analysis$Algorithm == "GLM"&
                     data_analysis$Model == "Interaction", ]

a$Collinearity_shift <- factor(a$Collinearity_shift, ordered = T)

Sys.time()
set.seed(1234)
JonckheereTerpstraTest(a$RMSE, a$Collinearity_shift,
                       nperm = 1000)
Sys.time()


# 4/7/2022 Built a function and simulate the high training collinearity and interaction terms to see how
# GLM generates skewed distributions of RMSE for collinearity shift = -0.9 using linear_model_1 function
Sys.time()
set.seed(1001)
high_train_coli_interaction_tracker <- replicate(1000, linear_model_1(vars = c("X1", "X2", "X3"),
                                                                      f = function(X) 25 + 2*X[,1] + 1.5*X[,2]- 3*X[,3],
                                                                      mu = c(0, 0, 0),
                                                                      sd = c(1, 1, 1),
                                                                      rho = c(0.9, 0.8, 0.7),
                                                                      response = "y",
                                                                      predictors = c("X1", "X2", "X3")))

Sys.time()

high_inter_tracker_rmse <- setDF(rbindlist(high_train_coli_interaction_tracker["RMSE", ]))
high_inter_tracker_CN_all <- setDF(rbindlist(high_train_coli_interaction_tracker["CN_all", ]))
high_inter_tracker_r_square <- unlist(high_train_coli_interaction_tracker["R_square", ])


high_inter_tracker_intercept <- unlist(high_train_coli_interaction_tracker["Coef_intercept", ])
high_inter_tracker_coef_X1 <- unlist(high_train_coli_interaction_tracker["Coef_X1", ])
high_inter_tracker_coef_X2 <- unlist(high_train_coli_interaction_tracker["Coef_X2", ])
high_inter_tracker_coef_X1_square <- unlist(high_train_coli_interaction_tracker["Coef_X1_square", ])
high_inter_tracker_coef_X3 <- unlist(high_train_coli_interaction_tracker["Coef_X3", ])
high_inter_tracker_coef_X1X2 <- unlist(high_train_coli_interaction_tracker["Coef_X1X2", ])

Sys.time()

# Find the worst model with the highest RMSE
# Compare the lowest RMSE in col shift = 0 with col shift = -0.9
(high_inter_tracker_rmse$test_same)[which.min(high_inter_tracker_rmse$test_same)]
(high_inter_tracker_rmse$`r=-0.9`)[which.min(high_inter_tracker_rmse$`r=-0.9`)]

# Compare the highest RMSE in col shift = 0 with col shift = -0.9
(high_inter_tracker_rmse$test_same)[which.max(high_inter_tracker_rmse$test_same)]
(high_inter_tracker_rmse$`r=-0.9`)[which.max(high_inter_tracker_rmse$`r=-0.9`)]

# The worst and the best
which.max(high_inter_tracker_rmse$`r=-0.9`)
which.min(high_inter_tracker_rmse$`r=-0.9`)

hist(high_inter_tracker_rmse$test_same)
hist(high_inter_tracker_rmse$`r=-0.9`)
(high_inter_tracker_rmse$`r=-0.9`)[15]

#function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2]- 1.5*(X[,1])*(X[,2])- 3*X[,3]

# Best
# f(best) ~ 25.052 + 1.89X1 - 1.98X1^2 + 1.43X2 - 1.42*X1X2 - 2.98X3
# f(worst) ~ 25.084 + 1.56X1 - 1.56X1^2 + 2.07X2 - 1.56*X1X2 - 2.92X3

# Find the data kept for the model with the highest RMSE
data_highest_RMSE_model <- high_train_coli_interaction_tracker["data_kept", 15] 



plot(data_highest_RMSE_model$data_kept$test_same$X1, data_highest_RMSE_model$data_kept$test_same$y)

reg_train <- lm(y ~ X1 + X2 + X3, data = data_highest_RMSE_model$data_kept$train)
reg_test <- lm(y ~ X1 + X2 + X3, data = data_highest_RMSE_model$data_kept$test_same)
reg_train

# Coefficients:
#   (Intercept)           X1           X2           X3  
# 25.1762       3.3157       0.4269      -3.3127

reg_test
# Coefficients:
#   (Intercept)           X1           X2           X3  
# 25.016        1.988        1.526       -3.023  

library(car)
#produce added variable plots
avPlots(reg_train, xlim = c(-1.5, 1.5), ylim = c(-10, 10))
avPlots(reg_test, xlim = c(-1.5, 1.5), ylim = c(-10, 10))


# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)

anova_model_training_2_left <- lm(RMSE ~ Collinearity_shift, 
                                  data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                         data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                         data_analysis$Algorithm == "GLM"&
                                                         data_analysis$Model == "Interaction", ])

anova(anova_model_training_2_left)
TukeyHSD(aov(anova_model_training_2_left))


# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_2_right <- lm(RMSE ~ Collinearity_shift, 
                                   data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                          data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                          data_analysis$Algorithm == "GLM"&
                                                          data_analysis$Model == "Interaction", ])

anova(anova_model_training_2_right)
TukeyHSD(aov(anova_model_training_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_model_training_3_left <- lm(RMSE ~ Collinearity_shift, 
                                  data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                         data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                         data_analysis$Algorithm == "GLM"&
                                                         data_analysis$Model == "Interaction", ])

anova(anova_model_training_3_left)
TukeyHSD(aov(anova_model_training_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_3_right <- lm(RMSE ~ Collinearity_shift, 
                                   data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                          data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                          data_analysis$Algorithm == "GLM"&
                                                          data_analysis$Model == "Interaction", ])

anova(anova_model_training_3_right)
TukeyHSD(aov(anova_model_training_3_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift

# Run three GLMs for different training collinearity
# Only consider quadratic model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_model_training_quad_1 <- lm(RMSE ~ Collinearity_shift, 
                                  data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                                         data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                         data_analysis$Algorithm == "GLM"&
                                                         data_analysis$Model == "Quadratic", ])
anova(anova_model_training_quad_1)
# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)

anova_model_training_quad_2_left <- lm(RMSE ~ Collinearity_shift, 
                                       data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                              data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                              data_analysis$Algorithm == "GLM"&
                                                              data_analysis$Model == "Quadratic", ])

anova(anova_model_training_quad_2_left)
TukeyHSD(aov(anova_model_training_quad_2_left))
# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_quad_2_right <- lm(RMSE ~ Collinearity_shift, 
                                        data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                               data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                               data_analysis$Algorithm == "GLM"&
                                                               data_analysis$Model == "Quadratic", ])

anova(anova_model_training_quad_2_right)
TukeyHSD(aov(anova_model_training_quad_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_model_training_quad_3_left <- lm(RMSE ~ Collinearity_shift, 
                                       data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                              data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                              data_analysis$Algorithm == "GLM"&
                                                              data_analysis$Model == "Quadratic", ])

anova(anova_model_training_quad_3_left)
TukeyHSD(aov(anova_model_training_quad_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_quad_3_right <- lm(RMSE ~ Collinearity_shift, 
                                        data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                               data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                               data_analysis$Algorithm == "GLM"&
                                                               data_analysis$Model == "Quadratic", ])

anova(anova_model_training_quad_3_right)
TukeyHSD(aov(anova_model_training_quad_3_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift

# Run three GLMs for different training collinearity
# Only consider linear model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_model_training_linear_1 <- lm(RMSE ~ Collinearity_shift, 
                                    data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                                           data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                           data_analysis$Algorithm == "GLM"&
                                                           data_analysis$Model == "Linear", ])
anova(anova_model_training_linear_1)
# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)

anova_model_training_linear_2_left <- lm(RMSE ~ Collinearity_shift, 
                                         data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == .4&
                                                                data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                                data_analysis$Algorithm == "GLM"&
                                                                data_analysis$Model == "Linear", ])

anova(anova_model_training_linear_2_left)
TukeyHSD(aov(anova_model_training_linear_2_left))
# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_linear_2_right <- lm(RMSE ~ Collinearity_shift, 
                                          data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                                 data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                                 data_analysis$Algorithm == "GLM"&
                                                                 data_analysis$Model == "Linear", ])

anova(anova_model_training_linear_2_right)
TukeyHSD(aov(anova_model_training_linear_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_model_training_linear_3_left <- lm(RMSE ~ Collinearity_shift, 
                                         data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                                data_analysis$Algorithm == "GLM"&
                                                                data_analysis$Model == "Linear", ])

anova(anova_model_training_linear_3_left)
TukeyHSD(aov(anova_model_training_linear_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_model_training_linear_3_right <- lm(RMSE ~ Collinearity_shift, 
                                          data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                 data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                                 data_analysis$Algorithm == "GLM"&
                                                                 data_analysis$Model == "Linear", ])

anova(anova_model_training_linear_3_right)
TukeyHSD(aov(anova_model_training_linear_3_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
# When Range shift is equal to 0
# The influence of left-side collinearity shift on RMSE was significant for all levels of training collinearity, all levels of model complexity.
# To note, High training collinearity only had left-side collinearity shift
# The influence of right-side collinearity shift on RMSE was dependent on the training collinearity at which the collinearity started to shift
# For low training collinearity that could shift more than mid training collinearity
# We found both-side collinearity shift was significant
#
# More generally, the training collinearity at which the collinearity started to shift matters because it decides how much the collinearity could shift
# If the magnitude of collinearity shift is enough, we will be able to find the significant difference between shift and no shift.
# Therefore, it will be worth plotting the change in collinearity vs RMSE for all significant comparison?
#
# Repeat all the ANOVAs for RF models
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################

# Run three RFs for different training collinearity
# Only consider interaction model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_RF_model_training_1 <- lm(RMSE ~ Collinearity_shift, 
                                data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == .6&
                                                       data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                       data_analysis$Algorithm == "RF"&
                                                       data_analysis$Model == "Interaction", ])

anova(anova_RF_model_training_1)

ggplot(data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                              data_analysis$Algorithm == "RF"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()


ggplot(data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                              data_analysis$Algorithm == "GLM"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()



# Try to plot RMSEs of RF and GLM in the same plot with different color
data_analysis_high_interaction_RF_GLM <- data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift %in% c(0)&
                                                         data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                         (data_analysis$Algorithm == "GLM")&
                                                         data_analysis$Model == "Interaction", ]

ggplot(data = data_analysis_high_interaction_RF_GLM, aes(Collinearity_shift, log(RMSE), color = Range_shift)) + geom_boxplot()
ggplot(data = data_analysis_high_interaction_RF_GLM, aes(Range_shift, log(RMSE), fill = Collinearity_shift)) + geom_boxplot(colour = "black", position=position_dodge(1)) + scale_fill_manual(values = getPalette(colourCount))

# Color palette interpolation
colourCount <- length(unique(data_analysis_high_interaction_RF_GLM$Collinearity_shift))
getPalette <- colorRampPalette(brewer.pal(11, "RdYlBu"))








hist(data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                     data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                     data_analysis$Algorithm == "GLM"&
                     data_analysis$Model == "Linear", "RMSE"])


ggplot(data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(-0.9)&
                              data_analysis$Algorithm == "GLM"&
                              data_analysis$Model == "Linear", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()








# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)


anova_RF_model_training_2_left <- lm(RMSE ~ Collinearity_shift, 
                                     data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                            data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                            data_analysis$Algorithm == "RF"&
                                                            data_analysis$Model == "Interaction", ])

ggplot(data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(seq(-.9, .9, by = .1))&
                              data_analysis$Algorithm == "RF"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()


anova(anova_RF_model_training_2_left)
TukeyHSD(aov(anova_RF_model_training_2_left))
# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_2_right <- lm(RMSE ~ Collinearity_shift, 
                                      data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                             data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                             data_analysis$Algorithm == "RF"&
                                                             data_analysis$Model == "Interaction", ])

anova(anova_RF_model_training_2_right)
TukeyHSD(aov(anova_RF_model_training_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_RF_model_training_3_left <- lm(RMSE ~ Collinearity_shift, 
                                     data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                            data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                            data_analysis$Algorithm == "RF"&
                                                            data_analysis$Model == "Interaction", ])


ggplot(data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                              data_analysis$Algorithm == "RF"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()

anova(anova_RF_model_training_3_left)
TukeyHSD(aov(anova_RF_model_training_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_3_right <- lm(RMSE ~ Collinearity_shift, 
                                      data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                             data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                             data_analysis$Algorithm == "RF"&
                                                             data_analysis$Model == "Interaction", ])

ggplot(data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                              data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                              data_analysis$Algorithm == "RF"&
                              data_analysis$Model == "Interaction", ], aes(Collinearity_shift, RMSE)) + geom_boxplot()

anova(anova_RF_model_training_3_right)
TukeyHSD(aov(anova_RF_model_training_3_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift

# Run three RFs for different training collinearity
# Only consider quadratic model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_RF_model_training_quad_1 <- lm(RMSE ~ Collinearity_shift, 
                                     data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                                            data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                            data_analysis$Algorithm == "RF"&
                                                            data_analysis$Model == "Quadratic", ])
anova(anova_RF_model_training_quad_1)
# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)

anova_RF_model_training_quad_2_left <- lm(RMSE ~ Collinearity_shift, 
                                          data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                                 data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                                 data_analysis$Algorithm == "RF"&
                                                                 data_analysis$Model == "Quadratic", ])

anova(anova_RF_model_training_quad_2_left)
TukeyHSD(aov(anova_RF_model_training_quad_2_left))
# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_quad_2_right <- lm(RMSE ~ Collinearity_shift, 
                                           data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                                  data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                                  data_analysis$Algorithm == "RF"&
                                                                  data_analysis$Model == "Quadratic", ])

anova(anova_RF_model_training_quad_2_right)
TukeyHSD(aov(anova_RF_model_training_quad_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_RF_model_training_quad_3_left <- lm(RMSE ~ Collinearity_shift, 
                                          data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                 data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                                 data_analysis$Algorithm == "RF"&
                                                                 data_analysis$Model == "Quadratic", ])

anova(anova_RF_model_training_quad_3_left)
TukeyHSD(aov(anova_RF_model_training_quad_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_quad_3_right <- lm(RMSE ~ Collinearity_shift, 
                                           data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                  data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                                  data_analysis$Algorithm == "RF"&
                                                                  data_analysis$Model == "Quadratic", ])

anova(anova_RF_model_training_quad_3_right)
TukeyHSD(aov(anova_RF_model_training_quad_3_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift

# Run three RFs for different training collinearity
# Only consider linear model complexity
# High training collinearity and all range shift (0, .2, .4, .6, 2, 4, 6)
# This will tell us which training collinearity has significant collinearity shift across all range shift gradients

anova_RF_model_training_linear_1 <- lm(RMSE ~ Collinearity_shift, 
                                       data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift == 0&
                                                              data_analysis$Collinearity_shift %in% c(seq(-.9, .8, by = .1), 0.9)&
                                                              data_analysis$Algorithm == "RF"&
                                                              data_analysis$Model == "Linear", ])
anova(anova_RF_model_training_linear_1)
# Significant left side for all range gradients (0, .2, .4, .6, 2, 4, 6)

anova_RF_model_training_linear_2_left <- lm(RMSE ~ Collinearity_shift, 
                                            data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == .4&
                                                                   data_analysis$Collinearity_shift %in% c(seq(-.9, .5, by = .1), 0.6)&
                                                                   data_analysis$Algorithm == "RF"&
                                                                   data_analysis$Model == "Linear", ])

anova(anova_RF_model_training_linear_2_left)
TukeyHSD(aov(anova_RF_model_training_linear_2_left))
# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_linear_2_right <- lm(RMSE ~ Collinearity_shift, 
                                             data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift == 0&
                                                                    data_analysis$Collinearity_shift %in% c(seq(.7, .9, by = .1), 0.6)&
                                                                    data_analysis$Algorithm == "RF"&
                                                                    data_analysis$Model == "Linear", ])

anova(anova_RF_model_training_linear_2_right)
TukeyHSD(aov(anova_RF_model_training_linear_2_right))

# Significant right side for range gradients (0)
# RMSE statistically decreased from no shift


anova_RF_model_training_linear_3_left <- lm(RMSE ~ Collinearity_shift, 
                                            data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                   data_analysis$Collinearity_shift %in% c(seq(-.9, .2, by = .1), 0.3)&
                                                                   data_analysis$Algorithm == "RF"&
                                                                   data_analysis$Model == "Linear", ])

anova(anova_RF_model_training_linear_3_left)
TukeyHSD(aov(anova_RF_model_training_linear_3_left))

# Significant left side for range gradients (0, .2, .4,)
# RMSE statistically increased from right to left

anova_RF_model_training_linear_3_right <- lm(RMSE ~ Collinearity_shift, 
                                             data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                                    data_analysis$Collinearity_shift %in% c(seq(.4, .9, by = .1), 0.3)&
                                                                    data_analysis$Algorithm == "RF"&
                                                                    data_analysis$Model == "Linear", ])

anova(anova_RF_model_training_linear_3_right)
TukeyHSD(aov(anova_RF_model_training_linear_3_right))


# Same as the results for GLMs
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################
###################################################################################################################################################################



#######################################################################################################################################################################
#######################################################################################################################################################################
# ANOVA
# Treat collinearity shift as a quantitative predictor
# The number of rows for the data is 2, 376, 000
#######################################################################################################################################################################
#######################################################################################################################################################################

# Small test
# Transpose all collinearity shift
a <- high_inter_CN_all[1:3, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
a$Collinearity_shift <- c(rep(0.9, 3),
                          rep(0.8, 3),
                          rep(0.7, 3),
                          rep(0.6, 3),
                          rep(0.5, 3),
                          rep(0.4, 3),
                          rep(0.3, 3),
                          rep(0.2, 3),
                          rep(0.1, 3),
                          rep(0, 3),
                          rep(0.-0.1, 3),
                          rep(0.-0.2, 3),
                          rep(0.-0.3, 3),
                          rep(0.-0.4, 3),
                          rep(0.-0.5, 3),
                          rep(0.-0.6, 3),
                          rep(0.-0.7, 3),
                          rep(0.-0.8, 3),
                          rep(0.-0.9, 3))
a$Collinearity_shift <- as.factor(a$Collinearity_shift)
a <- a[, -6]

# Transpose all range shift
b <- high_inter_CN_all[1:3, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, 3),
                   rep(0.4, 3),
                   rep(0.6, 3),
                   rep(2, 3),
                   rep(4, 3),
                   rep(6, 3))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]



# Transpose all interactions
c <- high_inter_CN_all[1:3, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, 3*19),
                   rep(0.4, 3*19),
                   rep(0.6, 3*19),
                   rep(2, 3*19),
                   rep(4, 3*19),
                   rep(6, 3*19))

c$Collinearity_shift <- rep(c(rep(0.9, 3),
                              rep(0.8, 3),
                              rep(0.7, 3),
                              rep(0.6, 3),
                              rep(0.5, 3),
                              rep(0.4, 3),
                              rep(0.3, 3),
                              rep(0.2, 3),
                              rep(0.1, 3),
                              rep(0, 3),
                              rep(0.-0.1, 3),
                              rep(0.-0.2, 3),
                              rep(0.-0.3, 3),
                              rep(0.-0.4, 3),
                              rep(0.-0.5, 3),
                              rep(0.-0.6, 3),
                              rep(0.-0.7, 3),
                              rep(0.-0.8, 3),
                              rep(0.-0.9, 3)), 6)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)
c <- c[, -6]

d_analysis <- rbind(a, b, c)

# Manipulate high_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(a, b, c)

# Manipulate high_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate high_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate high_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate high_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate high_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- high_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "High"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- high_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "High"
b$Range_shift <- 0
b$Collinearity_shift <- 0.9 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- high_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "High"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate mid_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate mid_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate mid_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate mid_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate mid_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate mid_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- mid_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "Mid"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- mid_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "Mid"
b$Range_shift <- 0
b$Collinearity_shift <- 0.6 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- mid_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "Mid"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate low_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Interaction"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Interaction"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Interaction"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate low_inter_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_inter_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Interaction"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_inter_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Interaction"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_inter_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Interaction"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate low_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Quadratic"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Quadratic"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Quadratic"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


# Manipulate low_quad_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_quad_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Quadratic"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_quad_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Quadratic"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_quad_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Quadratic"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate low_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "GLM"
a$Model <- "Linear"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "GLM"
b$Model <- "Linear"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "GLM"
c$Model <- "Linear"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)

# Manipulate low_linear_CN_all
# Generalize the data manipulation to n cases
# Transpose all collinearity shift
nrow_analysis <- 1000
a <- low_linear_CN_all[1:nrow_analysis, 3:21]

a$Algorithm <- "RF"
a$Model <- "Linear"
a$Predictor_coli <- "Low"
a$Range_shift <- 0
a$Collinearity_shift <- 0

a <- melt(a, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

a$Algorithm <- as.factor(a$Algorithm)
a$Model <- as.factor(a$Model)
a$Predictor_coli <- as.factor(a$Predictor_coli)
a$Range_shift <- as.factor(a$Range_shift)
a$Collinearity_shift <- as.factor(a$Collinearity_shift)

names(a)[7] <- "CN"
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
a <- a[, -6]

# Transpose all range shift
b <- low_linear_CN_all[1:nrow_analysis, 22:27]
b$Algorithm <- "RF"
b$Model <- "Linear"
b$Predictor_coli <- "Low"
b$Range_shift <- 0
b$Collinearity_shift <- 0.3 # This is the training collinearity

b <- melt(b, id.vars = c("Model", "Algorithm", "Predictor_coli", "Range_shift", "Collinearity_shift"))

b$Algorithm <- as.factor(b$Algorithm)
b$Model <- as.factor(b$Model)
b$Predictor_coli <- as.factor(b$Predictor_coli)
b$Range_shift <- as.factor(b$Range_shift)
b$Collinearity_shift <- as.factor(b$Collinearity_shift)

names(b)[7] <- "CN"

b$Range_shift <- c(rep(0.2, nrow_analysis),
                   rep(0.4, nrow_analysis),
                   rep(0.6, nrow_analysis),
                   rep(2, nrow_analysis),
                   rep(4, nrow_analysis),
                   rep(6, nrow_analysis))
b$Range_shift <- as.factor(b$Range_shift)
b <- b[, -6]

# Transpose all interactions
c <- low_linear_CN_all[1:nrow_analysis, 28:141]
c$Algorithm <- "RF"
c$Model <- "Linear"
c$Predictor_coli <- "Low"
c$Range_shift <- 0
c$Collinearity_shift <- 0

c <- melt(c, id.vars = c("Model", "Algorithm", "Predictor_coli", "Collinearity_shift", "Range_shift"))
c$Algorithm <- as.factor(c$Algorithm)
c$Model <- as.factor(c$Model)
c$Predictor_coli <- as.factor(c$Predictor_coli)
c$Range_shift <- as.factor(c$Range_shift)
c$Collinearity_shift <- as.factor(c$Collinearity_shift)

names(c)[7] <- c("CN")
c$Range_shift <- c(rep(0.2, nrow_analysis*19),
                   rep(0.4, nrow_analysis*19),
                   rep(0.6, nrow_analysis*19),
                   rep(2, nrow_analysis*19),
                   rep(4, nrow_analysis*19),
                   rep(6, nrow_analysis*19))

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
c <- c[, -6]

data_analysis_CN <- rbind(data_analysis_CN, a, b, c)


if(nrow(data_analysis) == nrow(data_analysis_CN)){
  data_analysis_combined <- cbind(data_analysis, as.data.frame(data_analysis_CN[, 6], drop = FALSE))
  names(data_analysis_combined)[7] <- "CN"
}

plot(data_analysis_combined$CN[data_analysis_combined$Collinearity_shift == -0.9&data_analysis_combined$Range_shift ==0], data_analysis_combined$RMSE[data_analysis_combined$Collinearity_shift == -0.9&data_analysis_combined$Range_shift ==0])

# Export data
write.csv(data_analysis_combined, "data_analysis_combined.csv", row.names = FALSE)


#########################################################################################################################################################################











































anova_model_training_3 <- lm(RMSE ~ Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift == 0&
                                                    data_analysis$Collinearity_shift %in% c(seq(-.9, .1, by = .1), 0.3)&
                                                    data_analysis$Algorithm == "GLM"&
                                                    data_analysis$Model == "Interaction", ])
anova(anova_model_training_3)

TukeyHSD(aov(anova_model_training_3))$Collinearity_shift







# High training collinearity and range shift (2, 4, 6)
anova_model_training_2 <- lm(RMSE ~ Model + Algorithm + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "High"&data_analysis$Range_shift %in% c(2, 4, 6), ])
anova(anova_model_training_2)

# Mid training collinearity and range shift (.2, .4, .6)
anova_model_training_3 <- lm(RMSE ~ Model + Algorithm + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift %in% c(0.2, 0.4, 0.6), ]) 
anova(anova_model_training_3)

# Mid training collinearity and range shift (2, 4, 6)
anova_model_training_4 <- lm(RMSE ~ Model + Algorithm + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "Mid"&data_analysis$Range_shift %in% c(2, 4, 6), ]) 
anova(anova_model_training_4)

# Low training collinearity and range shift (.2, .4, .6)
anova_model_training_5 <- lm(RMSE ~ Model + Algorithm + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift %in% c(0.2, 0.4, 0.6), ])
anova(anova_model_training_5)

# Low training collinearity and range shift (2, 4, 6)
anova_model_training_6 <- lm(RMSE ~ Model + Algorithm + Range_shift + Collinearity_shift + Range_shift:Collinearity_shift, 
                             data = data_analysis[data_analysis$Predictor_coli == "Low"&data_analysis$Range_shift %in% c(2, 4, 6), ])
anova(anova_model_training_6)



# Create five lists of loop variables
list_model <- levels(data_analysis$Model)
list_algorithm <- levels(data_analysis$Algorithm)
list_training_collinearity <- levels(data_analysis$Predictor_coli)
list_range_shift <- levels(data_analysis$Range_shift)
list_collinearity_shift <- levels(data_analysis$Collinearity_shift)













# Plot the heatmaps for the interaction
result_heatmap <- read.csv("result_heatmap_6.csv")
result_heatmap$Col <- as.factor(result_heatmap$Col)
result_heatmap$Model <- as.factor(result_heatmap$Model)

library(wesanderson)
library(ggplot2)
library(cowplot)

# Gradient color
pal <- wes_palette("Zissou1", 1000, type = "continuous")
ggplot(heatmap, aes(x = X2, y = X1, fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = pal) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal()

heatmap_high_linear <- ggplot(result_heatmap[result_heatmap$Range == "Range shift = 0.6"&result_heatmap$Model != "linear+quad+inter", ], 
                              mapping = aes(x = Col, y = Model, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of range\nshift") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2))) + 
  scale_y_discrete(breaks = c())+
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(.~Training_col, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "bottom")
heatmap_high_linear

heatmap_high_inter_range_3 <- ggplot(result_heatmap[result_heatmap$Range == "Range shift = 0.6"&result_heatmap$Model == "linear", ], 
                                     mapping = aes(x = Col, y = Model, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of training\ncollinearity") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2))) + 
  scale_y_discrete(breaks = c())+
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(Training_col~Range, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(axis.title=element_text(size = 10),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.key.size = unit(.4, 'cm'))
heatmap_high_inter_range_3

heatmap_high_inter_range_2 <- ggplot(result_heatmap[result_heatmap$Range == "Range shift = 0.4"&result_heatmap$Model == "linear", ], 
                                     mapping = aes(x = Col, y = Model, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of training\ncollinearity") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if bined color palette is preferred)
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2))) + 
  scale_y_discrete(breaks = c())+
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(Training_col~Range, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(axis.title=element_text(size = 10),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.key.size = unit(.4, 'cm'))
heatmap_high_inter_range_2

heatmap_high_inter_range_1 <- ggplot(result_heatmap[result_heatmap$Range == "Range shift = 0.2"&result_heatmap$Model == "linear", ], 
                                     mapping = aes(x = Col, y = Model, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of training\ncollinearity") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_x_discrete(breaks = c(seq(-0.9, 0.9, by = 0.2))) + 
  scale_y_discrete(breaks = c())+
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(Training_col~Range, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(axis.title=element_text(size = 10),
        axis.text.x = element_text(angle = 45),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.key.size = unit(.4, 'cm'))
heatmap_high_inter_range_1


heatmap_high_three <- plot_grid(NULL,
                                heatmap_high_inter_range_3,
                                NULL,
                                heatmap_high_inter_range_2,
                                NULL,
                                heatmap_high_inter_range_1,
                                nrow = 6,
                                labels = c("(a)", "", "(b)", "", "(c)", ""),
                                label_size = 10,
                                rel_heights = c(0.5, 1.8, 0.5, 1.8, 0.5, 1.8),
                                label_x = 0, label_y = 1)

heatmap_high_three


# Plot the histograms for R-square between high and low training collinearity
high_r_square <- as.data.frame(high_inter_r_square)
low_r_square <- as.data.frame(low_inter_r_square)
high_r_square$train_col <- "High training collinearity"
low_r_square$train_col <- "Low training collinearity"
names(high_r_square) <- c("Rsquare" , "train_col")
names(low_r_square) <- c("Rsquare" , "train_col")
hist_r_square <- rbind(high_r_square, low_r_square)

# change fill and outline color manually 
ggplot(hist_r_square, aes(x = Rsquare)) +
  geom_histogram(aes(color = train_col, fill = train_col), 
                 position = "identity", bins = 30, alpha = 0.4) +
  xlab("R-square") +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))

# Plot the histograms for coefficients between high and low training collinearity
high_se_X2 <- as.data.frame(high_inter_std_error_X2)
low_se_X2 <- as.data.frame(high_inter_std_error_X2)
high_se_X2$train_col <- "High training collinearity"
low_se_X2$train_col <- "Low training collinearity"
names(high_r_square) <- c("SE" , "train_col")
names(low_r_square) <- c("SE" , "train_col")
hist_se_X2 <- rbind(high_se_X2, low_se_X2)

# change fill and outline color manually 
ggplot(hist_r_square, aes(x = Rsquare)) +
  geom_histogram(aes(color = train_col, fill = train_col), 
                 position = "identity", bins = 30, alpha = 0.4) +
  xlab("Standard error of X2") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


# Visualize the RMSE vs condition number with error bars
high_inter_rmse_mean <- colMeans(high_inter_rmse)
library(resample)
high_inter_rmse_var <- colStdevs(high_inter_rmse)
high_inter_CN_essential_mean <- colMeans(high_inter_CN_essential)
high_inter_CN_non_essential_mean <- colMeans(high_inter_CN_non_essential)

scatter_rmse <- as.data.frame(high_inter_rmse_mean)
names(scatter_rmse) <- "rmse"
scatter_rmse$sd <- high_inter_rmse_var
scatter_rmse$CN <- colMeans(high_linear_CN_all)
scatter_rmse$rho <- c(rep(.7, 2), seq(.9, -.9 , by = -0.1)[-10], rep(.7, 3), rep(seq(.9, -.9 , by = - 0.1)[-10], 3))

ggplot(scatter_rmse[3:20, ], aes(x=rho, y=rmse)) + 
  geom_point()+
  geom_errorbar(aes(ymin=rmse-sd, ymax=rmse+sd), width=.05,
                position=position_dodge(0.05))


high_inter_rmse_mean <- colMeans(high_inter_rmse)
library(resample)
high_inter_rmse_var <- colStdevs(high_inter_rmse)
high_inter_CN_essential_mean <- colMeans(high_inter_CN_essential)
high_inter_CN_non_essential_mean <- colMeans(high_inter_CN_non_essential)

scatter_rmse <- as.data.frame(high_inter_rmse_mean)
names(scatter_rmse) <- "rmse"
scatter_rmse$sd <- high_inter_rmse_var
scatter_rmse$CN <- colMeans(high_linear_CN_all)
scatter_rmse$rho <- c(rep(.7, 2), seq(.9, -.9 , by = -0.1)[-10], rep(.7, 3), rep(seq(.9, -.9 , by = - 0.1)[-10], 3))

ggplot(scatter_rmse[2:20, ], aes(x=rho, y=rmse)) + 
  geom_point()+
  geom_errorbar(aes(ymin=rmse-sd, ymax=rmse+sd), width=.2,
                position=position_dodge(0.05))

low_inter_rmse_mean <- colMeans(low_inter_rmse)
library(resample)
low_inter_rmse_var <- colStdevs(low_inter_rmse)
low_inter_CN_essential_mean <- colMeans(low_inter_CN_essential)
low_inter_CN_non_essential_mean <- colMeans(low_inter_CN_non_essential)

scatter_rmse <- as.data.frame(low_inter_rmse_mean)
names(scatter_rmse) <- "rmse"
scatter_rmse$sd <- low_inter_rmse_var
scatter_rmse$CN <- colMeans(low_linear_CN_all)
scatter_rmse$rho <- c(rep(.7, 2), seq(.9, -.9 , by = -0.1)[-10], rep(.7, 3), rep(seq(.9, -.9 , by = - 0.1)[-10], 3))

ggplot(scatter_rmse[3:20, ], aes(x=rho, y=rmse)) + 
  geom_point()+
  geom_errorbar(aes(ymin=rmse-sd, ymax=rmse+sd), width=.05,
                position = position_dodge(0.05))


# Explore one RMSE
plot(seq(0.9, -0.9, by = -0.1)[-10], try_rmse_mean[3:20])
# Plot preliminary results to investigate the pattern
rmse_pos <- try_rmse_mean[3:11]
rmse_neg <- try_rmse_mean[12:20]

CN_essential_pos <- try_CN_essential_mean[3:11]
CN_all_pos <- try_CN_all[3:11]

CN_essential_neg <- try_CN_essential_mean[12:20]
CN_all_neg <- try_CN_all[12:20]

plot(x = CN_essential_pos, y = rmse_pos, main = "col = 0.3", xlim = c(-30, 30), ylim = c(3.1, 8.5))
points(x = -CN_essential_neg, y = rmse_neg, col = "red")
abline(v = 1.3)
abline(v = -1.3)

plot(x = log(CN_essential_pos), y = rmse_pos, main = "col = 0.3", xlim = c(-20,20), ylim = c(4, 6.3))
points(x = -log(CN_essential_neg), y = rmse_neg, col = "red")



write.table(try_rmse_mean, "try_rmse_mean.csv", sep = ",", row.names = T)

result_heatmap <- read.csv("result_heatmap_two_decays.csv")
result_heatmap$Col <- as.factor(result_heatmap$Col)

library(wesanderson)
# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(heatmap, aes(x = X2, y = X1, fill = value)) +
  geom_tile() + 
  scale_fill_gradientn(colours = pal) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal() 

# For predictor decay = 0.05, high predictor collinearity
heatmap_result_decay_0.05_positive <- ggplot(result_heatmap[result_heatmap$Model == "B", ], 
                                             mapping = aes(x = Col, y = Range, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of range\nshift") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  # facet_grid(.~Shift, 
  #            scales = "free_x", 
  #            space = "free_x") +
  geom_segment(aes(x = 16, y = .5, xend = 16, yend = 4.5), size = 1, linetype = "dashed") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
heatmap_result_decay_0.05_positive

heatmap_result_decay_0.05_negative <- ggplot(result_heatmap[result_heatmap$Model == "B"&result_heatmap$Shift == "Negative collinearity shift", ], 
                                             mapping = aes(x = rev(Col), y = Range, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity") +
  ylab(label = "Degree of range\nshift") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(.~Shift, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
heatmap_result_decay_0.05_negative

heatmap_two_dacays <- plot_grid(NULL,
                                heatmap_result_decay_0.05_positive,
                                NULL,
                                heatmap_result_decay_0.05_negative,
                                ncol = 4,
                                labels = c("(e)" ,"", "(f)", ""),
                                label_size = 10,
                                rel_widths = c(0.1, 1, 0.1, 1),
                                label_x = 0, label_y = 1)

heatmap_two_dacays


# For predictor decay = 0.5, low predictor collinearity
heatmap_result_decay_0.05_positive <- ggplot(result_heatmap[result_heatmap$Model == "A"&result_heatmap$Shift == "Positive collinearity", ], 
                                             mapping = aes(x = Col, y = Range, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity decay") +
  ylab(label = "Degree of range\nshift") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(.~Shift, 
             scales = "free_x", 
             space = "free_x") +
  geom_segment(aes(x = 8, y = .5, xend = 8, yend = 4.5), size = 1, linetype = "dashed") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
heatmap_result_decay_0.05_positive

heatmap_result_decay_0.05_negative <- ggplot(result_heatmap[result_heatmap$Model == "A"&result_heatmap$Shift == "Negative collinearity", ], 
                                             mapping = aes(x = rev(Col), y = Range, fill = RMSE)) +  
  geom_tile() +
  xlab(label = "Degree of collinearity decay") +
  ylab(label = "Degree of range\nshift") +
  #scale_fill_stepsn(n.breaks = 10, colours = pal) + #(if biined color palettee is preferred)
  scale_fill_gradientn(name = "\u0394RMSE", colors = pal) +
  facet_grid(.~Shift, 
             scales = "free_x", 
             space = "free_x") +
  #geom_vline(xintercept = 4.5, linetype = "dotted", size = 1.2) +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
heatmap_result_decay_0.05_negative

heatmap_two_dacays <- plot_grid(NULL,
                                heatmap_result_decay_0.05_positive,
                                NULL,
                                heatmap_result_decay_0.05_negative,
                                ncol = 4,
                                labels = c("(a)" ,"", "(b)", ""),
                                label_size = 10,
                                rel_widths = c(0.1, 1, 0.1, 1),
                                label_x = 0, label_y = 1
)

heatmap_two_dacays



# Random Forests simulator based on collineariser_4
linear_model_3 <- function(vars = c("X1", "X2", "X3"),
                           f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 + 1.5*X[,2] - 1.5*(X[,1])*(X[,2]) - 3*X[,3], 
                           mu = c(0, 0, 0),
                           sd = c(1, 1, 1),
                           rho = c(0.9, 0.3, 0.2),
                           response = "y",
                           predictors = c("X1", "X2", "I(X1^2)", "X1*X2", "X3")){
  
  s <- sample(1:100000, 1)
  # Create simulated data set
  simulated_data <- collineariser_4(FUN = f,
                                    varnames = vars,
                                    sd.training = sd,
                                    mu.training = mu,
                                    rho.training = rho,
                                    rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                                    range.vec = c(.2, .4, .6),
                                    seed = s,
                                    empirical = FALSE,
                                    y.noise = 0.5)
  
  # Fit a linear model for the simulated data
  formula <- as.formula(
    paste(response,
          paste(predictors, collapse = "+"),
          sep = "~")
  )
  
  fitted <- randomForest(formula, data = as.data.frame(simulated_data$train), mtry=3)
  
  # Calculate RMSE and Variance
  # RMSE
  rmse <- lapply(simulated_data, function(X) (sqrt(mean((predict(fitted, as.data.frame(X))-as.data.frame(X)$y)^2))))
  # Variance
  variance <- lapply(simulated_data, function(X) (mean((mean(predict(fitted, as.data.frame(X))) - predict(fitted, as.data.frame(X)))^2)))
  
  # Calculate the correlation between X1 and X2
  cor_pair <- lapply(simulated_data, function(X)(cor(X[,c(2,4)])))
  
  # Calculate the upper triangle of all predictors
  cor_uppertri <- lapply(simulated_data, function(X)(sum(cor(X[,2:4])[upper.tri(cor(X[,2:4]), diag = FALSE)])))
  
  model <- fitted
  #rmse_variance_list <- list("Rmse" = rmse, "Variance" = variance)
  
  return(list("RMSE" = rmse, "Variance" = variance, "Model" = model, "Cor_pair" = cor_pair, "Cor_uppertri" = cor_uppertri))
}


# Interaction simulation
set.seed(1001)
try <- replicate(10, linear_model_3(vars = c("X1", "X2", "X3"),
                                    f = function(X) 25 + 2*X[,1] - 2*X[,1]^2 - 2*X[,2]^2 + 1.5*X[,2]- 1.5*(X[,1])*(X[,2])- 3*X[,3] + 4*(X[,2])*(X[,3]),
                                    mu = c(0, 0, 0),
                                    sd = c(1, 1, 1),
                                    rho = c(0.7, 0.3, 0.2),
                                    response = "y",
                                    predictors = c("X1", "X2", "I(X1^2)", "I(X2^2)","X3")))

# Check results
try_rmse <- setDF(rbindlist(try["RMSE", ]))
try_cor_pair <- setDF(rbindlist(try["Cor_pair", ]))
try_cor_uppertri <- setDF(rbindlist(try["Cor_uppertri", ]))

try_rmse_mean <- colMeans(try_rmse)
try_cor_pair <- colMeans(try_cor_pair)
try_cor_uppertri<- colMeans(try_cor_uppertri)

try_rmse_mean

# Investigate partial dependence plot
library(ggplot2)
library(pdp)

# Fitted a RF model

library(faux)
simulated_data_rf <- collineariser_4(FUN = function(X) 25 + 2*X[,1] - 2*X[,1]^2 - 2*X[,2]^2 + 1.5*X[,2]- 1.5*(X[,1])*(X[,2])- 3*X[,3] + 4*(X[,2])*(X[,3]),
                                     varnames = c("X1", "X2", "X3"),
                                     sd.training = c(1, 1, 1),
                                     mu.training = c(0, 0, 0),
                                     rho.training = c(0.7, 0.3, 0.2),
                                     rho.test = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                                     range.vec = c(2, 3, 4),
                                     seed = 1001,
                                     empirical = FALSE)

formula_rf <- as.formula(paste("y",
                               paste(c("X1", "X2", "I(X1^2)", "I(X2^2)","X3"), collapse = "+"),
                               sep = "~")
)

fitted_rf <- randomForest(formula, data = as.data.frame(simulated_data_rf$train), mtry=3)

# Default partial dependence plot
fitted_rf_p1 <- partial(fitted_rf, pred.var = "X3", plot = TRUE, rug = TRUE)
plot(simulated_data_rf$train$X3, simulated_data_rf$train$y)

Sys.time()
pd_X1X2 <- partial(fitted_rf, pred.var = c("X1", "X2"))
Sys.time()

pdp2D_X1X2 <- plotPartial(pd_X1X2)

pdp3D_X1X2 <- plotPartial(pd_X1X2, levelplot = FALSE, zlab = "y", colorkey = TRUE, 
                          screen = list(z = -20, x = -60))






# ANOVA

data(ToothGrowth)

head(ToothGrowth, 10)