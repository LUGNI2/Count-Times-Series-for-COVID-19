## Application donnees covid 19 Senegal ####
############################################

rm(list = ls())
library(tidyverse)
library(MASS)
library(prettyR)
library(caret)
library(gridExtra)
library(cowplot)
library(maxLik)
library(Metrics)


Sys.setlocale("LC_TIME","English")
base <- read.table(file.choose(),sep = ";",header = T)
#base <- read.csv2("covid_south_africa.csv")
#base$date = as.Date(base$date,format="%d/%m/%Y")
base$date <- as.Date(base$date,format = "%Y-%m-%d")
describe(base,num.desc = c("mean","median","var","sd","valid.n"),xname = NA,horizontal = FALSE)

testingSet_size <- round(0.05*dim(base)[1])
trainingSet_size <- dim(base)[1] - testingSet_size
base_train <- head(base,trainingSet_size)  # training = 95%
base_test  <- tail(base,testingSet_size)  #validation = 5%


#d=param[1]     beta_0
#b=param[2]     beta_1
#a=param[3]     alpha_1
#base
v0 <- 1
y0 <- 0
v <- function(param,base, t){
  if (t == 1)
    return(param[1] + param[3]*v0 + param[2]*log(1 + y0))
  else
    return(param[1] + param[3]*v(param,base,t - 1) + param[2]*log(1 + base$cas[t - 1]))
}

# u=c(0.04586651,0.48555019,0.50466975) 
nu <- function(param,base){
  result <- c()
  for (i in 1:dim(base)[1])
    result[i] <- v(param,base,i)
  return(result)
}


###############################################################
# NEGATIVE BINOMIAL AR
###############################################################
#(d,b,a,phi)
#d=param[1], b=param[2], a=param[3], phi=param[4]
LogLikelihood <- function(param,base){
  m <- exp(nu(param[1:3],base))
  p <- param[4]/(m + param[4]) 
  return(sum(dnbinom(base$cas,size = param[4],prob = p,log = TRUE)))
}

#LogLik fix la base à base_train
LogLik <- function(param){
  return(LogLikelihood(param,base_train))
}

negLogLik <- function(param){
  -LogLik(param)
}

# start values
y <- base_train$cas
df <- data.frame(u = y[-1], beta1 = log(1 + y[-length(y)]))
init <- glm.nb(u~beta1, data = df)
phi <- init$theta
#init=glm(u~beta1, family = quasipoisson(link = "log"), data=df)
#summary(init)
#$phi= 1/69320.78 #valeur déterminée 
i <- c(coef(init)[1],coef(init)[2],alpha1 = 0,phi = phi)
constrA <- matrix(c(rep(0,7),c(-1,0,1,0,-1,1,0),c(0,-1,0,1,-1,1,0)
                  ,c(0,0,0,0,0,0,1)),ncol = 4,byrow = F) 

constrB <- c(rep(-1,6),0)

estimNB <- maxLik(LogLik,start = i,constraints = list(ineqA = constrA, ineqB = -constrB))
summary(estimNB)


#Prediction
# Goodness of Fit 
paramNB <- coef(estimNB)
predictedNB <- exp(nu(paramNB,base_test))
rmse <- sqrt(mean((base_test$cas - predictedNB)^2))
rmse
rmse(base_test$cas, predictedNB)
nb <- exp(nu(paramNB,base_train))
save(paramNB,predictedNB, estimNB,nb, file = "predictSA.Rdata")
AIC <- -2*LogLik(coef(estimNB)) + 8
BIC <- -2*LogLik(coef(estimNB)) + 8*log(length(y))
plot(predictedNB,(base_test$cas - predictedNB)^2)

trend <- function(){
  day <- c()
  alpha <- c()
  beta <- c()
  for (i in seq(30, 590, by = 40)) {
    z <- y[1:i]
    df <- data.frame(u = z[-1], beta1 = log(1 + z[-length(z)]))
    init <- glm.nb(u~beta1, data = df)
    phi <- init$theta
    start <- c(coef(init)[1],coef(init)[2],alpha1 = 0,phi = phi)
    esti <- maxLik(LogLik,start = start, 
            constraints = list(ineqA = constrA, ineqB = -constrB))
    beta[i] <- coef(estimNB)["beta1"]
    alpha[i] <- coef(estimNB)["alpha"]
  }
  
  don <- data.frame(day = seq(30, 590, by = 40), alpha = alpha, beta = beta)
  return(don)
}

don <- trend()

ggplot(data = don) +
  xlab("Day") +
  ylab("Estimated Parameters") +
  geom_line(mapping = aes(x = day,y = alpha, color = expression(alpha)), size = 1) +
  geom_line(mapping = aes(x = day,y = beta, color = expression(beta)), size = 1) +
  scale_color_manual(values = c(
    expression(alpha) = 'blue',
    expression(beta) = 'red',
    )) +
  labs(color = 'Parameter')