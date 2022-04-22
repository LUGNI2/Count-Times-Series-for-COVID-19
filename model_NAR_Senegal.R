###############################################################
# NEGATIVE BINOMIAL AR
###############################################################

############################################
## Application donnees covid 19 Senegal ####
############################################

rm(list = ls())
library(tidyverse)
library(MASS)
library(prettyR)
#library(GWRM)
library(caret)
library(gridExtra)
library(cowplot)
library(maxLik)



Sys.setlocale("LC_TIME","English")
base=read.table(file.choose(),sep=";",header=T)
base$date = as.Date(base$date,format="%d/%m/%Y")
describe(base,num.desc=c("mean","median","var","sd","valid.n"),xname=NA,horizontal=FALSE)


testingSet_size = round(0.05*dim(base)[1])
trainingSet_size = dim(base)[1]-testingSet_size
base_train = head(base,trainingSet_size)  # training = 95%
base_test  = tail(base,testingSet_size)  #validation = 5%


v0=1
y0=0
#param[1]     beta_0
#param[2]     beta_1
#param[3]     alpha_1
#param[4]     eta_1
#param[5]     eta_2
#param[6]     eta_3
v = function(param,base,t){
  if(t==1)
    return(param[1]+param[3]*v0+param[2]*log(1+y0)
           +param[4]*base$communautaire[1]+param[5]*base$importe[1]
           +param[6]*base$contact[1])
  else
    return(param[1]+param[3]*v(param,base,t-1)+param[2]*log(1+base$cas[t-1])
           +param[4]*base$communautaire[t]+param[5]*base$importe[t]
           +param[6]*base$contact[t])
}

# u=c(1.0356,0.7081,0.0021,0.0008,-0.0036,0.0042,5.1648)

nu = function(param,base){
  result = c()
  for(i in 1:dim(base)[1])
    result[i] = v(param,base,i)
  return(result)
}

#(beta0,beta1,alpha1,eta1,eta2.eta3,phi)
LogLikelihood <- function(param,base){
  m = exp(nu(param[1:6],base))
  p = param[7]/(m+param[7]) 
  return (sum(dnbinom(base$cas,size=param[7],prob=p,log=TRUE)))
}
#negative log likelihood
negLogLikelihood = function(param,base){
  return(-LogLikelihood(param,base))
}

#fix base à base_train
# log likelihood
LogLik = function(param){
  LogLikelihood(param,base_train)
}
#negative Loglik
negLogLik = function(param){
  -negLogLikelihood(param,base_train)
}


# start values
y = base_train$cas
communautaire = base_train$communautaire
importe = base_train$importe
contact = base_train$contact

df = data.frame(u=y[-1], beta1=log(1+y[-length(y)]),
                eta1=communautaire[-1],eta2=importe[-1],
                eta3=contact[-1])
init=glm.nb(u~beta1+eta1+eta2+eta3, data=df)

i=c(coef(init)[1],coef(init)[2],alpha1=0,coef(init)[3],
    coef(init)[4],coef(init)[5],phi=init$theta)

constrA= matrix(c(rep(0,7),c(-1,0,1,0,-1,1,0),c(0,-1,0,1,-1,1,0),
                  c(0,0,0,0,0,0,0),c(0,0,0,0,0,0,0),c(0,0,0,0,0,0,0)
                  ,c(0,0,0,0,0,0,1)),ncol=7,byrow=F) 

constrB=c(rep(-1,6),0)

#ui %*% theta+ci >= 0
solutionNB = maxLik(LogLik, start=i,
           constraints=list(ineqA=constrA, ineqB=-constrB))
summary(solutionNB)

#Autre méthode
#ui %*% theta - ci >= 0
#solution <- constrOptim(theta=i, f=negLogLik,grad=NULL,
#                        ui=constrA, ci=constrB)

## Predicted Values

# Goodness of Fit 
paramNB =coef(solutionNB)
#paramNB["eta2"]=0
#paramNB["alpha1"]=0
n <- length(base_test$cas)
predictedNB <- exp(nu(paramNB,base_test))
nb <- exp(nu(paramNB,base_train))
sqrt(mean((base_test$cas - predictedNB)^2))
save(paramNB,predictedNB, solutionNB,nb,file = "predict.Rdata")

###### Utilisation du package tscount 
###############################################################
regressors <- cbind(contact,communautaire, importe)
library(tscount)
par_model <- tsglm(ts = y, model = 
                    list(past_obs = 1,past_mean = 1),
                  xreg = regressors,link = "log")
summary(par_model)


nar_model <- tsglm(ts = y, model = list(past_obs = 1,past_mean = 1),
                  xreg = regressors, distr = "nbinom", link = "log")
summary(nar_model)

#comparaison de modèles
rbind(PAR = scoring(par_model),NAR = scoring(nar_model))
predict(par_model,10)
predict(nar_model,10)
#par(mfrow=c(2,2))
marcal(par_model,main = "PAR marginal calibration")
marcal(nar_model,main = "NAR marginal calibration")
pit(par_model,main = "PIT PAR")
pit(nar_model,main = "PIT NAR")

library(png)
library(grid)
marcalP <- readPNG("marcalP.png")
marcalNB <- readPNG("marcalNB.png")
pitPar <- readPNG("pitPar.png")
pitNar <- readPNG("pitNar.png")
grid.arrange(rasterGrob(marcalP),rasterGrob(marcalNB),
             rasterGrob(pitPar),rasterGrob(pitNar), ncol = 2)
