rm(list=ls())

#Author: J. N. (Yannick) Wiegers

#This script contains the R and Winbugs code to run the 
#integrated population model for the Mallard as described in Wiegers et al (2022).

##########################################
#Input
##########################################

#Set number of JAGS iterations. Parameters almost converge with ~50k iterations, which takes around 30 minutes to run.
#For the manuscript I ran 200k iterations to ensure total convergence. 

n.iter <- 200000

#Number of simulations for the LSA
LSA.sim <- 10000

##########################################
#Loading data
##########################################

#Load data
Duckling_Project <- read.csv("Duckling_Project.csv")
clutch.data <- read.csv("Clutch size.csv")
NS <- read.csv("Nest success.csv")
EHR.data <- read.csv("Hatching rate.csv")

#Renesting according to Hoekman et al. 2002
r <- 0.26
p <- c(0.968)

for(i in 2:5){
  p[i] <- 1-r*(i-1)
}
p
pN <- c(p[1],p[1]*p[2],p[1]*p[2]*p[3],p[1]*p[2]*p[3]*p[4],p[1]*p[2]*p[3]*p[4]*p[5])

#For Limosa
#pN[3:5] <- 0

data <- list(y=c(100.66, 90.3, 93, 98.91, 94, 94.12,92.44, 87.67, 79.3, 80.04,
                 81.62, 74.76, 76.82, 71.82, 72, 71.92,69.69), #Population index
             T1=18,
             T2=17, 
             Time = as.numeric(scale(1:18)), 
             #For adult survival
             m=structure(c(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 0, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 
                           0, 0, 2, 1, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 
                           2, 12, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 12, 0, 
                           3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 3, 5, 0, 0, 
                           0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 2, 2, 1, 5, 0, 0, 0, 0, 0, 
                           0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 12, 0, 0, 0, 0, 0, 1, 0, 
                           0, 0, 0, 1, 2, 0, 0, 0, 1, 6, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
                           1, 0, 2, 1, 1, 6, 2, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 
                           2, 2, 2, 1, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 7, 0, 2, 
                           6, 14, 37, 16, 9, 62, 27, 138, 29, 263, 101, 189, 153, 173, 250, 
                           167, 198, 311, 289), .Dim = 17:18), #The dead recovery m-array
             rel= c(45, 25, 13, 63, 29, 144, 36, 306, 109, 202, 162, 186, 283,175, 208, 326, 303), #No. released birds per year
             # for Duckling survival
             count = Duckling_Project$Count1, nextCount = Duckling_Project$Count2, interval = Duckling_Project$Interval, 
             age = Duckling_Project$Age, d.year = Duckling_Project$Year, 
             u.days = (c(1:56)-mean(Duckling_Project$Mean_age))/sd(Duckling_Project$Mean_age),
             # for Productivity
             pN=pN[1:4], nd=NS$Exposure.Days,nest_deaths = NS$Failed.nests, #For nest success
             c=clutch.data$C, c.year=clutch.data$Year, #For clutch size
             laid = EHR.data$C, hatched=EHR.data$Hatched) #For egg hatch rate

##########################################
#WinBugs model
##########################################

library(R2jags)

sink("IPM.bug")
cat("model{

##########################################
#Clutch size model
##########################################

#Priors
b0c ~ dnorm(0,0.01)
b1c ~ dnorm(0,0.01)
sdc ~ dunif(0,10)
tauc <- pow(sdc,-2)

for(t in 1:T1){
  ec[t] ~ dnorm(0,tauc)
  log(C[t]) <- b0c + ec[t] + b1c * Time[t]
  
} 

for(i in 1:length(c)){
  c[i] ~ dpois(C[c.year[i]])
}

##########################################
#Egg hatch rate model
##########################################

#Prior
E ~ dunif(0,1)

for(i in 1:length(hatched)){
  hatched[i] ~ dbin(E, laid[i])
}

##########################################
#Nest success model
##########################################

#Priors
b0n ~ dnorm(0,0.01)
b1n ~ dnorm(0,0.01)
sdn ~ dunif(0,10)
tau.n <- pow(sdn,-2)

for(t in 1:length(nd)){

  en[t] ~ dnorm(0,tau.n)
  logit(P[t]) <- b0n + en[t] + b1n * Time[t]
  
  nest_deaths[t] ~ dbin(P[t],nd[t])
  
  NS[t] <- (1-P[t])^37
  
  #Cumulative nest success (accounts for renesting)
  for(j in 1:4){
    ns.temp[t,j] <- NS[t]*(1-NS[t])^(j-1)*pN[j]
  }
  NS.cu[t] <- sum(ns.temp[t,])
  
}

##########################################
#State-space model - Code adapted from Brooks et al. 2004. 
##########################################

#Prior 
N[1] ~ dnorm(101, 1)I(0,) 

# Define population growth rate lambda[t] 
for(t in 2:(T1-1)){
  lambda[t] <- N[t+1]/N[t]
}

#Define the system process for the census/index data using the Normal/Binomial approximation
for(t in 2:T1){

mean1[t] <- rho[t-1]*N[t-1]*S[t-1]
meana[t] <- S[t-1]*N[t-1]


tau1[t] <- 1/(N[t-1]*rho[t-1]*S[t-1]) #Poisson approximation
taua[t] <- 1/(N[t-1]*S[t-1]*(1-S[t-1])) #Binomial approximation

N1[t] ~ dnorm(mean1[t],tau1[t])
Na[t] ~ dnorm(meana[t],taua[t])
N[t] <- N1[t]+Na[t]
}

sd.y ~ dunif(0,10)
tau.y <- pow(sd.y, -2)

# Define the observation process for the census/index data
for(t in 2:(T1-1)){
  y[t] ~ dnorm(N[t],tau.y)
}


## Productivity ##

for(t in 1:T1){
  rho[t] <- NS.cu[t] * C[t] * E * D[t] * 0.5
}

##########################################
#Adult survival model - Code adapted from Brooks et al. 2004
##########################################

#Priors

#Intercept and slope for logistic survival model
b0s ~ dnorm(0,0.01)
b1s ~ dnorm(0,0.01)

#Intercept and slope for logistic recovery rate model
b1r ~ dnorm(0,0.01)
b0r ~ dnorm(0,0.01)

#Survival SD
sds ~ dunif(0,10) 
tau.s <- pow(sds, -2)

#Recovery SD
sdr ~ dunif(0,10)
tau.r <- pow(sdr, -2)

#Logistic models for survival and recovery
for(t in 1:(T1-1)){

  es[t] ~ dnorm(0,tau.s) 
  esr[t] ~ dnorm(0,tau.r) 
  
  logit(S[t]) <- b0s + es[t] 
  logit(r[t]) <- b0r +  esr[t]
}

# Model the observed recoveries per cohort per year
# as a multinomial distribution with p[t] being the expected recovery probabilities.
# p[t] is filled in below. rel[t] is the number of released birds per cohort and is supplied
# in the data.  

#Define the likelihood
for(t in 1:(T1-1)){
    m[t, 1:(T2 + 1)] ~ dmulti(p[t, ], rel[t]) 
}
  
# Here we calculate the cell probabilities for the recovery table. In practice, we construct
# the Seber table from Program MARK - 'A Gentle Introduction', ch. 8 p. 27.


# Columns 1 to (T2-1) give the dead recovery probabilities. Column T1 gives the birds
# you never expect to see again. 

for(t1 in 1:(T2-1)){ 

 p[t1,t1] <- r[t1] * (1-S[t1]) #the diagonal

  p[t1,t1+1] <- r[t1+1] * S[t1]*(1-S[t1+1]) #the diagonal above that diagonal

    for(t2 in (t1+2):T2){

      for(t in (t1+1):(t2-1)){

          logS[t1,t2,t] <- log(S[t]) #this a a logarithm trick so you can use sum() to do a product sum below
      }

     p[t1,t2] <- r[t2] * S[t1]*(1-S[t2]) * exp(sum(logS[t1,t2,(t1+1):(t2-1)])) #and all the others

   }
   
 for(t2 in 1:(t1-1)){
    # Zero probabilities in lower triangle of table
    p[t1,t2] <- 0
  }
  
  # Probability of an animal never being seen again (1 - the other probabilities combined)
  p[t1,T2+1] <- 1 - sum(p[t1,1:T2])
  
  }
  
  p[T2,T2] <- r[T2]*(1-S[T2])

  for(t in 1:(T2-1)){
  p[T2,t] <- 0

  }
  
  p[T2,T2+1] <- 1 - p[T2,T2] #The ones not found dead are ones you'll technically never see again


##########################################
#Duckling survival model
##########################################

#Slope and intercept priors for duckling survival logistic model
for(i in 1:3){
  alpha.phiD[i] ~ dnorm(0,0.01)
  beta.phiD[i] ~ dnorm(0,0.01)
}

#Priors for the years where D is unknown

for(t in 1:(T1-3)){
 D[t] ~ dunif(0,1)
}

# Model daily duckling survival rate (DSRd)
for (i in 1:length(count)){

  logit(phiD[i]) <- alpha.phiD[d.year[i]] + beta.phiD[d.year[i]] * age[i]
  nextCount[i] ~  dpois(count[i] * (phiD[i]^interval[i])) 
}

#Calculate duckling survival to fledging age 
for (year in 1:3){
  for(d in 1:length(u.days)){ 
    DSRd[year,d] <- exp(alpha.phiD[year] + beta.phiD[year]*u.days[d])/(1+exp(alpha.phiD[year] + beta.phiD[year]*u.days[d]))
  }
  D[year+(T1-3)] <- exp(sum(log(DSRd[year,1:56])))
}

##########################################
#Extract mean vital rates
##########################################

S.mean <- sum(S)/length(S)
NS.mean <- sum(NS)/length(NS)
NS.cu.mean <- sum(NS.cu)/length(NS.cu)
D.mean <- sum(D)/length(D)
C.mean <- sum(C)/length(C)
rho.mean <- sum(rho)/length(rho)
l.mean <- sum(lambda[2:(T1-1)])/length(lambda[2:(T1-1)])
r.mean <- sum(r)/length(r)

}")
sink(NULL)

##########################################
#Run model
##########################################

#Parameters to be outputted
parameters <- c("S.mean","rho.mean","C.mean","D.mean", "NS.mean", "NS.cu.mean","E","r.mean",
                "sds","sdn","sdc","l.mean","sdr","S","rho","C","D","NS", "NS.cu","r","lambda","sd.y","N",
                "beta.phiD","b1n","b1c","DSRd")

#Run the model
model <- jags(data=data,parameters.to.save=parameters,model.file="IPM.bug",n.chains=3,n.iter=n.iter,n.burnin=n.iter/2,n.thin=n.iter/2000)

##########################################
#Retrieve model output
##########################################

options(max.print=10000)
print(model)

########### Life-stage sensitivity analysis ###########

#Get sd for egg hatch rate
E.df <- EHR.data %>%
  group_by(Year) %>%
  dplyr::summarize(C = sum(C), hatched = sum(Hatched))

E.sd <- sd(E.df$hatched/E.df$C)

NS.cu <- model[["BUGSoutput"]][["mean"]][["NS.cu"]] 
C <- model[["BUGSoutput"]][["mean"]][["C"]] 
E <- as.numeric(model[["BUGSoutput"]][["mean"]][["E"]])
D <- model[["BUGSoutput"]][["mean"]][["D"]] 
S <- model[["BUGSoutput"]][["mean"]][["S"]]


#Parameters to vary
params <- data.frame(Parameters = c("C","E", "NS.cu","D","S"),
                     Mean = c(mean(C),E,mean(NS.cu),mean(D),mean(S)),
                     sd = c(NA,E.sd,sd(NS.cu),sd(D),sd(S)),
                     Var = c(NA,E.sd,sd(NS.cu),sd(D),sd(S))^2)

#estimate alpha and beta for the beta distributions
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

params$alpha <- unlist(estBetaParams(params$Mean,params$Var)[1])
params$beta <- unlist(estBetaParams(params$Mean,params$Var)[2])

repParams <- matrix(nrow=5,ncol=LSA.sim) #10000 replicates
repParams[1,] <- runif(LSA.sim,min(C),max(C)) #Clutch size has a uniform distribution with as range the min and max annual C values

#The rest has beta distributions
for(i in 2:5){
  repParams[i,] <- rbeta(LSA.sim,params$alpha[i],params$beta[i])
}

#The model
#Parameters

sensMat <- matrix(nrow=LSA.sim,ncol=5)
elasMat <- matrix(nrow=LSA.sim,ncol=5)
lambda <- c()

for(i in 1:LSA.sim){
  
  c <- repParams[1,i] #clutch size
  e <- repParams[2,1]
  ns.cu <- repParams[3,i] 
  d <- repParams[4,i] 
  s <- repParams[5,i] 
  
  duck.vr <- list(c=c,e=e,ns.cu=ns.cu,d=d,s=s)
  
  duck.el <- expression(s*c*e*d*ns.cu*0.5, s*c*e*d*ns.cu*0.5,s,s)
  
  output <- vitalsens(duck.el,duck.vr)
  
  sensMat[i,] <- output$sensitivity
  elasMat[i,] <- output$elasticity
  
  Fec <-  s*c*0.5*e*d*ns.cu
  
  library(popbio)
  
  stages <- c("1st-year","adult")
  A <- matrix(c(Fec,s,Fec,s),ncol=2,nrow=2,dimnames=list(stages,stages))
  n <- c(5,5)
  p <- pop.projection(A,n,150)
  lambda[i] <- p$lambda
  
}

mean(lambda)
