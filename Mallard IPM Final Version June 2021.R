rm(list=ls())

setwd("~/Msc minor/IPM")

#POPULATION INDEX DATA 
#Census data for the Mallard population in the Netherlands.

library(readxl)
trend <- read_excel("trends_broedvogels_1990-2019_sovon.xlsx", 
                    sheet = "Landelijke trends 1990-2019")
trend <- trend[trend$Soort=="Wilde Eend",-c(1:4,35:36)]

############## NEST DATA
#Used to calculate egg hatch rate, clutch size, and nest success. 
library(readr)
nest.data <- read_delim("nestkaarten_eenden.csv", 
                        ";", escape_double = FALSE, trim_ws = TRUE)
nest.data$legbegin <- scale(nest.data$legbegin)
nest.data <- nest.data[nest.data$euring==1860,] #1860 = Mallard 

startYear <- 2000


#EGG HATCH RATE
EHR.data <- nest.data[nest.data$jaar>=startYear,]
EHR.data <- EHR.data[!is.na(EHR.data$legselgr),]
EHR.data <- EHR.data[!is.na(EHR.data$uitgekomen),]
EHR.data <- EHR.data[EHR.data$brd_succes>2 & EHR.data$brd_succes<12,] #only successful nests
EHR.data <- EHR.data[EHR.data$uitgekomen>0,] #at least one egg must have hatched
EHR.data <- EHR.data[EHR.data$legselgr>=EHR.data$uitgekomen,]
EHR.data <- na.omit(EHR.data)


### CLUTCH SIZE
clutch.data <- nest.data[nest.data$jaar>=startYear,]
clutch.data <- clutch.data[!is.na(clutch.data$legselgr),] #Need to know clutch size
clutch.data <- clutch.data[clutch.data$brd_succes>2 & clutch.data$brd_succes<12,] #only successful nests

clutch.data <- subset(clutch.data, jaar!="2014") #No successful nests in 2014

library(dplyr)
CS <- clutch.data %>%
  group_by(jaar) %>%
  dplyr::summarize(mean = mean(legselgr),n = n())

#NEST SUCCES
ns.data <- nest.data[nest.data$jaar>=startYear,]
ns.data <- ns.data[!is.na(ns.data$novl),]
ns.data <- ns.data[!is.na(ns.data$nestdagen),]
ns.data <- ns.data[ns.data$nestdagen!=0,]

NS <- ns.data %>%
  group_by(jaar) %>%
  dplyr::summarize(nestdagen = sum(nestdagen), novl = sum(novl))
NS$p <- 1- (NS$nestdagen-NS$novl)/NS$nestdagen
NS$NS <- NS$p^37

##RENESTING INTENSITY
#Here we calculate the chance of renesting p[x] for the xth nest attempt according to the formula in Hoekman et al. 2002. 

r <- 0.26
p <- c(0.968)
for(i in 2:3){
  p[i] <- 1-r*i
}
p <- c(p[1],p[1]*p[2],p[1]*p[2]*p[3])


### ADULT SURVIVAL ####
#Data wrangling in 'Ring recoveyr to marrays June 2021.R' 

#The m-array derived from the ring-recovery data. 

m_array <- read_csv("2k m-array.csv", col_names = FALSE)

release <- m_array[,1] #Number of birds released per year
m_array <- m_array[,-1]
m_array <- cbind(m_array,release-rowSums(m_array)) #Expected number of birds that are never seen again

### DUCKLING SURVIVAL
#Data wrangling in 'Duckling survival.R'

d <- readRDS("duckling.df.rData")

data <- list(y=unlist(trend[11:30]),T1=21,T2=20, Time = as.numeric(scale(1:21)), #Global variables
             m=m_array, rel=unlist(release), #For adult survival
             count = d$aantal, nextCount = d$nextAantal, interval = d$interval, age = as.numeric(scale(d$ageMid)), #for Duckling survival 
             u.days = (c(1:56)-mean(d$ageMid))/sd(d$ageMid), d.year = (d$jaar-2017), # for Duckling survival
             pN=p, nd=NS$nestdagen, nest_deaths = (NS$nestdagen-NS$novl),  #For nest success 
             c=clutch.data$legselgr,c.year=clutch.data$jaar-1999, #For clutch size
             laid = EHR.data$legselgr, hatched=EHR.data$uitgekomen) #For egg hatch rate

library(R2jags)


sink("IPM.bug")
cat("model{


########## Clutch size  ####################

############################################

b0c ~ dnorm(0,0.001)
b1c ~ dnorm(0,0.001)

sdc ~ dunif(0,10)
tauc <- pow(sdc,-2)


for(t in 1:21){
  ec[t] ~ dnorm(0,tauc)T(-10,10)
  log(C[t]) <- b0c + ec[t] + b1c * Time[t]
  
} 

for(i in 1:length(c)){
  c[i] ~ dpois(C[c.year[i]])
}

######## Egg hatch rate #####################

#############################################

E ~ dunif(0,1)

for(i in 1:length(hatched)){
  hatched[i] ~ dbin(E, laid[i])
}

####### Nest success ######################

############################################

b0n ~ dnorm(0,0.001)
b1n ~ dnorm(0,0.001)

sdn ~ dunif(0,5)
tau.n <- pow(sdn,-2)

for(t in 1:length(nd)){

  en[t] ~ dnorm(0,tau.n)
  logit(P[t]) <- b0n + en[t] + b1n * Time[t]
  
  nest_deaths[t] ~ dbin(P[t],nd[t])
  
  NS[t] <- (1-P[t])^37
  
  #Cumulative nest success (accounts for renesting)
  for(j in 1:3){
    ns.temp[t,j] <- NS[t]*(1-NS[t])^(j-1)*pN[j]
  }
  NS.cu[t] <- sum(ns.temp[t,])
  
}

############ State-space model ##############

#############################################

# Initial population prior
N[1] ~ dnorm(100.2,0.000001)

# Define r[t]
for(t in 2:(T1-1)){
  r[t] <- N[t+1]/N[t]
}

# Define the system process for the census/index data using the Normal/Binomial approximation
for(t in 2:T1){
mean1[t] <- rho[t-1]*N[t-1]*phi[t-1]
meana[t] <- phi[t-1]*N[t-1]


tau1[t] <- 1/(N[t-1]*rho[t-1]*phi[t-1]) #Poisson approximation
taua[t] <- 1/(N[t-1]*phi[t-1]*(1-phi[t-1])) #Binomial approximation

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


########### Adult survival  ################

############################################

#Priors

#Intercept and slope for logistic survival model
b0s ~ dnorm(0,0.001)
b1s ~ dnorm(0,0.001)

#Intercept and slope for logistic recovery rate model
b1l ~ dnorm(0,0.001)
b0l ~ dnorm(0,0.001)

#Survival SD
sds ~ dunif(0,5) 
tau.s <- pow(sds, -2)

#Recovery SD
sdl ~ dunif(0,5)
tau.l <- pow(sdl, -2)

#Logistic models for survival and recovery
for(t in 1:(T1-1)){

  es[t] ~ dnorm(0,tau.s)T(-10,10) 
  esl[t] ~ dnorm(0,tau.l)T(-10,10)
  
  logit(phi[t]) <- b0s + es[t] +  b1s * Time[t]
  logit(lambda[t]) <- b0l +  esl[t] +  b1l * Time[t]
}

# Model the observed recoveries per cohort per year
# as a multinomial distribution with p[t] being the expected recovery probabilities.
# p[t] is filled in below. rel[t] is the number of released birds per cohort and is supplied
# in the data.  

#m is for the adults
  
#Define the likelihood
for(t in 1:(T1-1)){
    m[t, 1:(T2 + 1)] ~ dmulti(p[t, ], rel[t]) 
}
  
# Here we calculate the cell probabilities for the recovery table. We basically construct
# the Seber table from Mark ch. 8 p. 27.


# Columns 1 to (T2-1) give the dead recovery probabilities. Column T1 gives the birds
# you never expect to see again. 

for(t1 in 1:(T2-1)){ 

 p[t1,t1] <- lambda[t1] * (1-phi[t1]) #the diagonal

  p[t1,t1+1] <- lambda[t1+1] * phi[t1]*(1-phi[t1+1]) #the diagonal above that diagonal

    for(t2 in (t1+2):T2){

      for(t in (t1+1):(t2-1)){

          lphi[t1,t2,t] <- log(phi[t]) #this a a log trick so you can use sum() to do a product sum below
      }

     p[t1,t2] <- lambda[t2] * phi[t1]*(1-phi[t2]) * exp(sum(lphi[t1,t2,(t1+1):(t2-1)])) #and all the others

   }
   
 for(t2 in 1:(t1-1)){
    # Zero probabilities in lower triangle of table
    p[t1,t2] <- 0
  }
  
  # Probability of an animal never being seen again (1 - the other p's combined)
  p[t1,T2+1] <- 1 - sum(p[t1,1:T2])
  
  }
  
  # Final ROW (not column you dipshit)
  p[T2,T2] <- lambda[T2]*(1-phi[T2])

  for(t in 1:(T2-1)){
  p[T2,t] <- 0

  }
  
  p[T2,T2+1] <- 1 - p[T2,T2] #The ones not found dead are ones you'll technically never see again


########## Duckling survival ###############

############################################

#Slope and intercept for duckling survival logistic model
for(i in 1:3){
  alphaS[i] ~ dnorm(0,0.01)
  betaS[i] ~ dnorm(0,0.01)
}

#Priors for the years where D is unknown
for(t in 1:18){
 D[t] ~ dunif(0,1)
}

# Model daily duckling survival
for (i in 1:length(count)){

  logit(S[i]) <- alphaS[d.year[i]] + betaS[d.year[i]] * age[i]
  nextCount[i] ~  dpois(count[i] * (S[i]^interval[i])) 
}

#Calculate duckling survival to fledging age 
for (year in 1:3){
  for(d in 1:length(u.days)){ 
    DSRd[year,d] <- exp(alphaS[year] + betaS[year]*u.days[d])/(1+exp(alphaS[year] + betaS[year]*u.days[d]))
  }
  D[year+18] <- exp(sum(log(DSRd[year,1:56])))
}

######### Mean vital rates ########################

###################################################

phi.mean <- sum(phi)/length(phi)
NS.mean <- sum(NS)/length(NS)
NS.cu.mean <- sum(NS.cu)/length(NS.cu)
D.mean <- sum(D)/length(D)
C.mean <- sum(C)/length(C)
rho.mean <- sum(rho)/length(rho)
r.mean <- sum(r[2:(T1-1)])/length(r[2:(T1-1)])
l.mean <- sum(lambda)/length(lambda)

}")
sink(NULL)

#Parameters to be outputted
parameters <- c("phi.mean","rho.mean","C.mean","D.mean", "NS.mean", "NS.cu.mean","E","r.mean",
                "sds","sdn","sdc","l.mean","sdl","phi","rho","C","D","NS", "NS.cu","r","lambda","sd.y")

#Mean parameters only (for the table)
parameters <- c("phiF.mean","phi.mean","rho.mean","C.mean","D.mean", "NS.mean", "NS.cu.mean","E","r.mean",
                "sds","sdsF","sdn","sdc","l.mean","sdl","sd.y")

#Run the model
model <- jags(data=data,parameters.to.save=parameters,model.file="IPM.bug",n.chains=3,n.iter=100000,n.burnin=30000,n.thin=10)

print(model)
