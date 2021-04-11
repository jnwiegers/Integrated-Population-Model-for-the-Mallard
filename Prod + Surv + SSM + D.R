rm(list=ls())

setwd("~/Msc minor/IPM")


#POPULATION INDEX DATA 

library(readxl)
trend <- read_excel("trends_broedvogels_1990-2019_sovon.xlsx", 
                    sheet = "Landelijke trends 1990-2019")
trend <- trend[trend$Soort=="Wilde Eend",-c(1:4,35:36)]


#CLIMATE COVARIATE DATA

KNMI <- read_excel("etmgeg_260.xlsx")
KNMI <- as.data.frame(KNMI)
colnames(KNMI)[1] <- "Day"

KNMI$Day <- as.Date(as.character(substr(KNMI$Day,5,12)),format="%Y%m%d")

library(lubridate)
KNMI <- KNMI[year(KNMI$Day)>1999,]
KNMI <- KNMI[year(KNMI$Day)<2021,]

KNMI <- KNMI[,colnames(KNMI) %in% c("Day","TX,","TN,")]
colnames(KNMI)[2:3] <- c("Tmin","Tmax")
KNMI$Tmax <- as.numeric(gsub(",","",KNMI$Tmax))/10
KNMI$Tmin <- as.numeric(gsub(",","",KNMI$Tmin))/10


library(dplyr)

#Days per year where Tmax and Tmin are below 0 degrees 
df.W <- KNMI %>%
  group_by(year(KNMI$Day)) %>%
  summarize(fdays = sum(Tmin<0),idays = sum(Tmax<0))

#Per year, the no. days in April below 10 degrees 

KNMI.april <- KNMI[month(KNMI$Day) == 4,]

df.A <- KNMI.april %>%
  group_by(year(KNMI.april$Day)) %>%
  summarize(cold.April.days = sum(Tmax<10))

weather.covar <- data.frame(jaar = c(2000:2020),
                            fdays = scale(df.W$fdays),
                            idays = scale(df.W$idays),
                            adays = scale(df.A$cold.April.days))


############## NEST DATA
library(readr)
nest.data <- read_delim("nestkaarten_eenden.csv", 
                        ";", escape_double = FALSE, trim_ws = TRUE)
nest.data$legbegin <- scale(nest.data$legbegin)
nest.data <- nest.data[nest.data$euring==1860,]

startYear <- 2000


#EGG HATCH RATE
EHR.data <- nest.data[!is.na(nest.data$legselgr),]
EHR.data <- EHR.data[!is.na(EHR.data$uitgekomen),]
EHR.data <- EHR.data[EHR.data$brd_succes>2 & EHR.data$brd_succes<12,] #only successful nests
EHR.data <- EHR.data[EHR.data$uitgekomen>0,] #nest success = at least one egg hatched. 
EHR.data <- EHR.data[EHR.data$legselgr>=EHR.data$uitgekomen,]
EHR.data <- na.omit(EHR.data)

### CLUTCH SIZE
clutch.data <- nest.data[nest.data$jaar>=startYear,]
clutch.data <- clutch.data[!is.na(clutch.data$legselgr),] #Need to know clutch size
clutch.data <- clutch.data[clutch.data$legselgr != 0,] #Females need to have laid at least one egg
clutch.data <- merge(clutch.data, weather.covar,by="jaar")

#CS <- clutch.data %>%
  #group_by(jaar) %>%
  #summarize(mean = mean(legselgr))

#NEST SUCCES
ns.data <- nest.data[nest.data$jaar>=startYear,]
ns.data <- ns.data[!is.na(ns.data$novl),]
ns.data <- ns.data[!is.na(ns.data$nestdagen),]
ns.data <- ns.data[ns.data$nestdagen!=0,]
ns.data <- merge(ns.data, weather.covar,by="jaar")
#ns.data <- ns.data[ns.data$jaar!=2014,] #Year with little monitoring; too low sample size


library(readr)
m_array <- read_csv("2k m-array.csv", 
                    col_names = FALSE)

release <- m_array[,1]
m_array <- m_array[,-1]
m_array <- cbind(m_array,release-rowSums(m_array))

### DUCKLING SURVIVAL
#Data wrangling in 'Duckling survival.R'

d <- readRDS("duckling.df.rData")


data <- list(y=unlist(trend[11:30]),T1=21,T2=20,m=m_array, rel=unlist(release), time = 2021-startYear,
             u.fdays=weather.covar$fdays, u.year=unique(as.factor(scale(ns.data$jaar))),c = clutch.data$legselgr, yearC = as.factor(scale(clutch.data$jaar)),
             EHR.laid = EHR.data$legselgr, EHR.hatch = EHR.data$uitgekomen, nobs = ns.data$nestdagen
             ,nsurv=ns.data$novl, yearN = as.factor(scale(ns.data$jaar)),
             fdaysC = clutch.data$fdays,fdaysN = ns.data$fdays, count = d$aantal, nextCount = d$nextAantal, interval = d$interval,
             age = as.numeric(scale(d$ageMid)), year = as.numeric(scale(d$jaar)),
             u.days = (c(1:56)-mean(d$ageMid))/sd(d$ageMid),
             adays=d$adays,doy = d$midDoyScaled,u.adays = unique(d$adays))


library(R2jags)





sink("IPM April.bug")
cat("model{

### Priors ###

#Egg hatch rate
E ~ dunif(0,1) 

#Clutch size
alphac ~ dnorm(0,0.01) 
betac1 ~ dnorm(0,0.01) #year
betac2 ~ dnorm(0,0.01) #days below 0

#Nest succes
alphaN ~ dnorm(0,0.01)
#betaN1 ~ dnorm(0,0.01)
betaN2 ~ dnorm(0,0.01)

#Observation error and initial population priors
tauy ~ dgamma(0.001,0.001) 
N[1] ~ dnorm(100.2,0.000001)

#Survival and recovery
phiC ~ dunif(0,1) 
lambdaC ~ dunif(0,1)


for(t in 1:(T1-1)){
  phi[t] <- phiC
  lambda[t] <- lambdaC
  #rho[t] ~ dunif(0,1)
}

#Duckling survival
alphaS ~ dnorm(0,0.01)
betaS ~ dnorm(0,0.01)
betaAdays ~ dnorm(0,0.01)

#Priors for the years where D is unknown
for(t in 1:18){
  D[t] ~ dunif(0,1)
}

#Model egg hatch rate
for(i in 1:length(EHR.laid)){
  EHR.hatch[i] ~ dbin(E, EHR.laid[i])
}

#Model clutch size
for(i in 1:length(c)){
  log(C[i]) <- alphac + betac1*yearC[i]+betac2*fdaysC[i]
  c[i] ~ dpois(C[i])
}

#Model nest success
for(i in 1:length(nobs)){
  logit(DSRn[i]) <- alphaN +betaN2*fdaysN[i] 
  nsurv[i] ~ dbin(DSRn[i],nobs[i])
}

#Estimate annual mean clutch size and nest success
for(t in 1:time){ 
  C.year[t] <- exp(alphac + betac1*u.year[t]+betac2*u.fdays[t])
  NS.year[t] <- (exp(alphaN + betaN2*u.fdays[t])/(1+exp(alphaN + betaN2*u.fdays[t])))^37
}

### State-space model ###

# Define r[t]
for(t in 2:(T1-1)){
r[t] <- N[t+1]/N[t]
}


# Define the system process for the census/index data using the Normal approximation
for(t in 2:T1){
mean1[t] <- rho[t-1]*phi[t-1]*N[t-1]
meana[t] <- phi[t-1]*N[t-1]
tau1[t] <- 1/(N[t-1]*rho[t-1]*phi[t-1])
taua[t] <- 1/(N[t-1]*phi[t-1]*(1-phi[t-1]))

N1[t] ~ dnorm(mean1[t],tau1[t])
Na[t] ~ dnorm(meana[t],taua[t])
N[t] <- N1[t]+Na[t]
}

# Define the observation process for the census/index data
for(t in 2:(T1-1)){
  y[t] ~ dnorm(N[t],tauy)
}


### Productivity ###
# 0.968 is mallard breeding propensity (Hoekman et al. 2002)
# The 1.57 factor accounts for failed nests replaced through renesting attempts (Arnold et al. 2010) 
# The 0.50 excludes the males

for(t in 1:T1){
  rho[t] <- 0.968 * NS.year[t]*1.57 * C.year[t] * E * D[t] * 0.50
}

### Adult survival with ring recovery ###

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

### Duckling survival

for (i in 1:length(count)){

  logit(S[i]) <- alphaS + betaS * age[i] + betaAdays * adays[i]
  nextCount[i] ~  dpois(count[i] * (S[i]^interval[i])) 
}

for (year in 1:3){
for(d in 1:length(u.days)){ 
    DSRd[year,d] <- exp(alphaS + betaS*u.days[d]+ betaAdays*u.adays[year])/(1+exp(alphaS + betaS*u.days[d] + betaAdays*u.adays[year]))
}
  D[year+18] <- exp(sum(log(DSRd[year,1:56])))
}

}")
sink(NULL)

parameters <- c("phi[1]","lambda[1]","N","E", "C.year","alphaN","betac1","NS.year","betac2","betaN2", "D","betaS","betaAdays","r")

model <- jags(data=data,parameters.to.save=parameters,model.file="IPM April.bug",n.chains=1,n.iter=5000,n.burnin=2000)

print(model)

prod(model[["BUGSoutput"]][["mean"]][["r"]][-1])



 
