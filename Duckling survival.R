# Here we calculate capture history of eack duckling as in Schmidt et al. 2010.
rm(list=ls())

### Data preparation

library(readr)
table <-  read_delim("~/Msc minor/waarnemingen 2016-2020 bewerkt.csv", 
                     ";", escape_double = FALSE, col_types = cols(ID = col_character(), 
                                                                  ID_gezin = col_character(), causeOfDeath = col_character(), 
                                                                  extraFeedings = col_character(), 
                                                                  gezinEerderGemeld = col_character(), 
                                                                  jaar = col_character(), numberOfDeaths = col_number(), 
                                                                  shore = col_character(), surface = col_character(), 
                                                                  tijd = col_time(format = "%H:%M:%S"), 
                                                                  valid = col_character(), water = col_character(), 
                                                                  zekerheid = col_character()), trim_ws = TRUE)

library(stringr)
library(lubridate)

#Remove invalid observations and months with sample size < 20
#NB: regular subsetting messes with NAs, so here I use %in%
table <- table[!(table$valid %in% "n") & !(table$valid %in% "nn"),]
table <- table[table$soort=="mal",]
#table <- table[table$mo >3 & table$mo < 9,]


table$maand <- month(table$datum)
y <- table[table$leeftijd==2,]
table(y$maand)



### Data visualisation

library(data.table)
table <- data.table(table)

# get the observations where a family was reported to have been observed at least twice.
repObs <- names(table(table$ID_gezin)[table(table$ID_gezin)>1])

d <- table[table$ID_gezin %in% repObs,]

d$datum <-  as.Date(d$datum,format="%d-%m-%Y")


d<-d[order(d$ID_gezin,d$datum),]

#Only interested in the pre-fledgeling period
d <- d[d$leeftijd_dagen<56,]


# determine interval length
d$interval<-d$nextAantal<-d$nextLeeftijd_dagen<-NA
d$midDay<-d$datum

#This loop adds to each row the no. ducklings of the next observation, and the age of those ducklings. 

for (i in 1:(nrow(d)-1)) { #Loop through the broods
  d$interval[i]<-as.numeric(d$datum[i+1]-d$datum[i]) 
  d$midDay[i]<-mean(d$datum[i:(i+1)])
  d$nextAantal[i]<-d$aantal[i+1]
  d$nextLeeftijd_dagen[i]<-d$leeftijd_dagen[i+1]
}
d$midDoy<-as.numeric(strftime(d$midDay,format="%j"))
d$midDoyScaled<-as.numeric(scale(d$midDoy))

# remove last observation (intervals are only for within)
d<-d[c(d$ID_gezin[1:(nrow(d)-1)]==d$ID_gezin[2:nrow(d)],FALSE),]


#Too long intervals are unreliable.
d <- d[d$interval<=30 & d$interval > 0,]


#These are the parameters used in Mayfield nest survival estimation (we use if for ducklings)
d$days <- d$aantal*d$interval #Weigh the number of days with the brood size
d$days_surv <- as.integer(d$days-(d$aantal-d$nextAantal)*d$interval/2) #Days survived by the ducklings (in case of intervals of unknowns, ducklings are assumed to have survived half of the interval)
d$days_notsurv <- d$aantal-d$nextAantal #deaths
d$ageMid <- (d$leeftijd_dagen+d$nextLeeftijd_dagen)/2 #average age 

#This removes observations that recorded more ducklings than there were at the previous observation 
d<-subset(d,days_notsurv>=0)

d <- d[d$days!=0,]

d$jaar = as.numeric(d$jaar)

d <- d[d$jaar %in% c(2018:2020),]


library(R2jags)

d <- merge(d, weather.covar,by="jaar")

setwd("~/Msc minor/IPM structures")
saveRDS(d, "duckling.df.rData")

data <- list(count = d$aantal, nextCount = d$nextAantal, interval = d$interval,
             age = as.numeric(scale(d$ageMid)), year = as.numeric(scale(d$jaar)),
             u.year = unique(as.numeric(scale(d$jaar))),u.days = (c(1:56)-mean(d$ageMid))/sd(d$ageMid),
             adays=d$adays,doy = d$midDoyScaled,u.adays = unique(d$adays))



sink("D April.bug")
cat("model{

  alphaS ~ dnorm(0,0.01)
  betaS ~ dnorm(0,0.01)
  betaAdays ~ dnorm(0,0.01)
  #betaDoy ~ dnorm(0,0.01)
  

for (i in 1:length(count)){

  logit(S[i]) <- alphaS + betaS * age[i] + betaAdays * adays[i] + 
  nextCount[i] ~  dpois(count[i] * (S[i]^interval[i])) 
}

for (year in 1:3){
for(d in 1:length(u.days)){ 
    DSR[year,d] <- exp(alphaS + betaS*u.days[d]+ betaAdays*u.adays[year])/(1+exp(alphaS + betaS*u.days[d] + betaAdays*u.adays[year]))
}
  D[year] <- exp(sum(log(DSR[year,1:56])))
}



}")
sink(NULL)

parameters <- c("D","betaAdays","betaS")

model <- jags(data=data,parameters.to.save=parameters,model.file="D April.bug",n.chains=1,n.iter=5000,n.burnin=2000)

print(model)
