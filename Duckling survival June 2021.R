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

table <- table[!(table$valid %in% "n") & !(table$valid %in% "nn"),] #Remove uncertain observations
table <- table[table$soort=="mal",] #Only mallards

library(lubridate)
table$maand <- month(table$datum)

### Data visualization

library(data.table)
table <- data.table(table)

# get the observations where a family was reported to have been observed at least twice.
repObs <- names(table(table$ID_gezin)[table(table$ID_gezin)>1])

d <- table[table$ID_gezin %in% repObs,]

d$datum <-  as.Date(d$datum,format="%d-%m-%Y")

d<-d[order(d$ID_gezin,d$datum),]

#Only interested in the pre-fledging period
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