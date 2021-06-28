rm(list=ls())

library(readr)

#Read in the ringing dataset

setwd("~/")
data <- read_delim("Msc minor//2020-04 wilde eend Erik Kleyheeg Sovon - kopie.csv", 
                   ";", escape_double = FALSE, col_types = cols(Datum = col_character(), 
                                                                InsertDate = col_character(), ReportingDate = col_character(), 
                                                                UpdateDate = col_character(), Ringleeftijd = col_character()), trim_ws = TRUE)

#We use data from March 2000 to March 2020
library(lubridate)
startYear <- 2000 #First year
startDate <- as.Date("2000-03-01") #Begin of first breeding season
endDate <- as.Date("2020-03-01") #End of last breeding season

#Remove inaccurate or missing data
data <- data[year(data$Datum)>=startYear,]
data <- data[data$Datum<endDate,]
data <- data[data$`Nauwkeurigheid Datum` < 8,]
data <- data[!is.na(data$RingID),]
data <- data[data$Conditie != 0 ,]

library(data.table)

#Sort the data to retrieve the date the birds were first captured
mydata <- data[with(data, do.call(order, list(RingID,Datum))), ]
DT <- data.table(mydata, key = "RingID")
DT <- DT[, head(.SD, 1), by = key(DT)] #first capture date of each bird
data <- merge(data,DT[,c(1,8)],by="RingID")
colnames(data)[48] <- "FirstCapDate"

#Make an event category
Event <- c(length=nrow(data)) 
Event[data$MetalringInformation %in% c(1,2,3)]<- "Ringed"
Event[data$MetalringInformation %in% c(4,5,6,7)]<- "Live recapture"
Event[data$Conditie %in% c(1,2,3)]<- "Dead recovery"
data$Event <- Event

#We only look at birds ringed in the breeding season (March - July)
data <- data[month(data$FirstCapDate) %in% c(3,4,5,6,7),]

#Birds ringed in 2020 don't matter so we remove those
data <- data[!(data$Event == "Ringed" & year(data$Datum.x) == "2020"),]

data <- data[data$Omstandigheden != 56,] #Remove deaths by botulism 

library(ggplot2)

#Loop to fill in the m-array based on the data. 

markVec <- c() #Vector of individual capture histories
mistakes <- c() #Keeps track of impossible observations (e.g. two deaths)

#Keeps track of sex and age if needed later
sex <- c()
age <- c()

for (i in 1:length(unique(data$RingID))){ #loop through the ducks
  indiv <-  unique(data$RingID)[i]
  
  
  set <- data[data$RingID == indiv,c('Conditie','Datum.x','Ringleeftijd','Geslacht', 'Event')]
  set <- set[order(set$Datum),]
  set <- set[!apply(is.na(set), 1, all),]
  
  
  seenVec <- rep("0",year(endDate)-startYear) #Vector of live recaptures (we only use the first one here, i.e. Seber format)
  recVec <- rep("0",year(endDate)-startYear) #Vector of dead recoveries
  
  encountered <- F 
  for(j in 1:nrow(set)){ #Loop through the individual history
    if(set$Conditie[j] > 3 && encountered == F){ #If alive
      
      seenVec[1+floor((as.Date(set$Datum[j])-startDate)/365.25)] <- 1
      encountered <- T #By setting this to T you only get the first live recapture (Seber format)
      
    }else if (set$Conditie[j] < 4){ #If dead
      recVec[1+floor((as.Date(set$Datum[j])-startDate)/365.25)] <- 1  
    }
  }
  
  #Get sex information
  if("M" %in% set$Geslacht && "F" %in% set$Geslacht){
    sex[i] <- "U"
  }else if("M" %in% set$Geslacht){
    sex[i] <- "M"
  }else if("F" %in% set$Geslacht){
    sex[i] <- "F"
  }else{
    sex[i] <- "U"
  }
  
  #Get age information
  if(set$Ringleeftijd[1] %in% c(1,3)){ #Ringed as first-year or pullus
    age[i] <- "first-year"
  }else if(set$Ringleeftijd[1] %in% c(0,2)){ #
    age[i] <- "unknown"
  }else{
    age[i] <- "second-year+"
  }
  
  #Paste the recaptures and recoveries together ('LDLD' format as used in MARK)
  tempVec<-paste(seenVec,recVec,sep="")
  
  #Removes observations that are incomplete/wrong (mostly birds ringed before 2000 that are not incorporated here)
  if(!(1 %in% seenVec)){ #No initial capture event
    mistakes <- c(mistakes,i)
  }else if(1 %in% recVec){ 
    if(sum(recVec %like% 1)>1){ #multiple reported death events
      mistakes <- c(mistakes,i)
    }else if(which(recVec==1) < max(which(seenVec==1))){ #capture event after reported death
      mistakes <- c(mistakes,i)
    }
  }
  
 
  
  markVec[i] <- paste(tempVec,collapse="")
}

mark.df <- data.frame(capHis = markVec,Sex=sex,Age=age)
mark.df <- mark.df[-mistakes,] 

mark.df$Count <- rep(1,nrow(mark.df))

#Discard sex and age information
mark.df <- mark.df[,-c(2,3)]

#Tally the unique capture histories 
mark.df %>% group_by(capHis) %>% 
  summarise_all(sum) %>%
  data.frame() -> no_group

#Write to a .inp file for MARK, where the m-arrays are generated
output<-c()
for(i in 1:nrow(no_group)){
  output[i] <- paste(no_group$capHis[i],no_group$Count[i], ";",sep=" ")
}

#Load this into MARK and run the basic model. You can then retrieve the m-array.
write.table(output,file="S2000_nogroup.inp",quote=F,row.names=F)


