# Integrated Population Model of the Dutch Mallard population

This repository contains the data and code required to run the Mallard IPM as described by Wiegers et al. 2022. The same data can be found on Dryad as well. 

This dataset and code are the basis of the integrated population model for the Mallard in the Netherlands as described in Wiegers et al (2022).

The data consists of four .csv files, the integrated population model in WinBugs, and the R code required to run this model and the life-stage sensitivity analysis.

The rows of ‘Clutch size.csv’ are individual observations of Mallard nests. The columns indicate the observation number; the clutch size C; and the observation year, where the years 2003-2020 are indicated by Year = 1-18.

The rows of ‘Egg hatch rate.csv’ are individual observations of Mallard nests for which the number of hatched eggs was known. Only successful nests are used to exclude incompletely observed nest. The columns indicate the observation number; the clutch size C; the number of egg hatched; and the year. 

The rows of ‘Nest success.csv’ depict annual data on nest survival of Mallards in the Netherlands. Nest success was calculated using the Mayfield method. The columns indicate the years, where 2003 is year 1; the total number of days that nests were exposed to potential nest failure; the total number of nests monitored; the total number of nests that failed to produce at least one egg; the corresponding Mayfield estimate of daily nest survival; and the total nest success rate of female Mallards, assuming an average incubation time of 37 days.

The rows of ‘Duckling_project.csv’ are repeated observations of Mallard broods. The columns indicate the observation number; the number of ducklings counted during the first observation; the number of ducklings counted during the next observation; the interval (in days) between these observations; the mean estimated age (in days) of the ducklings at these observations; the observation year, where 2003 is year 1; and finally the scaled estimated duckling age. 



