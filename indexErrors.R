## R Script to investigate the effect of adding ordinal levels on Regression Dilution
## Author: Paul Scherer
## Institute: University of Cambridge MRC Epidemiology Unit 
## Date: 15.06.2017

# Set up the Data

## We pull raw PAEE data from a normal distribution representing the amount of Physical Activity a person does
## A study attempts to record this with a device that has a set amount of random measurement error, and then
## puts this measurement into one of n indices of PA, this is done for sample of raw data.

## the PAEE data has a set true linear relationship with some continous condition C. We know that the backtransformation
## of index values to the level of the PAEE variable introduces regression dilution due to the random measurement error
## introduced at the measurement stage. We investigate, all things kept constant, whether changing the number of n
## indices introduces more regression dilution to the correclation coefficient between the backtransformed activity level.

library(ggplot2)
library(reshape2)
library(parallel)

# Raw Data Set Properties
