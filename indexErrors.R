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
library(stats)

# Raw Data Set Properties
trueBeta = 0.5
constant = 0

study_size = 20000
rawMean = 50
rawStdDev = 15

# Simple creation of the raw data and the related condition C without error
rawData = data.frame(paee = rnorm(n = study_size, mean = rawMean, sd = rawStdDev))
rawData$C = (trueBeta*rawData$paee) + constant

###############################################################################
########################### Functions #########################################
###############################################################################

createStudyData <- function(raw_data, measurement_error, number_of_indices){
    paee_errors = rnorm(n = study_size, mean = 0, sd = measurement_error)
    raw_data$paee_error = raw_data$paee + paee_errors
    raw_data$index = cut_number(x = raw_data$paee_error,n = number_of_indices, labels =FALSE)
    # raw_data$index = kmeans(x = raw_data$paee_error, centers = number_of_indices, iter.max = 100)
    raw_data = raw_data[with(raw_data, order(index)),]
    return(raw_data)
}

createValidationData <- function(val_size, measurement_error, number_of_indices) {
    validation_data = data.frame(paee =  rnorm(n = val_size, mean = rawMean, sd = rawStdDev))
    paee_errors = rnorm(n = val_size, mean = 0, sd = measurement_error)
    validation_data$paee_error = validation_data$paee + paee_errors
    validation_data$index = cut_number(x = validation_data$paee_error,n = number_of_indices, labels =FALSE)
    return (validation_data)
}

absDiff <- function(x,y){
    return (abs(x-y))
}

createMeansList <- function(data, number_of_indices) {
    meansList = vector(mode="list", length = number_of_indices)
    for (i in 1:number_of_indices) {
        meansList[i] = mean(unname(unlist((split(x=data$paee, f= as.factor(data$index)))[i])))
    }
    return(meansList)
}

probabilisticIndexer <- function(study_means_list, value){
    # This function provides an index for a value given a list of numbers. The euclidean distance between
    # the value and each member of the study_means_list is taken and then the value is given a weighted choice
    # by distance.
    distances_from_means = lapply(study_means_list, value, FUN = absDiff)
    distances_normed = distances_from_means/sum(distances_from_means)
    # to make the smaller distance gain the higher probabilities we will subtract 1 by the normed value
    probabilities = rep(1, length(distances_normed)) - distances_normed
    index_val = sample(1:length(study_means_list, size = 1, replace = TRUE, prob = probabilities)
    return(index_val)
}

###############################################################################
########################### DATA AND SETTINGS #################################
###############################################################################

measureError = 10
numLevels = 4

levels = c(2:120)

numTrials = 1000
betasData = data.frame(trial = c(1:numTrials))
stdErrorData = data.frame(trial = c(1:numTrials))

for (numLevels in levels) {
    betas <- vector("numeric", length = numTrials)
    std_errs <- vector("numeric", length = numTrials)
    for (i in 1:numTrials) {
        studyData = createStudyData(raw_data = rawData, measurement_error = measureError, number_of_indices = numLevels)
        validation_data = createValidationData(val_size = 200, measurement_error = measureError, number_of_indices = numLevels )

        validation_means = createMeansList(validation_data, numLevels)

        # find our estimated beta using the index_means
        studyData$ind_mean <- unlist(lapply(X=studyData$index, FUN=function(index_val){
            output =  validation_means[index_val]
            }))
        reg_ind_mean <- lm(formula=C~ind_mean, data=studyData)
        estimated_beta = reg_ind_mean$coefficients['ind_mean']
        standardError_reg_ind = summary(reg_ind_mean)$coefficients["ind_mean","Std. Error"]

        betas[i] = estimated_beta
        std_errs[i] = standardError_reg_ind
    }
    betasData = cbind(betasData, betas)
    stdErrorData = cbind(stdErrorData, std_errs)
}
columnnames = c("trial", levels)
colnames(betasData) = columnnames
colnames(stdErrorData) = columnnames

# To display the means of the betas per level
means_vector = lapply(X = betasData, FUN = mean)
means_df = data.frame(means_vector)
colnames(means_df) = columnnames
means_df = means_df[2:length(levels)+1]

# To graph the means of the betas per level
means_vector_y = means_vector[2:length(means_vector)]
xs = levels
plot(x = xs, y = means_vector_y)
lines(lowess(xs,means_vector_y), col="blue") # lowess line (x,y)


# To display the stderrors of the betas per level
stdErrorMeans = lapply(X = stdErrorData, FUN = mean)
stdError_df = data.frame(stdErrorMeans)
colnames(stdError_df) = columnnames
stdError_df = stdError_df[2:length(levels)+1]

# To graph the stderrors of the betas per level
stdErrorMeans_y = stdErrorMeans[2:length(stdErrorMeans)]
plot(x = xs, y = stdErrorMeans_y)
lines(lowess(xs,stdErrorMeans_y), col="blue") # lowess line (x,y)






