# Shrikant Pawar, Random Forest for methylation status analysis in renal cell cancer tissues
setwd("/Users/yalegenomecenter/Desktop/Brien Kidney Grant/Data/")
# Reading in the CSV file with methylation values
Cancer_Data <- read.csv(file = 'Cancer.csv', header=TRUE)
# Rename all columns with specific names
colnames(Cancer_Data) <- rep(c("id", "val"), 58)
# Make an empty matrix
m = matrix(, nrow = 58269240, ncol = 2)
# First rbind for stacking 1:4 columns in one
m <- rbind(Cancer_Data[,1:2], Cancer_Data[,3:4])
# For loop for stacking all columns
x <- seq(2,116,by=2)
y <- seq(1, by = 2, len = 58)
library(foreach)
foreach(val = y, val2 = x) %do% {
  m <- rbind(m, Cancer_Data[,val:val2])
}
# From the populated m matrix remove NA's 
Final <- m[,1:2][complete.cases(m[,1:2]),]
# Export Final object as CSV
write.csv(Final, file = "Final.csv", sep = ",", row.names = FALSE)
# Merging the tags for selected ids
module load R/3.6.1-foss-2018b

#sbatch bash_script.sh, Request --mem-per-cpu=99G and --ntasks=1 --nodes=10
#Following in myscript.R

setwd("/home/ps965/scratch60") 
Final <- read.csv(file = 'Final.csv', header=TRUE)                                                                                                                   
FinalID <- read.csv(file = 'Ids_Cancer.csv', header=FALSE)
colnames(FinalID) <- rep(c("id", "val"), 1)
Final_Merged <- merge(Final, FinalID, by="id", all = TRUE)
# Ids total and Ids with tags as disease: Same as number of patients
# Replaces NA's in third column with tag as Normal
library(forcats)
Final_Merged[,3] <- fct_explicit_na(Final_Merged$val.y, na_level = "Normal")
# Applying random forest on object final\
library(randomForest)
# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(100)
train <- sample(nrow(Final_Merged), 0.7*nrow(Final_Merged), replace = FALSE)
TrainSet <- Final_Merged[train,]
ValidSet <- Final_Merged[-train,]
summary(TrainSet)
summary(ValidSet)
# Creating new object with following structure
TrainSet <- na.omit(TrainSet)
TrainSet[!complete.cases(TrainSet),]
# Remove first column with ids as they are duplicate and R frame does not allow duplicate names, so can track back the ids latter.
TrainSet_New <- TrainSet[,2:3]
# Random forest model creation
# Memory allocation error so can't take all the data, data trimming needed
Disease <- TrainSet_New[TrainSet_New[,2]=="Disease",]
TrainSet_New_Trimmed <- rbind(TrainSet_New[1:10531,], Disease)
set.seed(100)
train <- sample(nrow(TrainSet_New_Trimmed), 0.99*nrow(TrainSet_New_Trimmed), replace = FALSE)
TrainSet <- TrainSet_New_Trimmed[train,]
ValidSet <- TrainSet_New_Trimmed[-train,]

model2 <- randomForest(val.y ~ ., data = TrainSet, na.action = na.pass, ntree = 500, mtry = 6, importance = TRUE)
print(model2)


# Predicting on Validation set
predValid <- predict(model2, ValidSet, type = "class")
# Checking classification accuracy
mean(predValid == ValidSet$val.y)                    
table(predValid,ValidSet$val.y)



####################Calculating sensitivity and specificity##########
The sensitivity is defined as the proportion of positive results out of the number of samples which were actually positive. When there are no positive results, sensitivity is not defined and a value of NA is returned. Similarly, when there are no negative results, specificity is not defined and a value of NA is returned. Similar statements are true for predictive values.
The positive predictive value is defined as the percent of predicted positives that are actually positive while the negative predictive value is defined as the percent of negative positives that are actually negative.
Suppose a 2x2 table with notation

Reference

Predicted
Event
No Event
Event
A
B
No Event
C
D

The formulas used here are:
Sensitivity = A/(A+C)
Specificity = D/(B+D)
Prevalence = (A+C)/(A+B+C+D)
PPV = (sensitivity * Prevalence)/((sensitivity*Prevalence) + ((1-specificity)*(1-Prevalence)))
NPV = (specificity * (1-Prevalence))/(((1-sensitivity)*Prevalence) + ((specificity)*(1-Prevalence)))
####################Calculating sensitivity and specificity##########