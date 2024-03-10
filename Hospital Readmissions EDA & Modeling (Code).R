#' Author: Pablo Diaz
#' Data: Mar 22,2023
#' Purpose: patients EDA
#' 

# libs
library(caret)
library(corrr)
library(DataExplorer)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggcorrplot)
library(ggplot2)
library(ggthemes)
library(lubridate)
library(naniar)
library(tidyr)
library(vtreat)
library(leaflet)
library(scales)
library(RColorBrewer)
library(GGally)
library(tm)
library(qdapRegex)
library(stringr)
library(wordcloud)
library(ranger)
library(MLmetrics)
library(ROSE)
library(Boruta)
library(pROC)
library(psych)
library(randomForest)

# Set WD
setwd("~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/PersonalFiles")

# Options & Functions
options(stringsAsFactors = FALSE)
Sys.setlocale('LC_ALL','C')

# Define functions for NLP
tryTolower <- function(x){
  y         = NA
  try_error = tryCatch(tolower(x), error = function(e) e)
  if (!inherits(try_error, 'error')){}
  y         = tolower(x)
  return(y)
}

cleanCorpus <- function(corpus, customStopwords){
  corpus <- tm_map(corpus, content_transformer(qdapRegex::rm_url))
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)
  corpus <- tm_map(corpus, removeNumbers)
  corpus <- tm_map(corpus, content_transformer(tryTolower))
  corpus <- tm_map(corpus, removeWords, customStopwords)
  return(corpus)
}

bigramTokens <-function(x){
  unlist(lapply(NLP::ngrams(words(x), 2), paste, collapse = " "), 
         use.names = FALSE)
}


# Create custom stop words
customStopwords <- c(stopwords("english"), "and")

# Import the Data 
patientsTrain <- read.csv('~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/PersonalFiles/patientsTrain.csv')
patientsTest <- read.csv('~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/PersonalFiles/patientsTest.csv')

# Combine the datasets and identify which dataset each row belongs to
patientsCombined <- bind_rows(patientsTrain %>% mutate(dataset = "train"),
                              patientsTest %>% mutate(dataset = "test"))

# Trandsform values to NA
patientsCombined <- patientsCombined %>% mutate_all(~ ifelse(. %in% c("?","Not Available","", "Not Mapped", "NA"), NA_character_, .))

# Loop through each variable in the dataset to see value frequency
for (var in names(patientsCombined)) {
  freq_table <- table(patientsCombined[[var]])
  freq_table <- as.data.frame(freq_table)
  freq_table <- freq_table %>% 
    rename(n = Freq) %>% 
    mutate(value = as.character(Var1)) %>% 
    select(value, n) %>% 
    arrange(desc(n))
  print(paste0("Variable: ", var))
  print(head(freq_table, n = 20))
  print("")
}

# Find the columns with more than 1 unique value
uniq_cols <- sapply(patientsCombined, function(x) length(unique(x))) > 1

# Find the columns where the most frequent unique value has a frequency of more than 95% to drop
freq_cols <- sapply(patientsCombined[, uniq_cols], function(x) {
  freq_table <- table(x)
  max_freq <- max(freq_table) / sum(freq_table)
  max_freq < 0.95
})

# Keep only the selected columns, excluding "diag_2_desc" and "diag_3_desc"
selected_cols <- names(patientsCombined[, uniq_cols][, freq_cols])
selected_cols <- selected_cols[!selected_cols %in% c("diag_2_desc", "diag_3_desc")]
patientsCombined <- patientsCombined[, selected_cols]

# Create a vector of rows to drop, those patients are not likely to be readmitted 
values_to_drop <- c("Expired", "Hospice / home", "Hospice / medical facility")

# Drop rows containing values_to_drop
patientsCombined <- patientsCombined[!patientsCombined$discharge_disposition %in% values_to_drop, ]

# Group the remaining columns for transformation
numericCols <- c("time_in_hospital",
                 "num_lab_procedures",
                 "num_procedures",
                 "num_medications",
                 "number_outpatient",
                 "number_emergency",
                 "number_inpatient",
                 "number_diagnoses",
                 "age",
                 "wgt"
)

factorCols <- c("admission_type_id",
                "discharge_disposition_id",
                "admission_source_id",
                "medical_specialty",
                "diag_1_desc",
                "max_glu_serum",
                "A1Cresult",
                "metformin",
                "glipizide",
                "glyburide",
                "pioglitazone",
                "rosiglitazone",
                "insulin",
                "change",
                "diabetesMed",
                "race",
                "gender",
                "payer_code"
)

# Change the datatype of the previously groupped columns
patientsCombined <- patientsCombined %>%
  mutate_if(names(.) %in% numericCols, as.numeric) %>%
  mutate_if(names(.) %in% factorCols, as.factor)

# Create a function to remove outliers based in the box plots
outliers_remover <- function(a){
  df <- a
  aa <- c()
  count <- 1
  for(i in 2:ncol(df)){
    if(is.integer(df[,i])){
      Q3 <- quantile(df[,i], 0.75, na.rm = TRUE)
      Q1 <- quantile(df[,i], 0.25, na.rm = TRUE) 
      IQR <- Q3 - Q1
      upper <- Q3 + 1.5 * IQR
      lower <- Q1 - 1.5 * IQR
      for(j in 1:nrow(df)){
        if(is.na(df[j,i]) == TRUE){
          next
        }
        else if(df[j,i] > upper | df[j,i] < lower){
          aa[count] <- j
          count <- count+1                  
        }
      }
    }
  }
  if (length(aa) > 0) {
    df <- df[-aa,]
  }
  return(df)
}

# Apply the function to the dataset patientsCombined
patientsCombined <- outliers_remover(patientsCombined)

# Build a volatile corpus
txtCorpus <- VCorpus(VectorSource(patientsTrain$diag_1_desc))
head(patientsTrain$diag_1_desc)

# Preprocess the corpus
txtCorpus <- cleanCorpus(txtCorpus, customStopwords)
head(patientsTrain$diag_1_desc)

# Make DocT Matrix
diagDTM  <- DocumentTermMatrix(txtCorpus)
diagDTMm <- as.matrix(diagDTM)
dim(diagDTMm)

diagFreq <- colSums(diagDTMm)
diagFreq <- data.frame(word=names(diagFreq),
                        frequency=diagFreq, 
                        row.names = NULL)

# Simple barplot; values greater than 50 
topWords      <- subset(diagFreq, diagFreq$frequency >= 100) 
topWords      <- topWords[order(topWords$frequency, decreasing=F),]

# Chg to factor for ggplot
topWords$word <- factor(topWords$word, 
                        levels=unique(as.character(topWords$word))) 

topWords %>%
  top_n(50, frequency) %>%
  ggplot(aes(x=word, y=frequency)) + 
  geom_bar(stat="identity", fill='darkred') + 
  coord_flip()+ theme_gdocs() +
  geom_text(aes(label=frequency), colour="white",hjust=1.25, size=5.0) 


# Number missing values
colSums(is.na(patientsCombined))

# Count the number of rows with at least one NA to see if we can drop them
n_missing_rows <- sum(!complete.cases(patientsCombined))
n_missing_rows

# Scale the numeric columns to log (in case it improved the models), did not.
#patientsCombined[numericCols] <- log(patientsCombined[numericCols])

# Compute correlation matrix
cor_matrix <- cor(patientsCombined[numericCols], use="pairwise.complete.obs")

# Plot correlation matrix with ggcorrplot
ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 2.5, method="circle", 
           colors = c("#6D9EC1", "white", "#E46726"),
           title = "Correlation Plot")


# Split train and test data
set.seed(2021)
trainData   <- subset(patientsCombined, patientsCombined[, ncol(patientsCombined)] == "train")
testData    <- subset(patientsCombined, patientsCombined[, ncol(patientsCombined)] == "test")

# sample from train data to use in preparation plan
idxPrep     <- sample(1:nrow(trainData),.1*nrow(trainData))
prepData    <- trainData[idxPrep,]
nonPrepData <- trainData[-idxPrep,]

# select target and information variables
targetVar       <- names(prepData)[30]
informativeVars <- names(prepData)[1:29]

# Design a "C" Categorical variable plan 
plan <- designTreatmentsC(prepData, 
                          informativeVars,
                          targetVar,'TRUE')

# Partition to avoid overfitting
set.seed(2021)
idx <- sample(1:nrow(nonPrepData), 0.8 * nrow(nonPrepData))
train <- nonPrepData[idx, ]
validation <- nonPrepData[-idx, ]

# Now apply the variable treatment plan
treatedTrain <- prepare(plan, train)
treatedVal <- prepare(plan, validation)
treatedTest <- prepare(plan, testData)

table(treatedTrain$readmitted_y)
table(treatedVal$readmitted_y)
table(treatedTest$readmitted_y)

# ensure results are repeatable
set.seed(2021)

## Run the algorythm
#boruta <- Boruta(readmitted_y ~., data = treatedTrain, doTrace = 2)
#plot(boruta, las = 2, cex.axis = 0.5)
#attStats(boruta)
#boruta

#Tentative Fix
#bor <- TentativeRoughFix(boruta)
#print(bor)

## Save bor object as an .rds file in the specified directory
#saveRDS(bor, file = "~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/personalFiles/boruta_results.rds")

# Load bor object from saved file
bor <- readRDS(file = "~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/personalFiles/boruta_results.rds")

# Get the important variable names from bor object
importantVars <- getSelectedAttributes(bor, withTentative = FALSE)

# Subset the treatedTrain data to include only the important variables and response variable
importantTrainwoid <- treatedTrain[c(importantVars, "readmitted_y")]
importantTrain <- treatedTrain[c(importantVars, "readmitted_y","tmpID")]
importantVal   <- treatedVal[c(importantVars, "readmitted_y", "tmpID")]
importantTest  <- treatedTest[c(importantVars, "readmitted_y", "tmpID")]

# Use ROSE to rebalance the data (Just in case it improved the models) did not
#trainingRose <- ROSE(readmitted_y ~., importantTrain)$data
#treatedVal_rose <- ROSE(readmitted_y~ ., data = importantVal)$data
#treatedTest_rose <- ROSE(readmitted_y~ ., data = importantTest)$data

# 10 folds cross validation
trControl <- trainControl(method = "CV", number = 10)

# Model
# Fit the glm model with cross-validation
logitMod2_CV <- glm(readmitted_y ~., 
                    importantTrainwoid, 
                    family = 'binomial')

summary(logitMod2_CV)

# Make predictions on training and validation data
trainPreds <- predict(logitMod2_CV, importantTrain, type = 'response')
valPreds  <- predict(logitMod2_CV, importantVal, type = 'response')

# Classify 
cutoff       <- 0.5
trainClasses <- ifelse(trainPreds >= cutoff, 1, 0)

# Organize with actual
trainResults <- data.frame(importantTrain$tmpID,
                           actual = importantTrain$readmitted_y,
                           probablity = trainPreds,
                           classes = trainClasses)
head(trainResults)

# Homogenize the outcome for the confusion matrix
trainResults$actual <- ifelse(trainResults$actual == TRUE, "yes", "no")
trainResults$classes <- ifelse(trainResults$classes == 1, "yes", "no")

# Change columns as factor
trainResults$actual <- as.factor(trainResults$actual)
trainResults$classes <- as.factor(trainResults$classes)

# Create the confusion matrix
confMat <- confusionMatrix(trainResults$classes, trainResults$actual)
confMat

# Calculate accuracy
Accuracy(trainResults$classes, trainResults$actual)

# Perform backward stepwise regression to eliminate noise
set.seed(2021)
#backFit <- step(logitMod2_CV, direction = 'backward', trace = 5)
#saveRDS(backFit,'~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/personalFiles/backFit.rds')
#summary(backFit)

# Load the saved backFit object
backFit <- readRDS('~/Desktop/Hult/Hult_Visualizing-Analyzing-Data-with-R/personalFiles/backFit.rds')
backFitVars <- names(backFit$data)

# Train a new model with the selected variables
logitModSelected <- train(as.factor(readmitted_y) ~., 
                          data = importantTrainwoid[, backFitVars], 
                          method = "glm", 
                          family = "binomial", 
                          trControl = trControl) 
                      
# Make predictions on training and validation data
valPredsSelected  <- predict(logitModSelected, importantVal[, backFitVars])
valProbSelected <- predict(logitModSelected, importantVal[, backFitVars], type = 'prob')

# Classify 
cutoff               <- 0.6
ValReadmittedSelected <- ifelse(valProbSelected >= cutoff, 1, 0)

# Make predictions on the test set
trainResultsSelected <- data.frame(importantVal$tmpID,
                           actual  =importantVal$readmitted_y,
                           probablity = valProbSelected,
                           pred = ValReadmittedSelected)
head(trainResultsSelected)

# Homogenize the outcome for the confusion matrix
trainResultsSelected$actual <- ifelse(trainResultsSelected$actual == TRUE, "yes", "no")
trainResultsSelected$pred <- ifelse(trainResultsSelected$pred.TRUE == 1, "yes", "no")

# Change columns as factor
trainResultsSelected$actual <- as.factor(trainResultsSelected$actual)
trainResultsSelected$pred <- as.factor(trainResultsSelected$pred)

# Create the confusion matrix
confMatSelected <- confusionMatrix(trainResultsSelected$pred, trainResultsSelected$actual)
confMatSelected

# Calculate and compare the accuracy between base model and backFitted
Accuracy(trainResults$pred, trainResults$actual)
Accuracy(trainResultsSelected$pred, trainResultsSelected$actual)


# Fit decision tree model
set.seed(2021)
DTMod_CV <- train(as.factor(readmitted_y) ~ ., 
                  data = importantTrainwoid[, backFitVars], 
                  method = "rpart",
                  trControl = trControl)

# Predict with the validation data
DT_pred_CV <- predict(DTMod_CV, newdata = importantVal)

# Create a confusion matrix
DTConfMat <- confusionMatrix(as.factor(DT_pred_CV), as.factor(importantVal$readmitted_y))
DTConfMat

# Fit naive Bayes model
set.seed(2021)
NBMod_CV <- train(as.factor(readmitted_y) ~ ., 
                  data = importantTrainwoid[, backFitVars], 
                  method = "naive_bayes",
                  trControl = trControl)

# Predict with the validation data
NB_pred_CV <- predict(NBMod_CV, newdata = importantVal)

# Create a confusion matrix
NBConfMat <- confusionMatrix(as.factor(NB_pred_CV), as.factor(importantVal$readmitted_y))
NBConfMat

# Fit random forest model
set.seed(2021)
RFMod_CV <- train(as.factor(readmitted_y) ~ ., 
                  data = importantTrainwoid[, backFitVars], 
                  method = "rf",
                  trControl = trControl,
                  verbose = FALSE,
                  ntree = 500,
                  tuneGrid = data.frame(mtry = c(1,2)))

# Predict with the validation data
RF_pred_CV <- predict(RFMod_CV, newdata = importantVal)

# Create a confusion matrix
RFConfMat <- confusionMatrix(as.factor(RF_pred_CV), as.factor(importantVal$readmitted_y))
RFConfMat

# Fit neural networks model
set.seed(2021)
NNMod_CV <- train(as.factor(readmitted_y) ~ ., 
                  data = importantTrainwoid[, backFitVars], 
                  method = "nnet",
                  trControl = trControl)

# Predict with the validation data
NN_pred_CV <- predict(NNMod_CV, newdata = importantVal)

# Create a confusion matrix
NNConfMat <- confusionMatrix(as.factor(NN_pred_CV), as.factor(importantVal$readmitted_y))
NNConfMat

# Fit Gradient Boosting model
set.seed(2021)
GBMod_CV <- train(as.factor(readmitted_y) ~ ., 
                  data = importantTrainwoid[, backFitVars], 
                  method = "gbm",
                  trControl = trControl)

# Predict with the validation data
GBM_pred_CV <- predict(GBMod_CV, newdata = importantVal)

# Create a confusion matrix
GBMConfMat <- confusionMatrix(as.factor(GBM_pred_CV), as.factor(importantVal$readmitted_y))
GBMConfMat


# Combine the predictions of all models
combined_preds <- data.frame(LogisticRegression = valPredsSelected,
                             DecisionTree = DT_pred_CV,
                             NaiveBayes = NB_pred_CV,
                             RandomForest = RF_pred_CV,
                             NeuralNetworks = NN_pred_CV,
                             GradientBoosting = GBM_pred_CV)

# Create a list of models and predictions
model_list <- list(LR = logitModSelected, DT = DTMod_CV, RF = RFMod_CV, 
                   NB = NBMod_CV, NN = NNMod_CV, GB = GBMod_CV)

# generating resampled estimates of the performance of each of the models in the list
res <- resamples(model_list)
summary(res)

# Ploting the ROC Curve for all models
roc.curve(importantVal$readmitted_y, valPredsSelected, plotit = T, col = "blue")
roc.curve(importantVal$readmitted_y, DT_pred_CV, plotit = T, add.roc = T, col = "darkred")
roc.curve(importantVal$readmitted_y, NB_pred_CV, plotit = T, add.roc = T, col = "darkorange")
roc.curve(importantVal$readmitted_y, RF_pred_CV, plotit = T, add.roc = T, col = "darkgreen")
roc.curve(importantVal$readmitted_y, NN_pred_CV, plotit = T, add.roc = T, col = "gold")
roc.curve(importantVal$readmitted_y, GBM_pred_CV, plotit = T, add.roc = T, col = "green")

legend(.8, .7, legend = c("LG", "DT", "NB", "RF","NN", "GBM"),
       col = c("blue", "darkred", "darkorange", "darkgreen", "gold", "green"),
      lty = c(1,2,3,4,5,6), ncol = 1)

# box-and-whisker plot for all models (median and interquartile)
bwplot(res)

# Plotting the variable importance for the LR model
ggplot(varImp(logitModSelected), top = 20) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  ggtitle("Variable Importance Plot") +
  xlab("Variable") +
  ylab("Importance")

# # Plotting the variable importance for the RF model
ggplot(varImp(RFMod_CV), top = 20) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  ggtitle("Variable Importance Plot") +
  xlab("Variable") +
  ylab("Importance")


# predict based on the Test Data
RFTestpred_CV <- predict(RFMod_CV, newdata = importantTest)
RFTestprob_CV <- predict(RFMod_CV, newdata = importantTest, type = 'prob')

# Create a DataFrame with the results
RFResultsSelected <- data.frame(importantTest$tmpID,
                                actual = importantTest$readmitted_y,
                                probablity = RFTestprob_CV,
                                pred = RFTestpred_CV)
# Sort the values by probability
RFProbSelected_sorted <- RFResultsSelected[order(RFResultsSelected[,4], decreasing = TRUE),]
head(RFProbSelected_sorted)

# Select the highest 100
top100 <- head(RFProbSelected_sorted, 100)
tail(top100)
# Change the name of the column to merge
colnames(top100)[which(colnames(top100) == "importantTest.tmpID")] <- "tmpID"

# Merge with the base DataFrame to inspect the patients
importantTest_top100 <- merge(patientsCombined, top100, by = "tmpID", all = FALSE)

# Create plots with the different variables
ggplot(importantTest_top100,aes(x=num_lab_procedures,fill="lab proc"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=5)+theme_bw()

ggplot(importantTest_top100,aes(x=num_lab_procedures,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=5)+theme_bw()

ggplot(importantTest_top100,aes(x=age,fill="age"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=5)+theme_bw()

ggplot(importantTest_top100,aes(x=age,group=readmitted_y,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=5)+theme_bw()

ggplot(importantTest_top100,aes(x=num_medications,fill="medications"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=2)+theme_bw()

ggplot(importantTest_top100,aes(x=num_medications,group=readmitted_y,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=2)+theme_bw()

ggplot(importantTest_top100,aes(x=number_diagnoses,fill="diagnoses"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

ggplot(importantTest_top100,aes(x=number_diagnoses,group=readmitted_y,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

ggplot(importantTest_top100,aes(x=number_inpatient,fill="inpatient"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

ggplot(importantTest_top100,aes(x=number_inpatient,group=readmitted_y,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

ggplot(importantTest_top100,aes(x=time_in_hospital,fill="time hosp"))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

ggplot(importantTest_top100,aes(x=time_in_hospital,group=readmitted_y,fill=readmitted_y))+
  geom_histogram(position="identity",alpha=0.8,binwidth=1)+theme_bw()

# Build a volatile corpus
txtCorpus <- VCorpus(VectorSource(importantTest_top100$diag_1_desc))
head(importantTest_top100$diag_1_desc)
# Preprocess the corpus
txtCorpus <- cleanCorpus(txtCorpus, customStopwords)
head(patientsTrain$diag_1_desc)

# Make TDM
diagDTM  <- DocumentTermMatrix(txtCorpus)
diagDTMm <- as.matrix(diagDTM)

diagFreq <- colSums(diagDTMm)
diagFreq <- data.frame(word=names(diagFreq),
                       frequency=diagFreq, 
                       row.names = NULL)

# Simple barplot; values greater than 50 
topWords      <- subset(diagFreq, diagFreq$frequency >= 6) 
topWords      <- topWords[order(topWords$frequency, decreasing=F),]

# Chg to factor for ggplot
topWords$word <- factor(topWords$word, 
                        levels=unique(as.character(topWords$word))) 

topWords %>%
  top_n(50, frequency) %>%
  ggplot(aes(x=word, y=frequency)) + 
  geom_bar(stat="identity", fill='steelblue') + 
  coord_flip()+ theme_gdocs() +
  geom_text(aes(label=frequency), colour="white",hjust=1.25, size=5.0) 

# See a bi-gram
idx <- grep('failure', colnames(diagDTMm))
diagDTMm[,idx]

# Get Row Sums & organize
diagTDMmVec <- sort(colSums(diagDTMm), decreasing = TRUE)
wordFreqDF   <- data.frame(word      = names(diagTDMmVec), 
                           freq      = diagTDMmVec, 
                           row.names = NULL)

# Choose a color & drop light ones
pal <- brewer.pal(8, "Blues")
pal <- pal[-(1:2)]

# Make simple word cloud
# Reminder to expand device pane
wordcloud(wordFreqDF$word,
          wordFreqDF$freq,
          max.words=50,
          random.order=FALSE,
          colors=pal)

top100
write.csv(top100[, c("tmpID", "probablity.TRUE")], "top100PatientsAndProbs.csv", row.names = FALSE)
