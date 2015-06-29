library(protr)
library(readxl)
library(caret)
library(Biostrings)
library(RWeka)
library(dplyr)
library(e1071)
library(randomForest)
library(nnet)
library(Rcpi)
library(RCurl)

setwd("/Volumes/SAW SIMEON/PPI")

x <- getURL("https://github.com/Rnewbie/PPI/blob/master/data.xlsx")

data <- read_excel("data.xlsx")
## Subset the data with only reported pKd
data <- subset(data, pKd != "NA")
seq <- data$Protein1_Sequence
names(seq) <- data$Protein1_PDB
x <- AAStringSet(seq)
writeXStringSet(x, file = "protein1.fasta", width = 80)

seq2 <- data$Protein2_Sequence
names(seq2) <- data$Protein2_PDB
y <- AAStringSet(seq2)
writeXStringSet(y, file = "protein2.fasta", width = 80)

#### EMAIL From BOSS Last AUGUST 19, 2014 


####Dear Champ,
###
##This is a data containing protein and protein interaction. Please generate "Composition, Transition, Distribution" descriptors for each protein and compute all possible cross-terms (13 models like our previous GFP project).
###
###Best regards,
###Chanin
###
###Saw,
###
###I re-formatted your excel file into the attached data set.
###
###Best regards,
##Chanin
##Associate Professor Chanin Nantasenamat, Ph.D.
##Head, Center of Data Mining and Biomedical Informatics
##Faculty of Medical Technology, Mahidol University
##999 Phutthamonthon 4 Road, Salaya, Nakhon Pathom 73170
### OK this will be real now
### 

## extractCTDC() Composition  extractCTDT(x) Transition extractCTDD(x) Distribution

## bug fixed by koefoed for calculation of protien descriptors (Distribution) for protr R package
extractCTDD = function (x) {
  
  if (protcheck(x) == FALSE) stop('x has unrecognized amino acid type')
  
  group1 = list(hydrophobicity  = c('R', 'K', 'E', 'D', 'Q', 'N'),
                normwaalsvolume = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
                polarity        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
                polarizability  = c('G', 'A', 'S', 'D', 'T'),
                charge          = c('K', 'R'),
                secondarystruct = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
                solventaccess   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
  
  group2 = list(hydrophobicity  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
                normwaalsvolume = c('N', 'V', 'E', 'Q', 'I', 'L'),
                polarity        = c('P', 'A', 'T', 'G', 'S'),
                polarizability  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
                charge          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 
                                    'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
                secondarystruct = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
                solventaccess   = c('R', 'K', 'Q', 'E', 'N', 'D'))
  
  group3 = list(hydrophobicity  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
                normwaalsvolume = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
                polarity        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
                polarizability  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
                charge          = c('D', 'E'),
                secondarystruct = c('G', 'N', 'P', 'S', 'D'),
                solventaccess   = c('M', 'S', 'P', 'T', 'H', 'Y'))
  
  xSplitted = strsplit(x, split = '')[[1]]
  n  = nchar(x)
  
  G = vector('list', 7)
  for (i in 1:7) G[[i]] = rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:7) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- 'G1')
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- 'G2')
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- 'G3')
  }
  
  # Compute Distribution
  
  D = vector('list', 7)
  for (i in 1:7) D[[i]] = matrix(ncol = 5, nrow = 3)
  
  for (i in 1:7) {
    inds = which(G[[i]] == 'G1')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][1, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G2')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][2, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
    inds = which(G[[i]] == 'G3')
    quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
    D[[i]][3, ] = ifelse(length(inds)>0,
                         (inds[c(1, ifelse(quartiles>0, quartiles, 1), length(inds))])*100/n,
                         0)
  }
  
  D = do.call(rbind, D)
  D = as.vector(t(D))
  
  names(D) = paste(rep(paste('prop', 1:7, sep = ''), each = 15),
                   rep(rep(c('.G1', '.G2', '.G3'), each = 5), times = 7),
                   rep(paste('.residue', c('0', '25', '50', '75', '100'), 
                             sep = ''), times = 21), sep = '')
  
  return(D)
  
}

protein_A <- readFASTA("protein1.fasta")
protein_B <- readFASTA("protein2.fasta")
a <- protein_A[(sapply(protein_A, protcheck))]
b <- protein_B[(sapply(protein_B, protcheck))]
bossrequestedDes <- function(x) {
  c(extractCTDC(x), ### composition
    extractCTDT(x), ### transition
    extractCTDD(x)) ### distribution
}
des_A <- t(sapply(a, bossrequestedDes))
protein_A_column_names <- colnames(des_A)
label_names_protein_A <- paste("Protein_A", protein_A_column_names)
colnames(des_A) <- label_names_protein_A
des_B <- t(sapply(b, bossrequestedDes))
protein_B_column_names <- colnames(des_B)
label_names_protein_B <- paste("Protein_B", protein_B_column_names)
colnames(des_B) <- label_names_protein_B
rownames(des_A) <- NULL
rownames(des_B) <- NULL

affinity <- data$Affinity
affinity <- as.factor(affinity)
set.seed(299)
protein_a_cor <- cor(des_A, use="complete.obs")
dim(protein_a_cor)
cor_protein_a <- protein_a_cor[1:147, 1:147]
high_Corr_protein_a <- findCorrelation(cor_protein_a, cutoff = 0.90)
length(high_Corr_protein_a)
highCorrRemoveProtein_A <- des_A[, -high_Corr_protein_a]
input_protein_a <- highCorrRemoveProtein_A
input_partition_protein_a <- cbind(affinity, data.frame(input_protein_a))
set.seed(300)
protein_b_cor <- cor(des_B, use = "complete.obs")
dim(protein_b_cor)
cor_protein_b <- protein_b_cor[1:147, 1:147]
high_Corr_protein_b <- findCorrelation(cor_protein_b, cutoff = 0.90)
length(high_Corr_protein_b)
highCorrRemoveProtein_B <- des_B[, -high_Corr_protein_b]
input_protein_b <- highCorrRemoveProtein_B
inputpartitionprotein_b <- cbind(affinity, data.frame(input_protein_b))
### preparing for strafied sampling
protein_A_high <- subset(input_partition_protein_a, affinity == "High")
protein_A_low <- subset(input_partition_protein_a, affinity =="Low")
protein_B_high <- subset(inputpartitionprotein_b, affinity == "High")
protein_B_low <- subset(inputpartitionprotein_b, affinity == "Low")
### protein alone data
protein_A_data <- rbind(protein_A_high, protein_A_low)
protein_B_data <- rbind(protein_B_high, protein_B_low)
protein_A_data <-  protein_A_data[,! apply(protein_A_data, 2, function(x) any(is.na(x)))]
protein_B_data <- protein_B_data[,! apply(protein_B_data, 2, function(x) any(is.na(x)))]
affinity <- protein_A_data$affinity
### cross terms data
p_A_data <- protein_A_data[-1]
p_B_data <- protein_B_data[-1]
crossTerms <- getCPI(p_A_data, p_B_data, type = "tensorprod")
crossTerms <- as.data.frame(crossTerms)
df_protein_A <- names(data.frame(protein_A_data[,2:49]))
df_protein_B <- names(data.frame(protein_B_data[,2:53]))
protein_A_Namecross <- rep(df_protein_A, each=52)
Protein_B_Namecross <- rep(df_protein_B,times=48)
label <- paste(protein_A_Namecross, Protein_B_Namecross, sep="_")
colnames(crossTerms) <- label
cross_terms_data <- cbind(affinity, crossTerms)
### Protein A self cross terms
selfcross_protein_A <- getCPI(p_A_data, p_A_data, type = "tensorprod")
df <- names(data.frame(protein_A_data[, 2:49]))
protein_A_name2 <- rep(df, times = 48)
protein_A_name1 <- rep(df, each = 48)
label <- paste(protein_A_name1, protein_A_name2, sep = "_")
colnames(selfcross_protein_A) <- label
selfcross_protein_A <- data.frame(selfcross_protein_A)
index <- seq(1, 2304, by = 49)
selfcross_protein_A <- selfcross_protein_A[, -index]
transposedIndex_protein_A <- t(selfcross_protein_A)
index1 <- which(duplicated(transposedIndex_protein_A))
removed_protein_A_train <- transposedIndex_protein_A[-index1, ]
protein_A_finalselfcrosstrain <- t(removed_protein_A_train)
protein_A_finalselfcrosstrain <- data.frame(protein_A_finalselfcrosstrain)
AxA_data <- cbind(affinity, protein_A_finalselfcrosstrain)
## selfcross B self cross terms
selfcross_protein_B <- getCPI(p_B_data, p_B_data, type = "tensorprod")
df <- names(data.frame(protein_B_data[, 2:53]))
protein_B_name2 <- rep(df, times = 52)
protein_B_name1 <- rep(df, each = 52)
label <- paste(protein_B_name1, protein_B_name2, sep = "_")
colnames(selfcross_protein_B) <- label
selfcross_protein_B <- data.frame(selfcross_protein_B)
index2 <- seq(1, 2704, by = 53)
selfcross_protein_B <- selfcross_protein_B[, -index2]
transposedIndex_protein_B <- t(selfcross_protein_B)
index3 <- which(duplicated(transposedIndex_protein_B))
removed_protein_B_train <- transposedIndex_protein_B[-index3, ]
protein_B_finalselfcrosstrain <- t(removed_protein_B_train)
protein_B_finalselfcrosstrain <- data.frame(protein_B_finalselfcrosstrain)
BxB_data <- cbind(affinity, protein_B_finalselfcrosstrain)
## data input

protein_A_data
protein_B_data
AxA_data
BxB_data
A_B_data <- cbind(affinity, p_A_data, p_B_data)
A_B_AxB_data_block_scale <- cbind(p_A_data, p_B_data, crossTerms) * (1/sqrt(length(p_A_data)+length(p_B_data)+length(crossTerms)))
A_B_AxB_data <- cbind(affinity, A_B_AxB_data_block_scale)
A_B_AxA_data_block_scale <- cbind(p_A_data, p_B_data, 
                                  protein_A_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+length(protein_A_finalselfcrosstrain)))
A_B_AxA_data <- cbind(affinity, A_B_AxA_data_block_scale)
A_B_BxB_data_block_scale <- cbind(p_A_data, p_B_data,
                                  protein_B_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+length(protein_B_finalselfcrosstrain)))
A_B_BxB_data <- cbind(affinity, A_B_BxB_data_block_scale)
A_B_AxB_AxA_data_block_scale <- cbind(p_A_data, p_B_data, crossTerms,
                                      protein_A_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+length(crossTerms)+length(protein_A_finalselfcrosstrain)))
A_B_AxB_AxA_data <- cbind(affinity, A_B_AxB_AxA_data_block_scale)
A_B_AxB_BxB_data_block_scale <- cbind(p_A_data, p_B_data, crossTerms,
                                      protein_B_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+length(crossTerms)+length(protein_B_finalselfcrosstrain)))
A_B_AxB_BxB_data <- cbind(affinity, A_B_AxB_BxB_data_block_scale)
A_B_AxA_BxB_data_block_scale <- cbind(p_A_data, p_B_data, protein_A_finalselfcrosstrain,
                                      protein_B_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+
                                                                                 length(protein_A_finalselfcrosstrain)+length(protein_B_finalselfcrosstrain)))
A_B_AxA_BxB_data <- cbind(affinity, A_B_AxA_BxB_data_block_scale)
A_B_AxB_AxA_BxB_data_block_scale <- cbind(p_A_data, p_B_data, crossTerms,
                                          protein_A_finalselfcrosstrain,
                                          protein_B_finalselfcrosstrain) * (1/sqrt(length(p_A_data)+length(p_B_data)+
                                                                                     length(protein_A_finalselfcrosstrain)+
                                                                                     length(protein_B_finalselfcrosstrain)))
A_B_AxB_AxA_BxB_data <- cbind(affinity, A_B_AxB_AxA_BxB_data_block_scale)

### data set name
#protein_A_data
#protein_B_data
#cross_terms_data
#BxB_data
#A_B_data 
#A_B_AxB_data 
#A_B_AxA_data
#A_B_BxB_data
#A_B_AxB_AxA_data
#A_B_AxB_BxB_data 
#A_B_AxA_BxB_data
#A_B_AxB_AxA_BxB_data

### function for training J48 training


J48_training <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

### mean and SD value

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

J48_train <- function(x) {
  ok <- J48_training(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

#protein_A_data
#protein_B_data
#cross_terms_data
#AxA_data
#BxB_data
#A_B_data 
#A_B_AxB_data 
#A_B_AxA_data
#A_B_BxB_data
#A_B_AxB_AxA_data
#A_B_AxB_BxB_data 
#A_B_AxA_BxB_data
#A_B_AxB_AxA_BxB_data

#### results for training J48
protein_A_training <- J48_train(protein_A_data)
protein_B_training <- J48_train(protein_B_data)
protein_cross_terms_training <- J48_train(cross_terms_data)
protein_AxA_training <- J48_train(AxA_data)
protein_BxB_training <- J48_train(BxB_data)
protein_A_B_training <- J48_train(A_B_data)
protein_A_B_AxB_training <- J48_train(A_B_AxB_data)
protein_A_B_AxA_training <- J48_train(A_B_AxA_data)
protein_A_B_BxB_training <- J48_train(A_B_BxB_data)
protein_A_B_AxB_AxA_training <- J48_train(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_training <- J48_train(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_training <- J48_train(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_training <- J48_train(A_B_AxB_AxA_BxB_data)

### function for 10 fold cross validation

J48_10fold <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    eval_j48 <- evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_cross_validation <- function(x) {
  ok <- J48_10fold(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

#### results for 10 fold cross validation

protein_A_cross_validation <- J48_cross_validation(protein_A_data)
protein_B_cross_validation <- J48_cross_validation(protein_B_data)
protein_crossterms_cross_validation <- J48_cross_validation(cross_terms_data)
protein_AxA_cross_validation <- J48_cross_validation(AxA_data)
protein_BxB_cross_validation <- J48_cross_validation(BxB_data)
protein_A_B_cross_validation <- J48_cross_validation(A_B_data)
protein_A_B_AxB_cross_validation <- J48_cross_validation(A_B_AxB_data)
protein_A_B_AxA_cross_validation <- J48_cross_validation(A_B_AxA_data)
protein_A_B_BxB_cross_validation <- J48_cross_validation(A_B_BxB_data)
protein_A_B_AxB_AxA_cross_validation <- J48_cross_validation(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_cross_validation <- J48_cross_validation(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_cross_validation <- J48_cross_validation(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_cross_validation <- J48_cross_validation(A_B_AxB_AxA_BxB_data)



#### J48 modeling testing results
J48_testing <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    model_train <- J48(affinity~., data = train)
    eval_external <- evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_external <- function(x) {
  ok <- J48_testing(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

### results for testing set


protein_A_testing <- J48_external(protein_A_data)
protein_B_testing <- J48_external(protein_B_data)
protein_crossterms_testing <- J48_external(cross_terms_data)
protein_AxA_testing <- J48_external(AxA_data)
protein_BxB_testing <- J48_external(BxB_data)
protein_A_B_testing <- J48_external(A_B_data)
protein_A_B_AxB_testing <- J48_external(A_B_AxB_data)
protein_A_B_AxA_testing <- J48_external(A_B_AxA_data)
protein_A_B_BxB_testing <- J48_external(A_B_BxB_data)
protein_A_B_AxB_AxA_testing <- J48_external(A_B_AxB_AxA_data)
protein_A_B_AxB_BxB_testing <- J48_external(A_B_AxB_BxB_data)
protein_A_B_AxA_BxB_testing <- J48_external(A_B_AxA_BxB_data)
protein_A_B_AxB_AxA_BxB_testing <- J48_external(A_B_AxB_AxA_BxB_data)




RF_training <- function(x) {
  high <- subset(x, affinity == "High")
  low <- subset(x, affinity =="Low")
  results <- list(100)
  for (i in 1:100) {
    train_high <- sample_n(high, size = 54)
    test_high <- sample_n(high, size = 14)
    train_low <- sample_n(low, size = 51)
    test_low <- sample_n(low, size = 13)
    train <- rbind(train_high, train_low)
    test <- rbind(test_high, test_low)
    ctrl <- trainControl(method = "cv",  number=5,savePredictions=FALSE)
    rf <- train(affinity~., data = train,  method = "rf", trControl = ctrl)
    model_train <- randomForest(affinity~., data = train, mtry = rf$bestTune[[1]], ntree = 500)
    confusionmatrix <- model_train$confusion[, 1:2]
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

### mean and SD value

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

RF_train <- function(x) {
  ok <- RF_training(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}

random_forest_protein_A_training <- RF_train(protein_A_data)
random_forest_protein_B_training <- RF_train(protein_B_data)
random_forest_protein_cross_terms_training <- RF_train(cross_terms_data)
random_forest_protein_AxA_training <- RF_train(AxA_data)
random_forest_protein_BxB_training <- RF_train(BxB_data)
random_forest_protein_A_B_training <- RF_train(A_B_data)
random_forest_protein_A_B_AxB_training <- RF_train(A_B_AxB_data)
random_forest_protein_A_B_AxA_training <- RF_train(A_B_AxA_data)
random_forest_protein_A_B_BxB_training <- RF_train(A_B_BxB_data)
random_forest_protein_A_B_AxB_AxA_training <- RF_train(A_B_AxB_AxA_data)
random_forest_protein_A_B_AxB_BxB_training <- RF_train(A_B_AxB_BxB_data)
random_forest_protein_A_B_AxA_BxB_training <- RF_train(A_B_AxA_BxB_data)
random_forest_protein_A_B_AxB_AxA_BxB_training <- RF_train(A_B_AxB_AxA_BxB_data)


