# There are four samples of yeast: wild type, a strain with the
# deletion of the gat1 gene, a strain with the deletion of the 
# gcn4 gene, and a strain with the deletion of the leu3 gene. 

#The objective of this program is to predict the expression of 
# particular genes in the gat1 deletion sample given the gene's 
# expression levels in the other three samples. 

# This program uses kmeans clustering to assign centers to the
# gene expression values in each sample except for the gat1 deletion
# sample. These centers are then used to create a random forest
# and make a numerical prediction for the gene expression levels
# in the gat1 deletion strain. Once these expression levels have
# been predicted, they are ranked. The Spearman correlation 
# determines how accurate the rankings are compared to the Gold
# Standard rankings provided on the challenge website. 

# Choose a time at which the samples were tested. There are 8 
# possible times: 0, 10, 20, 30, 45, 60, 90, or 120 seconds.
# You can also choose how many centers the kmeans algorithm uses.
# The default number of centers is 6.

# Load the dataset
setwd("~/Desktop/Bio/Binf/DREAM3 Gene Expression Prediction/DREAM3geneexpressionprediction")
expdata <- read.table("DREAM3_GeneExpressionChallenge_ExpressionData_UPDATED.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
library(dplyr)
expdata <- tbl_df(expdata)


# Let's set the seed to a specific number, so that whenever we 
# have to use this number, we can just use the variable "seed" 
# instead.
seed = 13

# MAKE IT INTO A FUNCTION
tm = 120 #default is 0
tmcoldefault <- c(4, 12, 20, 28) # these are the columns of the dataset to use if the default value of time is chosen
goldcoldefault <- 3 # this is the column of the gold standard rankings table to use if the default value of time is chosen
# the following conditional statement alters the columns of the dataframe to use, the column of the gold standard to use, and the number of centers for kmeans (if necessary) based on the time of the sample given
if (tm == 0) {
  tmcol <- tmcoldefault
  goldcol <- goldcoldefault
} else if (tm == 10) {
  tmcol <- tmcoldefault + 1
  goldcol <- goldcoldefault + 1
} else if (tm == 20) {
  tmcol <- tmcoldefault + 2
  goldcol <- goldcoldefault + 2
} else if (tm == 30) {
  tmcol <- tmcoldefault + 3
  goldcol <- goldcoldefault + 3
} else if (tm == 45) {
  tmcol <- tmcoldefault + 4
  goldcol <- goldcoldefault + 4
} else if (tm == 60) {
  tmcol <- tmcoldefault + 5
  goldcol <- goldcoldefault + 5
} else if (tm == 90) {
  tmcol <- tmcoldefault + 6
  goldcol <- goldcoldefault + 6
} else if (tm == 120) {
  tmcol <- tmcoldefault + 7
  goldcol <- goldcoldefault + 7
} else {
  print("Error: time not found")
}


# Make a dataset for the time specified. Columns = geneID, gene, wt, gat1, gcn4, and leu3
df <- expdata[ , c(1, 2, tmcol)]
names(df) <- c("ProbeID", "geneName", "wt", "gat1", "gcn4", "leu3")

# Divide dataset into train and test data. These datasets will contain numerical expression values for each gene that has been observed. 
test <- subset(df, gat1 == "PREDICT")
train <- subset(df, gat1 != "PREDICT")
train$gat1 <- as.numeric(train$gat1)

# Perform PCA on the training dataset each gene.
mat <- as.matrix(train[ , c(3:6)])
rownames(mat) <- train$ProbeID
pca <- prcomp(t(mat), scale = TRUE) 

# Select the highest scoring genes to include in the final training dataset
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
top_genes <- names(gene_score_ranked[1:4500]) # gives the top 4500 genes


# Use rbind to combine the rows of the output of PCA to the test dataset
top <- subset(train, ProbeID %in% top_genes)
t10 <- rbind(top, test)

# Use kmeans clustering to assign cluster numbers to each expression level for each gene in each sample
# use NbClust package to determine the optimal number of clusters for expression levels in each sample
library(NbClust)
if (tm == 0) {
  set.seed(seed)
  kmwt <- kmeans(as.matrix(t10$wt), centers = 1) # the wild type expression data at time = 0 seconds all have gene expression values of 0, so they cannot have multiple centers
} else {
  nbwt <- NbClust(as.matrix(t10$wt), 
                  distance = "euclidean", 
                  min.nc = 1, 
                  max.nc = 10, 
                  method = "complete", 
                  index = "ch") # this function determines the optimal number of clusters for the gene expression data of one sample (in this case, it's the wild type sample)
  nclustwt <- nbwt$Best.nc[1] # the optimal number of clusters is assigned to this variable, then used in the kmeans calculation for this sample
  set.seed(seed)
  kmwt <- kmeans(as.matrix(t10$wt), centers = nclustwt)
}
# each expression value is replaced with an assignment to a cluster from the kmeans algorithm
clusterwt <- kmwt$cluster


# this process is repeated for the other two samples that are not being predicted.
# None of the other samples have all 0's as their expression values, so a conditional (like the one used for the wild type data) is not necessary.
nbgcn4 <- NbClust(as.matrix(t10$gcn4), 
                distance = "euclidean", 
                min.nc = 1, 
                max.nc = 10, 
                method = "complete", 
                index = "ch")
nclustgcn4 <- nbgcn4$Best.nc[1]
set.seed(seed)
kmgcn4 <- kmeans(as.matrix(t10$gcn4), centers = nclustgcn4)
clustergcn4 <- kmgcn4$cluster



nbleu3 <- NbClust(as.matrix(t10$leu3), 
                  distance = "euclidean", 
                  min.nc = 1, 
                  max.nc = 10, 
                  method = "complete", 
                  index = "ch")
nclustleu3 <- nbleu3$Best.nc[1]
set.seed(seed)
kmleu3 <- kmeans(as.matrix(t10$leu3), centers = nclustleu3)
clusterleu3 <- kmleu3$cluster


# A dataframe is created that contains the Probe IDs, Gene names, cluster numbers for the predictor variables, and numerical expression values for the gat1 sample that is to be predicted.
clust <- tbl_df(t10[ , 1:2])
clust <- mutate(clust, 
                     wt = clusterwt, 
                     gcn4 = clustergcn4, 
                     leu3 = clusterleu3, 
                     gat1 = t10$gat1)

# The cluster numbers are converted into factors
clust[ , 3:5] <- lapply(clust[ , 3:5], as.factor)

# The dataset containing the cluster numbers as factors is divided into a testing and training dataset. 
clustertest <- subset(clust, gat1 == "PREDICT")
clustertrain <- subset(clust, gat1 != "PREDICT")
clustertrain$gat1 <- as.numeric(clustertrain$gat1)
clustertrain[ , 3:5] <- lapply(clustertrain[ , 3:5], droplevels)

# Large unnecessary data is removed from the environment
rm(expdata)
rm(kmgcn4)
rm(clustergcn4)
rm(kmleu3)
rm(clusterleu3)
rm(kmwt)
rm(clusterwt)
rm(clust)
rm(t10)
rm(mat)
rm(pca)

# The next step is to create a random forest using the training dataset and create a regression model that can predict expression values for various genes in the gat1 deletion sample based on the expression values of the same genes in the wild type, gcn4, and leu3 samples.
# The following packages are needed to construct the random forest:
library(caret)
library(randomForest)
library(ranger)


# In order to maximize the effectiveness of the random forest, a tuning grid search will be performed to determine the optimal parameters of the randome forest.
# The following tuning grid tests the optimal values for the number of trees, variables to consider at each split, the target size of the node, and the optimal sample size. The best combination of parameters will be the one with the lowest RMSE value. 
rfgrid <- expand.grid(numtree = c(2:20)*100, 
                      mtry = c(1, 2, 3), 
                      node_size = c(6:12)*10,
                      sample_size = c(0.55, 0.632, 0.7, 0.8), 
                      RMSE = 0)
# The tuning grid is tested by creating random forests for each combination of parameters. 
for (i in 1:nrow(rfgrid)) {
  set.seed(seed)
  rftune <- ranger(
    formula = gat1 ~ wt + gcn4 + leu3, 
    data = clustertrain, 
    num.trees = rfgrid$numtree[i],
    min.node.size = rfgrid$node_size[i], 
    mtry = rfgrid$mtry[i], 
    sample.fraction = rfgrid$sample_size[i]
  )
  rfgrid$RMSE[i] <- sqrt(rftune$prediction.error)
}

# The resulting tuning grid containing the RMSE values for each combination of parameters is sorted from least to greatest. The combination of parameters that minimized the RMSE is the combination that will be used to create the random forest for the predictive regression model. 
rfgrid <- rfgrid %>%
  arrange(RMSE)

# Create the model based on optimal paramters from the tuning grid
set.seed(seed = seed)
rfmodel <- ranger(formula = gat1 ~ wt + gcn4 + leu3, 
               data = clustertrain, 
               num.trees = rfgrid$numtree[1], 
               mtry = rfgrid$mtry[1], 
               min.node.size = rfgrid$node_size[1], 
               sample.fraction = rfgrid$sample_size[1]
               )


# Make predictions for the expression levels of the genes in the gat1 deletion sample
final <- clustertest[ , 2:6]
pred <- predict(rfmodel, data = clustertest)
final$gat1 <- pred$predictions


# Arrange the final dataset that contains the gene expression predictions to be in alphabetical order by gene_name, so that it may be compared to the Gold Standard rankings. 
final <- arrange(final, geneName)

# Load the Gold Standard rankings in the form of a dataframe 
goldstandard <- read.table("DREAM3GoldStandard_ExpressionChallenge.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
goldstandard <- tbl_df(goldstandard)

# There are 3 genes from the original dataset that are absent in the updated dataset, so the ranks from the Golden Standard will be adjusted.
# The absent genes are BAT1 (#7), DSE2 (#12), SPL2 (#39)
nogene <- c(7, 12, 39)
gs <- goldstandard [-nogene, c(1, 2, goldcol)] #this creates a dataframe that deletes the rows of absent genes and restricts the columns to the Probe IDs, gene names, and the specific time we are looking at. 
names(gs) <- c("probeID", "gene_name", "ranktime") # the third column, now named "ranktime", has been converted into a generic name so that this column can be manipulated no matter what time these rankings are for
gs <- arrange(gs, ranktime) # the genes are arranged by their rank
gs$ranktime <- 1:nrow(gs) # these ranks are reordered to reflect the changes that have been made by deleteing three genes. 
gs <- arrange(gs, gene_name) # the genes are then arranged by alphabetical order again so that they can be compared

# Predicted ranks of gene expression for this particular time in the gat1 deletion sample are compared to the Gold Standard rankings using Spearman correlation.
cor.test(gs$ranktime, order(final$gat1, decreasing = TRUE), method = "spearman")

# The predicted rankings are then saved in a vector based on which time these predictions represent
if (tm == 0) {
  time0 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 10) {
  time10 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 20) {
  time20 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 30) {
  time30 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 45) {
  time45 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 60) {
  time60 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 90) {
  time90 <- order(final$gat1, decreasing = TRUE)
} else if (tm == 120) {
  time120 <- order(final$gat1, decreasing = TRUE)
}

# Large unnecessary data is removed from the environment
rm(tmcol)
rm(goldcol)
rm(final)

# This is repeated 8 times: one for each time the gene expressions in the samples were measured







# Create a dataframe of the expression level rankings for each time tested
exprank <- as.data.frame(gs$gene_name) %>%
  mutate(time0 = time0, 
         time10 = time10, 
         time20 = time20, 
         time30 = time30, 
         time45 = time45, 
         time60 = time60, 
         time90 = time90, 
         time120 = time120) %>%
  tbl_df()

goldstandard <- goldstandard[ , 2:10] # the gold standard dataframe is now restricted to show the gene_name column and the expression level rankings for the gat1 deletion sample at each time it was tested.


# Create a new Gold Standard dataframe that modifies the rankings. This is because three genes were absent in the testing dataset.
gold <- data.frame(gene_name = gs$gene_name)
for (i in 1:8) {
  goldonly <- as.data.frame(goldstandard[ -nogene, c(1, i + 1)])
  names(goldonly) <- c("gene_name", "ranktime")
  goldonly <- arrange(goldonly, ranktime)
  goldonly$ranktime <- 1:nrow(goldonly)
  goldonly <- arrange(goldonly, gene_name)
  gold[ , i + 1] <- goldonly[ , 2]
  names(gold)[i + 1] <- colnames(goldstandard)[i + 1]
}
# The dataset gold will be used from now on as the Gold Standard dataset

# Calculate the overall gene-profile p-values for each time
pvalcol <- numeric(length = 8) # create an empty vector where the p-values will be stored. There are 8 times, meaning 8 p-values, so the vector should be numeric and of length 8.
exprankonly <- as.data.frame(exprank[ , 2:9])
for (i in 1:8) {
  pvalcol[i] <- cor.test(as.vector(gold[ , i + 1]), 
                         as.vector(exprankonly[ , i]), 
                         method = "spearman")$p.value
} # this for loop takes every column with predicted rankings and compares them to the gold standard rankings for all genes within one time
library(psych)
pg <- geometric.mean(pvalcol) # this variable stores the overall gene-profile p-value
print(paste("The overall gene-profile p-value is", geometric.mean(pvalcol)))


# Calculate overall time profile p values for each gene
pvalrow <- numeric(length = 47) # create an empty vector with 47 spaces, that will be filled with 47 p-values for the rankings of each gene
exprankonly <- as.data.frame(exprank)
for (i in 1:47) {
  pvalrow[i] <- cor.test(as.vector(unlist(gold[i , 2:9])), 
                         as.vector(unlist(exprank[i , 2:9])), 
                         method = "spearman")$p.value
} # this for loop takes every row (gene) with predicted rankings and compares them to the Gold Standard rankings of the same genes using Spearman correlation
library(psych)

py <- geometric.mean(pvalrow) # this variable stores the overall time profile p-value
print(paste("The overall time-profile p-value is", geometric.mean(pvalrow)))


# Calculate the overall score (a log transformed average of the overall gene-profile p-value and the overall time profile p-value)
-0.5 * log10(pg * py)
print(paste("The overall score is", -0.5 * log10(pg * py)))

