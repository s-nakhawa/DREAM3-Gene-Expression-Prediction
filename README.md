# DREAM3 Challenge: Gene Expression Prediction
An algorithm to predict gene expression in a particular sample of yeast.


## Background
GAT1, GCN4, and LEU3 are yeast transcription factors. Each of these transcription factors has something to do with controlling genes involved in nitrogen or amino acid metabolism. The genes are not essential because strains that have perfect deletions of any of these genes are viable. In this challenge, we provide gene expression data from four strains: (i) a strain that is wild-type for all three transcription factors (wt, or parental), (ii) a strain that is identical to the parental strain except that it has a deletion of the GAT1 gene (gat1Δ), (iii) a strain that is identical to the parental strain except that it has a deletion of the GCN4 gene (gcn4Δ), and (iv) a strain that is identical to the parental strain except that it has a deletion of the the LEU3 gene (leu3Δ).

Expression levels were assayed separately in all four strains following the addition of 3-aminotriazole (3AT). 3AT is an inhibitor of an enzyme in the histidine biosynthesis pathway and, in the appropriate media (which is the case in these experiments) inhibition of the histidine biosynthetic pathway has the effect of starving the cells for this essential amino acid.

Data from eight time points was obtained from 0 to 120 minutes. Time t=0 means the absence of 3AT.


## The Challenge
Predict, for a set of 50 genes, the expression levels in the gat1Δ strain in the absence of 3-aminotriazole (t=0) and at 7 time points ( t=10, 20, 30, 45, 60, 90 and 120 minutes) following the addition of 3AT. Absolute expression levels are not required or desired; instead, the fifty genes should be ranked according to relative induction or repression relative to the expression levels observed in the wild-type parental strain in the absence of 3AT.


## The Data
The file DREAM3_GeneExpressionChallenge_TargetList.txt is a tab-delimited file that lists the target genes whose relative induction/repression are to be predicted. The first column lists the Affymetrix probeset IDs. The second column lists the corresponding commonly-used gene names, as extracted from files obtained from Affymetrix. This file should also be used as a template for submission of predictions. Consequently, there are headings for eight additional columns (see section on Format of Predictions).

The file DREAM3_GeneExpressionChallenge_ExpressionData.txt is a tab-delimited file that provides the relevant expression data. Columns are labeled, and are summarized here as well. The first column gives the Affymetrix probeset ID. The second column lists the commonly used gene name if there is one for that probeset. The third column represents the absolute expression level (in arbitrary units) for the probeset in the parental strain at time t=0. The next set of 8 columns contains the time course data for the wild-type strain, the following set of 8 columns contains the time course data for the gat1Δ strain, the next set of 8 columns contains the time course data for the gcn4Δ strain, and final set of 8 columns contains the time course data for the leu3Δ strain. Within each set of columns, the time points are t=0, 10, 20, 30, 45, 60, 90 and 120 minutes. The values in all of these columns express transcript levels as the log (base 2) of the ratio of expression in the indicated strain and time point to the expression level in the parental strain at time t=0. Thus, positive values indicate higher levels of expression than is observed for that probeset in the parental strain at time t=0, and negative values indicate lower expression. Data is provided for all probesets and in all strains, and at all time points, except for the 50 probesets (genes) whose expression is to be predicted (DREAM3_GeneExpressionChallenge_TargetList.txt). For those genes, the text "PREDICT" was inserted in the corresponding entries in the columns that correspond to the gat1Δ data in the file DREAM3_GeneExpressionChallenge_ExpressionData.txt.

PLEASE NOTE. The data that is being provided initially is derived from two technical replicates, using a single biological replicate. An additional biological replicate will be obtained soon, and a new version of the DREAM3_GeneExpressionChallenge_ExpressionData.txt file will be provided.
UPDATE NOTE (July 15, 2008)

As noted in the original posting of this challenge, the data set that was provided initially DREAM3_GeneExpressionChallenge_ExpressionData.txt, was based on a single biological replicate, with two technical replicates. We noted that the data file was going to be updated as additional data were obtained. Challenge participants are hereby notified that the original data file has now been superseded by the file

DREAM3_GeneExpressionChallenge_ExpressionData_UPDATED.txt.

The values in this file are based on the original data, plus a new biological replicate. All array data been reprocessed using the RMA algorithm within the commercial program GeneSpring. Probeset hybridization values were median normalized within arrays prior to the calculation of fold-change. This is the dataset that will be used in the evaluation of challenge predictions.


## Scoring
Predictions will be assessed based on rank order metrics such as Spearman’s rank correlation coefficient, and its corresponding p-value under the null hypothesis that the ranks are randomly distributed.


This challenge can be found on the DREAM3 website: https://www.synapse.org/#!Synapse:syn3033083/wiki/74369
