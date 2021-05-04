#-----------------------------------------------------------------
# Compare performance of different methods for data imputation
#
# Method 1 - MICE, absolute amino acid level
# Method 2 - random forest (missForest), absolute amino acid level
# Method 3 - MICE, amino acid level scaled to protein level (AA/(Protein+1))
# Method 4 - random forest, amino acid level scaled to protein level (AA/(Protein+1))
#-----------------------------------------------------------------

#Load libraries
library(missForest)
library(mice)
library(Hmisc)

#Load and filter data
AANames <- c("Leucine","Phenylalanine","Isoleucine","Alanine","Valine","Threonine",
             "Tyrosine","Lysine","Methionine","Histidine","Serine","Aspartic.acid",
             "Glutamic.acid","Arginine","Tryptophan","Glycine","Cystine","Proline")
x <- read.csv("input_data/FoodMatrix.csv")
x_filtered <- x[,3:180]
missing_ratio <- colSums(is.na(x_filtered))/8788
x_filtered <- x_filtered[,missing_ratio<0.6]  #only keep nutrients with missing rate < 60%
x_scaled <- x_filtered
x_scaled[,AANames] <- x_filtered[,AANames]/(1+x_filtered[,"Protein"]) #Transform absolute AA level to AA/(1+Protein)

#Generate labels for cross-validation
cv_labels <- matrix(sample.int(5,prod(dim(x_filtered)),replace = TRUE),dim(x_filtered)[1],dim(x_filtered)[2])
#colname(cv_labels) <- colnames(x_filtered)

#missForest, no transformation
x_addNA <- x_filtered
x_impute <- x_filtered
for(i in 1:5)
{
  x_addNA[cv_labels==i] <- NA
  imp_cv <- missForest(x_addNA,verbose = TRUE)
  x_impute[cv_labels==i] <- imp_cv$ximp[cv_labels==i]
  x_addNA <- x_scaled
}
ximp_missForest <- x_impute
ximp_missForest[ximp_missForest<0] <- 0

#mice, no transformation
for(i in 1:5)
{
  x_addNA[cv_labels==i] <- NA
  imp_cv <- complete(mice(x_addNA,m = 1))
  x_impute[cv_labels==i] <- imp_cv[cv_labels==i]
  x_addNA <- x_scaled
}
ximp_mice <- x_impute
ximp_mice[ximp_mice<0] <- 0

#missForest, with transformation
x_addNA <- x_scaled
x_impute <- x_scaled
for(i in 1:5)
{
  x_addNA[cv_labels==i] <- NA
  imp_cv <- missForest(x_addNA,verbose = TRUE)
  x_impute[cv_labels==i] <- imp_cv$ximp[cv_labels==i]
  x_addNA <- x_scaled
}
ximp_missForest_transf <- x_impute
ximp_missForest_transf[ximp_missForest<0] <- 0
ximp_missForest_transf[,AANames] <- ximp_missForest[,AANames]*(x_filtered[,"Protein"]+1)

#mice, with transformation
for(i in 1:5)
{
  x_addNA[cv_labels==i] <- NA
  imp_cv <- complete(mice(x_addNA,m = 1))
  x_impute[cv_labels==i] <- imp_cv[cv_labels==i]
  x_addNA <- x_scaled
}
ximp_mice_transf <- x_impute
ximp_mice_transf[ximp_mice<0] <- 0
ximp_mice_transf[,AANames] <- ximp_mice[,AANames]*(x_filtered[,"Protein"]+1)

c_missForest <- diag(cor(ximp_missForest,x_filtered, use = "pairwise.complete.obs",method = "spearman"))
c_mice <- diag(cor(ximp_mice,x_filtered, use = "pairwise.complete.obs",method = "spearman"))
c_missForest_transf <- diag(cor(ximp_missForest_transf,x_filtered, use = "pairwise.complete.obs",method = "spearman"))
c_mice_transf <- diag(cor(ximp_mice_transf,x_filtered, use = "pairwise.complete.obs",method = "spearman"))
c_compare <- cbind(as.matrix(c_missForest),as.matrix(c_mice),as.matrix(c_missForest_transf),as.matrix(c_mice_transf))
colnames(c_compare) <- c("missForest","mice","missForest_transformation","mice_transformation")

#Compare ratios between amino acids before and after imputation
aaimp_missForest <- ximp_missForest[,AANames]
aaimp_mice <- ximp_mice[,AANames]
aaimp_missForest_transf <- ximp_missForest_transf[,AANames]
aaimp_mice_transf <- ximp_mice_transf[,AANames]
aa_filtered <- x_filtered[,AANames]

count <- 0
c_aaratio <- matrix(ncol = 2,nrow = 153) #correlation between pairwise amino acid ratios before/after imputation
colnames(c_aaratio) <- c("missForest","mice","missForest_transformation","mice_transformation")
for(i in 1:17)
{
  for(j in (i+1):18)
  {
    count <- count+1
    a <- aa_filtered[,i]/aa_filtered[,j]
    b1 <- aaimp_missForest[,i]/aaimp_missForest[,j]
    b2 <- aaimp_mice[,i]/aaimp_mice[,j]
    b3 <- aaimp_missForest_transf[,i]/aaimp_missForest_transf[,j]
    b4 <- aaimp_mice_transf[,i]/aaimp_mice_transf[,j]
    c_aaratio[count,1] <- cor(a,b1,method = "spearman", use = "na.or.complete")
    c_aaratio[count,2] <- cor(a,b2,method = "spearman", use = "na.or.complete")
    c_aaratio[count,3] <- cor(a,b3,method = "spearman", use = "na.or.complete")
    c_aaratio[count,4] <- cor(a,b4,method = "spearman", use = "na.or.complete")
  }
}
