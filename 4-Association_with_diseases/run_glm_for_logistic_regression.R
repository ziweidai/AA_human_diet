library(glmnet)
library(ROCR)
library(dplyr)
library(ggplot2)
library(broom)

df_dep_var <- read.csv('disease_variables_for_lr.csv')
df_indep_var <- read.csv('independent_variables_for_lr.csv')
df_indep_var <- log(df_indep_var+0.01) #log transform the independent variables
df_weights <- read.csv('survey_weights_for_lr.csv')
weights_scaled <- df_weights/max(df_weights)
df_for_lr <- cbind(df_indep_var,df_dep_var)

df_nutrient_categories <- read.csv('NHANES_Nutrients_Annotation.csv', row.names = 1)
df_nutrients_nhanes <- df_indep_var[,17:55]
df_aas <- df_indep_var[,56:73]
auc_1se_mean <- matrix(nrow = 4, ncol = 8)
auc_min_mean <- matrix(nrow = 4, ncol = 8)
auc_1se_std <- matrix(nrow = 4, ncol = 8)
auc_min_std <- matrix(nrow = 4, ncol = 8)
nrepeats <- 10
auc_1se_vec <- vector(mode = "numeric", length = nrepeats)
auc_min_vec <- vector(mode = "numeric", length = nrepeats)

for (i in 1:4)
{
  for (j in 1:8)
  {
    if (j<=6)
    {
      x <- cbind(df_indep_var[,1:16],df_nutrients_nhanes[,df_nutrient_categories[,j]==1])
      y <- df_dep_var[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), i]
      weights <- weights_scaled[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), ]
      x <- x[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])),]
    }
    else if (j==7)
    {
      x <- cbind(df_indep_var[,1:16],df_aas)
      y <- df_dep_var[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), i]
      weights <- weights_scaled[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), ]
      x <- x[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])),]
    }
    else
    {
      x <- df_indep_var
      y <- df_dep_var[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), i]
      weights <- weights_scaled[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), ]
      x <- x[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])),]
    }
    
    nsamples <- length(weights)
    
    for (k in 1:nrepeats)
    {
      idx_permut <- sample(nsamples)
      ntrain <- as.integer(nsamples*0.7)
      ntest <- nsamples - ntrain
      xtrain <- head(x[idx_permut,], n = ntrain)
      ytrain <- head(y[idx_permut], n = ntrain)
      xtest <- tail(x[idx_permut,], n = ntest)
      ytest <- tail(y[idx_permut], n = ntest)
      
      lr_model_cv <- cv.glmnet(as.matrix(xtrain), ytrain, family = "binomial",
                               weights = head(weights[idx_permut], n = ntrain),
                               alpha = 0.5, nfolds = 5)
      
      plot(lr_model_cv)
      lr_model_min <- glmnet(as.matrix(xtrain), ytrain, family = "binomial",
                             weights = head(weights[idx_permut], n = ntrain),
                             alpha = 0.5, lambda = lr_model_cv$lambda.min)
      lr_model_1se <- glmnet(as.matrix(xtrain), ytrain, family = "binomial",
                             weights = head(weights[idx_permut], n = ntrain),
                             alpha = 0.5, lambda = lr_model_cv$lambda.1se)
      prob_min <- predict(lr_model_min, newx = as.matrix(xtest))
      prob_1se <- predict(lr_model_1se, newx = as.matrix(xtest))
      pred_min <- prediction(prob_min, ytest)
      pred_1se <- prediction(prob_1se, ytest)
      perf_min <- performance(pred_min, "tpr", "fpr")
      perf_1se <- performance(pred_1se, "tpr", "fpr")
      auc_min_vec[k] <- performance(pred_min, "auc")@y.values[[1]]
      auc_1se_vec[k] <- performance(pred_1se, "auc")@y.values[[1]]
    }
    auc_min_mean[i,j] <- mean(auc_min_vec)
    auc_min_std[i,j] <- sd(auc_min_vec)
    auc_1se_mean[i,j] <- mean(auc_1se_vec)
    auc_1se_std[i,j] <- sd(auc_1se_vec)
  }
}

model_coefs_min <- matrix(nrow = ncol(df_indep_var), ncol = 4)
model_coefs_1se <- matrix(nrow = ncol(df_indep_var), ncol = 4)
var_sd <- sapply(df_indep_var, sd)
for (i in 1:4)
{
  x <- df_indep_var
  y <- df_dep_var[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), i]
  weights <- weights_scaled[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), ]
  x <- x[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])),]
  lr_model_cv <- cv.glmnet(as.matrix(x), y, family = "binomial",
                           weights = weights,
                           alpha = 0.5, nfolds = 5)
  lr_model_min <- glmnet(as.matrix(x), y, family = "binomial",
                         weights = weights,
                         alpha = 0.5, lambda = lr_model_cv$lambda.min)
  lr_model_1se <- glmnet(as.matrix(x), y, family = "binomial",
                         weights = weights,
                         alpha = 0.5, lambda = lr_model_cv$lambda.1se)
  model_coefs_1se[,i] <- as.matrix(lr_model_1se$beta)
  model_coefs_min[,i] <- as.matrix(lr_model_min$beta)
}
rownames(model_coefs_min) <- colnames(df_indep_var)
rownames(model_coefs_1se) <- colnames(df_indep_var)
colnames(model_coefs_min) <- colnames(df_dep_var)
colnames(model_coefs_1se) <- colnames(df_dep_var)

model_coefs_1se_std <- model_coefs_1se*replicate(4,var_sd)
model_coefs_min_std <- model_coefs_min*replicate(4,var_sd)
write.csv(model_coefs_min_std, file = "logistic_regression_feature_importance.csv")

for (i in 1:4)
{
  x <- df_indep_var
  y <- df_dep_var[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), i]
  weights <- weights_scaled[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])), ]
  x <- x[!is.nan(rowSums(x)) & (!is.nan(df_dep_var[, i])),]
  nsamples <- length(weights)
  idx_permut <- sample(nsamples)
  ntrain <- as.integer(nsamples*0.7)
  ntest <- nsamples - ntrain
  xtrain <- head(x[idx_permut,], n = ntrain)
  ytrain <- head(y[idx_permut], n = ntrain)
  xtest <- tail(x[idx_permut,], n = ntest)
  ytest <- tail(y[idx_permut], n = ntest)
  lr_model_cv <- cv.glmnet(as.matrix(xtrain), ytrain, family = "binomial",
                           weights = head(weights[idx_permut], n = ntrain),
                           alpha = 0.5, nfolds = 5)
  lr_model_min <- glmnet(as.matrix(xtrain), ytrain, family = "binomial",
                         weights = head(weights[idx_permut], n = ntrain),
                         alpha = 0.5, lambda = lr_model_cv$lambda.min)
  prob_min <- predict(lr_model_min, newx = as.matrix(xtest))
  pred_min <- prediction(prob_min, ytest)
  perf_min <- performance(pred_min, "tpr", "fpr")
  auc_min <- performance(pred_min, "auc")@y.values[[1]]
  
  df_roc <- data.frame(tpr = perf_min@y.values[[1]], fpr = perf_min@x.values[[1]])
  dev.new()
  plot0 <- ggplot() + geom_line(data = df_roc, aes(x = fpr, y = tpr), size = 1) +
    geom_line(aes(x = c(0,1), y = c(0,1)), size = 1, linetype = "dashed") +
    theme(title = element_text(size = 22), text = element_text(size = 22)) +
    annotate("text", x = 0.75, y = 0.25, size = 10, label = paste("AUC = ", round(auc_min,2))) +
    scale_x_continuous(name = "False positive rate") +
    scale_y_continuous(name = "True positive rate") + 
    ggtitle(colnames(df_dep_var)[i])
  print(plot0)
}