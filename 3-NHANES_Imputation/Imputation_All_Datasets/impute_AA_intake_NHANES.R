#------------------------------------------------------------------------
# Process datasets from USDA SR, FNDDS and NHANES
# Final output is daily nutrient intake profiles for participants
# in NHANES with calculated or imputed amino acid intake data
#------------------------------------------------------------------------
#------------------------------------------------------------------------
# Load libraries and define necessary function
#------------------------------------------------------------------------
library(Matrix)
library(Hmisc)
library(missForest)

shared_elements <- function(a,b) #Function for extracting shared elements
{
  return(b[match(a,b)][!is.na(b[match(a,b)])])
}

merge_matrix <- function(a,b) #Function for merging two matrice
{
  cols_all <- union(colnames(a),colnames(b))
  rows_all <- union(rownames(a),rownames(b))
  ab <- matrix(NA,ncol = length(cols_all),nrow = length(rows_all))
  rownames(ab) <- rows_all
  colnames(ab) <- cols_all
  ab[rownames(a),colnames(a)] <- a
  ab[rownames(b),colnames(b)] <- b
  return(ab)
}

#------------------------------------------------------------------------
# Set the folders
#------------------------------------------------------------------------
Folders <- c("2007-2008","2009-2010","2011-2012","2013-2014")

#------------------------------------------------------------------------
# Read nutrient retention factors
#------------------------------------------------------------------------
Retention <- as.matrix(read.csv(file = "input_data/Datasets_dietary/Nutrient_retention_factors.csv",row.names = 1))
colnames(Retention) <- gsub("X","",colnames(Retention))

#------------------------------------------------------------------------
# Go through the 4 folders
#------------------------------------------------------------------------
for(nf in 1:4)
{
  directory <- paste('input_data/Datasets_dietary/',Folders[nf],sep = "")
  setwd(directory)
  print(directory)
  #------------------------------------------------------------------------
  # Read and process USDA SR data files, perform data imputation
  # using missForest
  #------------------------------------------------------------------------
  print("Processing USDA SR data files...")
  
  #### Section
  nut_data <- read.csv(file = 'NUT_DATA.csv')
  foods <- unique(nut_data[,1])
  nuts <- unique(nut_data[,2])
  foods_idx <- match(nut_data[,1],foods)
  nuts_idx <- match(nut_data[,2],nuts)
  
  food_nut_matrix <- matrix(NA,ncol = length(nuts),nrow = length(foods))
  for(i in 1:dim(nut_data)[1])
  {
    i1 <- foods_idx[i]
    i2 <- nuts_idx[i]
    food_nut_matrix[i1,i2] <- nut_data[i,3]
  }
  
  rownames(food_nut_matrix) <- as.character(foods)
  colnames(food_nut_matrix) <- as.character(nuts)
  
  missing_ratio <- colSums(is.na(food_nut_matrix))/dim(food_nut_matrix)[1]
  
  sum(missing_ratio<0.6)
  ####
  
  
  #Keep nutrients with <60% missing values
  x_filtered <- food_nut_matrix[,missing_ratio<0.6]
  #Transform absolute AA level to AA/(1+Protein)
  x_filtered[,as.character(501:518)] <- x_filtered[,as.character(501:518)]/(1+x_filtered[,"203"]) 
  x_filtered[is.na(x_filtered[,"501"]) & x_filtered[,"203"]==0,as.character(501:518)] <- 0
  
  if(nf>1)
  {
    SR_matrix_old <- food_nut_matrix_imputed #Keep the nutrient data from the previous version of database
  }
  
  x_imputed <- missForest(x_filtered,verbose = TRUE)
  food_nut_matrix_imputed <- x_imputed$ximp
  food_nut_matrix_imputed[,as.character(501:518)] <- 
    food_nut_matrix_imputed[,as.character(501:518)]*(1+x_filtered[,"203"])
  write.csv(food_nut_matrix_imputed,file = "USDA_SR_Imputed.csv")
  
  #Add foods from the previous SR versions to avoid bugs due to discontinued food codes
  if(nf>1)
  {
    food_nut_matrix_imputed <- merge_matrix(food_nut_matrix_imputed,SR_matrix_old)
  }
  
  #-----------------------------------------------------------------------
  # Calculate nutrient values for foods in FNDDS
  #-----------------------------------------------------------------------
  print("Calculating nutrient values for foods in FNDDS...")
  
  Retention_full <- matrix(1,ncol = dim(food_nut_matrix_imputed)[2], nrow = dim(Retention)[1])
  colnames(Retention_full) <- colnames(food_nut_matrix_imputed)
  rownames(Retention_full) <- rownames(Retention)
  Retention_full[rownames(Retention),colnames(Retention)] <- Retention/100
  
  FNDDS_SR_links <- read.csv("FNDDSSRLinks.csv")
  FNDDS_list <- as.character(unique(FNDDS_SR_links[,"Food.code"])) #Extract list of foods in FNDDS
  SR_list <- rownames(food_nut_matrix_imputed)
  
  
  #Calculated values are stored in the matrix FNDDS_nut_matrix
  FNDDS_nut_matrix <- matrix(NA,nrow = length(FNDDS_list), ncol = dim(food_nut_matrix_imputed)[2])
  colnames(FNDDS_nut_matrix) <- colnames(food_nut_matrix_imputed)
  rownames(FNDDS_nut_matrix) <- FNDDS_list
  
  
  #Read adjustments of moisture and fats
  mf_adj <- read.csv(file = "MoistNFatAdjust.csv")
  fat_FNDDS_list <- setdiff(unique(mf_adj[,"Type.of.fat"]),union(SR_list,"0"))
  
  #Update information for fats not included in USDA SR
  for(f in fat_FNDDS_list)
  {
    i <- which(f %in% fat_FNDDS_list)
    components <- FNDDS_SR_links[FNDDS_SR_links[,"Food.code"]==f,"SR.code"]
    weights <- as.matrix(FNDDS_SR_links[FNDDS_SR_links[,"Food.code"]==f,"Weight"])
    components_SR <- intersect(components,SR_list)
    components_FNDDS <- intersect(components,FNDDS_list)
    components_unknown <- setdiff(components,union(components_SR,components_FNDDS))
    if(length(components_unknown)==0)
    {
      t1 <- 0
      t2 <- 0
      pos_SR <- which(!is.na(match(components,SR_list)))
      pos_FNDDS <- which(!is.na(match(components,FNDDS_list)))
      if(length(components_FNDDS)>0)
      {
        t1 <- t(weights[pos_FNDDS]) %*% 
          FNDDS_nut_matrix[shared_elements(components,FNDDS_list),]
      }
      if(length(components_SR)>0)
      {
        t2 <- t(weights[pos_SR]) %*% 
          food_nut_matrix_imputed[shared_elements(components,SR_list),]
      }
      FNDDS_nut_matrix[as.character(fat_FNDDS_list[i]),] <- (t1 + t2)/sum(weights)
    }
  }
  
  count_na_new <- sum(is.na(FNDDS_nut_matrix[,1]))
  count_na_old <- 0
  while(count_na_old != count_na_new)
  {
    count_na_old <- count_na_new
    for(i in 1:length(FNDDS_list))
    {
      if(is.na(FNDDS_nut_matrix[i,1]))
      {
        f <- FNDDS_list[i]
        components <- FNDDS_SR_links[FNDDS_SR_links[,"Food.code"]==f,"SR.code"]
        weights <- as.matrix(FNDDS_SR_links[FNDDS_SR_links[,"Food.code"]==f,"Weight"])
        components_SR <- intersect(components,SR_list)
        components_FNDDS <- intersect(components,FNDDS_list)
        components_unknown <- setdiff(components,union(components_SR,components_FNDDS))
        rc <- FNDDS_SR_links[FNDDS_SR_links[,"Food.code"]==f,"Retention.code"]
        retention_sub <- Retention_full[as.character(rc),,drop = FALSE]
        colnames(retention_sub) <- colnames(FNDDS_nut_matrix)
        if(length(components_unknown)==0) #all components of this food are known
        {
          t1 <- 0
          t2 <- 0
          pos_SR <- which(!is.na(match(components,SR_list)))
          pos_FNDDS <- which(!is.na(match(components,FNDDS_list)))
          if(length(components_FNDDS)>0)
          {
            retention_sub1 <- retention_sub[pos_FNDDS,]
            t1 <- t(weights[pos_FNDDS]) %*%
              (FNDDS_nut_matrix[shared_elements(components,FNDDS_list),]*retention_sub1)
          }
          if(length(components_SR)>0)
          {
            retention_sub2 <- retention_sub[pos_SR,]
            t2 <- t(weights[pos_SR]) %*% 
              (food_nut_matrix_imputed[shared_elements(components,SR_list),]*retention_sub2)
          }
          FNDDS_nut_matrix[i,] <- (t1 + t2)/sum(weights)
          nut_old <- FNDDS_nut_matrix[i,]
          #Adjust for moisture and fat changes
          weight_new <- 100 + mf_adj[mf_adj[,1]==f,"Moisture.change"] + mf_adj[mf_adj[,1]==f,"Fat.change"]
          nut_change1 <- 0*nut_old
          nut_change1["255"] <- mf_adj[mf_adj[,1]==f,"Moisture.change"] #Change of water
          nut_change2 <- 0*nut_old
          fat_Change <- mf_adj[mf_adj[,1]==f,"Fat.change"]
          if(fat_Change != 0)
          {
            fat_ID <- mf_adj[mf_adj[,1]==f,"Type.of.fat"]
            if(fat_ID %in% FNDDS_list)
            {
              nut_change2 <- FNDDS_nut_matrix[as.character(fat_ID),]*fat_Change/100
            }
            else
            {
              nut_change2 <- food_nut_matrix_imputed[as.character(fat_ID),]*fat_Change/100
            }
          }
          FNDDS_nut_matrix[i,] <- (nut_old+nut_change1+nut_change2)*100/weight_new
        }
      }
    }
    count_na_new <- sum(is.na(FNDDS_nut_matrix[,1]))
  }
  FNDDS_nut_matrix[FNDDS_nut_matrix<0] <- 0 #Correct for negative values
  
  #Do additional imputation of missing values in FNDDS
  fori <- read.csv(file = "FNDDSNutVal.csv") #Original nutritional values in FNDDS, for validation only
  forimat <- t(matrix(fori[,"Nutrient.value"],nrow = length(unique(fori[,2])), ncol = length(unique(fori[,1]))))
  colnames(forimat) <- unique(fori[,2])
  rownames(forimat) <- unique(fori[,1])
  
  Nut_mat_comb <- FNDDS_nut_matrix  #Replace NAs in the imputed FNDDS_nut_matrix with values from FNDDS
  Nut_mat_comb[is.na(Nut_mat_comb[,1]),colnames(forimat)] <- forimat[is.na(Nut_mat_comb[,1]),]
  Nut_mat_comb[,as.character(501:518)] <- Nut_mat_comb[,as.character(501:518)]/(1+Nut_mat_comb[,"203"])
  Nut_mat_comb[is.na(Nut_mat_comb[,"501"]) & Nut_mat_comb[,"203"]==0,as.character(501:518)] <- 0
  x_imputed <- missForest(Nut_mat_comb,verbose = TRUE) #Impute for the remaining missing values in Nut_mat_comb
  FNDDS_nut_imp <- x_imputed$ximp

  FNDDS_nut_imp[FNDDS_nut_imp<0] <- 0
  FNDDS_nut_imp[,as.character(501:518)] <- FNDDS_nut_imp[,as.character(501:518)]*(1+FNDDS_nut_imp[,"203"])
  write.csv(FNDDS_nut_imp, file = "FNDDS_Imputed.csv")
  
  #-----------------------------------------------------------------------------
  # Calculate nutrient intake profiles for individuals in NHANES
  #-----------------------------------------------------------------------------
  print("Calculating nutrient intakes for individuals in NHANES...")
  
  #Process the food frequency recall data from NHANES and calculate nutritional intakes for the individuals
  FF_Files <- c("DR1IFF.XPT","DR2IFF.XPT")
  f_d1 <- sasxport.get(FF_Files[1])
  f_d2 <- sasxport.get(FF_Files[2])
  
  seqn_d1 <- unique(f_d1[,"seqn"])
  seqn_d2 <- unique(f_d2[,"seqn"])
  Nut_D1 <- matrix(NA,ncol = dim(FNDDS_nut_imp)[2], nrow = length(seqn_d1))
  Nut_D2 <- matrix(NA,ncol = dim(FNDDS_nut_imp)[2], nrow = length(seqn_d2))
  rownames(Nut_D1) <- seqn_d1
  colnames(Nut_D1) <- colnames(FNDDS_nut_imp)
  rownames(Nut_D2) <- seqn_d2
  colnames(Nut_D2) <- colnames(FNDDS_nut_imp)
  for(i in 1:length(seqn_d1))
  {
    pos_ind <- which(f_d1[,"seqn"]==seqn_d1[i])
    food_codes <- as.character(f_d1[pos_ind,"dr1ifdcd"])
    weights <- f_d1[pos_ind,"dr1igrms"]
    nut_values <- t(weights) %*% FNDDS_nut_imp[food_codes,]/100
    Nut_D1[i,] <- nut_values
  }
  for(i in 1:length(seqn_d2))
  {
    pos_ind <- which(f_d2[,"seqn"]==seqn_d2[i])
    food_codes <- as.character(f_d2[pos_ind,"dr2ifdcd"])
    weights <- f_d2[pos_ind,"dr2igrms"]
    nut_values <- t(weights) %*% FNDDS_nut_imp[food_codes,]/100
    Nut_D2[i,] <- nut_values
  }
  write.csv(Nut_D1, file = "NHANES_D1.csv")
  write.csv(Nut_D2, file = "NHANES_D2.csv")
  setwd("../../..")
}

#-----------------------------------------------------------------------
# Merge nutrition intake profiles from different 2-year rounds
#-----------------------------------------------------------------------
for(nf in 1:4)
{
  setwd(paste('input_data/Datasets_dietary/',Folders[nf],sep = ""))
  data_d1 <- as.matrix(read.csv(file = "NHANES_D1.csv",header = TRUE, row.names = 1))
  data_d2 <- as.matrix(read.csv(file = "NHANES_D2.csv",header = TRUE, row.names = 1))
  if(nf == 1)
  {
    merged_d1 <- data_d1
    merged_d2 <- data_d2
  }
  else
  {
    merged_d1 <- merge_matrix(merged_d1,data_d1)
    merged_d2 <- merge_matrix(merged_d2,data_d2)
  }
  setwd("../..")
}
colnames(merged_d1) <- gsub("X","",colnames(merged_d1))
colnames(merged_d2) <- gsub("X","",colnames(merged_d2))

#---------------------------------------------------------------------------------------
# Calculate nutrient intakes averaged over day 1 and day 2
#---------------------------------------------------------------------------------------
ID_mapping <- read.csv(file = "Datasets_dietary/Nut_ID_Mapping.csv")
shared_seqn12 <- intersect(rownames(merged_d1),rownames(merged_d2))
merged_ave12 <- merged_d1[shared_seqn12,]/2 + merged_d2[shared_seqn12,]/2
colnames(merged_ave12) <- ID_mapping[match(colnames(merged_ave12),ID_mapping[,1]),4]

write.csv(merged_ave12, file = "../output_data/NHANES_reconstructed_D12.csv")

