#---------------------------------------------------------------------
# Process .XPT files obtained from NHANES database, write to csv files
#---------------------------------------------------------------------

library(Hmisc)
library(ggfortify)

filenames <- c("2007-2008/DR1TOT.XPT","2007-2008/DR2TOT.XPT",
               "2009-2010/DR1TOT.XPT","2009-2010/DR2TOT.XPT","2011-2012/DR1TOT.XPT",
               "2011-2012/DR2TOT.XPT","2013-2014/DR1TOT.XPT","2013-2014/DR2TOT.XPT")
dataset_tags_day <- rep(c("1","2"),4)
dataset_tags_year <- rep(c("2007-2008","2009-2010","2011-2012","2013-2014"),rep(2,4))

#Keep nutrients that have values available in USDA SR
NutIDs <- c("tzinc","tp182","tff","tmfat","tcaff","tp183","tvk","ts160","tcarb","ttheo","tfibe",
            "tniac","ttfat","tvb2","tiron","tcopp","tfdfe","tphos","tcalc","tchl","tsugr","tvc",
            "tpfat","tvara","tvb12","tatoc","tvb1","tvb6","tm181","tsodi","tsele","tmagn","tprot",
            "tkcal","ts180","tsfat","tpota","tvd","tfa")

seqn <- 0
for (i in 1:length(filenames))
{
  fn <- paste('input_data/Datasets_dietary/',filenames[i],sep = "")
  nt <- sasxport.get(fn)
  if (i %% 2 == 1)
  {
    colnames(nt) <- gsub("dr1","",colnames(nt))
  }
  else
  {
    colnames(nt) <- gsub("dr2","",colnames(nt))
  }
  rownames(nt) <- nt[,"seqn"]
  if (i==1)
  {
    seqn <- rownames(nt)
  }
  else
  {
    seqn <- c(seqn,rownames(nt))
  }
  if (i==1)
  {
    mat_combine <- nt
  }
  else
  {
    SharedIDs <- intersect(colnames(mat_combine),colnames(nt))
    mat_combine <- rbind(mat_combine[,SharedIDs],nt[,SharedIDs])
  }
}

#Write data to file
GoodPos <- which(mat_combine[,"tkcal"]!=0)
mat_combine <- mat_combine[GoodPos,]
write.csv(mat_combine,file = "output_data/NHANES_Nutrients.csv")

#Read all other .XPT files downloaded from NHANES
#Names of datasets that are available for all four 2-yr cycles are stored in
#Datasets/CompleteFileList.txt

DirList <- c("input_data/Datasets_clinical/2007-2008/","input_data/Datasets_clinical/2009-2010/",
             "input_data/Datasets_clinical/2011-2012/","input_data/Datasets_clinical/2013-2014/")
SuffixList <- c("_E.XPT","_F.XPT","_G.XPT","_H.XPT")
TypeList <- c("Complete","Demo","DietSupp")
OutputFileNameList <- c("output_data/ClinVars.csv","output_data/Demo.csv",
                        "output_data/DietSupp.csv")

for (it in 1:3)
{
  namelist_file <- paste('input_data/Datasets_clinical/',TypeList[it],'FileList.txt',sep = "")
  NameList <- read.table(namelist_file,stringsAsFactors = F)
  
  n_dataset <- dim(NameList)[1]
  
  DatasetList <- list()
  for(i in 1:4)
  {
    y <- data.frame() #data frame for merging all datasets
    for(j in 1:n_dataset)
    {
      x <- sasxport.get(paste(DirList[i],'XPT/',NameList[j,1],SuffixList[i],sep = ""))
      rns <- x[,"seqn"]
      cns <- colnames(x)[-1]
      x <- as.data.frame(x[,-1])
      rownames(x) <- rns
      colnames(x) <- cns
      if(j==1)
      {
        y <- x
      }
      else
      {
        y <- merge(y,x,all = TRUE, by = "row.names")
        rownames(y) <- y[,"Row.names"]
        y <- y[,-1]
        a <- dim(y)
        print(a)
      }
    }
    DatasetList[[i]] <- y
  }
  
  SharedVar <- intersect(intersect(colnames(DatasetList[[1]]),colnames(DatasetList[[2]])),
                         intersect(colnames(DatasetList[[3]]),colnames(DatasetList[[4]])))
  Dataset_Comb <- rbind(DatasetList[[1]][,SharedVar],DatasetList[[2]][,SharedVar],
                        DatasetList[[3]][,SharedVar],DatasetList[[4]][,SharedVar])
  #Dataset_match_nut <- Dataset_Comb[seqn,]
  
  AllVar <- union(union(colnames(DatasetList[[1]]),colnames(DatasetList[[2]])),
                  union(colnames(DatasetList[[3]]),colnames(DatasetList[[4]])))
  nvar_all <- length(AllVar)
  mat_ext <- matrix(ncol=nvar_all,nrow = dim(DatasetList[[1]])[1] +
                      dim(DatasetList[[2]])[1]+dim(DatasetList[[3]])[1]+dim(DatasetList[[4]])[1])
  rownames(mat_ext) <- c(rownames(DatasetList[[1]]),rownames(DatasetList[[2]]),rownames(DatasetList[[3]]),rownames(DatasetList[[4]]))
  colnames(mat_ext) <- AllVar
  for(i in 1:4)
  {
    a <- DatasetList[[i]]
    mat_ext[rownames(a),colnames(a)] <- as.matrix(a)
  }
  write.csv(Dataset_Comb,file = OutputFileNameList[it])
}
