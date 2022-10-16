# Reconstruction of amino acid landscape in human foods and diets
This project includes computer code for reconstruction of amino acid profiles in human foods, dietary patterns, and dietary records. Scripts were implemented and tested under `MATLAB R2020a`, `MATLAB R2021a`, `MATLAB R2022a`, and `R version 4.0.4`. `Mosek Optimization Tools Version 9.2` is required for solving the linear programming problems.

## 1-Analyze_AA_in_foods
In this part, nutritional profiles of 8,788 human foods obtained from the [USDA SR](https://fdc.nal.usda.gov/download-datasets.html) database are used to analyze the landscape of amino acid content in human foods. 

To execute the MATLAB code, make sure that you have all required packages installed (for this part, Mosek is not required so the only software needed is MATLAB) and are working under the directory `1-Analyze_AA_in_foods`, then use the command below:

``` bash
matlab -r analyze_AA_landscape_in_foods
```

## 2-Analyze_AA_in_dietary_patterns
In this part, ten dietary patterns are defined, and then absolute and relative amino acid content in these dietary patterns are computed using linear programming and hit-and-run sampling.

To execute the MATLAB code, make sure that you: (1) have executed analyze_AA_landscape_in_foods.m for at least once; (2) have installed and configured the Mosek Optimization Tools Version 9.2; (3) are working under the directory `2-Analyze_AA_in_dietary_patterns`. Use the command below:

``` bash
matlab -r analyze_AA_landscape_in_diets
```

## 3-Reconstruct_AA_in_dietary_records
In this part, daily food consumption records of human subjects in the [NHANES](https://www.cdc.gov/nchs/nhanes/index.htm) 2007-2008, 2009-2010, 2011-2012, and 2013-2014 datasets are combined with food nutritional profiles in the [USDA SR](https://fdc.nal.usda.gov/download-datasets.html) database and the [FNDDS](https://data.nal.usda.gov/dataset/food-and-nutrient-database-dietary-studies-fndds) database to reconstruct dietary amino acid intake in these food consumption records. This part includes two subsections, one comparing the performance of several missing data imputation algorithms, the other one performing the imputation of missing data in the USDA SR and FNDDS datasets, mapping of foods and dietary records in the three databases, estimation of nutrient retention during food preparation, and reconstruction of amino acid intake profiles in the NHANES dietary records.

Code for the analysis in this part requires the R packages [missForest](https://www.rdocumentation.org/packages/missForest/versions/1.4), [Hmisc](https://www.rdocumentation.org/packages/Hmisc/versions/4.5-0), [mice](https://www.rdocumentation.org/packages/mice/versions/3.13.0) and [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html). To install these packages, use the command below:

``` bash
Rscript install_R_packages.R
```

### Compare_Imputation_Methods

This part compares four schemes of data imputation, i.e. random forest (RF) and predictive mean matching (PMM) with/without data transformation. Make sure that you have installed the required packages by running the Rscript install_R_packages.R and are working under the directory `Compare_Imputation_Methods`, and then execute the command below:

``` bash
Rscript benchmark_imputation_methods.R
```

### Imputation_All_Datasets

This part performs the imputation of missing data and reconstruction of amino acid intake profiiles in NHANES. Make sure that you are working under the directory `Imputation_All_Datasets`, and then execute the commands below:

First, read and process the raw .xpt files in input_data/:

``` bash
Rscript convert_xpt2csv.R
```

Then the dietary amino acid intake values can be reconstructed by the command below:

``` bash
Rscript impute_AA_intake_NHANES.R
```

Downstream analysis of the reconstructed dietary amino acid intake profiles can be done using the command below:

``` bash
matlab -r analyze_AA_landscape_in_NHANES
```

## 4-Association_with_diseases

This part of analysis links the human dietary amino acid records reconstructed in the previous section to health outcomes based on the clinical records available in NHANES. Make sure that you are working under the directory `4-Association_with_diseases` and have executed all scripts in the previous sections, then run the command below:

``` bash
matlab -r correlate_dietary_AA_diseases
```

## 5-Pareto_analysis_for_obesity

This part performs an in-depth analysis of the association between dietary amino acid intake and obesity incidence, and searches for diets that best balance the needs of optimization the intake of different amino acids. Make sure that you are working under the directory `5-Pareto_analysis_for_obesity` and have executed all scripts in the previous sections, then run the command below:

``` bash
matlab -r Pareto_analysis_for_obesity
```

## 6-Database

This folder includes Microsoft Excel and Microsoft Access database files that contain the amino acid levels in human foods, dietary patterns, and dietary records reconstructed in this study.

* `AA_absolute_NHANES.xlsx`: Absolute intake (g/day) of amino acids in human subjects included in the NHANES datasets
* `AA_ratio_NHANES.xlsx`: Relative intake (g AA/g protein) of amino acids in human subjects included in the NHANES datasets
* `AA_diets_max.xlsx`: Maximal absolute intake of each amino acid (g/day) in human dietary patterns such as Mediterranean diet, vegetarian diet, ketogenic diet, and so on
* `AA_diets_min.xlsx`: Minimal absolute intake of each amino acid in ten human dietary patterns such as Mediterranean diet, vegetarian diet, ketogenic diet, etc
* `AA_diets_ratio_mean.xlsx`: Average values of relative amino acid intake (g AA/g protein) in the ten human dietary patterns
* `AA_foods.xlsx`: Content of amino acids in foods included in the USDA SR database
* `Amino_acids.accdb`: Microsoft Access database file containing all information in the .xlsx files listed above

# Contact
### Ziwei Dai
* GitHub: [https://github.com/ziweidai](https://github.com/ziweidai)
* Email: [daizw@sustech.edu.cn](mailto:daizw@sustech.edu.cn)
* Website: [https://drziweidai.com](https://drziweidai.com)
