Summary of Mendelian Randomization (MR) Code
This document provides a summary of the code used for Mendelian Randomization (MR) analysis. 

#File list: 
Package_pre.txt
DTI_fastGWA_to_csv.txt                    # LD clumping command line code
Add_missing_data.txt                      # Run MR from image derived phenotypes to clinical outcomes
IV_read_to_clump.txt                    # Run MR from clinical outcomes to image derived phenotypes
Forward_MR.txt       # modify the GRAPPLE co
Find_confounder_influence.txt                  # prepare for submit jobs
Reverse_MR.txt
MVMR_IV_and_Outcome_pre.txt                 # prepare for submit jobs
MVMR_analysis.txt                 # submit job to cluster 
MR_plot.txt


----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
Step 1: Preparing the Data
Reading and Preprocessing GWAS Summary Statistics

----------------------Some important points of R----------------------
# Load necessary libraries  
library(data.table)  
library(dplyr)  
library(stringr)  
  
# Set working directory  
setwd("working_directory_path")  
  
# Read exposure (IV) data  
exposure_data <- fread("exposure_data.csv")  
  
# Read outcome data  
outcome_data <- fread("outcome_data.csv")  
  
# Standardize column names for consistency  
colnames(exposure_data) <- c("CHR", "POS", "SNP", "A1", "A2", "EAF", "BETA", "SE", "PVAL")  
colnames(outcome_data) <- c("CHR", "POS", "SNP", "A1", "A2", "EAF", "BETA", "SE", "PVAL")
harmonised_data <- harmonise_data(exposure_data, outcome_data)
--------------------------------------



----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
Step 2: Performing MR Analysis
Univariate MR Analysis
MR Analysis: Perform MR analysis using different methods (e.g., Wald ratio, Egger regression, weighted median, and inverse variance weighted) to estimate the causal effect.

----------------------Some important points of R----------------------
# Load MR analysis library  
library(TwoSampleMR)  
  
# Perform MR analysis  
mr_results <- mr(harmonised_data,   
                 method_list = c("wald_ratio", "MR_egger", "weighted_median", "weighted_mode", "inverse_variance_weighted"))

--------------------------------------

Multivariable MR (MVMR) Analysis
MVMR Analysis: Perform MVMR analysis using multiple methods to account for potential pleiotropy and horizontal pleiotropy.

----------------------Some important points of R----------------------
# Load MVMR analysis libraries  
library(MendelianRandomization)  
library(MVMR)  
  
# Convert harmonised data to MVMR input format  
mv_input <- mv_dat_to_MV_Input(harmonised_data)  
  
# Perform MVMR analysis  
mv_results <- mvmr(mv_input, mv_method = c("mv_ivw", "mv_egger", "mv_lasso", "mv_median"))
----------------------



----------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------
Step 3: Post-Analysis and Interpretation
Outlier Detection and Removal
Detect Outliers: Use MR-PRESSO to detect outliers that may distort the MR estimates.
Remove Outliers: Remove detected outliers from the data.
Re-Analysis: Re-perform the MR analysis using the cleaned data to obtain more robust estimates.
Visualization and Reporting
----------------------Some important points of R----------------------
# Load MR-PRESSO for outlier detection  
library(MRPRESSO)  
  
# Detect outliers  
mrpresso_results <- mr_presso(harmonised_data, test_outliers = TRUE, test_distortion = TRUE)  
  
# Remove outliers  
outlier_free_data <- remove_outliers(harmonised_data, mrpresso_results)  
  
# Re-perform MR analysis  
mr_results_clean <- mr(outlier_free_data, method_list = c("wald_ratio", "egger_regression", "weighted_median", "inverse_variance_weighted"))

# Load libraries for visualization  
library(ggplot2)  
library(flextable)  
  
# Create forest plot  
forest_plot <- forest_plot(mr_results)  
  
# Print forest plot  
print(forest_plot)  
  
# Create a flextable report  
report <- flextable(mr_results)  
  
# Save report  
write.csv(report, "mr_analysis_report.csv")
----------------------
