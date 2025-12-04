# Load required packages
library(readxl)
library(dplyr)
library(discSurv)
library(survminer)
library(survHE)
library(tidyverse)
library(mgcv)
library(writexl)

options(warn = -1)

# Read example IPD data
data <- read_xls("Example IPD.xls", sheet = "Sheet1", col_names = FALSE)

# Load main functions
source("2_Extrapolation function.R")
source("3_Generate results.R")

# Run survival extrapolation for intervention and control arms
res_intervention <- Surv_analysis1(data)
res_control      <- Surv_analysis2(data)

# Generate within-arm evaluation results
Arm_res_intervention <- Arm_res(res_intervention)
Arm_res_control      <- Arm_res(res_control)

# Placeholder: between-arm comparison (to be added later)
# Arm_compare <- Compare_arms(res_intervention, res_control)

# Export output results
write_xlsx(res_intervention, "res_intervention.xlsx")
write_xlsx(res_control, "res_control.xlsx")