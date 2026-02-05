library(OTrecod)   # load existing package
library(transport)
library(dplyr)
library(purrr)
library(proxy)
library(ompr)
library(ompr.roi)
library(ROI)
library(ROI.plugin.glpk)
library(magrittr)
library(dml)
if(!require(lattice)) install.packages("lattice")
if(!require(viridis)) install.packages("viridis")

library(lattice)
library(viridis)

source("/home/fnkamsug/Téléchargements/OTrecod_0.1.2/OTrecod/R/OT_joint_unbalanced_robust.R")  # add the new function

data(ncds_14); data(ncds_5)

summary(ncds_14); summary(ncds_5)

merged_tab <- merge_dbs(ncds_14, ncds_5,
                        row_ID1 = 1, row_ID2 = 1,
                        NAME_Y = "GO90", NAME_Z = "RG91",
                        ordinal_DB1 = 3, ordinal_DB2 = 4,
                        impute = "MICE", R_MICE = 2,
                        seed_choice = 3023)

merged_fin = merged_tab$DB_READY[, -4]

# Implementation of the Unbalanced logic for your specific call
res = OT_Metric_Robust_KL_joint(merged_fin,
                              nominal = 1:4, ordinal = 5:6,
                              tau_A = 1.2,         # Stricter marginal for DB A
  tau_B = 0.9,         # Relaxed marginal for DB B (Unbalanced)
  epsilon = 0.05,      # Regularization strength
  k_huber = 1.345,     # M-estimator robustness threshold
  which.DB = "B")

res$DATA1_OT  # Predictions for A
res$DATA2_OT  # Predictions for B
