Data, Matlab & R Code corresponding to the paper "Double-wavelet transform for multi-subject resting state functional magnetic resonance imaging data"


Data

data_analysis/data/: 1000.nii ~ 2005.nii 12 subjects' resting state fMRI data (preprocessed)
data_analysis/data/AAL.nii: AAL atlas data file
data_analysis/data/ROI_pairs_aal.csv: ROI pairs information

Original data will be available from the authors upon request.

#####################################################################################

Code for data analysis and simulation study

data_analysis/code/resting_analysis_submitted.m: this Matlab code performs our double-wavelet approach and AVG-FC approach
data_analysis/code/Summarize_result.R: this R code integrates multiple-subject's results and compute p-values

simulation/data_generation/gen_exp.m: generate two ROIs data using exponential covariance function for spatial correlation
simulation/data_generation/gen_gau0.m: generate two ROIs data when all voxels are uncorrelated
simulation/data_generation/gen_gau100.m: generate two ROIs data when all voxels are identical

simulation/simulate_analysis/wave.m: Using double wavelet transform to analyze simulated two ROI's data
simulation/simulate_analysis/mean.m: Using AV-GLM to analyze simulated two ROI's data

