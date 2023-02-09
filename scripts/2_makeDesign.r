################################

# Create design matrix for random effects imaging analysis 
# Diana Smith
# March 2022
# Last Updated August 2022

###############################################################
# path to makeDesign function 
source('/home/d9smith/github/cmig_tools_internal/cmig_tools_utils/r/makeDesign.R')

###############################################################
# Load the data from mini RDS file
ndafile <- '/space/syn50/1/data/ABCD/d9smith/random_effects/data/nda4.0_offrel.RDS'
nda <- readRDS(ndafile)

# only include subjects which pass QC
idx_dmri_inc <- which(nda$imgincl_dmri_include==1)
nda_dmri_inc <- nda[idx_dmri_inc,]

idx_t1w_inc <- which(nda$imgincl_t1w_include==1)
nda_t1w_inc <- nda[idx_t1w_inc,]

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/designMat'

###############################################################
# Design Matrix 1 - dmri_inc sample, age + sex + scanner + software

# Define the name of your design matrix file 
fname <- 'designMat01_dmri_AgeSexScanSoft.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 2 - t1w_inc sample, age + sex + scanner + software

# Define the name of your design matrix file 
fname <- 'designMat02_t1w_AgeSexScanSoft.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 3 - dmri_inc sample, age + sex + scanner + software + educ + income + PCs

# Define the name of your design matrix file 
fname <- 'designMat03_dmri_AgeSexScanSoftEducIncPcs.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 4 - t1w_inc sample, age + sex + scanner + software + educ + income + PCs

# Define the name of your design matrix file 
fname <- 'designMat04_t1w_AgeSexScanSoftEducIncPcs.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 5 - t1w_inc sample, age + sex + scanner + software + mean cortical thickness

# Define the name of your design matrix file 
fname <- 'designMat05_t1w_AgeSexScanSoftThick.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'smri_thick_cdk_mean')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 6 - t1w_inc sample, age + sex + scanner + software + total surface area

# Define the name of your design matrix file 
fname <- 'designMat06_t1w_AgeSexScanSoftArea.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'smri_area_cdk_total')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 7 - t1w_inc sample, age + sex + scanner + software + intracranial volume

# Define the name of your design matrix file 
fname <- 'designMat07_t1w_AgeSexScanSoftIcv.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'smri_vol_scs_intracranialv')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 8 - t1w_inc sample, age + sex + scanner + software + educ + income + PCs + mean cortical thickness

# Define the name of your design matrix file 
fname <- 'designMat08_t1w_AgeSexScanSoftEducIncPcsThick.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'smri_thick_cdk_mean')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 9 - t1w_inc sample, age + sex + scanner + software + educ + income + PCs + total surface area

# Define the name of your design matrix file 
fname <- 'designMat09_t1w_AgeSexScanSoftEducIncPcsArea.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'smri_area_cdk_total')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 10 - t1w_inc sample, age + sex + scanner + software + educ + income + PCs + intracranial volume 

# Define the name of your design matrix file 
fname <- 'designMat10_t1w_AgeSexScanSoftEducIncPcsIcv.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8','PC9', 'PC10', 'smri_vol_scs_intracranialv')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'high.educ', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_t1w_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)

###############################################################
# Design Matrix 11 - dmri_inc sample, age + sex + scanner + software + gestational age

# Define the name of your design matrix file 
fname <- 'designMat11_dmri_AgeSexScanSoftGest.txt'
# and it's full path to your fave output directory
outfile <- paste0(outpath, '/', fname) 

# makeDesign encodes continuous and categorical variables differently, therefore they are 
# specified using different flags. Continuous variables are added using the "contvar" flag. 
# Here we include age and the genetic PCs as continuous variables.
contvar <- c('interview_age', 'gest_age')

# Categorical variables are dummy coded and added using the "catvar" flag. makeDesign automatically
# includes an intercept. For each categorical variable one category is defined as the reference category 
# and that column is dropped. makeDesign also checks whether your matrix is rank deficient and if so
# will automatically drop further categories to avoid this. Here we include sex, scanner info and 
# SES demographics as categorical.
catvar <- c('sex', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

# The time points for which we wish to extract data are specified. You do not have to specify multiple 
# time points, however if you do, specify them in chronoligical order ( this will be important for 
# longitudinal modelling)
time <- c('baseline_year_1_arm_1', '2_year_follow_up_y_arm_1') # order matters! start with baseline!

# You can specify whether you wish to demean your continuous variables or not. the default is set to 
# demean=TRUE. Other flags include "delta" for longitudinal modelling, "interact" to include interactions, 
# "subjs" if you wish to filter by subject and "quadratic" if you wish to include a quadratic term 
# (more details on these in the following sections). The defaults to these are set to null. 

# Now run makeDesign! 
makeDesign(nda_dmri_inc, outfile, time, contvar=contvar, catvar=catvar, delta=NULL, interact=NULL, subjs=NULL, demean=TRUE, quadratic=NULL)