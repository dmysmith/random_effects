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
ndafile <- '/space/syn50/1/data/ABCD/d9smith/random_effects/imaging/nda4.0_offrel.RDS'
nda <- readRDS(ndafile)

# It is STRONGLY recommended that you only include subjects which pass QC 
# for the imaging variable you are interested in testing.  For example this step
# filters using the imgincl_dmri_include variable from the abcd_imgincl01.txt instrument. 
idx_dmri_inc <- which(nda$imgincl_dmri_include==1)
nda_dmri_inc <- nda[idx_dmri_inc,]

# Define the path to the directory where you would like to save out your design matrix 
outpath <- '/space/syn50/1/data/ABCD/d9smith/random_effects/imaging/designMat'

###############################################################
# Design Matrix 1 - all covariates

# Define the name of your design matrix file 
fname <- 'designMat1_allcovs.txt'
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
catvar <- c('sex', 'high.educ', 'hisp', 'household.income', 'mri_info_deviceserialnumber', 'mri_info_softwareversion')

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

############################################################
# Design Matrix 2 - age and sex + scanner and software
# Define the name of your design matrix file 
fname <- 'designMat2_agesexscansoft.txt'
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

if (0) {

############################################################
# Design Matrix 0 - age and sex only

# Define the name of your design matrix file 
fname <- 'designMat0_agesex.txt'
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
catvar <- c('sex')

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

}
