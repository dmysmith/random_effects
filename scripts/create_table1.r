############################################################
# Create table 1 for flux poster (demographics)
# Diana Smith
# Sept 2022
############################################################

library(tableone)
library(plyr)

designmat_file = "/space/syn50/1/data/ABCD/d9smith/random_effects/imaging/designMat/designMat1_allcovs.txt"

vxl_idfile = "/space/syn50/1/data/ABCD/d9smith/random_effects/imaging/vxl_id_list.txt"

designmat = read.table(designmat_file, sep = "\t", header = TRUE)

table(designmat$eventname)

# vertexwise data
myvars = c("age", "interview_age", 
"sexM", "high.educHS.Diploma.GED", "high.educSome.College", 
"high.educBachelor", "high.educPost.Graduate.Degree", "hispYes", 
"household.income...50K....100K.", "household.income...100K.")

catvars = c("sexM", "high.educHS.Diploma.GED", "high.educSome.College", 
"high.educBachelor", "high.educPost.Graduate.Degree", "hispYes", 
"household.income...50K....100K.", "household.income...100K.")

CreateTableOne(data = designmat, vars = myvars, factorVars = catvars, strata = "eventname")

vxl_ids = read.table(vxl_idfile, sep = ",", header = TRUE)
names(vxl_ids) = c("eventname","src_subject_id")
designmat_vxl = join(vxl_ids, designmat)

CreateTableOne(data = designmat_vxl, vars = myvars, factorVars = catvars, strata = "eventname")
