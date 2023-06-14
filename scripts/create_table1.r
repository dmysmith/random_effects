############################################################
# Create table 1 for flux poster (demographics)
# Diana Smith
# Sept 2022
############################################################

library(tableone)
library(plyr)

obs_file = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/designMat02_t1w_AgeSexScanSoft/FASE/obs.txt'
nda_file = '/space/syn50/1/data/ABCD/d9smith/random_effects/data/nda4.0_offrel.RDS'

outfile = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/table1.csv' 

obs = read.table(obs_file, sep = ",", header = TRUE)
colnames(obs) = c("src_subject_id","eventname")
nda = readRDS(nda_file)

df <- join(obs, nda, by = c("src_subject_id","eventname"))

table(df$eventname)

# vertexwise data
myvars = c("interview_age", 
"sex", "household.income","high.educ","married","race.4level", 
"race.6level","hisp")

catvars = c("sex", "household.income","high.educ","married","race.4level", 
"race.6level","hisp")

table <- CreateTableOne(data = df, vars = myvars, factorVars = catvars, strata = "eventname")

## Then prepare table for export
table_p <- print(table, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

## Save to a CSV file
write.csv(table_p, file = outfile)
