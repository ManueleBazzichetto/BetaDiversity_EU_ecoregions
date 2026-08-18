
#Test fitting GDMs and computing deviance partitions using the gdm R package

library(gdm)

rm(list = ls())

gc()

#load example of tables formatted for GDMs

load('/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/Alps_cmf_grass.RData')

#the object loaded in the environment is named 'tmp_list'

class(tmp_list) #list
length(tmp_list) #2
names(tmp_list) #"Period1" "Period2"

sapply(tmp_list, nrow)

head(tmp_list$Period1)

# ---------------- fit gdm

#period 1
test_run_prd1 <- gdm::gdm(data = tmp_list$Period1, geo = T)

gc()

#period 2
test_run_prd2 <- gdm::gdm(data = tmp_list$Period2, geo = T)

gc()

#check summary
summary(test_run_prd1)
summary(test_run_prd2)

#expl. deviance
test_run_prd1$explained
(1 - (test_run_prd1$gdmdeviance/test_run_prd1$nulldeviance))*100

#is geo included as a predictor?
test_run_prd1$geo

#test_run_prd1$sample
#test_run_prd1$sumCoeff


# ---------------- compute deviance partitions

#create list of variables' groups to be used for deviance partitioning
#note that if env = c('Roughness'), or any other fourth predictor, is included, the function will return the following error: Cannot partition more than three variables sets
#for this reason, I'm only including climatic variables, the human modification index and geographic distance (setting partSpace = T)
test_part_vars <- list(climate = c('Prcp', 'Tavg'), human = c('Hmi_value')) #notice that geo dist is included by setting the partSpace = TRUE

test_period_devpart <- gdm::gdm.partition.deviance(sitePairTable = tmp_list$Period1, varSets = test_part_vars, partSpace = TRUE)

test_period_devpart














