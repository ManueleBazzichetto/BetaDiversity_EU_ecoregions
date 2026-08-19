
#This script includes the R code to fit GDMs and compute deviance partitions for grassland and forest ecoregions

#GDMs are fitted using the gdm function and including the following predictors:
#Temperature, Precipitation, Human Modification Index, Topographic Roughness, Geographic Distance and (difference) in Plot Size

#Arguments for setting number of splines basis functions and knots are left at their default values
#(3 basis functions and 3 knots, min, median and max of the predictor)

#Explained deviance of the fitted GDMs is also partitioned among the following components:
#climate (temp and prcp), human (hmi), and space (geographic distance)
#Notice that it is not possible to use more than 3 partitions, and than space is considered for the partitioning setting partSpace = T

#version 1.6-0-7
library(gdm)

rm(list = ls())

#these components will be used to partition deviance of GDMs of both grasslands and forests
part_vars <- list(climate = c('Prcp', 'Tavg'), human = c('Hmi_value'))


# ------------------------------ grasslands

#create object including path to tables formatted for GDMs
path_to_grass_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/'

#retrieve object names of tables formatted for GDMs
ecor_grass_obj <- list.files(path = path_to_grass_tables, pattern = 'grass.RData', full.names = F)

#extract ecoregion names
ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '.', fixed = T), function(i) i[1])

#create vector including period names matching names of elements in tmp_list (object's name of tables formatted for GDMs)
#names(tmp_list) #"Period1" "Period2"
period_nms <- c("Period1", "Period2")

#the loop runs across ecoregion names and computes GDMs and dev partitions separately for each period.
#outputs are then stored in period-specific element of already existing tmp_list, so that tables formatted for GDMs
#are replaced by outputs.

for(nm in ecor_grass_nm) {
  
  #load tmp_list
  load(paste0(path_to_grass_tables, nm, '.RData'))
  
  #fit GDMs separately for period
  for(prd_nm in period_nms) {
    
    #fit the gdm
    gdm_output <- gdm::gdm(data = tmp_list[[prd_nm]], geo = TRUE)
    
    #call gc to free memory
    gc()
    
    #compute dev partitions
    dev_part <- gdm::gdm.partition.deviance(sitePairTable = tmp_list[[prd_nm]], varSets = part_vars, partSpace = TRUE)
    
    #save output (a list) in period-specific position of tmp_list
    tmp_list[[prd_nm]] <- list(GDM_out = gdm_output, Dev_out = dev_part)
    
    rm(gdm_output, dev_part)
    
    }
  
  #save tmp_list to disk
  save(tmp_list, file = paste0('GDM_output_grassland/', nm, '_gdm_out.RData'))
  
  #message
  message(paste('Done with', nm))
  
  #rm tmp_list
  rm(tmp_list)
  
  #call gc
  gc()
  
}

rm(nm, prd_nm)



# ------------------------------ forests

rm(list = ls())

#these components will be used to partition deviance of GDMs of both grasslands and forests
part_vars <- list(climate = c('Prcp', 'Tavg'), human = c('Hmi_value'))

# --

#create object including path to tables formatted for GDMs
path_to_for_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_forest/'

#retrieve object names of tables formatted for GDMs
ecor_for_obj <- list.files(path = path_to_for_tables, pattern = 'forest.RData', full.names = F)

#extract ecoregion names
ecor_for_nm <- sapply(strsplit(x = ecor_for_obj, split = '.', fixed = T), function(i) i[1])

#create vector including period names matching names of elements in tmp_list (object's name of tables formatted for GDMs)
#names(tmp_list) #"Period1" "Period2"
period_nms <- c("Period1", "Period2")


for(nm in ecor_for_nm) {
  
  #load tmp_list
  load(paste0(path_to_for_tables, nm, '.RData'))
  
  #fit GDMs separately for period
  for(prd_nm in period_nms) {
    
    #fit the gdm
    gdm_output <- gdm::gdm(data = tmp_list[[prd_nm]], geo = TRUE)
    
    #call gc to free memory
    gc()
    
    #compute dev partitions
    dev_part <- gdm::gdm.partition.deviance(sitePairTable = tmp_list[[prd_nm]], varSets = part_vars, partSpace = TRUE)
    
    #save output (a list) in period-specific position of tmp_list
    tmp_list[[prd_nm]] <- list(GDM_out = gdm_output, Dev_out = dev_part)
    
    rm(gdm_output, dev_part)
    
  }
  
  #save tmp_list to disk
  save(tmp_list, file = paste0('GDM_output_forest/', nm, '_gdm_out.RData'))
  
  #message
  message(paste('Done with', nm))
  
  #rm tmp_list
  rm(tmp_list)
  
  #call gc
  gc()
  
}

rm(nm, prd_nm)

















