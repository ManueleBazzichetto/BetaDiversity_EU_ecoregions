


library(gdm)

part_vars <- list(climate = c('Prcp', 'Tavg'), human = c('Hmi_value'))



# ------------------------------ grasslands

ecor_grass_obj <- list.files(path = '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/',
                              pattern = 'grass.RData', full.names = F)

ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '.', fixed = T), function(i) i[1])

path_to_grass_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/'

period_nms <- c('Period1', 'Period2')

for(nm in ecor_grass_nm) {
  
  #load tmp_list
  load(paste0(path_to_grass_tables, ecor_grass_nm, '.RData'))
  
  
  #fit GDMs
  for(prd_nm in period_nms) {
    
    gdm_output <- gdm::gdm(data = tmp_list[[prd_nm]], geo = TRUE)
    
    gc()
    
    dev_part <- gdm::gdm.partition.deviance(sitePairTable = tmp_list[[prd_nm]], varSets = part_vars, partSpace = TRUE)
    
    tmp_list[[prd_nm]] <- list(GDM_out = gdm_output, Dev_out = dev_part)
    
    }
  
  #save tmp_list to disk
  
}













