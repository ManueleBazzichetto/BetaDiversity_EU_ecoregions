
#This code is simply to save number of dissimilarities for each ecoregion and period after excluding spatial duplicates
#and setting the cup to 100 million pairs

#The number of (dis)similarities is the same for tables used for fitting GDMs and for fitting distance-decay models

# -- grasslands

exists('path_to_grass_tables'); exists('ecor_grass_obj'); exists('ecor_grass_nm') #FALSE*3

#create object including path to tables formatted for GDMs
path_to_grass_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/'

#retrieve object names of tables formatted for GDMs
ecor_grass_obj <- list.files(path = path_to_grass_tables, pattern = 'grass.RData', full.names = F)

#extract ecoregion names
ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '.', fixed = T), function(i) i[1])

#check to del
#exists('tmp_list') #F

#load(paste0(path_to_grass_tables, "ItaScl_sdf_grass", '.RData'))

#head(tmp_list$Period1); nrow(tmp_list$Period1)
#head(tmp_list$Period2); nrow(tmp_list$Period2)

#sapply(tmp_list, nrow)

#rm(tmp_list)
#--

#create empty matrix to store result
exists('grass_smp_dis') #F
grass_smp_dis <- matrix(data = rep(0, times = length(ecor_grass_nm)*2), nrow = length(ecor_grass_nm), ncol = 2,
                        dimnames = list(ecor_grass_nm, c('Period1', 'Period2')))

for(nm in ecor_grass_nm) {
  
  load(paste0(path_to_grass_tables, nm, '.RData'))
  
  dis_smp_sz <- sapply(tmp_list, nrow)
  
  grass_smp_dis[nm, 'Period1'] <- dis_smp_sz[['Period1']]
  grass_smp_dis[nm, 'Period2'] <- dis_smp_sz[['Period2']]
  
  rm(tmp_list, dis_smp_sz)
  
}

rm(nm)

min(as.vector(grass_smp_dis)) #646,625
median(as.vector(grass_smp_dis)) #13,919,373


# -- forests

exists('path_to_for_tables'); exists('ecor_for_obj'); exists('ecor_for_nm') #FALSE*3

#create object including path to tables formatted for GDMs
path_to_for_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_forest/'

#retrieve object names of tables formatted for GDMs
ecor_for_obj <- list.files(path = path_to_for_tables, pattern = 'forest.RData', full.names = F)

#extract ecoregion names
ecor_for_nm <- sapply(strsplit(x = ecor_for_obj, split = '.', fixed = T), function(i) i[1])


#create empty matrix to store result
exists('for_smp_dis') #F

for_smp_dis <- matrix(0, nrow = length(ecor_for_nm), ncol = 2, dimnames = list(ecor_for_nm, c('Period1', 'Period2')))

exists('tmp_list') #F

for(nm in ecor_for_nm) {
  
  load(paste0(path_to_for_tables, nm, '.RData'))
  
  dis_smp_sz <- sapply(tmp_list, nrow)
  
  for_smp_dis[nm, 'Period1'] <- dis_smp_sz[['Period1']]
  for_smp_dis[nm, 'Period2'] <- dis_smp_sz[['Period2']]
  
  rm(tmp_list, dis_smp_sz)
  
}


rm(nm)

min(as.vector(for_smp_dis)) #765,171
median(as.vector(for_smp_dis)) #6,081,323

