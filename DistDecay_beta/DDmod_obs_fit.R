
#This script includes code for running distance-decay models of similarity.
#In short, the idea is to test whether homogenization or differentiation are occurring at given geographic distances between plots.
#To test this, I am fitting a distance-decay model which includes the statistical interaction between distance and Period,
#which is a categorical variable with two levels: Period1 and Period2.
#In addition, I am controlling for the difference between plot size, which affects similarity.

#The model is formulated as sim_i,j ~ geographic distance*Period + abs difference plot size,
#where sim_i,j is the similarity computed comparing composition in plot i and j,
#geographic distance is the distance between the plots, Period essentially differentiate between
#plot-pairs belonging to either Period1 or Period2, and abs difference plot size is the absolute
#of the difference in plot size (sqrt-transformed in the model, because most plots have similari size and, therefore,
#the distribution of differences is right-skewed).

#The model is fitted as a binomial with log link following Millar et al. (2011). The log link is not canonical,
#but it allows to recover the distance-decay model for similarity, while also handling zeros and ones.
#Zeros and ones are difficult to treat when fitting the log-log formulation of the distance-decay model.

#Also, the log-link model proposed by Millar et al. (2011) is equivalent to the GDM formulation
#for similarity. The GDM model uses dis_i,j = 1 - exp(-1*lin predictor_i,j). If sim_i,j = 1 - dis_i,j
#then sim_i,j = 1 - (1 - dis_i,j) = 1 - (1 - dis_i,j) = 1 - (1 - exp(-1*lin predictor_i,j)) =
#exp(-1*lin predictor_i,j) -> using log(sim_i,j) = -1*lin predictor_i,j (the distance-decay model)

#To assess homog/diff at given geographic distances, I am computing the partial effect of Period at certain distances between plots
#In case of interaction, the mean change in similarity between periods (difference in expected similarity between periods)
#is computed as the combination of the parameter for Period2 (which is the shift in intercept between Period2 at 0 distance)
#plus the product of the following quantities: the parameter of the interaction term times the distance between plots.

#The idea is to derive confidence intervals for these quantities. However, the non-independence between similarities (each plot
#is ideally used in N-1 comparisons) makes it necessary to use bootstrap. The bootstrap procedure will be implemented on
#a computer cluster.


library(data.table) #quick computations on large tables for fitting dd-models
library(glm2) #glm algorithm for improving convergence in case of non-canonical model formulations


#save model formula, the same formula will be used for both grasslands and forests
mod_formula <- as.formula('~ euc_dist*Period + Abs_diff_plot_size_sqrt')


# -- loop across ecoregions and run the dd-model, extract coefficients and other quantities of interest


# ------------- grasslands

exists('path_to_grass_tables'); exists('ecor_grass_obj'); exists('ecor_grass_nm') #FALSE*3

#create object including path to tables formatted for dd-analysis
path_to_grass_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_ddmodels_grassland/'

#retrieve object names of tables formatted for dd-analysis
ecor_grass_obj <- list.files(path = path_to_grass_tables, pattern = '_ddmod_grass.RData', full.names = F)

#extract ecoregion names
ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

exists('tmp_list')

#create empty list to store results
exists('grass_obs_ddmod_out'); exists('grass_obs_ddmod_dist')

grass_obs_ddmod_out <- setNames(vector(mode = 'list', length = length(ecor_grass_nm)), nm = ecor_grass_nm)
grass_obs_ddmod_dist <- grass_obs_ddmod_out


for(eco_nm in ecor_grass_nm) {
  
  #load tmp_list including tables for Period1 and Period2
  load(paste0(path_to_grass_tables, eco_nm, '_ddmod_grass.RData'))
  
  #process tmp_list to select only columns relevant for analyses, transform quantities and add Period column
  table_for_ddmod <- lapply(names(tmp_list), function(prd_nm) {
    
    #extract period-specific table
    dtf <- tmp_list[[prd_nm]]
    
    #transform to data.table for speeding computations up - also exclude fields that won't be used for fitting the models
    dtf <- as.data.table(dtf[c('bray', 'euc_dist', 'Abs_diff_plot_size')])
    
    #transform dissimilarity (bray) into similarity
    dtf[, bray := (1 - bray)]
    
    #sqrt-transform Abs_diff_plot_size
    dtf[, Abs_diff_plot_size_sqrt := sqrt(Abs_diff_plot_size)]
    
    #drop untransformed Abs_diff_plot_size
    dtf[, Abs_diff_plot_size := NULL]
    
    #scale euc_dist (Euclidean distance in meters) to km
    dtf[, euc_dist := euc_dist/1000]
    
    #add Period column
    dtf[, Period := prd_nm]
    
    #return dtf
    return(dtf)
    
  })
  
  #rm tmp_list to free memory
  rm(tmp_list)
  
  gc()
  
  #rbindlist and transform to data.frame
  table_for_ddmod <- data.frame(rbindlist(table_for_ddmod))
  
  #order Period levels
  table_for_ddmod$Period <- factor(table_for_ddmod$Period, levels = c('Period1', 'Period2'))
  
  #compute distances at which evaluate partial effect of Period
  med_dist <- round(mean(x = tapply(table_for_ddmod$euc_dist, INDEX = table_for_ddmod$Period, median)), digits = 2)
  max_dist <- round(min(tapply(table_for_ddmod$euc_dist, INDEX = table_for_ddmod$Period, max)), digits = 2)
  
  #save distances
  grass_obs_ddmod_dist[[eco_nm]] <-  c(Dist_mdn = med_dist, Dist_max = max_dist)
  
  #create design matrix
  desmat_for_ddmod <- model.matrix(object = mod_formula, data = table_for_ddmod)
  
  #save response vector
  sim_resp <- table_for_ddmod$bray
  
  #rm table_for_ddmod to free memory
  rm(table_for_ddmod)
  
  gc()
  
  #run the model
  mod_obj <- glm.fit2(x = desmat_for_ddmod, y = sim_resp,
                          family = binomial(link = 'log'))
  
  #extract model coefs
  mod_coefs <- coef(mod_obj)
  
  #extract quantities of interest
  grass_obs_ddmod_out[[eco_nm]] <- list(Coef = mod_coefs,
                                        Coef_combo = c(Dist1 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']], #distance 1 km - so no need to multiply
                                                       Dist2 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']]*med_dist, #median distance
                                                       Dist3 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']]*max_dist), #max distance
                                        Iter = mod_obj[['iter']],
                                        Conv = mod_obj[['converged']],
                                        Dev_expl = (1 - (deviance(mod_obj)/mod_obj[['null.deviance']]))*100)
  
  #rm objects to free memory
  rm(med_dist, max_dist, desmat_for_ddmod, sim_resp, mod_obj, mod_coefs)
  
  #message
  message(paste0('Done with: ', eco_nm))
  
  gc()
  
  }

rm(eco_nm)


#check content
grass_obs_ddmod_dist
grass_obs_ddmod_out



# ------------- forests

exists('path_to_for_tables'); exists('ecor_for_obj'); exists('ecor_for_nm') #FALSE*3

#create object including path to tables formatted for dd-analysis
path_to_for_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_ddmodels_forest/'

#retrieve object names of tables formatted for dd-analysis
ecor_for_obj <- list.files(path = path_to_for_tables, pattern = '_ddmod_forest.RData', full.names = F)

#extract ecoregion names
ecor_for_nm <- sapply(strsplit(x = ecor_for_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

exists('tmp_list') #FALSE

#create empty list to store results
exists('for_obs_ddmod_out'); exists('for_obs_ddmod_dist') #FALSE*2

for_obs_ddmod_out <- setNames(vector(mode = 'list', length = length(ecor_for_nm)), nm = ecor_for_nm)
for_obs_ddmod_dist <- for_obs_ddmod_out


for(eco_nm in ecor_for_nm) {
  
  #load tmp_list
  load(paste0(path_to_for_tables, eco_nm, '_ddmod_forest.RData'))
  
  #process tmp_list
  table_for_ddmod <- lapply(names(tmp_list), function(prd_nm) {
    
    #extract period-specific table
    dtf <- tmp_list[[prd_nm]]
    
    #transform to data.table for speeding computations up - also exclude fields that won't be used for fitting the models
    dtf <- as.data.table(dtf[c('bray', 'euc_dist', 'Abs_diff_plot_size')])
    
    #transform dissimilarity (bray) into similarity
    dtf[, bray := (1 - bray)]
    
    #sqrt-transform Abs_diff_plot_size
    dtf[, Abs_diff_plot_size_sqrt := sqrt(Abs_diff_plot_size)]
    
    #drop untransformed Abs_diff_plot_size
    dtf[, Abs_diff_plot_size := NULL]
    
    #scale euc_dist (Euclidean distance in meters) to km
    dtf[, euc_dist := euc_dist/1000]
    
    #add Period column
    dtf[, Period := prd_nm]
    
    #return dtf
    return(dtf)
    
    })
  
  #rm tmp_list to free memory
  rm(tmp_list)
  
  gc()
  
  #rbindlist and transform to data.frame
  table_for_ddmod <- data.frame(rbindlist(table_for_ddmod))
  
  #order Period levels
  table_for_ddmod$Period <- factor(table_for_ddmod$Period, levels = c('Period1', 'Period2'))
  
  #compute distances at which evaluate partial effect of Period
  med_dist <- round(mean(x = tapply(table_for_ddmod$euc_dist, INDEX = table_for_ddmod$Period, median)), digits = 2)
  max_dist <- round(min(tapply(table_for_ddmod$euc_dist, INDEX = table_for_ddmod$Period, max)), digits = 2)
  
  #save distances
  for_obs_ddmod_dist[[eco_nm]] <- c(Dist_mdn = med_dist, Dist_max = max_dist)
  
  #create design matrix
  desmat_for_ddmod <- model.matrix(object = mod_formula, data = table_for_ddmod)
  
  #save response vector
  sim_resp <- table_for_ddmod$bray
  
  #rm table_for_ddmod to free memory
  rm(table_for_ddmod)
  
  gc()
  
  #run the model
  mod_obj <- glm.fit2(x = desmat_for_ddmod, y = sim_resp,
                      family = binomial(link = 'log'))
  
  #extract model coefs
  mod_coefs <- coef(mod_obj)
  
  #extract quantities of interest
  for_obs_ddmod_out[[eco_nm]] <- list(Coef = mod_coefs,
                                        Coef_combo = c(Dist1 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']], #distance 1 km - so no need to multiply
                                                       Dist2 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']]*med_dist, #median distance
                                                       Dist3 = mod_coefs[['PeriodPeriod2']] + mod_coefs[['euc_dist:PeriodPeriod2']]*max_dist), #max distance
                                        Iter = mod_obj[['iter']],
                                        Conv = mod_obj[['converged']],
                                        Dev_expl = (1 - (deviance(mod_obj)/mod_obj[['null.deviance']]))*100)
  
  #rm objects to free memory
  rm(med_dist, max_dist, desmat_for_ddmod, sim_resp, mod_obj, mod_coefs)
  
  #message
  message(paste0('Done with: ', eco_nm))
  
  gc()
  
  }

rm(eco_nm)

#check content
for_obs_ddmod_dist
for_obs_ddmod_out





