
#TO BE UPDATED!
library(gdm)

library(mapview)
library(sf)
library(ggplot2)
library(ggpubr)
library(ggtext) #installed on Jan, 28, 2026
library(vegan)
library(data.table)
library(car) #for VIF


#it is not possible to run GDMs for all ecoregions simultaneously, as objects will completely fill the memory

#the comparisons (dissimilarities) between sites having the same coordinates should be dropped out from
#the table formatted as input data for GDMs (sitePairTable)
#to this aim, I'm using the EVA_duply list to create all pairs between PlotID of duplicates
#and I'm dropping these pairs out from the table that will be used as input data for GDMs

#at the same time, the table used as input for the GDMs should not include plots with the same coordinates
#for this reason, I am adding a small noise to the coordinates of spatial duplicates
#this create fake coordinates, that are not used for any practical reason, except for computing the geographic distance
#between a plots (the noise enters into play when at least one of the plot pairs is a spatial duplicate)
#the amount of noise is very small, to not significantly affect the computation of geographic distance between plots

#in case geo = TRUE in the gdms, check effect of distance on intercept - if duplicate pairs are excluded (see below)
#the effect should be reduced


#-----load data

#load Matched datasets for grasslands
load(file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Matched_dt_grass.RData')

#load Matched datasets for forests
load(file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Matched_dt_forest.RData')

#load vegetation data
load(file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/EVA_veg_data.RData')

#load EVA_duply - a list including, for each ecoregion and period, the groups of duplicated plots
load(file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/EVA_duply_list.RData')

#load ecoregion names
load(file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Selected_ecor_names.RData')

#check number of groups of duplicates for ecoregion and period
do.call(rbind, lapply(EVA_duply, function(i) sapply(i, length)))


#1) prepare data for gdm::formatsitepair()

#2) data prepared for gdm::formatsitepair() should be saved separately so that they can be used in a for loop

#3) The idea is to loop across ecoregions and run a series of functions on single ecoregions and save outputs [ideally on Prague's or any other available server]


#--------------------------------------------------------------operations to execute before splitting grasslands and forests

#check that there are no duplicated PlotIDs in Matched_datasets_grass and Matched_datasets_forest
sum(duplicated(unname(unlist(lapply(Matched_datasets_grass, function(eco) unlist(lapply(eco, function(prd) prd[['PlotID']]))))))) #0
sum(duplicated(unname(unlist(lapply(Matched_datasets_forest, function(eco) unlist(lapply(eco, function(prd) prd[['PlotID']]))))))) #0

#check hab type is unique
unique(do.call(rbind, lapply(Matched_datasets_grass, function(eco) do.call(rbind, eco)))[['EunisVerbose_lev1']]) #grasslands
unique(do.call(rbind, lapply(Matched_datasets_forest, function(eco) do.call(rbind, eco)))[['EunisVerbose_lev1']]) #forests

#subset EVA_duply to keep only duplicates of selected ecoregions

all(names(sel_ecor_names) %in% unique(c(names(Matched_datasets_grass), names(Matched_datasets_forest)))) #TRUE

all(sel_ecor_names %in% names(EVA_duply)) #TRUE

EVA_duply <- EVA_duply[sel_ecor_names]

#rename EVA_duply elements from lower-case period to upper-case Period

#https://stackoverflow.com/questions/18509527/first-letter-to-upper-case
#test_to_del <- c('period1', 'period2')
#substr(test_to_del, 1, 1) <- toupper(substr(test_to_del, 1, 1))
#test_to_del

for(i in seq_along(EVA_duply)) {
  
  prd_names <- names(EVA_duply[[i]])
  
  substr(prd_names, start = 1, stop = 1) <- toupper(substr(prd_names, start = 1, stop = 1))
  
  names(EVA_duply[[i]]) <- prd_names
  
  }

rm(prd_names)

unique(unlist(lapply(EVA_duply, names))) #"Period1" "Period2"

#rename EVA_duply elements: from long eco names to shortened versions
identical(names(EVA_duply), unname(sel_ecor_names))

names(EVA_duply) <- names(sel_ecor_names)

#!!REMOVE IF NOT NEEDED!!
#divide ecoregion names for grasslands and forests
#grass_sel_econm <- sel_ecor_names[names(sel_ecor_names) %in% names(Matched_datasets_grass)]
#for_sel_econm <- sel_ecor_names[names(sel_ecor_names) %in% names(Matched_datasets_forest)]

#--lists of PlotIDs of spatial duplicates

#get all possible pairs of spatial duplicates' PlotIDs (separately for each ecoregion and period)
#and an integer vector including unique PlotIDs of the spatial duplicates (again, separately for each ecoregion and period)

#the possible pairs will be used to filter out dissimilarities between spatial duplicates in the table generated using gdm::formatsitepair
#while the integer vector of PlotIDs will be used to add a small noise to the coordinates of the spatial duplicates (before generating the sitepair table)
#adding this noise is important because some functions (e.g., the function for computing variable importance)
#of the GDM R package assume sites are unique in terms of spatial location (see Fitzpatrick's email)

#--pairs of duplicates

#check number of groups of duplicates
do.call(rbind, lapply(EVA_duply, function(eco) sapply(eco, length)))

#create EVA_duply_pairs, which will include all possible pairs of spatial duplicates' PlotIDs
#run combo_of_plotIDs on all ecoregions to obtain all possible pairs of PlotIDs

EVA_duply_pairs <- lapply(EVA_duply, function(eco) {
  
  res <- lapply(eco, function(prd) unlist(lapply(prd, combo_of_plotIDs)))
  
  return(res)
  
})

#--obtain lists of PlotIDs having at least one duplicate

EVA_duply_unique <- lapply(EVA_duply, function(eco) {
  
  res <- lapply(eco, function(prd) {
    
    prd_int <- unlist(strsplit(prd, split = '-'))
    
    prd_int <- as.integer(prd_int)
    
    if(sum(duplicated(prd_int)) > 0) stop('Found duplicated PlotIDs')
    
    return(prd_int)
    
  })
  
  return(res)
  
  })


#--subset EVA_veg to keep only PlotIDs in Matched_datasets_grass and Matched_datasets_forest

anyNA(EVA_veg$Cover_perc) #FALSE

#--grass

EVA_veg_grass <- lapply(Matched_datasets_grass, function(eco) {
  
  #loop over periods
  res <- lapply(eco, function(prd) {
    
    #get PlotIDs for the period
    plot_ids_vec <- prd[['PlotID']]
    
    #extract vegetation data for these plots
    veg_data <- as.data.frame(EVA_veg[PlotID %in% plot_ids_vec])
    
    #count number of zero-cover observations
    zero_cov <- sum(veg_data[['Cover_perc']] == 0)
    
    #print the number
    message(paste('Dropping', zero_cov, 'rows', sep = ' '))
    
    #drop zero-cover observations
    veg_data <- veg_data[veg_data[['Cover_perc']] > 0, ]
    
    #order PlotIDs
    veg_data <- veg_data[order(veg_data[['PlotID']], decreasing = FALSE), ]
    
    return(veg_data)
    
    })
  
  #return result
  return(res)
  
  })


#--for

EVA_veg_forest <- lapply(Matched_datasets_forest, function(eco) {
  
  #loop over periods
  res <- lapply(eco, function(prd) {
    
    #get PlotIDs for the period
    plot_ids_vec <- prd[['PlotID']]
    
    #extract vegetation data for these plots
    veg_data <- as.data.frame(EVA_veg[PlotID %in% plot_ids_vec])
    
    #count number of zero-cover observations
    zero_cov <- sum(veg_data[['Cover_perc']] == 0)
    
    #print the number
    message(paste('Dropping', zero_cov, 'rows', sep = ' '))
    
    #drop zero-cover observations
    veg_data <- veg_data[veg_data[['Cover_perc']] > 0, ]
    
    #order PlotIDs
    veg_data <- veg_data[order(veg_data[['PlotID']], decreasing = FALSE), ]
    
    return(veg_data)
    
  })
  
  #return result
  return(res)
  
})

#rm EVA_veg to free memory
rm(EVA_veg)

gc()

#sanity check on remaining zero-cover rows
sum(unlist(lapply(EVA_veg_grass, function(eco) sapply(eco, function(prd) sum(prd[['Cover_perc']] == 0)))))

#check class of EVA_veg_grass and forest
unique(sapply(EVA_veg_grass, function(eco) unique(sapply(eco, class)))) #data.frame
unique(sapply(EVA_veg_forest, function(eco) unique(sapply(eco, class)))) #data.frame

#--create vector with colnames to use in the formatsite table

#dropping Slope because this won't be used in the GDMs
col_to_keep <- c('PlotID', 'Prcp', 'Tavg', 'Elevation', 'Roughness', 'Hmi_value', 'X_laea', 'Y_laea')

#--------------------------------------------------------------grasslands

#order datasets included in Matched_datasets_grass by PlotID and add a PlotID_cov column
#the PlotID_cov col, which is a copy of PlotID, will be used to keep track of the combinations of PlotIDs in the
#formatsite table for fitting GDMs
#finally, select columns that will be used for deriving the formatsite table

Matched_datasets_grass <- lapply(Matched_datasets_grass, function(eco) {
  
  #loop across periods
  res <- lapply(eco, function(prd) {
    
    #order by PlotID and select cols
    prd_ord <- prd[order(prd[['PlotID']], decreasing = FALSE), col_to_keep]
    #add PlotID_cov col
    prd_ord$PlotID_cov <- prd_ord[['PlotID']]
    #return period-specific result
    return(prd_ord)
    
  })
  #return result
  return(res)
  
})

#check order of PlotIDs in Matched_datasets_grass and EVA_veg_grass actually matches

all(mapply(function(x, y) {
  
  #extract PlotIDs in Matched_dataset_* for each eco
  plot_in_matched <- unlist(lapply(x, function(prd) prd[['PlotID']]))
  #extract PlotIDs in EVA_veg_*
  plot_in_veg <- unlist(lapply(y, function(prd) unique(prd[['PlotID']])))
  
  res <- identical(plot_in_matched, plot_in_veg)
  
  return(res)
  
}, x = Matched_datasets_grass, y = EVA_veg_grass)) #TRUE

#add a small amount of noise to coordinates of duplicates in Matched_datasets_grass

#check reproducibility
#set.seed(48)
#to_del <- add_noise_coords(mthc_lst = Matched_datasets_grass, duply_ids = EVA_duply_unique, eps_sd = 1, crs_proj = 3035)

set.seed(48)
Matched_datasets_grass <- add_noise_coords(mtch_lst = Matched_datasets_grass, duply_ids = EVA_duply_unique, eps_sd = 1, crs_proj = 3035)

#check reproducibility
#identical(to_del, Matched_datasets_grass) #TRUE
#rm(to_del)


#reformat EVA_veg data from long to wide format and add plot coordinates

EVA_veg_grass <- veg_longtowide(veg_lst = EVA_veg_grass, mtch_lst = Matched_datasets_grass)


#check potential NAs remaining in the list
any(sapply(Matched_datasets_grass, function(eco) sapply(eco, anyNA)))
any(sapply(EVA_veg_grass, function(eco) sapply(eco, anyNA)))

#check correlation (vif) among predictors

#probably not so useful to set a hard threshold (e.g., 2) since it's hard to predict effect of vif on coef's variance and CIs..
#simply report max vif observed among ecoregions (and periods)

#drop Elevation and keep roughness in
#note that Elevation has a low vif in some ecoregions, however I want to use the same predictors for all ecoregions - this is why I'm dropping only Elevation
check_multicoll(mtch_lst = Matched_datasets_grass, vars = c('Prcp', 'Tavg', 'Elevation', 'Roughness', 'Hmi_value'), vif_thr = 2)
check_multicoll(mtch_lst = Matched_datasets_grass, vars = c('Prcp', 'Tavg', 'Roughness', 'Hmi_value'), vif_thr = 2)
#check_multicoll(mtch_lst = Matched_datasets_grass, vars = c('Prcp', 'Tavg', 'Hmi_value'), vif_thr = 2)

#drop Elevation

Matched_datasets_grass <- lapply(Matched_datasets_grass, function(eco) {
  
  res_prd <- lapply(eco, function(prd) {
    
    prd <- prd[setdiff(colnames(prd), 'Elevation')]
    
    return(prd)
    
    
  })
  
  return(res_prd)
  
})


#--check sample size of matched datasets

Smp_size_datasets_gr <- as.data.frame(do.call(rbind, lapply(Matched_datasets_grass, function(eco) {
  
  res <- sapply(eco, nrow)
  
  return(res)
  
  })))

#add column with total sample size

Smp_size_datasets_gr$Total <- with(Smp_size_datasets_gr, Period1 + Period2)


#--create table formatted as input data for GDMs, drop dissimilarities between spatial duplicates and save data in local

#save ecoregion specific lists including data for each period

#check eco and periods match among matched datasets, veg data and EVA_duply_pairs
identical(names(Matched_datasets_grass), names(EVA_veg_grass)) #TRUE
identical(names(Matched_datasets_grass), names(EVA_duply_pairs)) #FALSE (this is correct, as EVA_duply_pairs include forest ecoregions too)
identical(lapply(Matched_datasets_grass, names), lapply(EVA_veg_grass, names)) #TRUE
identical(lapply(Matched_datasets_grass, names), lapply(EVA_duply_pairs[names(Matched_datasets_grass)], names)) #TRUE

#check PlotID match between matched datasets and veg data
all(sapply(names(Matched_datasets_grass), function(nm) {
  
  all(sapply(names(Matched_datasets_grass[[nm]]), function(prd_nm) {
    
    identical(Matched_datasets_grass[[c(nm, prd_nm, 'PlotID')]], EVA_veg_grass[[c(nm, prd_nm, 'PlotID')]])
    
  }))
  
  })) #TRUE

#process data and save objects

grass_names <- names(Matched_datasets_grass)

prd_names <- names(Matched_datasets_grass[[1]])

tmp_list <- setNames(vector(mode = 'list', length = length(prd_names)), nm = prd_names)

for(nm in grass_names) {
  
  for(prd in prd_names) {
    
    tmp_list[[prd]] <- gdm::formatsitepair(bioData = EVA_veg_grass[[nm]][[prd]], bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                   XColumn = 'X_laea', YColumn = 'Y_laea', predData = Matched_datasets_grass[[nm]][[prd]])
    
    tmp_list[[prd]] <- drop_unwanted_combos(x = tmp_list[[prd]], combos = EVA_duply_pairs[[nm]][[prd]], col1 = 's1.PlotID_cov', col2 = 's2.PlotID_cov')
    
    }
  
  save(tmp_list, file = paste('/Temporary_proj_run_GDM/tmp_obj_for_gdm_grass/', nm, '_grass.RData', sep = ''))
  
  tmp_list <- setNames(vector(mode = 'list', length = length(prd_names)), nm = prd_names)
  
  }

rm(nm, prd, tmp_list)

#samples size of formatted tables for GDMs before excluding dissimilarities among spatial duplicates
#I am creating this vector to evaluate range of sample sizes and set proportion of dissimilarities to use
#when computing variable importance

diss_size_grass <- do.call(rbind, lapply(Matched_datasets_grass, function(eco) {
  
  res <- sapply(eco, nrow)
  
  res <- (res*(res - 1))/2
  
  return(res)
  
  }))


quantile(as.vector(diss_size_grass), probs = c(.25, .5, .75))



#--------------------------------------------------------------forests


#order datasets included in Matched_datasets_forest by PlotID and add a PlotID_cov column
#the PlotID_cov col, which is a copy of PlotID, will be used to keep track of the combinations of PlotIDs in the
#formatsite table for fitting GDMs
#finally, select columns that will be used for deriving the formatsite table

Matched_datasets_forest <- lapply(Matched_datasets_forest, function(eco) {
  
  #loop across periods
  res <- lapply(eco, function(prd) {
    
    #order by PlotID and select cols
    prd_ord <- prd[order(prd[['PlotID']], decreasing = FALSE), col_to_keep]
    #add PlotID_cov col
    prd_ord$PlotID_cov <- prd_ord[['PlotID']]
    #return period-specific result
    return(prd_ord)
    
  })
  #return result
  return(res)
  
})


#check order of PlotIDs in Matched_datasets_forest and EVA_veg_forest actually matches

all(mapply(function(x, y) {
  
  #extract PlotIDs in Matched_dataset_* for each eco
  plot_in_matched <- unlist(lapply(x, function(prd) prd[['PlotID']]))
  #extract PlotIDs in EVA_veg_*
  plot_in_veg <- unlist(lapply(y, function(prd) unique(prd[['PlotID']])))
  
  res <- identical(plot_in_matched, plot_in_veg)
  
  return(res)
  
}, x = Matched_datasets_forest, y = EVA_veg_forest)) #TRUE


#add a small amount of noise to coordinates of duplicates in Matched_datasets_forest

set.seed(58)
Matched_datasets_forest <- add_noise_coords(mtch_lst = Matched_datasets_forest, duply_ids = EVA_duply_unique, eps_sd = 1, crs_proj = 3035)


#reformat EVA_veg data from long to wide format and add plot coordinates

EVA_veg_forest <- veg_longtowide(veg_lst = EVA_veg_forest, mtch_lst = Matched_datasets_forest)


#check potential NAs remaining in the list
any(sapply(Matched_datasets_forest, function(eco) sapply(eco, anyNA)))
any(sapply(EVA_veg_forest, function(eco) sapply(eco, anyNA)))


#check correlation (vif) among predictors

#probably not so useful to set a hard threshold (e.g., 2) since it's hard to predict effect of vif on coef's variance and CIs..
#simply report max vif observed among ecoregions (and periods)

#dropping only Elevation for same reasons as reported for grasslands
check_multicoll(mtch_lst = Matched_datasets_forest, vars = c('Prcp', 'Tavg', 'Elevation', 'Roughness', 'Hmi_value'), vif_thr = 2)
check_multicoll(mtch_lst = Matched_datasets_forest, vars = c('Prcp', 'Tavg', 'Roughness', 'Hmi_value'), vif_thr = 2)
#check_multicoll(mtch_lst = Matched_datasets_forest, vars = c('Prcp', 'Tavg', 'Hmi_value'), vif_thr = 2)

#drop Elevation

Matched_datasets_forest <- lapply(Matched_datasets_forest, function(eco) {
  
  res_prd <- lapply(eco, function(prd) {
    
    prd <- prd[setdiff(colnames(prd), 'Elevation')]
    
    return(prd)
    
    
  })
  
  return(res_prd)
  
})


#--check sample size of matched datasets

Smp_size_datasets_for <- as.data.frame(do.call(rbind, lapply(Matched_datasets_forest, function(eco) {
  
  res <- sapply(eco, nrow)
  
  return(res)
  
})))

#add column with total sample size

Smp_size_datasets_for$Total <- with(Smp_size_datasets_for, Period1 + Period2)


#--create table formatted as input data for GDMs, drop dissimilarities between spatial duplicates and save data in local

#save ecoregion specific lists including data for each period

#check eco and periods match among matched datasets, veg data and EVA_duply_pairs
identical(names(Matched_datasets_forest), names(EVA_veg_forest)) #TRUE
identical(names(Matched_datasets_forest), names(EVA_duply_pairs)) #FALSE (this is correct, as EVA_duply_pairs include grass ecoregions too)
identical(lapply(Matched_datasets_forest, names), lapply(EVA_veg_forest, names)) #TRUE
identical(lapply(Matched_datasets_forest, names), lapply(EVA_duply_pairs[names(Matched_datasets_forest)], names)) #TRUE

#check PlotID match between matched datasets and veg data
all(sapply(names(Matched_datasets_forest), function(nm) {
  
  all(sapply(names(Matched_datasets_forest[[nm]]), function(prd_nm) {
    
    identical(Matched_datasets_forest[[c(nm, prd_nm, 'PlotID')]], EVA_veg_forest[[c(nm, prd_nm, 'PlotID')]])
    
  }))
  
})) #TRUE


#process data and save objects

exists(x = 'nm'); exists(x = 'prd'); exists(x = 'tmp_list') #FALSE*3


forest_names <- names(Matched_datasets_forest)

#using prd_names created for grasslands
#prd_names <- names(Matched_datasets_forest[[1]])

tmp_list <- setNames(vector(mode = 'list', length = length(prd_names)), nm = prd_names)

for(nm in forest_names) {
  
  for(prd in prd_names) {
    
    tmp_list[[prd]] <- gdm::formatsitepair(bioData = EVA_veg_forest[[nm]][[prd]], bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                           XColumn = 'X_laea', YColumn = 'Y_laea', predData = Matched_datasets_forest[[nm]][[prd]])
    
    tmp_list[[prd]] <- drop_unwanted_combos(x = tmp_list[[prd]], combos = EVA_duply_pairs[[nm]][[prd]], col1 = 's1.PlotID_cov', col2 = 's2.PlotID_cov')
    
  }
  
  save(tmp_list, file = paste('/Temporary_proj_run_GDM/tmp_obj_for_gdm_forest/', nm, '_forest.RData', sep = ''))
  
  tmp_list <- setNames(vector(mode = 'list', length = length(prd_names)), nm = prd_names)
  
}

rm(nm, prd, tmp_list)


#samples size of formatted tables for GDMs before excluding dissimilarities among spatial duplicates
#I am creating this vector to evaluate range of sample sizes and set proportion of dissimilarities to use
#when computing variable importance

diss_size_for <- do.call(rbind, lapply(Matched_datasets_forest, function(eco) {
  
  res <- sapply(eco, nrow)
  
  res <- (res*(res - 1))/2
  
  return(res)
  
}))


quantile(as.vector(diss_size_for), probs = c(.25, .5, .75))


#save EVA_veg_* datasets to be used in another project to assess how beta diversity changes along geographical distance
#save EVA_duply_pairs for the same reason

save(EVA_veg_grass, EVA_veg_forest, file = '/Temporary_proj_beta_dist/EVA_veg_datasets.RData')
save(EVA_duply_pairs, file = '/Temporary_proj_beta_dist/EVA_duply_pairs_list.RData')







