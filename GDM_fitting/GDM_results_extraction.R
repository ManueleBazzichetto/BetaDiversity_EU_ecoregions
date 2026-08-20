
#This script includes the R code to extract and process GDMs output.

library(gdm)

#clean up the environment

rm(list = ls())

#loop across objects containing GDMs outputs and extract results (deviance partitions and isplines)


# ------------------------------ grasslands


#create object including path to GDMs outputs
path_to_grass_out <- '/MOTIVATE/GDM_EuropeanEcoregions/GDM_fitting/GDM_output_grassland/'

#retrieve object names of GDMs outputs
ecor_grass_obj <- list.files(path = path_to_grass_out, pattern = 'grass_gdm_out.RData', full.names = F)

#extract ecoregion names
ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

#create vector including period names matching names of elements in tmp_list
period_nms <- c("Period1", "Period2")

#create empty matrix to store values of explained deviance
expl_dev_grass <- matrix(0, nrow = length(ecor_grass_nm), ncol = length(period_nms),
                         dimnames = list(ecor_grass_nm, period_nms))

#expl_dev_grass['Alps_cmf', 'Period1']

#create empty data.frame for storing dev part outputs
dev_part_grass <- data.frame(VARIABLE_SET = character(0), DEVIANCE = numeric(0),
                             DEVIANCE_scaled = numeric(0), ECO_NM = character(0),
                             Period = character(0)) 

#create empty list for storing isplines
isplines_grass <- vector(mode = 'list', length = length(ecor_grass_nm))

names(isplines_grass) <- ecor_grass_nm

#!!The reformat_ispl_obj function is in the Test_code script
exists('reformat_ispl_obj') #TRUE

#run the for loop to extract results

for(nm in ecor_grass_nm) {
  
  load(paste0(path_to_grass_out, nm, '_grass_gdm_out.RData'))
  
  for(prd_nm in period_nms) {
    
    # ----- DEVIANCE PARTITIONING
    
    #extract total explained deviance
    tot_expl_dev <- tmp_list[[c(prd_nm, 'GDM_out', 'explained')]]
    
    #save expl dev value
    expl_dev_grass[nm, prd_nm] <- tot_expl_dev
    
    #extract dev partitions
    dev_part_out <- tmp_list[[c(prd_nm, 'Dev_out')]]
    
    #add DEVIANCE_scaled (deviance explained by each individual or combined component scaled by total explained deviance)
    dev_part_out$DEVIANCE_scaled <- round((dev_part_out[['DEVIANCE']]/tot_expl_dev)*100, digits = 2)
    
    #add ECO_NM
    dev_part_out$ECO_NM <- nm
    
    #add Period
    dev_part_out$Period <- prd_nm
    
    #drop row including UNEXPLAINED
    dev_part_out <- dev_part_out[dev_part_out[['VARIABLE_SET']] != 'UNEXPLAINED', ]
    
    #save output
    dev_part_grass <- rbind(dev_part_grass, dev_part_out)
    
    # ----- EXTRACT ISPLINES
    
    isplines_grass[[nm]][[prd_nm]] <- reformat_ispl_obj(x = tmp_list[[c(prd_nm, 'GDM_out')]], eco_nm = nm, prd_nm = prd_nm)
    
    rm(tot_expl_dev, dev_part_out)
    
  }
  
  isplines_grass[[nm]] <- do.call(rbind, isplines_grass[[nm]])
  
  rm(tmp_list)
  
  gc()
  
  }

rm(nm, prd_nm)


# ------------------------------ forests

exists('tmp_list') #FALSE
exists('tot_expl_dev') #FALSE
exists('dev_part_out') #FALSE


#create object including path to GDMs outputs
path_to_for_out <- '/MOTIVATE/GDM_EuropeanEcoregions/GDM_fitting/GDM_output_forest/'

#retrieve object names of GDMs outputs
ecor_for_obj <- list.files(path = path_to_for_out, pattern = 'forest_gdm_out.RData', full.names = F)

#extract ecoregion names
ecor_for_nm <- sapply(strsplit(x = ecor_for_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

#create vector including period names matching names of elements in tmp_list
#period_nms <- c("Period1", "Period2")

#create empty matrix to store values of explained deviance
expl_dev_for <- matrix(0, nrow = length(ecor_for_nm), ncol = length(period_nms),
                         dimnames = list(ecor_for_nm, period_nms))

#create empty data.frame for storing dev part outputs
dev_part_for <- data.frame(VARIABLE_SET = character(0), DEVIANCE = numeric(0),
                             DEVIANCE_scaled = numeric(0), ECO_NM = character(0),
                             Period = character(0)) 

#create empty list for storing isplines
isplines_for <- vector(mode = 'list', length = length(ecor_for_nm))

names(isplines_for) <- ecor_for_nm


#run the for loop to extract results

for(nm in ecor_for_nm) {
  
  load(paste0(path_to_for_out, nm, '_forest_gdm_out.RData'))
  
  for(prd_nm in period_nms) {
    
    # ----- DEVIANCE PARTITIONING
    
    #extract total explained deviance
    tot_expl_dev <- tmp_list[[c(prd_nm, 'GDM_out', 'explained')]]
    
    #save expl dev value
    expl_dev_for[nm, prd_nm] <- tot_expl_dev
    
    #extract dev partitions
    dev_part_out <- tmp_list[[c(prd_nm, 'Dev_out')]]
    
    #add DEVIANCE_scaled (deviance explained by each individual or combined component scaled by total explained deviance)
    dev_part_out$DEVIANCE_scaled <- round((dev_part_out[['DEVIANCE']]/tot_expl_dev)*100, digits = 2)
    
    #add ECO_NM
    dev_part_out$ECO_NM <- nm
    
    #add Period
    dev_part_out$Period <- prd_nm
    
    #drop row including UNEXPLAINED
    dev_part_out <- dev_part_out[dev_part_out[['VARIABLE_SET']] != 'UNEXPLAINED', ]
    
    #save output
    dev_part_for <- rbind(dev_part_for, dev_part_out)
    
    # ----- EXTRACT ISPLINES
    
    isplines_for[[nm]][[prd_nm]] <- reformat_ispl_obj(x = tmp_list[[c(prd_nm, 'GDM_out')]], eco_nm = nm, prd_nm = prd_nm)
    
    rm(tot_expl_dev, dev_part_out)
    
  }
  
  isplines_for[[nm]] <- do.call(rbind, isplines_for[[nm]])
  
  rm(tmp_list)
  
  gc()
  
}

rm(nm, prd_nm)




