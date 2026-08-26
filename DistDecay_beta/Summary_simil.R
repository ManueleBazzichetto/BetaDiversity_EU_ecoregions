
#This script includes code for analysing how quickly beta diversity (here intended as similarity among communities) decreases
#with increasing geographic distance, while also controlling for differences in plot size.
#The analysis consists in computing summary stats of similarity within periods (ignoring geographic distance), within bins of geographic distance,
#and also estimating mean change in similarity using distance decay models.
#The idea is to compare quantities estimated for each period and assess whether homogenization (increase in similarity) occurred
#at short geographic distance.

#The code relies on tables (datasets) computed in the Data_for_analysis project. The tables include the values of three different
#beta diversity indices (Horn-Morisita, Jaccard and Bray-Curtis), the geographic distance among sites (expressed in meters), and the
#absolute difference in plot size.

library(data.table)


# ------------- grasslands

#empty env
rm(list = ls())

#create object including path to tables formatted for dd-analysis
path_to_grass_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_ddmodels_grassland/'

#retrieve object names of tables formatted for dd-analysis
ecor_grass_obj <- list.files(path = path_to_grass_tables, pattern = '_ddmod_grass.RData', full.names = F)

#extract ecoregion names
ecor_grass_nm <- sapply(strsplit(x = ecor_grass_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

#create an empty list that will store summaries computed across different settings
#also name the list according to object names
summary_grass <- vector(mode = 'list', length = length(ecor_grass_nm))

names(summary_grass) <- ecor_grass_nm


for(ecor_nm in ecor_grass_nm) {
  
  #load table formatted for dd-analysis - the object will be named tmp_list
  load(paste0(path_to_grass_tables, ecor_nm, '_ddmod_grass.RData'))
  
  #compute number of plots - this is a check to assess that all plots are considered in the analysis
  #after spatial duplicates and dissimilarities exceeding cap (> 100 milion plots) were excluded
  #these numbers will be compared with the number of plots from the matched datasets
  final_num_plots <- sapply(tmp_list, function(dtf) {
    
    n_pl <- length(union(unique(dtf[['s1.PlotID_cov']]), unique(dtf[['s2.PlotID_cov']])))
    
    return(n_pl)
    
  })
  
  #rbind period-specific tables after adding Period column and coercing them to data.table
  table_for_summ <- rbindlist(lapply(names(tmp_list), function(prd_nm) {
    
    dtf <- tmp_list[[prd_nm]]
    
    dtf[['Period']] <- prd_nm
    
    dtf <- as.data.table(dtf)
    
    return(dtf)
    
  }))
  
  #rm tmp_list to free memory
  rm(tmp_list)
  
  #drop PlotID and Abs_diff_plot_size cols
  table_for_summ[, c('s1.PlotID_cov', 's2.PlotID_cov', 'Abs_diff_plot_size') := NULL]
  
  #replace dissimiliraties with similarities (similarity = 1 - dissimilarity)
  table_for_summ[, names(.SD) := lapply(.SD, function(i) 1 - i), .SDcols = c('bray', 'horn', 'jaccard')]
  
  #scale distance (in meters) so that it is expressed in km
  table_for_summ[, euc_dist := euc_dist/1000]
  
  # -- create bins of geographic distance
  
  #extract max distance
  max_dist <- table_for_summ[, max(euc_dist)]
  
  #create intervals - by is the resolution
  bin_int <- seq(0, max_dist, by = 100)
  
  #create names of bin_int
  bin_int_nm <- c(paste(bin_int[-length(bin_int)], bin_int[-1], sep = '-'), paste0('>', max(bin_int)))
  
  #create column with bins - modify by reference
  table_for_summ[, Geo_bins := findInterval(x = euc_dist, vec = bin_int)]
  
  #rename values of Geo_bins according to bin_int_nm
  table_for_summ[, Geo_bins := bin_int_nm[Geo_bins]]
  
  #transform Geo_bins to a factor
  table_for_summ[, Geo_bins := factor(Geo_bins, levels = bin_int_nm)]
  
  # -- compute summaries of similiraty
  
  #period-specific summary
  prd_summary <- table_for_summ[, {
    
    rbindlist(lapply(.SD, function(x) {
      .(Mean = mean(x),
        Median = median(x),
        fst_qrt = quantile(x, probs = 0.25),
        trd_qrt = quantile(x, probs = 0.75))
      
      }), idcol = 'Index')
    
    }, by = Period, .SDcols = c('bray', 'horn', 'jaccard')]
  
  #distribution of geographic distances
  geo_distr_summary <- table_for_summ[, .(Mean = mean(euc_dist),
                                          Median = median(euc_dist),
                                          fst_qrt = quantile(euc_dist, probs = 0.25),
                                          trd_qrt = quantile(euc_dist, probs = 0.75)), by = Period]
  
  #period-specific summaries within geo bins
  prd_dist_summary <- table_for_summ[, {
    
    rbindlist(lapply(.SD, function(x) {
      
      .(Mean = mean(x),
        Median = median(x),
        fst_qrt = quantile(x, probs = 0.25),
        trd_qrt = quantile(x, probs = 0.75))
      
      }), idcol = 'Index')
    
    }, by = .(Period, Geo_bins), .SDcols = c('bray', 'horn', 'jaccard')]
  
  
  #store results
  summary_grass[[ecor_nm]] <- list(N_plots = final_num_plots,
                                   Prd_summ = as.data.frame(prd_summary),
                                   Dist_summ = as.data.frame(geo_distr_summary),
                                   Prd_dist_summ = as.data.frame(prd_dist_summary))
  
  message(paste('Done with'), ecor_nm)
  
  rm(final_num_plots, table_for_summ, max_dist, bin_int, bin_int_nm, prd_summary, geo_distr_summary, prd_dist_summary)
  
  gc()
  
  }


#check number of plots - this looks fine
Num_plots_grass <- do.call(rbind, lapply(summary_grass, function(ecor_out) ecor_out[['N_plots']]))



# ------------- forests

#rm all objects in the env, but summary_grass

obj_in_env <- ls()

obj_in_env <- c(obj_in_env[obj_in_env != 'summary_grass'], 'obj_in_env')

rm(list = obj_in_env)


#create object including path to tables formatted for dd-analysis
path_to_for_tables <- '/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_ddmodels_forest/'

#retrieve object names of tables formatted for dd-analysis
ecor_for_obj <- list.files(path = path_to_for_tables, pattern = '_ddmod_forest.RData', full.names = F)

#extract ecoregion names
ecor_for_nm <- sapply(strsplit(x = ecor_for_obj, split = '_', fixed = T), function(i) paste(i[1], i[2], sep = '_'))

#create an empty list that will store summaries computed across different settings
#also name the list according to object names
summary_for <- vector(mode = 'list', length = length(ecor_for_nm))

names(summary_for) <- ecor_for_nm


for(ecor_nm in ecor_for_nm) {
  
  #load table formatted for dd-analysis - the object will be named tmp_list
  load(paste0(path_to_for_tables, ecor_nm, '_ddmod_forest.RData'))
  
  #compute number of plots - this is a check to assess that all plots are considered in the analysis
  #after spatial duplicates were excluded
  #these numbers will be compared with the number of plots from the matched datasets
  final_num_plots <- sapply(tmp_list, function(dtf) {
    
    n_pl <- length(union(unique(dtf[['s1.PlotID_cov']]), unique(dtf[['s2.PlotID_cov']])))
    
    return(n_pl)
    
  })
  
  #rbind period-specific tables after adding Period column and coercing them to data.table
  table_for_summ <- rbindlist(lapply(names(tmp_list), function(prd_nm) {
    
    dtf <- tmp_list[[prd_nm]]
    
    dtf[['Period']] <- prd_nm
    
    dtf <- as.data.table(dtf)
    
    return(dtf)
    
  }))
  
  #rm tmp_list to free memory
  rm(tmp_list)
  
  #drop PlotID and Abs_diff_plot_size cols
  table_for_summ[, c('s1.PlotID_cov', 's2.PlotID_cov', 'Abs_diff_plot_size') := NULL]
  
  #replace dissimiliraties with similarities (similarity = 1 - dissimilarity)
  table_for_summ[, names(.SD) := lapply(.SD, function(i) 1 - i), .SDcols = c('bray', 'horn', 'jaccard')]
  
  #scale distance (in meters) so that it is expressed in km
  table_for_summ[, euc_dist := euc_dist/1000]
  
  # -- create bins of geographic distance
  
  #extract max distance
  max_dist <- table_for_summ[, max(euc_dist)]
  
  #create intervals - by is the resolution
  bin_int <- seq(0, max_dist, by = 100)
  
  #create names of bin_int
  bin_int_nm <- c(paste(bin_int[-length(bin_int)], bin_int[-1], sep = '-'), paste0('>', max(bin_int)))
  
  #create column with bins - modify by reference
  table_for_summ[, Geo_bins := findInterval(x = euc_dist, vec = bin_int)]
  
  #rename values of Geo_bins according to bin_int_nm
  table_for_summ[, Geo_bins := bin_int_nm[Geo_bins]]
  
  #transform Geo_bins to a factor
  table_for_summ[, Geo_bins := factor(Geo_bins, levels = bin_int_nm)]
  
  # -- compute summaries of similiraty
  
  #period-specific summary
  prd_summary <- table_for_summ[, {
    
    rbindlist(lapply(.SD, function(x) {
      .(Mean = mean(x),
        Median = median(x),
        fst_qrt = quantile(x, probs = 0.25),
        trd_qrt = quantile(x, probs = 0.75))
      
    }), idcol = 'Index')
    
  }, by = Period, .SDcols = c('bray', 'horn', 'jaccard')]
  
  #distribution of geographic distances
  geo_distr_summary <- table_for_summ[, .(Mean = mean(euc_dist),
                                          Median = median(euc_dist),
                                          fst_qrt = quantile(euc_dist, probs = 0.25),
                                          trd_qrt = quantile(euc_dist, probs = 0.75)), by = Period]
  
  #period-specific summaries within geo bins
  prd_dist_summary <- table_for_summ[, {
    
    rbindlist(lapply(.SD, function(x) {
      
      .(Mean = mean(x),
        Median = median(x),
        fst_qrt = quantile(x, probs = 0.25),
        trd_qrt = quantile(x, probs = 0.75))
      
    }), idcol = 'Index')
    
  }, by = .(Period, Geo_bins), .SDcols = c('bray', 'horn', 'jaccard')]
  
  
  #store results
  summary_for[[ecor_nm]] <- list(N_plots = final_num_plots,
                                   Prd_summ = as.data.frame(prd_summary),
                                   Dist_summ = as.data.frame(geo_distr_summary),
                                   Prd_dist_summ = as.data.frame(prd_dist_summary))
  
  message(paste('Done with '), ecor_nm)
  
  rm(final_num_plots, table_for_summ, max_dist, bin_int, bin_int_nm, prd_summary, geo_distr_summary, prd_dist_summary)
  
  gc()
  
}


rm(ecor_nm, ecor_for_nm, ecor_for_obj, path_to_for_tables)


#check number of plots - this looks fine
Num_plots_for <- do.call(rbind, lapply(summary_for, function(ecor_out) ecor_out[['N_plots']]))


rm(Num_plots_for)





