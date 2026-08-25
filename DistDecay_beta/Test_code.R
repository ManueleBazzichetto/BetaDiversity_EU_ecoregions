
#This script includes test code for the analysis of the relationship between similarity and geographic distance.

library(data.table)
library(ggplot2)

#import a table formatted for dd analysis: this appears as tmp_list in the environment
load('/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_ddmodels_forest/TyrAdr_smf_ddmod_forest.RData')

length(tmp_list) #2
names(tmp_list) #"Period1" "Period2"
str(tmp_list$Period1)
#the distribution of plot sizes will likely be asymmetric for most ecoregions because most of the plots
#belong to one or a couple of size classes
hist(tmp_list$Period1$Abs_diff_plot_size) 

#extract period-specific table
#to_del <- tmp_list$Period1

#count number of plots
#this is done mostly to check that no plots were excluded when dropping dissimilarities exceeding cap 
num_plots <- sapply(tmp_list, function(dtf) {
  
  n_pl <- length(union(unique(dtf[['s1.PlotID_cov']]), unique(dtf[['s2.PlotID_cov']])))
  
  return(n_pl)
  
})

#add Period column and rbind
#this step will be implemented differently when processing the data - in this script I still need to compare outcomes of
#operations implemented using the dtf (data.frame) and the data.table
#in the loop, the period-specific datasets will be coerced to data.table in the lapply and then rbindlist will
#be used inplace of do.call + rbind
to_del <- do.call(rbind, lapply(names(tmp_list), function(prd_nm) {
  
  dtf <- tmp_list[[prd_nm]]
  
  dtf[['Period']] <- prd_nm
  
  return(dtf)
  
}))


#coerce period-specific dtf to a data.table
to_del_dtb <- as.data.table(to_del)

#drop PlotID cols
to_del_dtb[, c('s1.PlotID_cov', 's2.PlotID_cov') := NULL]

#replace dissimiliraties with similarities
to_del_dtb[, names(.SD) := lapply(.SD, function(i) 1 - i), .SDcols = c('bray', 'horn', 'jaccard')]

#check equality
all.equal(1-to_del$bray, to_del_dtb[['bray']]) #TRUE
all.equal(1-to_del$horn, to_del_dtb[['horn']]) #TRUE
all.equal(1-to_del$jaccard, to_del_dtb[['jaccard']]) #TRUE

#scale distance so that it is expressed in km
to_del_dtb[, euc_dist := euc_dist/1000]

all.equal(to_del_dtb$euc_dist, to_del$euc_dist/1000) #TRUE

# -- create column including distance bins

#extract max distance
max_dist <- to_del_dtb[, max(euc_dist)]

max_dist2 <- max(to_del$euc_dist/1000)
all.equal(max_dist, max_dist2) #TRUE

#create intervals - by is the resolution
bin_int <- seq(0, max_dist, by = 100)

#create names of bin_int
bin_int_nm <- c(paste(bin_int[-length(bin_int)], bin_int[-1], sep = '-'), paste0('>', max(bin_int)))

#create column with bins - modify by reference
to_del_dtb[, Geo_bins := findInterval(x = euc_dist, vec = bin_int)]

#rename values of Geo_bins according to bin_int_nm
to_del_dtb[, Geo_bins := bin_int_nm[Geo_bins]]

#transform Geo_bins to a factor
to_del_dtb[, Geo_bins := factor(Geo_bins, levels = bin_int_nm)]

#summary stats per period
period_summary <- to_del_dtb[, {
  
  rbindlist(lapply(.SD, function(x) {
    .(Mean = mean(x),
      Median = median(x),
      fst_qrt = quantile(x, probs = 0.25),
      trd_qrt = quantile(x, probs = 0.75))
    
    }), idcol = 'Index')
  
  }, by = Period, .SDcols = c('bray', 'horn', 'jaccard')]

#tapply(1 - to_del$bray, INDEX = list(to_del$Period), mean)
#tapply(1 - to_del$bray, INDEX = list(to_del$Period), median)
#tapply(1 - to_del$bray, INDEX = list(to_del$Period), quantile, probs = 0.25)

#check distribution of geographic distances
geo_distr <- to_del_dtb[, .(Mean = mean(euc_dist),
                          Median = median(euc_dist),
                          fst_qrt = quantile(euc_dist, probs = 0.25),
                          trd_qrt = quantile(euc_dist, probs = 0.75)), by = Period]

tapply(to_del$euc_dist/1000, INDEX = list(to_del$Period), mean)
tapply(to_del$euc_dist/1000, INDEX = list(to_del$Period), median)
tapply(to_del$euc_dist/1000, INDEX = list(to_del$Period), quantile, probs = 0.25)

#summary stats per period per dist bins
period_dist_summary <- to_del_dtb[, {
  
  rbindlist(lapply(.SD, function(x) {
    
    .(Mean = mean(x),
      Median = median(x),
      fst_qrt = quantile(x, probs = 0.25),
      trd_qrt = quantile(x, probs = 0.75))
    
    }), idcol = 'Index')
  
  }, by = .(Period, Geo_bins), .SDcols = c('bray', 'horn', 'jaccard')]


#plot
ggplot(period_dist_summary, aes(x = Geo_bins, y = Mean, group = Period, col = Period)) +
  geom_line() +
  facet_wrap(~ Index)









