
#This script includes code for importing results of the distance-decay models fitted on the Czech HPC.

library(ggplot2)
library(ggpubr)

##!!!!!CHECK SIGNS OF COEFS!!!!!


#rm(list_out_fin, distances_used)

#this loads 2 objects: list_out_fin (with model outputs)
#and distances_used (saves distances at which expected change in similarity was evaluated)

#load('DDmod_output/DDmod_grass_out/Celtic_bf_bt_out_RAM_optim.RData')

#check convergence
#all(sapply(list_out_fin, function(i) i[['Conv']])) #T
#check number of iterations
#table(sapply(list_out_fin, function(i) i[['Iter']])) #6

#plot time of execution
#plot(1:1000, sapply(list_out_fin, function(i) i[['Time_repl']]))

#extract Coef_combo
#coef_combo_to_del <- do.call(rbind, lapply(list_out_fin, function(i) i[['Coef_combo']]))
#coef_combo_to_del <- as.data.frame(coef_combo_to_del)

#compute percentiles
#apply(coef_combo_to_del, 2, quantile, probs = c(0.025))
#apply(coef_combo_to_del, 2, quantile, probs = c(0.975))

#rm(list_out_fin, distances_used, coef_combo_to_del)

# ----------------------------- grasslands

# -- create a list including results for the ecoregions

#get object names
grass_ddmod_nm <- list.files(path = 'DDmod_output/DDmod_grass_out/', pattern = 'optim.RData', full.names = F)

#extract names of ecoregions
grass_ddmod_nm <- sapply(strsplit(grass_ddmod_nm, split = '_'), function(i) paste(i[1], i[2], sep = '_'))

#create empty lists to store mod outputs and distances
grass_ddmod_out <- setNames(vector(mode = 'list', length = length(grass_ddmod_nm)), nm = grass_ddmod_nm)
grass_ddmod_dist <- grass_ddmod_out

#load mod out and distances in the env, and store them in lists

for(eco_nm in names(grass_ddmod_out)) {
  
  load(paste0('DDmod_output/DDmod_grass_out/', eco_nm, '_bt_out_RAM_optim.RData'))
  
  grass_ddmod_out[[eco_nm]] <- list_out_fin
  grass_ddmod_dist[[eco_nm]] <- distances_used
  
  rm(list_out_fin, distances_used)
  
  }

rm(eco_nm)

#check list structure
length(grass_ddmod_out$Alps_cmf) #1000 (1 element for each bt replicate)
length(grass_ddmod_out$Alps_cmf[[1]]) #6 elements within each bt replicate 
names(grass_ddmod_out$Alps_cmf[[1]]) #names of objects associated with each replicate

# -- check convergence

#check convergence
grass_conv <- sapply(grass_ddmod_out, function(i) all(sapply(i, function(.) .[['Conv']])))
all(grass_conv)

#check number of iterations
grass_iter <- lapply(grass_ddmod_out, function(i) table(sapply(i, function(.) .[['Iter']])))

# -- check time of execution

#extract Time_repl, which measures the time (in seconds) each replicate takes to be completed
grass_time_repl <- stack(lapply(grass_ddmod_out, function(i) sapply(i, function(.) .[['Time_repl']])))

#add col with replicate id
grass_time_repl$id_repl <- rep(1:1000, times = length(unique(grass_time_repl$ind)))
#grass_time_repl[999:1002, ]

#rename cols
colnames(grass_time_repl) <- c('Time_repl', 'ECO_NM', 'Id_repl')

ggplot(data = grass_time_repl, aes(x = Id_repl, y = Time_repl, group = ECO_NM)) +
  geom_line() +
  facet_wrap(~ ECO_NM, scales = 'free_y') +
  theme_pubr()

# -- deviance

#extract deviance - stack allows creating a data.frame where each column includes bt deviance values for each ecoregion
grass_mod_dev <- stack(lapply(grass_ddmod_out, function(i) sapply(i, function(.) .[['Dev_expl']])))

#rename columns
colnames(grass_mod_dev) <- c('Deviance', 'ECO_NM')

ggplot(data = grass_mod_dev, aes(x = Deviance)) +
  geom_histogram() +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()


# -- check coefficients


# -- extract 'Coef_combo', which includes the partial effect of period at 1 km, median and max dist between plots

grass_coef_combo <- lapply(names(grass_ddmod_out), function(eco_nm) {
  
  #extract bt realizations of coef combos -> rbind them in a data.frame
  mat_coef_combo <- data.frame(do.call(rbind, lapply(grass_ddmod_out[[eco_nm]], function(.) .[['Coef_combo']])))
  
  #add ecoregion column
  mat_coef_combo$ECO_NM <- eco_nm
  
  #rename distance columns
  colnames(mat_coef_combo)[c(1, 2, 3)] <- paste0('Dist_', c('1_km', 'median', 'max'))
  
  return(mat_coef_combo)
  
  })

#check distribution of coef combos - look for outlying realizations
grass_coef_combo_dtf <- do.call(rbind, grass_coef_combo)

ggplot(grass_coef_combo_dtf, aes(x = Dist_1_km)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('1km distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()

ggplot(grass_coef_combo_dtf, aes(x = Dist_median)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Median distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()

ggplot(grass_coef_combo_dtf, aes(x = Dist_max)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Max. distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()


# -- compute quantiles of bt replicates of coef combos

#compute 2.5% percentile of coef combos
grass_combo_q2.5 <- data.frame(rbind(with(grass_coef_combo_dtf, tapply(Dist_1_km, INDEX = ECO_NM, quantile, probs = c(0.025))),
                                     with(grass_coef_combo_dtf, tapply(Dist_median, INDEX = ECO_NM, quantile, probs = c(0.025))),
                                     with(grass_coef_combo_dtf, tapply(Dist_max, INDEX = ECO_NM, quantile, probs = c(0.025)))),
                               Dist = c('1km', 'Median', 'Max.'),
                               Q = '2.5')

#compute 97.5% percentile of coef combos
grass_combo_q97.5 <- data.frame(rbind(with(grass_coef_combo_dtf, tapply(Dist_1_km, INDEX = ECO_NM, quantile, probs = c(0.975))),
                                      with(grass_coef_combo_dtf, tapply(Dist_median, INDEX = ECO_NM, quantile, probs = c(0.975))),
                                      with(grass_coef_combo_dtf, tapply(Dist_max, INDEX = ECO_NM, quantile, probs = c(0.975)))),
                                Dist = c('1km', 'Median', 'Max.'),
                                Q = '97.5')

#bind data.frames
grass_combo_qnts <- rbind(grass_combo_q2.5, grass_combo_q97.5)

#make the data.frame long for plotting results
grass_combo_qnts <- data.frame(tidyr::pivot_longer(grass_combo_qnts, cols = Alps_cmf:Sar_mf,
                                                   names_to = 'ECO_NM',
                                                   values_to = 'Percentile'))

#split percentile columns into Q2.5% and Q97.5%
grass_combo_qnts <- data.frame(tidyr::pivot_wider(grass_combo_qnts, names_from = Q,
                                                  names_prefix = 'Q_',
                                                  values_from = Percentile))

#order levels of Dist columns
grass_combo_qnts$Dist <- factor(grass_combo_qnts$Dist, levels = c('1km', 'Median', 'Max.'))

ggplot(grass_combo_qnts, aes(x = Dist, group = ECO_NM)) +
  geom_hline(yintercept = 0, col = 'red') +
  geom_errorbar(aes(ymin = Q_2.5, ymax = Q_97.5), width = 0.1) +
  facet_wrap(~ ECO_NM) +
  theme_pubr()


# ----------------------------- forests

# -- create a list including results for the ecoregions

for_ddmod_nm <- list.files(path = 'DDmod_output/DDmod_for_out/', pattern = 'optim.RData', full.names = F)

for_ddmod_nm <- sapply(strsplit(for_ddmod_nm, split = '_'), function(i) paste(i[1], i[2], sep = '_'))

#empty lists
for_ddmod_out <- setNames(vector(mode = 'list', length = length(for_ddmod_nm)), nm = for_ddmod_nm)
for_ddmod_dist <- for_ddmod_out

#load mod out and distances in the env, and store them in lists

for(eco_nm in names(for_ddmod_out)) {
  
  load(paste0('DDmod_output/DDmod_for_out/', eco_nm, '_bt_out_RAM_optim.RData'))
  
  for_ddmod_out[[eco_nm]] <- list_out_fin
  for_ddmod_dist[[eco_nm]] <- distances_used
  
  rm(list_out_fin, distances_used)
  
  }

rm(eco_nm)

# -- check convergence

#check convergence
for_conv <- sapply(for_ddmod_out, function(i) all(sapply(i, function(.) .[['Conv']])))
all(for_conv)

#check number of iterations
for_iter <- lapply(for_ddmod_out, function(i) table(sapply(i, function(.) .[['Iter']])))

# -- check time of execution

#extract Time_repl, which measures the time (in seconds) each replicate takes to be completed
for_time_repl <- stack(lapply(for_ddmod_out, function(i) sapply(i, function(.) .[['Time_repl']])))

#add col with replicate id
for_time_repl$id_repl <- rep(1:1000, times = length(unique(for_time_repl$ind)))
#for_time_repl[999:1002, ]

#rename cols
colnames(for_time_repl) <- c('Time_repl', 'ECO_NM', 'Id_repl')

ggplot(data = for_time_repl, aes(x = Id_repl, y = Time_repl, group = ECO_NM)) +
  geom_line() +
  facet_wrap(~ ECO_NM, scales = 'free_y') +
  theme_pubr()

# -- deviance

#extract deviance - stack allows creating a data.frame where each column includes bt deviance values for each ecoregion
for_mod_dev <- stack(lapply(for_ddmod_out, function(i) sapply(i, function(.) .[['Dev_expl']])))

#rename columns
colnames(for_mod_dev) <- c('Deviance', 'ECO_NM')

ggplot(data = for_mod_dev, aes(x = Deviance)) +
  geom_histogram() +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()


# -- check coefficients


# -- extract 'Coef_combo', which includes the partial effect of period at 1 km, median and max dist between plots

for_coef_combo <- lapply(names(for_ddmod_out), function(eco_nm) {
  
  mat_coef_combo <- data.frame(do.call(rbind, lapply(for_ddmod_out[[eco_nm]], function(i) i[['Coef_combo']])))
  
  mat_coef_combo$ECO_NM <- eco_nm
  
  #rename distance columns
  colnames(mat_coef_combo)[c(1, 2, 3)] <- paste0('Dist_', c('1_km', 'median', 'max'))
  
  return(mat_coef_combo)
  
})

#check distribution of coef combos - look for outlying realizations
for_coef_combo_dtf <- do.call(rbind, for_coef_combo)

ggplot(for_coef_combo_dtf, aes(x = Dist_1_km)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('1km distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()

ggplot(for_coef_combo_dtf, aes(x = Dist_median)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Median distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()

ggplot(for_coef_combo_dtf, aes(x = Dist_max)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Max. distance') +
  facet_wrap(~ ECO_NM, scales = 'free') +
  theme_pubr()


# -- compute quantiles of bt replicates of coef combos

#compute 2.5% percentile of coef combos
for_combo_q2.5 <- data.frame(rbind(with(for_coef_combo_dtf, tapply(Dist_1_km, INDEX = ECO_NM, quantile, probs = c(0.025))),
                                   with(for_coef_combo_dtf, tapply(Dist_median, INDEX = ECO_NM, quantile, probs = c(0.025))),
                                   with(for_coef_combo_dtf, tapply(Dist_max, INDEX = ECO_NM, quantile, probs = c(0.025)))),
                             Dist = c('1km', 'Median', 'Max.'),
                             Q = '2.5')

#compute 97.5% percentile of coef combos
for_combo_q97.5 <- data.frame(rbind(with(for_coef_combo_dtf, tapply(Dist_1_km, INDEX = ECO_NM, quantile, probs = c(0.975))),
                                    with(for_coef_combo_dtf, tapply(Dist_median, INDEX = ECO_NM, quantile, probs = c(0.975))),
                                    with(for_coef_combo_dtf, tapply(Dist_max, INDEX = ECO_NM, quantile, probs = c(0.975)))),
                              Dist = c('1km', 'Median', 'Max.'),
                              Q = '97.5')

#bind data.frames
for_combo_qnts <- rbind(for_combo_q2.5, for_combo_q97.5)

#make the data.frame long for plotting results
for_combo_qnts <- data.frame(tidyr::pivot_longer(for_combo_qnts, cols = Alps_cmf:TyrAdr_smf, 
                                                 names_to = 'ECO_NM', values_to = 'Percentile'))

#split percentile columns into Q2.5% and Q97.5%
for_combo_qnts <- data.frame(tidyr::pivot_wider(for_combo_qnts, names_from = Q, values_from = Percentile,
                                                names_prefix = 'Q_'))

#order levels of Dist columns
for_combo_qnts$Dist <- factor(for_combo_qnts$Dist, levels = c('1km', 'Median', 'Max.'))

#use the same y-ranges as for grasslands - TBD
ggplot(for_combo_qnts, aes(x = Dist, group = ECO_NM)) +
  geom_hline(yintercept = 0, col = 'red') +
  geom_errorbar(aes(ymin = Q_2.5, ymax = Q_97.5), width = 0.1) +
  facet_wrap(~ ECO_NM) +
  theme_pubr()


