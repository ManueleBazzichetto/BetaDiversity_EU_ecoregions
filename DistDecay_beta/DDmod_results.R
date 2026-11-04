
#This script includes code for importing results of the distance-decay model bootstrap fitted on the
#Czech High-Performance Computing cluster.

library(ggplot2)
library(ggpubr)

#rm(list_out_fin, distances_used)

#this loads 2 objects: list_out_fin (with model outputs - 6 output for bt replicate)
#and distances_used (saved distances at which expected change in similarity was evaluated)

#load('DDmod_output/DDmod_grass_out/drive-download-20261103T080718Z-1-001/CenEu_mf_bt_out_RAM_optim.RData')

#check number of replicates
#length(list_out_fin) #5000

#check convergence
#all(sapply(list_out_fin, function(i) i[['Conv']])) #T
#check number of iterations
#table(sapply(list_out_fin, function(i) i[['Iter']])) #6

#plot time of execution
#plot(1:5000, sapply(list_out_fin, function(i) i[['Time_repl']]))

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
grass_ddmod_nm <- list.files(path = 'DDmod_output/DDmod_grass_out/drive-download-20261103T080718Z-1-001/',
                             pattern = 'optim.RData', full.names = F)

#extract names of ecoregions
grass_ddmod_nm <- sapply(strsplit(grass_ddmod_nm, split = '_'), function(i) paste(i[1], i[2], sep = '_'))

#create empty lists to store mod outputs and distances
grass_ddmod_out <- setNames(vector(mode = 'list', length = length(grass_ddmod_nm)), nm = grass_ddmod_nm)
grass_ddmod_dist <- grass_ddmod_out

#load mod out and distances in the env, and store them in lists

for(eco_nm in names(grass_ddmod_out)) {
  
  load(paste0('DDmod_output/DDmod_grass_out/drive-download-20261103T080718Z-1-001/', eco_nm, '_bt_out_RAM_optim.RData'))
  
  grass_ddmod_out[[eco_nm]] <- list_out_fin
  grass_ddmod_dist[[eco_nm]] <- distances_used
  
  rm(list_out_fin, distances_used)
  
  }

rm(eco_nm)

#check list structure
length(grass_ddmod_out$Alps_cmf) #5000 (1 element for each bt replicate)
length(grass_ddmod_out$Alps_cmf[[1]]) #6 elements within each bt replicate 
names(grass_ddmod_out$Alps_cmf[[1]]) #names of objects associated with each replicate

#Coef: coefficients estimated in each bt replicate
#Coef_combo: partial effect of period at given distances between plots
#Iter: number of iterations IWLS
#Conv: logical indicating whether model converged or not
#Dev_expl: explained deviance for each bt model replicate
#Time_repl: time of execution of each bt replicate (useful to identify potential weird pattern - weird replicates - can depend on many other things, though)

# -- check number of bt replicates

lengths(grass_ddmod_out) #all 5,000

# -- check convergence

#check convergence
grass_conv <- sapply(grass_ddmod_out, function(i) all(sapply(i, function(.) .[['Conv']])))
all(grass_conv) #TRUE

#check number of iterations
grass_iter <- lapply(grass_ddmod_out, function(i) table(sapply(i, function(.) .[['Iter']]))) #always 6, except for ItaScl_sdf (7)

# -- check time of execution

#extract Time_repl, which measures the time (in seconds) each replicate takes to be completed
grass_time_repl <- stack(lapply(grass_ddmod_out, function(i) sapply(i, function(.) .[['Time_repl']])))

#add col with replicate id
grass_time_repl$id_repl <- rep(1:5000, times = length(unique(grass_time_repl$ind)))
grass_time_repl[4999:5002, ]

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


# -- compare distances used for computing partial effect of period obs vs bt model

identical(grass_ddmod_dist, grass_obs_ddmod_dist) #TRUE

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
grass_combo_qnts <- data.frame(tidyr::pivot_longer(grass_combo_qnts, cols = Alps_cmf:WesEu_bf,
                                                   names_to = 'ECO_NM',
                                                   values_to = 'Percentile'))

#split percentile columns into Q2.5% and Q97.5%
grass_combo_qnts <- data.frame(tidyr::pivot_wider(grass_combo_qnts, names_from = Q,
                                                  names_prefix = 'Q_',
                                                  values_from = Percentile))

#order levels of Dist columns
grass_combo_qnts$Dist <- factor(grass_combo_qnts$Dist, levels = c('1km', 'Median', 'Max.'))

#add column with estimated coef combo value from the obs model
grass_combo_qnts_est <- dplyr::left_join(x = grass_combo_qnts, y = obs_coef_combo_grass, by = c('Dist', 'ECO_NM'))

grass_qnts_est_plot <- ggplot(grass_combo_qnts_est[grass_combo_qnts_est$Dist != 'Max.', ], aes(x = Dist, y = Est, group = ECO_NM)) +
  geom_errorbar(aes(ymin = Q_2.5, ymax = Q_97.5), width = 0.1, lwd = 1.2, col = 'black') +
  geom_point(size = 2, col = 'green') +
  geom_hline(yintercept = 0, col = 'grey', lwd = 1) +
  ylab('Change in mean similarity between periods') + xlab(NULL) +
  scale_x_discrete(labels = c('1km' = 'Short distance',
                              'Median' = 'Median distance')) +
  ggtitle('Grassland') +
  facet_wrap(~ ECO_NM) +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 25, hjust = .5, vjust = .6))

ggsave(filename = 'Results_figs/DDmod_qnts_est_grass.jpeg', plot = grass_qnts_est_plot, device = 'jpeg',
       width = 30, height = 30, units = 'cm', dpi = 300)


# -- check on coef (not combos)

#check sign of coefs for Geographic distance (euc_dist) and Abs_diff_plot_size_sqrt
#these should be negative (or anyway super close to 0), otherwise there may be problems in fitting the models

exists('grass_coef') #F

grass_coef <- do.call(rbind, lapply(names(grass_ddmod_out), function(eco_nm) {
  
  coef_mat <- data.frame(do.call(rbind, lapply(grass_ddmod_out[[eco_nm]], function(i) i[['Coef']])))
  
  coef_mat$ECO_NM <- eco_nm
  
  return(coef_mat)
  
  }))


ggplot(grass_coef, aes(x = euc_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Geographic distance') +
  facet_wrap(~ ECO_NM, scales = 'free_x') +
  theme_pubr()

ggplot(grass_coef, aes(x = Abs_diff_plot_size_sqrt)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Absolute difference plot size (sqrt)') +
  facet_wrap(~ ECO_NM, scales = 'free_x') +
  theme_pubr()


# ----------------------------- forests

# -- create a list including results for the ecoregions

for_ddmod_nm <- list.files(path = 'DDmod_output/DDmod_for_out/drive-download-20261103T165405Z-1-001/', pattern = 'optim.RData', full.names = F)

for_ddmod_nm <- sapply(strsplit(for_ddmod_nm, split = '_'), function(i) paste(i[1], i[2], sep = '_'))

#empty lists
for_ddmod_out <- setNames(vector(mode = 'list', length = length(for_ddmod_nm)), nm = for_ddmod_nm)
for_ddmod_dist <- for_ddmod_out

#load mod out and distances in the env, and store them in lists

for(eco_nm in names(for_ddmod_out)) {
  
  load(paste0('DDmod_output/DDmod_for_out/drive-download-20261103T165405Z-1-001/', eco_nm, '_bt_out_RAM_optim.RData'))
  
  for_ddmod_out[[eco_nm]] <- list_out_fin
  for_ddmod_dist[[eco_nm]] <- distances_used
  
  rm(list_out_fin, distances_used)
  
  }

rm(eco_nm)

#check structure
length(for_ddmod_out$Alps_cmf) #5,000
length(for_ddmod_out$Alps_cmf[[1]]) #6
names(for_ddmod_out$Alps_cmf[[1]])

# -- check number of bt replicates

lengths(for_ddmod_out) #all 5,000

# -- check convergence

#check convergence
for_conv <- sapply(for_ddmod_out, function(i) all(sapply(i, function(.) .[['Conv']])))
all(for_conv)

#check number of iterations
for_iter <- lapply(for_ddmod_out, function(i) table(sapply(i, function(.) .[['Iter']]))) #mostly 5, but also cases of 6

# -- check time of execution

#extract Time_repl, which measures the time (in seconds) each replicate takes to be completed
for_time_repl <- stack(lapply(for_ddmod_out, function(i) sapply(i, function(.) .[['Time_repl']])))

#add col with replicate id
for_time_repl$id_repl <- rep(1:5000, times = length(unique(for_time_repl$ind)))
#for_time_repl[4999:5002, ]

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

# -- compare distances used for computing partial effect of period obs vs bt model

identical(for_obs_ddmod_dist, for_ddmod_dist) #TRUE

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
for_combo_qnts <- data.frame(tidyr::pivot_longer(for_combo_qnts, cols = Alps_cmf:WesEu_bf, 
                                                 names_to = 'ECO_NM', values_to = 'Percentile'))

#split percentile columns into Q2.5% and Q97.5%
for_combo_qnts <- data.frame(tidyr::pivot_wider(for_combo_qnts, names_from = Q, values_from = Percentile,
                                                names_prefix = 'Q_'))

#order levels of Dist columns
for_combo_qnts$Dist <- factor(for_combo_qnts$Dist, levels = c('1km', 'Median', 'Max.'))

#add column with estimates from obs model
for_combo_qnts_est <- dplyr::left_join(x = for_combo_qnts, y = obs_coef_combo_for, by = c('Dist', 'ECO_NM'))

#use the same y-ranges as for grasslands -> or show results from grasslands and forests in the same plot
for_qnts_est_plot <- ggplot(for_combo_qnts_est[for_combo_qnts_est$Dist != 'Max.', ], aes(x = Dist, y = Est, group = ECO_NM)) +
  geom_errorbar(aes(ymin = Q_2.5, ymax = Q_97.5), width = 0.1, lwd = 1.2, col = 'black') +
  geom_point(size = 2, col = 'darkgreen') +
  geom_hline(yintercept = 0, col = 'grey', lwd = 1) +
  ylab('Change in mean similarity between periods') + xlab(NULL) +
  scale_x_discrete(labels = c('1km' = 'Short distance', 'Median' = 'Median distance')) +
  ggtitle('Forest') +
  facet_wrap(~ ECO_NM) +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 25, hjust = .5, vjust = .6))

ggsave(filename = 'Results_figs/DDmod_qnts_est_for.jpeg', plot = for_qnts_est_plot, device = 'jpeg',
       width = 30, height = 30, units = 'cm', dpi = 300)

# -- check on coef (not combos)

#check sign of coefs for Geographic distance (euc_dist) and Abs_diff_plot_size_sqrt
#these should be negative (or anyway super close to 0), otherwise there may be problems in fitting the models

exists('for_coef') #F

for_coef <- do.call(rbind, lapply(names(for_ddmod_out), function(eco_nm) {
  
  coef_mat <- data.frame(do.call(rbind, lapply(for_ddmod_out[[eco_nm]], function(i) i[['Coef']])))
  
  coef_mat$ECO_NM <- eco_nm
  
  return(coef_mat)
  
}))


ggplot(for_coef, aes(x = euc_dist)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Geographic distance') +
  facet_wrap(~ ECO_NM, scales = 'free_x') +
  theme_pubr()

ggplot(for_coef, aes(x = Abs_diff_plot_size_sqrt)) +
  geom_histogram() +
  geom_vline(xintercept = 0, col = 'red') +
  xlab('Absolute difference plot size (sqrt)') +
  facet_wrap(~ ECO_NM, scales = 'free_x') +
  theme_pubr()



# ------ single plot showing dd-model results for grassland and forest

identical(colnames(grass_combo_qnts_est), colnames(for_combo_qnts_est)) #TRUE

grass_for_qnts_est_dtf <- rbind(data.frame(grass_combo_qnts_est, Habitat = 'Grassland'),
                                data.frame(for_combo_qnts_est, Habitat = 'Forest'))

grass_for_qnts_est_dtf$Habitat <- factor(grass_for_qnts_est_dtf$Habitat, levels = c('Grassland', 'Forest'))

#ggplot(grass_for_qnts_est_dtf[grass_for_qnts_est_dtf$Dist != 'Max.', ],
#       aes(x = Dist, y = Est, group = ECO_NM, col = Habitat)) +
#  geom_errorbar(aes(ymin = Q_2.5, ymax = Q_97.5), width = 0.1, lwd = 1.2, col = 'black') +
#  geom_point(size = 2) +
#  geom_hline(yintercept = 0, col = 'grey', lwd = 1) +
#  ylab('Change in mean similarity between periods') + xlab(NULL) +
#  scale_color_manual(values = c('Grassland' = 'green', 'Forest' = 'blue')) +
#  scale_x_discrete(labels = c('1km' = 'Short distance', 'Median' = 'Median distance')) +
#  facet_wrap(~ ECO_NM) +
#  theme_pubr() +
#  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
#        axis.text.x.bottom = element_text(size = 14, angle = 25, hjust = .5, vjust = .6))

#reverse order of ecoregion labels
grass_for_qnts_est_dtf$ECO_NM <- factor(grass_for_qnts_est_dtf$ECO_NM, 
                                        levels = sort(unique(grass_for_qnts_est_dtf$ECO_NM), decreasing = T))


grass_for_qnts_est_unq_plot <- ggplot(grass_for_qnts_est_dtf[grass_for_qnts_est_dtf$Dist != 'Max.', ],
                                      aes(x = Est, y = ECO_NM, color = Habitat, group = Habitat)) +
  geom_vline(xintercept = 0, col = 'grey', linewidth = 1) +
  geom_errorbar(aes(xmin = Q_2.5, xmax = Q_97.5), width = 0, linewidth = 2, color = 'black',
                position = position_dodge(width = 0.3)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('Grassland' = "#009E73", 'Forest' = "#A6761D")) +
  geom_segment(x = 0, xend = 0.4, y = .6, yend = .6, col = 'black', linewidth = 1.3,
               arrow = arrow(length = unit(0.5, 'cm'))) +
  annotate(geom = 'text', x = 0.6, y = .7, label = 'homog.', size = 10) +
  geom_segment(x = 0, xend = -0.4, y = .6, yend = .6, col = 'black', linewidth = 1.3, 
               arrow = arrow(length = unit(0.5, 'cm'))) +
  annotate(geom = 'text', x = -0.5, y = .7, label = 'diff.', size = 10) +
  facet_wrap(~ Dist, labeller = as_labeller(c('1km' = 'Short distance', 'Median' = 'Median distance'))) +
  ylab(NULL) + xlab('Change in mean similarity between periods') +
  theme_pubr() +
  theme(strip.text = element_text(size = 22), axis.text.y.left = element_text(size = 22),
        axis.title.x.bottom = element_text(size = 22), axis.text.x.bottom = element_text(size = 18), 
        legend.title = element_blank(), legend.text = element_text(size = 22))


ggsave(filename = 'Results_figs/DDmod_qnts_est_grass_for_unq.jpeg', plot = grass_for_qnts_est_unq_plot, device = 'jpeg',
       width = 45, height = 24, units = 'cm', dpi = 300)


