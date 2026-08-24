

#This script includes the R code to process and plot GDMs output.

library(ggplot2)
library(ggpubr)
library(ggrepel)


# ---- objects used in plots

#get vector of ecoregions shared by grasslands and forests
ecor_union_nm <- intersect(ecor_grass_nm, ecor_for_nm)

#get vector of ecoregions unique to grasslands
ecor_unq_grass_nm <- setdiff(ecor_grass_nm, ecor_for_nm)

#get vector of ecoregions unique to forests
ecor_unq_for_nm <- setdiff(ecor_for_nm, ecor_grass_nm)

#load data.frame with information on ecoregion position along elevation, long and lat gradients
load('/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Lon_lat_alt_ecoregions.RData')

#add ECO_NM col to lon_lat_alt_* data.frames
lon_lat_alt_grass$ECO_NM <- row.names(lon_lat_alt_grass)
lon_lat_alt_for$ECO_NM <- row.names(lon_lat_alt_for)

#check ecor names match
setdiff(x = lon_lat_alt_grass$ECO_NM, y = unique(dev_part_grass$ECO_NM)) #chr(0)
setdiff(x = unique(dev_part_grass$ECO_NM), y = lon_lat_alt_grass$ECO_NM) #chr(0)

setdiff(x = lon_lat_alt_for$ECO_NM, y = unique(dev_part_for$ECO_NM)) #chr(0)
setdiff(x = unique(dev_part_for$ECO_NM), y = lon_lat_alt_for$ECO_NM) #chr(0)


# ------------------------------ grasslands


# ------ deviance explained

expl_dev_grass_long <- data.frame(Expl_dev = as.vector(expl_dev_grass),
                                  ECO_NM = rownames(expl_dev_grass), #this gets recycled
                                  Period = rep(colnames(expl_dev_grass), each = nrow(expl_dev_grass))) 


#order ecoregions so that those included in both grass and for sets appear first
expl_dev_grass_long$ECO_NM <- factor(expl_dev_grass_long$ECO_NM, levels = c(sort(ecor_union_nm), sort(ecor_unq_grass_nm)))


expl_dev_grass_plot <- ggplot(data = expl_dev_grass_long, aes(x = ECO_NM, y = Expl_dev, fill = Period)) +
  geom_col(position = position_dodge2(width = .5)) +
  scale_fill_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  geom_hline(yintercept = 5, colour = 'orange', linetype = 'dashed', lwd = 1.2) +
  ylab('Explained deviance (%)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y.left = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16), legend.position = 'bottom',
        title = element_text(size = 18))

# ------ deviance partitions

class(dev_part_grass$VARIABLE_SET) #chr
unique(dev_part_grass$VARIABLE_SET) #

#consider unique contribution of climate and human, meaning that the part of deviance shared by the two components is excluded
#this should be the quantity more similar to the variable imporance computed by the gdm package
cl_hum_part_grass <- dev_part_grass[dev_part_grass$VARIABLE_SET %in% c('climate alone', 'human alone'), ]

#reformat data.frame to have separate fields for deviance in Period1 and Period2

cl_hum_part_grass_prd1 <- cl_hum_part_grass[cl_hum_part_grass$Period == 'Period1', ]
cl_hum_part_grass_prd2 <- cl_hum_part_grass[cl_hum_part_grass$Period == 'Period2', ]

#check fields actually match
identical(cl_hum_part_grass_prd1$ECO_NM, cl_hum_part_grass_prd2$ECO_NM) #T
identical(cl_hum_part_grass_prd1$VARIABLE_SET, cl_hum_part_grass_prd2$VARIABLE_SET) #T

#modify colnames for period-specific quantities
colnames(cl_hum_part_grass_prd1)[c(2, 3)] <- paste0(colnames(cl_hum_part_grass_prd1)[c(2, 3)], '_Prd1')
colnames(cl_hum_part_grass_prd2)[c(2, 3)] <- paste0(colnames(cl_hum_part_grass_prd2)[c(2, 3)], '_Prd2')

#cbind data.frames
cl_hum_part_grass <- cbind(cl_hum_part_grass_prd1[, c('ECO_NM', 'VARIABLE_SET', 'DEVIANCE_Prd1', 'DEVIANCE_scaled_Prd1')],
                           cl_hum_part_grass_prd2[, c('DEVIANCE_Prd2', 'DEVIANCE_scaled_Prd2')])


#rm period-specific datasets
rm(cl_hum_part_grass_prd1, cl_hum_part_grass_prd2)

#add data on ecor position along elev, long and lat gradients
cl_hum_part_grass <- dplyr::left_join(x = cl_hum_part_grass, y = lon_lat_alt_grass, by = 'ECO_NM')

#exclude ecoregions with explained deviance never equal to or larger than 5
which(rowSums((expl_dev_grass >= 5)*1) == 0) #EuAtl_mf

#exclude EuAtl_mf from cl_hum_part_grass
cl_hum_part_grass <- cl_hum_part_grass[cl_hum_part_grass$ECO_NM != 'EuAtl_mf', ]

#re-code columns with ecoregion position along elev, long and lat gradients
cl_hum_part_grass$X_ord <- as.numeric(as.factor(cl_hum_part_grass$X_ord))
cl_hum_part_grass$Y_ord <- as.numeric(as.factor(cl_hum_part_grass$Y_ord))
cl_hum_part_grass$Alt_ord <- as.numeric(as.factor(cl_hum_part_grass$Alt_ord))


cl_hum_grass_alt_plot <- ggplot(cl_hum_part_grass, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = Alt_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Elevation',
                        breaks = c(min(cl_hum_part_grass$Alt_ord), max(cl_hum_part_grass$Alt_ord)),
                        labels = c('Low elevation', 'High elevation'), option = 'plasma') +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Grassland - Elevation') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')


cl_hum_grass_long_plot <- ggplot(cl_hum_part_grass, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = X_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Longitude',
                        breaks = c(min(cl_hum_part_grass$X_ord), max(cl_hum_part_grass$X_ord)),
                        labels = c('Westward', 'Eastward')) +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Grassland - Longitude') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')


cl_hum_grass_lat_plot <- ggplot(cl_hum_part_grass, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = Y_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Latitude',
                        breaks = c(min(cl_hum_part_grass$Y_ord), max(cl_hum_part_grass$Y_ord)),
                        labels = c('Southward', 'Northward'), option = 'mako') +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Grassland - Latitude') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')


# ------ warping functions

#for each variable, scale estimated splines value by their corresponding maximum

isplines_grass_dtf <- do.call(rbind, lapply(isplines_grass, function(dtf) {
  
  #compute max value of estimated spline for each variable
  max_ispl <- tapply(dtf[['Values_y']], INDEX = list(dtf[['Variable_x']]), max)
  
  #add column with max estimated value of isplines
  dtf[['Max_y_value']] <- as.double(unname(max_ispl[dtf[['Variable_x']]]))
  
  #add column with scaled values of estimated splines
  dtf[['Scaled_y_values']] <- dtf[['Values_y']]/dtf[['Max_y_value']]
  
  #return result
  return(dtf)
  
}))

#modify row names
row.names(isplines_grass_dtf) <- as.character(seq_len(nrow(isplines_grass_dtf)))

#scale geographic distance to express it in km (rather than meters)
isplines_grass_dtf[isplines_grass_dtf$Variable_x == 'Geographic', 'Values_x'] <- isplines_grass_dtf[isplines_grass_dtf$Variable_x == 'Geographic', 'Values_x']/1000 

#create two plots: one including climate and hmi and another for the other predictors

# -- climate and hmi
isplines_grass_cl_hmi_dtf <- isplines_grass_dtf[isplines_grass_dtf$Variable_x %in% c('Tavg', 'Prcp', 'Hmi_value'), ]

#re-order levels of Variable_x
isplines_grass_cl_hmi_dtf$Variable_x <- factor(isplines_grass_cl_hmi_dtf$Variable_x, levels = c('Tavg', 'Prcp', 'Hmi_value'))

#Tavg, Prcp, Hmi_value
wfuncs_cl_hmi_grass_plot <- ggplot(isplines_grass_cl_hmi_dtf, aes(x = Values_x, y = Scaled_y_values, group = Period, col = Period)) +
  geom_line(lwd = 1.5) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_grid(ECO_NM ~ Variable_x, scales = 'free',
             labeller = labeller(Variable_x = c('Tavg' = 'Temperature (C°)',
                                                'Prcp' = 'Precipitation (mm)',
                                                'Hmi_value' = 'Land use'))) +
  ylab('Partial ecological distance (scaled)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubr() +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 16), strip.text.y.right = element_text(size = 14),
        strip.text.x.top = element_text(size = 16), legend.position = 'bottom', plot.title = element_text(size = 20))

ggsave(plot = wfuncs_cl_hmi_grass_plot, filename = 'Results_figs/Warping_functions_cl_hmi_grass.jpeg', device = 'jpeg',
       width = 26, height = 30, units = 'cm', dpi = 300)

# -- remaining predictors
isplines_grass_rem_dtf <- isplines_grass_dtf[!isplines_grass_dtf$Variable_x %in% c('Tavg', 'Prcp', 'Hmi_value'), ]

wfuncs_rem_grass_plot <- ggplot(isplines_grass_rem_dtf, aes(x = Values_x, y = Scaled_y_values, group = Period, col = Period)) +
  geom_line(lwd = 1.5) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_grid(ECO_NM ~ Variable_x, scales = 'free',
             labeller = labeller(Variable_x = c('Geographic' = 'Geo. distance (km)',
                                                'Releve_area_m2' = 'Plot size (m^2)',
                                                'Roughness' = 'Topo. roughness'))) +
  ylab('Partial ecological distance (scaled)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubr() +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 16), strip.text.y.right = element_text(size = 14),
        strip.text.x.top = element_text(size = 16), legend.position = 'bottom', plot.title = element_text(size = 20))

ggsave(plot = wfuncs_rem_grass_plot, filename = 'Results_figs/Warping_functions_rem_grass.jpeg', device = 'jpeg',
       width = 26, height = 30, units = 'cm', dpi = 300)

# ------------------------------ forests

#deviance explained

expl_dev_for_long <- data.frame(Expl_dev = as.vector(expl_dev_for),
                                ECO_NM = rownames(expl_dev_for),
                                Period = rep(colnames(expl_dev_for), each = nrow(expl_dev_for)))

#order ecoregions so that those included in both grass and for sets appear first
expl_dev_for_long$ECO_NM <- factor(expl_dev_for_long$ECO_NM, levels = c(sort(ecor_union_nm), sort(ecor_unq_for_nm)))

expl_dev_for_plot <- ggplot(data = expl_dev_for_long, aes(x = ECO_NM, y = Expl_dev, fill = Period)) +
  geom_col(position = position_dodge2(width = .5)) +
  scale_fill_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  geom_hline(yintercept = 5, colour = 'orange', linetype = 'dashed', lwd = 1.2) +
  ylab('Explained deviance (%)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubr() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
        axis.text.y.left = element_text(size = 14), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16), legend.position = 'bottom',
        title = element_text(size = 18))


# ------ deviance partitions

class(dev_part_for$VARIABLE_SET) #chr
unique(dev_part_for$VARIABLE_SET)

#consider unique contribution of climate and human
cl_hum_part_for <- dev_part_for[dev_part_for$VARIABLE_SET %in% c('climate alone', 'human alone'), ]

#reformat to have separate fields for Period1 and Period2
cl_hum_part_for_prd1 <- cl_hum_part_for[cl_hum_part_for$Period == 'Period1', ]
cl_hum_part_for_prd2 <- cl_hum_part_for[cl_hum_part_for$Period == 'Period2', ]

#check fields match
identical(cl_hum_part_for_prd1$ECO_NM, cl_hum_part_for_prd2$ECO_NM) #T
identical(cl_hum_part_for_prd1$VARIABLE_SET, cl_hum_part_for_prd2$VARIABLE_SET) #T

#modify columns for period-specific quantities
colnames(cl_hum_part_for_prd1)[c(2, 3)] <- paste(colnames(cl_hum_part_for_prd1)[c(2, 3)], 'Prd1', sep = '_')
colnames(cl_hum_part_for_prd2)[c(2, 3)] <- paste(colnames(cl_hum_part_for_prd2)[c(2, 3)], 'Prd2', sep = '_')

#cbind data.frames
cl_hum_part_for <- cbind(cl_hum_part_for_prd1[c('ECO_NM', 'VARIABLE_SET', 'DEVIANCE_Prd1', 'DEVIANCE_scaled_Prd1')],
                         cl_hum_part_for_prd2[c('DEVIANCE_Prd2', 'DEVIANCE_scaled_Prd2')])

#rm period-specific datasets
rm(cl_hum_part_for_prd1, cl_hum_part_for_prd2)

#add data on ecor position along elev, long and lat gradients
cl_hum_part_for <- dplyr::left_join(x = cl_hum_part_for, y = lon_lat_alt_for, by = 'ECO_NM')

#exclude ecor with dev expl never equal to or greater than 5%
which(rowSums((expl_dev_for >= 5)*1) == 0) #EuAtl_mf

#exclude EuAtl_mf from cl_hum_part_for
cl_hum_part_for <- cl_hum_part_for[cl_hum_part_for$ECO_NM != 'EuAtl_mf', ]

#re-code columns with ecoregion position along elev, long and lat gradients
cl_hum_part_for$X_ord <- as.numeric(as.factor(cl_hum_part_for$X_ord))
cl_hum_part_for$Y_ord <- as.numeric(as.factor(cl_hum_part_for$Y_ord))
cl_hum_part_for$Alt_ord <- as.numeric(as.factor(cl_hum_part_for$Alt_ord))


cl_hum_for_alt_plot <- ggplot(cl_hum_part_for, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = Alt_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Elevation',
                        breaks = c(min(cl_hum_part_for$Alt_ord), max(cl_hum_part_for$Alt_ord)),
                        labels = c('Low elevation', 'High elevation'), option = 'plasma') +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Forest - Elevation') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')


cl_hum_for_long_plot <- ggplot(cl_hum_part_for, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = X_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Longitude',
                        breaks = c(min(cl_hum_part_for$X_ord), max(cl_hum_part_for$X_ord)),
                        labels = c('Westward', 'Eastward')) +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Forest - Longitude') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')



cl_hum_for_lat_plot <- ggplot(cl_hum_part_for, aes(x = DEVIANCE_scaled_Prd1, y = DEVIANCE_scaled_Prd2)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey', lty = 'dashed') +
  geom_point(aes(colour = Y_ord), size = 8, alpha = .6) +
  geom_text_repel(aes(label = ECO_NM, size = 2), max.overlaps = Inf,
                  box.padding = .8, show.legend = FALSE, alpha = .8, segment.alpha = 0.6) +
  scale_color_viridis_c(name = 'Latitude',
                        breaks = c(min(cl_hum_part_for$Y_ord), max(cl_hum_part_for$Y_ord)),
                        labels = c('Southward', 'Northward'), option = 'mako') +
  xlab('Explained deviance - Period1 (%)') + ylab('Explained deviance - Period2 (%)') +
  ggtitle('Forest - Latitude') +
  facet_wrap(~ VARIABLE_SET, labeller = as_labeller(c('climate alone' = 'Climate', 'human alone' = 'Land use'))) +
  theme_pubr() +
  theme(plot.title = element_text(size = 18), legend.title = element_blank(), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16), axis.title = element_text(size = 14), legend.position = 'right')


# ------ warping functions

#scale estimated splines' value by their max
isplines_for_dtf <- do.call(rbind, lapply(isplines_for, function(dtf) {
  
  #compute max estimated ispline value
  max_ispl <- tapply(dtf[['Values_y']], INDEX = list(dtf[['Variable_x']]), max)
  
  #add field with max estimated value
  dtf[['Max_y_value']] <- as.double(unname(max_ispl[dtf[['Variable_x']]]))
  
  #add column with scaled values of estimated splines
  dtf[['Scaled_y_values']] <- dtf[['Values_y']]/dtf[['Max_y_value']]
  
  #return result
  return(dtf)
  
  }))


#modify row.names
row.names(isplines_for_dtf) <- as.character(seq_len(nrow(isplines_for_dtf)))

#scale geographic distance to express it in km
isplines_for_dtf[isplines_for_dtf$Variable_x == 'Geographic', 'Values_x'] <- isplines_for_dtf[isplines_for_dtf$Variable_x == 'Geographic', 'Values_x']/1000


#create two plots: one including climate and hmi and another for the remaining predictors

# -- climate and hmi

isplines_for_cl_hmi_dtf <- isplines_for_dtf[isplines_for_dtf$Variable_x %in% c('Tavg', 'Prcp', 'Hmi_value'), ]

#re-order levels of Variable_x
isplines_for_cl_hmi_dtf$Variable_x <- factor(isplines_for_cl_hmi_dtf$Variable_x, levels = c('Tavg', 'Prcp', 'Hmi_value'))

#Tavg, Prcp, Hmi_value
wfuncs_cl_hmi_for_plot <- ggplot(isplines_for_cl_hmi_dtf, aes(x = Values_x, y = Scaled_y_values, group = Period, col = Period)) +
  geom_line(lwd = 1.5) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_grid(ECO_NM ~ Variable_x, scales = 'free', labeller = labeller(Variable_x = c('Tavg' = 'Temperature (C°)',
                                                                                      'Prcp' = 'Precipitation (mm)',
                                                                                      'Hmi_value' = 'Land use'))) +
  ylab('Partial ecological distance (scaled)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubr() +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.title.y = element_text(size = 16),
        strip.text.y.right = element_text(size = 14), strip.text.x.top = element_text(size = 16),
        legend.position = 'bottom', plot.title = element_text(size = 20))


ggsave(plot = wfuncs_cl_hmi_for_plot, filename = 'Results_figs/Warping_functions_cl_hmi_for.jpeg', device = 'jpeg',
       width = 26, height = 30, units = 'cm', dpi = 300)


# -- remaining predictors

isplines_for_rem_dtf <- isplines_for_dtf[!isplines_for_dtf$Variable_x %in% c('Tavg', 'Prcp', 'Hmi_value'), ]

wfuncs_rem_for_plot <- ggplot(isplines_for_rem_dtf, aes(x = Values_x, y = Scaled_y_values, group = Period, col = Period)) +
  geom_line(lwd = 1.5) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_grid(ECO_NM ~ Variable_x, scales = 'free', labeller = labeller(Variable_x = c('Geographic' = 'Geo. distance (km)',
                                                                                      'Releve_area_m2' = 'Plot size (m^2)',
                                                                                      'Roughness' = 'Topo. roughness'))) +
  ylab('Partial ecological distance (scaled)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubr() +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.title.y = element_text(size = 16),
        strip.text.y.right = element_text(size = 14), strip.text.x.top = element_text(size = 16),
        legend.position = 'bottom', plot.title = element_text(size = 20))

ggsave(plot = wfuncs_rem_for_plot, filename = 'Results_figs/Warping_functions_rem_for.jpeg', device = 'jpeg',
       width = 26, height = 30, units = 'cm', dpi = 300)


# ------------------------------ combined grasslands and forests plots


expl_dev_combined <- ggarrange(expl_dev_grass_plot, expl_dev_for_plot, nrow = 2, common.legend = T, legend = 'right')

ggsave('Results_figs/Explained_dev_for_grass.jpeg', plot = expl_dev_combined, device = 'jpeg', dpi = 300,
       units = 'cm', width = 20, height = 15)


lon_lat_elev_grass_plot <- ggarrange(cl_hum_grass_long_plot, cl_hum_grass_lat_plot, cl_hum_grass_alt_plot, nrow = 3, ncol = 1)

ggsave('Results_figs/Dev_part_cl_hmi_grass.jpeg', plot = lon_lat_elev_grass_plot, device = 'jpeg', dpi = 300,
       units = 'cm', width = 26, height = 34)

lon_lat_elev_for_plot <- ggarrange(cl_hum_for_long_plot, cl_hum_for_lat_plot, cl_hum_for_alt_plot, nrow = 3, ncol = 1)

ggsave('Results_figs/Dev_part_cl_hmi_for.jpeg', plot = lon_lat_elev_for_plot, device = 'jpeg', dpi = 300,
       units = 'cm', width = 26, height = 34)

