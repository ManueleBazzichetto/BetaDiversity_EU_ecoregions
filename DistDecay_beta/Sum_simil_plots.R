
#This script includes code for plotting summary stats on similarities computed in the Summary_simil script.

library(ggplot2)
library(ggpubr)


# -- objects used in plots

#get vector of ecoregions shared by grasslands and forests
ecor_union_nm <- intersect(names(summary_grass), names(summary_for))

#get vector of ecoregions unique to grasslands
ecor_unq_grass_nm <- setdiff(names(summary_grass), names(summary_for))

#get vector of ecoregions unique to forests
ecor_unq_for_nm <- setdiff(names(summary_for), names(summary_grass))

#levels for ECO_NM

# -- grassland
eco_nm_lev_grass <- c(sort(ecor_union_nm), sort(ecor_unq_grass_nm))

eco_nm_lev_for <- c(sort(ecor_union_nm), sort(ecor_unq_for_nm))

# ------------------ grasslands


# ---- Period summary

#rbind ecoregion-specific outputs
Period_summ_grass <- do.call(rbind, lapply(names(summary_grass), function(ecor_nm) {
  
  dtf <- summary_grass[[c(ecor_nm, 'Prd_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
  }))


Period_summ_grass$Period <- factor(Period_summ_grass$Period, levels = c('Period1', 'Period2'))
Period_summ_grass$Index <- factor(Period_summ_grass$Index, levels = c('bray', 'horn', 'jaccard'),
                                  labels = c('Bray-Curtis', 'Horn-Morisita', 'Jaccard'))
Period_summ_grass$ECO_NM <- factor(Period_summ_grass$ECO_NM, levels = eco_nm_lev_grass)


period_summ_grass_plot <- ggplot(Period_summ_grass, aes(x = Index, y = Median, col = Period)) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2, lwd = 1.2, position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_wrap(~ ECO_NM) +
  ylab('Median similarity (first, third quartiles)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 45, vjust = 1, hjust = 1))



# ---- Distribution of geographic distances

Dist_summ_grass <- do.call(rbind, lapply(names(summary_grass), function(ecor_nm) {
  
  dtf <- summary_grass[[c(ecor_nm, 'Dist_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
  }))


Dist_summ_grass$Period <- factor(Dist_summ_grass$Period, levels = c('Period1', 'Period2'))
Dist_summ_grass$ECO_NM <- factor(Dist_summ_grass$ECO_NM, levels = eco_nm_lev_grass)


dist_summ_grass_plot <- ggplot(Dist_summ_grass, aes(x = Period, y = Median, col = Period)) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2, lwd = 1.2, position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_wrap(~ ECO_NM) +
  ylab('Median geo. distance (first, third quartiles)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        legend.position = 'bottom')


# ---- Period summary along distances


Period_dist_summ_grass <- do.call(rbind, lapply(names(summary_grass), function(ecor_nm) {
  
  dtf <- summary_grass[[c(ecor_nm, 'Prd_dist_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
}))

Period_dist_summ_grass$Index <- factor(Period_dist_summ_grass$Index,levels = c('bray', 'horn', 'jaccard'),
                                       labels = c('Bray-Curtis', 'Horn-Morisita', 'Jaccard'))

Period_dist_summ_grass$Period <- factor(Period_dist_summ_grass$Period, levels = c('Period1', 'Period2'))
Period_dist_summ_grass$ECO_NM <- factor(Period_dist_summ_grass$ECO_NM, levels = eco_nm_lev_grass)


period_dist_summ_grass_plot <- ggplot(Period_dist_summ_grass, aes(x = Geo_bins, y = Median, col = Period, group = Period)) +
  geom_line(lwd = 1.2) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2) +
  geom_point(size = 2) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) + #plot every second label on the x-axis ticks
  facet_grid(Index ~ ECO_NM, scales = 'free_x') +
  ylab('Median similarity (first, third quartiles)') + xlab(NULL) + ggtitle('Grassland') +
  theme_pubclean() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_blank())


ggsave(filename = 'Results_figs/Period_dist_summary_grass.jpeg', plot = period_dist_summ_grass_plot, device = 'jpeg',
       width = 34, height = 28, units = 'cm', dpi = 300)


# ------------------ forests


# ---- Period summary

Period_summ_for <- do.call(rbind, lapply(names(summary_for), function(ecor_nm) {
  
  dtf <- summary_for[[c(ecor_nm, 'Prd_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
  }))


Period_summ_for$Period <- factor(Period_summ_for$Period, levels = c('Period1', 'Period2'))
Period_summ_for$Index <- factor(Period_summ_for$Index, levels = c('bray', 'horn', 'jaccard'),
                                labels = c('Bray-Curtis', 'Horn-Morisita', 'Jaccard'))

Period_summ_for$ECO_NM <- factor(Period_summ_for$ECO_NM, levels = eco_nm_lev_for)


period_summ_for_plot <- ggplot(Period_summ_for, aes(x = Index, y = Median, col = Period)) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2, lwd = 1.2, position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_wrap(~ ECO_NM) +
  ylab('Median similarity (first, third quartiles)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 45, vjust = 1, hjust = 1))


# ---- Distribution of geographic distances

Dist_summ_for <- do.call(rbind, lapply(names(summary_for), function(ecor_nm) {
  
  dtf <- summary_for[[c(ecor_nm, 'Dist_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
  }))


Dist_summ_for$Period <- factor(Dist_summ_for$Period, levels = c('Period1', 'Period2'))
Dist_summ_for$ECO_NM <- factor(Dist_summ_for$ECO_NM, levels = eco_nm_lev_for)


dist_summ_for_plot <- ggplot(Dist_summ_for, aes(x = Period, y = Median, col = Period)) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2, lwd = 1.2, position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  facet_wrap(~ ECO_NM) +
  ylab('Median geo. distance (first, third quartiles)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubr() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
        legend.position = 'bottom')



# ---- Period summary along distances

Period_dist_summ_for <- do.call(rbind, lapply(names(summary_for), function(ecor_nm) {
  
  dtf <- summary_for[[c(ecor_nm, 'Prd_dist_summ')]]
  
  dtf[['ECO_NM']] <- ecor_nm
  
  return(dtf)
  
  }))


Period_dist_summ_for$Index <- factor(Period_dist_summ_for$Index, levels = c('bray', 'horn', 'jaccard'),
                                labels = c('Bray-Curtis', 'Horn-Morisita', 'Jaccard'))
Period_dist_summ_for$Period <- factor(Period_dist_summ_for$Period, levels = c('Period1', 'Period2'))
Period_dist_summ_for$ECO_NM <- factor(Period_dist_summ_for$ECO_NM, levels = eco_nm_lev_for)


period_dist_summ_for_plot <- ggplot(Period_dist_summ_for, aes(x = Geo_bins, y = Median, col = Period, group = Period)) +
  geom_line(lwd = 1.2) +
  geom_errorbar(aes(ymin = fst_qrt, ymax = trd_qrt), width = .2) +
  geom_point(size = 2) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) + #plot every second label on the x-axis ticks
  facet_grid(Index ~ ECO_NM, scales = 'free_x') +
  ylab('Median similarity (first, third quartiles)') + xlab(NULL) + ggtitle('Forest') +
  theme_pubclean() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_blank())


ggsave(filename = 'Results_figs/Period_dist_summary_for.jpeg', plot = period_dist_summ_for_plot, device = 'jpeg',
       width = 34, height = 28, units = 'cm', dpi = 300)


# ------------------ combine grassland and forest plots


period_summ_comb_plot <- ggarrange(period_summ_grass_plot, period_summ_for_plot, nrow = 2, ncol = 1,
                                   common.legend = T, legend = 'right')


ggsave(filename = 'Results_figs/Period_summary_comb.jpeg', plot = period_summ_comb_plot, device = 'jpeg',
       width = 28, height = 32, units = 'cm', dpi = 300)


dist_summ_comb_plot <- ggarrange(dist_summ_grass_plot, dist_summ_for_plot, nrow = 2, ncol = 1,
                                 common.legend = T, legend = 'right')


ggsave(filename = 'Results_figs/Dist_summary_comb.jpeg', plot = dist_summ_comb_plot, device = 'jpeg',
       width = 28, height = 32, units = 'cm', dpi = 300)




