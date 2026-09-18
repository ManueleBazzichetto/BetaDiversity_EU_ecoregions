
#This script includes code for plotting summary stats on similarities computed in the Summary_simil script.

library(ggplot2)
library(ggpubr)
#packages to save tables
library(flextable) #installed on 17 Sept 2026
#library(officer) #installed on 17 Sept 2026

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


#create and save table of summary stats

#exclude column including mean values
#transform the table to wide format (needed for creating period-specific fields)
tab_prd_summ_grass <- as.data.frame(tidyr::pivot_wider(data = Period_summ_grass[setdiff(colnames(Period_summ_grass), 'Mean')], names_from = Period,
                                         values_from = c(Median, fst_qrt, trd_qrt), id_cols = c(ECO_NM, Index)))

#rename quartiles' columns
colnames(tab_prd_summ_grass)[c(1, 5, 6, 7, 8)] <- c('Ecoregion', 'Q1_Period1', 'Q1_Period2', 'Q3_Period1', 'Q3_Period2')
#transform to a flextable obj
tab_prd_summ_grass <- flextable(tab_prd_summ_grass)
#round values to 2nd digit
tab_prd_summ_grass <- colformat_double(tab_prd_summ_grass, digits = 2)
#create sub-header for Period1 and Period2
tab_prd_summ_grass <- separate_header(tab_prd_summ_grass)
#bold header
tab_prd_summ_grass <- bold(tab_prd_summ_grass, bold = TRUE, part = 'header')
#re-size header
tab_prd_summ_grass <- fontsize(tab_prd_summ_grass, part = 'header', size = 12)
#cell merging
tab_prd_summ_grass <- merge_v(tab_prd_summ_grass, j = ~ Ecoregion)
#add bottom borders to each row of the body
tab_prd_summ_grass <- border_inner_h(tab_prd_summ_grass, border = fp_border(color="gray", width = 1), part = 'body')
#autofit for a nicer layout
tab_prd_summ_grass <- set_table_properties(x = tab_prd_summ_grass, layout = 'autofit', width = 1)
#save output
save_as_docx(tab_prd_summ_grass, path = 'Results_tabs/Period_sum_grass_tab.docx')


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


#create and save table of summary stats

tab_prd_summ_for <- as.data.frame(tidyr::pivot_wider(data = Period_summ_for[setdiff(colnames(Period_summ_for), 'Mean')],
                                                     id_cols = c(ECO_NM, Index), names_from = Period,
                                                     values_from = c(Median, fst_qrt, trd_qrt)))
#rename cols
colnames(tab_prd_summ_for)[c(1, 5, 6, 7, 8)] <- c('Ecoregion', 'Q1_Period1', 'Q1_Period2', 'Q3_Period1', 'Q3_Period2')
#trnsf to flextable obj
tab_prd_summ_for <- flextable(tab_prd_summ_for)
#round values to 2nd digit
tab_prd_summ_for <- colformat_double(tab_prd_summ_for, digits = 2)
#create sub-header for Period1 and 2
tab_prd_summ_for <- separate_header(tab_prd_summ_for)
#bold header
tab_prd_summ_for <- bold(tab_prd_summ_for, bold = TRUE, part = 'header')
#re-size font
tab_prd_summ_for <- fontsize(tab_prd_summ_for, size = 12, part = 'header')
#merge cells
tab_prd_summ_for <- merge_v(tab_prd_summ_for, j = ~ Ecoregion)
#add bottom borders to body
tab_prd_summ_for <- border_inner_h(tab_prd_summ_for, border = fp_border(color="gray", width = 1), part = 'body')
#autofit for a nicer layout
tab_prd_summ_for <- set_table_properties(tab_prd_summ_for, layout = 'autofit', width = 1)
#save output
save_as_docx(tab_prd_summ_for, path = 'Results_tabs/Period_sum_for_tab.docx')



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


#create another plot that uses the same range of the y-axis than the plot for grassland data

#extract y range of grassland plot
round(layer_scales(period_dist_summ_grass_plot)$y$range$range, digits = 2) #0 - 0.31

period_dist_summ_for_plot_zoom <- ggplot(Period_dist_summ_for, aes(x = Geo_bins, y = Median, col = Period, group = Period)) +
  geom_line(lwd = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c('Period1' = 'grey', 'Period2' = 'purple')) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) + #plot every second label on the x-axis ticks
  facet_grid(Index ~ ECO_NM, scales = 'free_x') +
  ylab('Median similarity') + xlab(NULL) + ggtitle('Forest') +
  ylim(c(0, 0.31)) +
  theme_pubclean() +
  theme(title = element_text(size = 18), strip.text = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_blank())

ggsave(filename = 'Results_figs/Period_dist_summary_for_zoom.jpeg', plot = period_dist_summ_for_plot_zoom, device = 'jpeg',
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




