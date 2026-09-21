#This script includes R code to export summary tables of best performing matching approaches

library(MatchIt)
library(flextable)
library(officer)

#example

exists('ex_table') #FALSE

#extract summary tables
ex_table <- summary(Alps_cmf_genetic1_ratio1_cal1_gr, interactions = T)
ex_table <- as.data.frame(ex_table$sum.matched)

#add Covariate column
ex_table$Covariate <- row.names(ex_table)

#re-order columns
ex_table <- ex_table[c('Covariate', setdiff(colnames(ex_table), 'Covariate'))]

#transform to felxtable
ex_table <- flextable(ex_table)
#round to 2 digit
ex_table <- colformat_double(ex_table, digits = 2)
#bold header text
ex_table <- bold(ex_table, bold = TRUE, part = 'header')
#set fontsize
ex_table <- fontsize(ex_table, size = 9, part = 'all')
#add inner borders
ex_table <- border_inner_h(ex_table, border = fp_border(color="gray", width = 1), part = 'body')
#set autofit with width = 1
ex_table <- set_table_properties(ex_table, layout = 'autofit', width = 1)
#save output
#save_as_docx('Alps_cmf grassland' = ex_table, path = '/Users/Manuele/Desktop/ex_table.docx')


#write a function that does the steps above
format_match_summ <- function(x) {
  require(MatchIt)
  require(officer)

  #extract summary table
  x <- as.data.frame(summary(x, interactions = T)$sum.matched)
  
  #add Covariate col and re-order cols
  x$Covariate <- row.names(x)
  x <- x[c('Covariate', setdiff(colnames(x), 'Covariate'))]
  
  #transform to flextable and modify appearence
  x <- flextable(x)
  x <- colformat_double(x, digits = 2)
  x <- bold(x, bold = TRUE, part = 'header')
  x <- fontsize(x, size = 9, part = 'all')
  x <- border_inner_h(x, border = fp_border(color="gray", width = 1), part = 'body')
  x <- set_table_properties(x, layout = 'autofit', width = 1)
  
  return(x)
  
  }

identical(ex_table, format_match_summ(x = Alps_cmf_genetic1_ratio1_cal1_gr)) #T

rm(ex_table)


# ------------- grasslands

# -- save each table separately

exists('Match_summary_tab_grass') #FALSE

Match_summary_tab_grass <- list(Alps_cmf = Alps_cmf_genetic1_ratio1_cal1_gr, 
                                Baltic_mf = Baltic_mf_genetic1_ratio1_cal1_gr,
                                Carp_mf = Carp_mf_genetic1_cal1_gr,
                                Celtic_bf = Celtic_bf_genetic1_ratio1_cal1_gr,
                                CenEu_mf = CenEu_mf_mahala1_ord1_gr,
                                EuAtl_mf = EuAtl_mf_mahala1_ord1_ratio1_cal1_gr,
                                ItaScl_sdf = ItaScl_sdf_mahala1_ord1_gr,
                                Pan_mf = Pan_mf_rob_mahala1_ord1_cal1_gr,
                                Sar_mf = Sar_mf_genetic1_ratio1_cal1_gr,
                                WesEu_bf = WesEu_bf_mahala1_ord1_cal1_gr)

Match_summary_tab_grass <- lapply(Match_summary_tab_grass, format_match_summ)

for(eco_nm in names(Match_summary_tab_grass)) {
  
  save_as_docx(Match_summary_tab_grass[[eco_nm]], path = paste0('Summary_tables_matching/Grasslands/', eco_nm, '_summ_tab_grass.docx'))
  
}

rm(eco_nm)

# -- create and save a unique table

exists('Match_unique_sumtab_grass') #FALSE

Match_unique_sumtab_grass <- list(Alps_cmf = Alps_cmf_genetic1_ratio1_cal1_gr, 
                                  Baltic_mf = Baltic_mf_genetic1_ratio1_cal1_gr,
                                  Carp_mf = Carp_mf_genetic1_cal1_gr,
                                  Celtic_bf = Celtic_bf_genetic1_ratio1_cal1_gr,
                                  CenEu_mf = CenEu_mf_mahala1_ord1_gr,
                                  EuAtl_mf = EuAtl_mf_mahala1_ord1_ratio1_cal1_gr,
                                  ItaScl_sdf = ItaScl_sdf_mahala1_ord1_gr,
                                  Pan_mf = Pan_mf_rob_mahala1_ord1_cal1_gr,
                                  Sar_mf = Sar_mf_genetic1_ratio1_cal1_gr,
                                  WesEu_bf = WesEu_bf_mahala1_ord1_cal1_gr)


Match_unique_sumtab_grass <- do.call(rbind, lapply(names(Match_unique_sumtab_grass), function(eco_nm) {
  
  #extract summary table
  x <- as.data.frame(summary(Match_unique_sumtab_grass[[eco_nm]], interactions = T)$sum.matched)
  
  #add Covariate and Ecoregion cols and re-order cols
  x$Covariate <- row.names(x)
  x$Ecoregion <- eco_nm
  x <- x[c('Ecoregion', 'Covariate', setdiff(colnames(x), c('Ecoregion', 'Covariate')))]
  
  return(x)
  
  }))


#transform to felxtable
Match_unique_sumtab_grass <- flextable(Match_unique_sumtab_grass)
#round to 2 digit
Match_unique_sumtab_grass <- colformat_double(Match_unique_sumtab_grass, digits = 2)
#bold header text
Match_unique_sumtab_grass <- bold(Match_unique_sumtab_grass, bold = TRUE, part = 'header')
#set fontsize
Match_unique_sumtab_grass <- fontsize(Match_unique_sumtab_grass, size = 10, part = 'header')
Match_unique_sumtab_grass <- fontsize(Match_unique_sumtab_grass, size = 9, part = 'body')
#add inner borders
Match_unique_sumtab_grass <- border_inner_h(Match_unique_sumtab_grass, border = fp_border(color="gray", width = 1), part = 'body')
#merge Ecoregion col
Match_unique_sumtab_grass <- merge_v(Match_unique_sumtab_grass, j = ~ Ecoregion)
#set autofit with width = 1
Match_unique_sumtab_grass <- set_table_properties(Match_unique_sumtab_grass, layout = 'autofit', width = 1)

save_as_docx('Grassland' = Match_unique_sumtab_grass, path = paste0('Summary_tables_matching/Grasslands/Unique_summtab_grass.docx'))

# ------------- forests

# -- save each table separately

exists('Match_summary_tab_for') #FALSE

Match_summary_tab_for <- list(Alps_cmf = Alps_cmf_mahala1_ord1_ratio1_for,
                              Carp_mf = Carp_mf_mahala1_ord1_for,
                              CenEu_mf = CenEu_mf_pscore1_ord1_cal1_for,
                              DinMon_mf = DinMon_mf_pscore1_ord1_for,
                              EuAtl_mf = EuAtl_mf_genetic1_ratio1_cal1_for,
                              ItaScl_sdf = ItaScl_sdf_genetic1_cal1_for,
                              Pan_mf = Pan_mf_genetic1_for,
                              Sar_mf = Sar_mf_mahala1_ord1_ratio1_cal1_for,
                              TyrAdr_smf = TyrAdr_smf_pscore1_ord1_cal1_for,
                              WesEu_bf = WesEu_bf_mahala1_ord1_for)


Match_summary_tab_for <- lapply(Match_summary_tab_for, format_match_summ)

for(eco_nm in names(Match_summary_tab_for)) {
  
  save_as_docx(Match_summary_tab_for[[eco_nm]], path = paste0('Summary_tables_matching/Forests/', eco_nm, '_summ_tab_for.docx'))
  
}

rm(eco_nm)

# -- create and save a unique table

exists('Match_unique_sumtab_for') #FALSE

Match_unique_sumtab_for <- list(Alps_cmf = Alps_cmf_mahala1_ord1_ratio1_for,
                                Carp_mf = Carp_mf_mahala1_ord1_for,
                                CenEu_mf = CenEu_mf_pscore1_ord1_cal1_for,
                                DinMon_mf = DinMon_mf_pscore1_ord1_for,
                                EuAtl_mf = EuAtl_mf_genetic1_ratio1_cal1_for,
                                ItaScl_sdf = ItaScl_sdf_genetic1_cal1_for,
                                Pan_mf = Pan_mf_genetic1_for,
                                Sar_mf = Sar_mf_mahala1_ord1_ratio1_cal1_for,
                                TyrAdr_smf = TyrAdr_smf_pscore1_ord1_cal1_for,
                                WesEu_bf = WesEu_bf_mahala1_ord1_for)


Match_unique_sumtab_for <- do.call(rbind, lapply(names(Match_unique_sumtab_for), function(eco_nm) {
  
  x <- as.data.frame(summary(Match_unique_sumtab_for[[eco_nm]], interactions = T)$sum.matched)
  
  #add Covariate and Ecoregion cols and re-order cols
  x$Covariate <- row.names(x)
  x$Ecoregion <- eco_nm
  x <- x[c('Ecoregion', 'Covariate', setdiff(colnames(x), c('Ecoregion', 'Covariate')))]
  
  return(x)
  
  }))


#transform to felxtable
Match_unique_sumtab_for <- flextable(Match_unique_sumtab_for)
#round to 2 digit
Match_unique_sumtab_for <- colformat_double(Match_unique_sumtab_for, digits = 2)
#bold header text
Match_unique_sumtab_for <- bold(Match_unique_sumtab_for, bold = TRUE, part = 'header')
#set fontsize
Match_unique_sumtab_for <- fontsize(Match_unique_sumtab_for, size = 10, part = 'header')
Match_unique_sumtab_for <- fontsize(Match_unique_sumtab_for, size = 9, part = 'body')
#add inner borders
Match_unique_sumtab_for <- border_inner_h(Match_unique_sumtab_for, border = fp_border(color="gray", width = 1), part = 'body')
#merge Ecoregion col
Match_unique_sumtab_for <- merge_v(Match_unique_sumtab_for, j = ~ Ecoregion)
#set autofit with width = 1
Match_unique_sumtab_for <- set_table_properties(Match_unique_sumtab_for, layout = 'autofit', width = 1)

save_as_docx('Forest' = Match_unique_sumtab_for, path = paste0('Summary_tables_matching/Forests/Unique_summtab_for.docx'))
