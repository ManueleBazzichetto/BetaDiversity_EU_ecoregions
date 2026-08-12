

# --------------------------------------- 1) Functions used for preparing tables formatted for GDMs

#-- function to obtain all possible pairs of duplicates

#test example
ex_to_del <- EVA_duply$Alps_cmf$Period1[7]
#split chr in elements and coerce to integer
ex_to_del <- unlist(strsplit(ex_to_del, split = '-'))
ex_to_del <- as.integer(ex_to_del)
#use the outer function to derive an N*N matrix of combinations
ex_to_del_mat <- outer(X = ex_to_del, Y = ex_to_del, FUN = function(x, y) paste(x, y, sep = '-'))
#check dimension
dim(ex_to_del_mat) #4 * 4 - this means I should have (4*3)/2 in the end
#extract upper diagonal
#!!Important: when looking for combos in the table formatted for GDMs
#both combinations of PlotID order must be looked for, that is Plot1-Plot2 and Plot2-Plot1
#the point is that the upper triangle contains only one possible combination of plot ids
#while both combos could be found in the formatted table
ex_to_del_mat_vec <- ex_to_del_mat[upper.tri(ex_to_del_mat)]
length(ex_to_del_mat_vec) #6

#write a function that replicates what achieved in the test example

combo_of_plotIDs <- function(x) {
  
  #from group of plot ids to integer
  x <- unlist(strsplit(x, split = '-'))
  
  #check duplicates
  if((sum(duplicated(x)) > 0)) stop('Found duplicated PlotIDs')
  
  #coerce to integer
  x <- as.integer(x)
  
  #derive combos
  x <- outer(X = x, Y = x, FUN = function(a, b) paste(a, b, sep = '-'))
  
  #extract vector of combos
  x <- x[upper.tri(x)]
  
  #return result
  return(x)
  
}

#check
identical(combo_of_plotIDs(EVA_duply$Alps_cmf$Period1[7]), 
          ex_to_del_mat_vec) #TRUE

rm(ex_to_del, ex_to_del_mat, ex_to_del_mat_vec)


#-- function for adding small noise to coordinates of spatial duplicates

#!remember to set.seed for reproducibility

#mthc_lst is Matched_datasets_*
#duply_ids is EVA_duply_unique (a list with PlotIDs of the spatial duplicates for ecoregions in both grass and forest sets)
#eps_sd is the standard deviation of the gaussian error (epsilon) added to coordinates 
#(default value of 1 will move coordinates by +-2 meters in long and lat 95% of cases)
#crs_proj is the crs of the coordinates

add_noise_coords <- function(mtch_lst, duply_ids, eps_sd = 1, crs_proj = 3035) {
  
  require(sf)
  
  #get names of ecoregions included in mtch_lst
  eco_nms <- names(mtch_lst)
  #if all names in eco_nms are present in duply_ids, subset duply_ids to only keep ecoregion in mtch_lst
  if(all(eco_nms %in% names(duply_ids))) duply_ids <- duply_ids[eco_nms] else stop('Eco names do not match')
  #check if lengths of the objects match
  if(!all(lengths(mtch_lst) == lengths(duply_ids))) stop('Lengths of mtch_lst and duply_ids are different')
  
  #loop over mtch_lst and duply_ids to add noise to coords in period-specific datasets
  res <- Map(function(a, b) {
    
    #loop over periods
    res_prd <- lapply(seq_along(a), function(i) {
      
      #PlotIDs of spatial duplicates for i-th period
      plot_ids <- b[[i]]
      
      #add noise to coordinates
      dtf <- a[[i]]
      
      #positions of spatial duplicates in dtf
      sp_dup_pos <- which(dtf[['PlotID']] %in% plot_ids)
      
      #create vectors representing small noise
      noise_x <- rnorm(n = length(sp_dup_pos), mean = 0, sd = eps_sd)
      noise_y <- rnorm(n = length(sp_dup_pos), mean = 0, sd = eps_sd)
      
      dtf[sp_dup_pos, 'X_laea'] <- dtf[sp_dup_pos, 'X_laea'] + noise_x
      
      dtf[sp_dup_pos, 'Y_laea'] <- dtf[sp_dup_pos, 'Y_laea'] + noise_y
      
      #check if there are new spatial duplicates due to adding noise
      
      #coerce dtf to a spatial obj
      dtf_sp <- st_as_sf(x = dtf[c('PlotID', 'X_laea', 'Y_laea')], coords = c('X_laea', 'Y_laea'), crs = crs_proj)
      
      new_duply <- st_equals(dtf_sp)
      
      new_duply <- sum(lengths(new_duply) > 1)
      
      message(paste('There are', new_duply, 'new duplicates', sep = ' '))
      
      return(dtf)
      
    })
    
    #rename res_prd
    names(res_prd) <- names(a)
    
    return(res_prd)
    
  }, a = mtch_lst, b = duply_ids)
  
  return(res)
  
  }


#test the function

test_to_del <- Matched_datasets_forest[1]

set.seed(3958)
test_res_del <- add_noise_coords(mtch_lst = test_to_del, duply_ids = EVA_duply_unique, eps_sd = 1, crs_proj = 3035)

pos_dup_prd1 <- which(test_to_del$Alps_cmf$Period1$PlotID %in% EVA_duply_unique$Alps_cmf$Period1)
pos_dup_prd2 <- which(test_to_del$Alps_cmf$Period2$PlotID %in% EVA_duply_unique$Alps_cmf$Period2)

set.seed(3958)
noise_x_to_del_prd1 <- rnorm(n = length(pos_dup_prd1), mean = 0, sd = 1)
noise_y_to_del_prd1 <- rnorm(n = length(pos_dup_prd1), mean = 0, sd = 1)

noise_x_to_del_prd2 <- rnorm(n = length(pos_dup_prd2), mean = 0, sd = 1)
noise_y_to_del_prd2 <- rnorm(n = length(pos_dup_prd2), mean = 0, sd = 1)

test_to_del$Alps_cmf$Period1$X_laea[pos_dup_prd1] <- test_to_del$Alps_cmf$Period1$X_laea[pos_dup_prd1] + noise_x_to_del_prd1
test_to_del$Alps_cmf$Period1$Y_laea[pos_dup_prd1] <- test_to_del$Alps_cmf$Period1$Y_laea[pos_dup_prd1] + noise_y_to_del_prd1

test_to_del$Alps_cmf$Period2$X_laea[pos_dup_prd2] <- test_to_del$Alps_cmf$Period2$X_laea[pos_dup_prd2] + noise_x_to_del_prd2
test_to_del$Alps_cmf$Period2$Y_laea[pos_dup_prd2] <- test_to_del$Alps_cmf$Period2$Y_laea[pos_dup_prd2] + noise_y_to_del_prd2

identical(test_res_del, test_to_del) #TRUE

identical(Matched_datasets_forest$Alps_cmf$Period1[-pos_dup_prd1, ], test_to_del$Alps_cmf$Period1[-pos_dup_prd1, ]) #TRUE
identical(Matched_datasets_forest$Alps_cmf$Period2[-pos_dup_prd2, ], test_to_del$Alps_cmf$Period2[-pos_dup_prd2, ]) #TRUE

rm(test_to_del, test_res_del, pos_dup_prd1, pos_dup_prd2, noise_x_to_del_prd1, noise_y_to_del_prd1, noise_x_to_del_prd2, noise_y_to_del_prd2)


#-- function for transforming veg data from long to wide format

veg_longtowide <- function(veg_lst, mtch_lst) {
  
  #check names are identical
  if(!identical(names(veg_lst), names(mtch_lst))) stop('Names do not match')
  
  res <- Map(function(x, y) {
    
    #loop over periods
    res_prd <- lapply(seq_along(x), function(i) {
      
      veg_prd <- x[[i]]
      mtch_dt <- y[[i]]
      
      
      veg_to_wd <- as.data.frame(tidyr::pivot_wider(data = veg_prd, names_from = 'Species_name', values_from = 'Cover_perc', values_fill = 0))
      
      if(identical(veg_to_wd[['PlotID']], mtch_dt[['PlotID']])) {
        
        veg_to_wd <- data.frame(veg_to_wd, mtch_dt[c('X_laea', 'Y_laea')])
        
      } else {
        
        stop("PlotID of veg and matched data do not match")
        
      }
      
      return(veg_to_wd)
      
      })
    
    names(res_prd) <- names(x)
    
    return(res_prd)
    
  }, x = veg_lst, y = mtch_lst)
  
  return(res)
  
  }

#check

test_to_del <- EVA_veg_grass['Sar_mf']
test_to_del_2 <- Matched_datasets_grass['Sar_mf']

test_res_del <- veg_longtowide(veg_lst = test_to_del, mtch_lst = test_to_del_2)

test_to_del$Sar_mf$Period1 <- as.data.frame(tidyr::pivot_wider(data = test_to_del$Sar_mf$Period1, names_from = 'Species_name', values_from = 'Cover_perc', values_fill = 0))
unique(sapply(test_to_del$Sar_mf$Period1, class)) #integer numeric
test_to_del$Sar_mf$Period2 <- as.data.frame(tidyr::pivot_wider(data = test_to_del$Sar_mf$Period2, names_from = 'Species_name', values_from = 'Cover_perc', values_fill = 0))

test_to_del$Sar_mf$Period1 <- data.frame(test_to_del$Sar_mf$Period1, test_to_del_2$Sar_mf$Period1[c('X_laea', 'Y_laea')])
test_to_del$Sar_mf$Period2 <- data.frame(test_to_del$Sar_mf$Period2, test_to_del_2$Sar_mf$Period2[c('X_laea', 'Y_laea')])

identical(test_to_del, test_res_del) #TRUE

rm(test_to_del, test_to_del_2, test_res_del)


#-- function for computing vif of predictors in Matched_datasets_*

#vif(lm(PlotID ~ Prcp + Tavg + Elevation + Roughness + Hmi_value + Releve_area_m2, data = Matched_datasets_grass$Alps_cmf$Period1))
#vif(lm(X_laea ~ Prcp + Tavg + Elevation + Roughness + Hmi_value + Releve_area_m2, data = Matched_datasets_grass$Alps_cmf$Period1))

check_multicoll <- function(mtch_lst, vars, vif_thr = 2) {
  
  require(car)
  
  if(!all(vars %in% colnames(mtch_lst[[c(1, 1)]]))) stop('Some of vars are not included in matched datasets')
  
  res <- do.call(rbind, lapply(names(mtch_lst), function(eco_nm) {
    
    
    res_prd <- do.call(rbind, lapply(names(mtch_lst[[eco_nm]]), function(prd_nm) {
      
      dtf <- mtch_lst[[c(eco_nm, prd_nm)]]
      
      prd_form <- as.formula(paste('PlotID', paste(vars, collapse = '+'), sep = '~'))
      
      prd_vif <- as.list(vif(lm(formula = prd_form, data = dtf)))
      
      pred_over_thr <- sum(sapply(prd_vif, function(i) i > vif_thr))
      
      prd_vif[['Period']] <- prd_nm
      
      prd_vif[['PredOverThr']] <- pred_over_thr
      
      prd_vif <- as.data.frame(prd_vif)
      
      return(prd_vif)
      
      }))
    
    
    res_prd[['ECO_NAME']] <- eco_nm
    
    return(res_prd)
    
    
  }))
  
  return(res)
  
  }


check_multicoll(mtch_lst = Matched_datasets_grass, vars = c('Prcp', 'Tavg', 'Elevation', 'Roughness', 'Hmi_value', 'Releve_area_m2'))


#-- function to drop dissimilarities (rows) between spatial duplicates from input table for GDMs 

#Important!! Check both Plot1-Plot2 and Plot2-Plot1 against duply_pairs
#Specifically, all possible orders of plot ids in *PlotID_cov columns
#this is because EVA_duply_pairs_* contains only one possible order
#(the upper triangle of the matrix obtained using outer() )
#however, it is possible that combo Plot1-Plot2 includes all pairs in the formatted table
#because PlotID were ordered in the Matched_datasets_* and plot combos were extracted as the upper triangle of the matrix computed as the 'outer product' of the PlotIDs

#function to drop pairs of duplicates

drop_unwanted_combos <- function(x, combos, col1 = 's1.PlotID_cov', col2 = 's2.PlotID_cov') {
  
  #paste columns col1-col2
  plot_combo1 <- paste(x[[col1]], x[[col2]], sep = '-')
  plot_combo2 <- paste(x[[col2]], x[[col1]], sep = '-')
  
  #get positions of combos to drop
  plot_to_exc <- which((plot_combo1 %in% combos) | (plot_combo2 %in% combos))
  
  x <- x[-(plot_to_exc), ]
  
  #print number of combos excluded
  message(paste('Removing', length(plot_to_exc), 'combos', sep = ' '))
  
  #rm objects that are no longer used
  rm(plot_combo1, plot_combo2)
  
  #drop combos columns -> this is done in the loop after removing plot pairs to reach cap
  #x[[col1]] <- NULL
  #x[[col2]] <- NULL
  
  #return result
  return(x)
  
}

#check
#Note: I originally tested WesEu_bf$Period1 to assess whether gdm could manage to handle the ecoregion with the largest number of plots
test_to_del <- Matched_datasets_grass$Sar_mf$Period1
test_to_del_2 <- EVA_veg_grass$Sar_mf$Period1
row.names(test_to_del_2) <- test_to_del_2$PlotID

test_to_del <- gdm::formatsitepair(bioData = test_to_del_2, bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                   XColumn = 'X_laea', YColumn = 'Y_laea', predData = test_to_del)

head(test_to_del)

#check structure of the SitePair table
str(test_to_del) #Releve_area_m2 is num; PlotID_cov is int


sum(paste(test_to_del$s1.PlotID_cov, test_to_del$s2.PlotID_cov, sep = '-') %in% EVA_duply_pairs$Sar_mf$Period1) #3,745 (out of 137,315 for both grass and for)
sum(paste(test_to_del$s2.PlotID_cov, test_to_del$s1.PlotID_cov, sep = '-') %in% EVA_duply_pairs$Sar_mf$Period1) #0

to_drop_del <- drop_unwanted_combos(x = test_to_del, combos = EVA_duply_pairs$Sar_mf$Period1)

any(EVA_duply_pairs$Sar_mf$Period1 %in% paste(to_drop_del$s1.PlotID_cov, to_drop_del$s2.PlotID_cov, sep = '-')) #F
any(EVA_duply_pairs$Sar_mf$Period1 %in% paste(to_drop_del$s2.PlotID_cov, to_drop_del$s1.PlotID_cov, sep = '-')) #F

any(paste(to_drop_del$s1.PlotID_cov, to_drop_del$s2.PlotID_cov, sep = '-') %in% EVA_duply_pairs$Sar_mf$Period1) #F
any(paste(to_drop_del$s2.PlotID_cov, to_drop_del$s1.PlotID_cov, sep = '-') %in% EVA_duply_pairs$Sar_mf$Period1) #F

#check equality with BR computed using vegdist

test_to_del_3 <- as.matrix(vegan::vegdist(x = test_to_del_2[setdiff(colnames(test_to_del_2), c('PlotID', 'X_laea', 'Y_laea'))]))

test_to_del_3[1:3, 1:3]

test_to_del_3 <- test_to_del_3[upper.tri(test_to_del_3)]

test_to_del_4 <- outer(X = test_to_del_2$PlotID, Y = test_to_del_2$PlotID, FUN = paste, sep = '-')

test_to_del_4 <- test_to_del_4[upper.tri(test_to_del_4)]

test_to_del_3 <- data.frame('distance' = test_to_del_3, 'plot_combo' = test_to_del_4)

test_to_del_3$plot1 <- sapply(strsplit(test_to_del_3$plot_combo, split = '-'), function(i) i[1])
test_to_del_3$plot2 <- sapply(strsplit(test_to_del_3$plot_combo, split = '-'), function(i) i[2])

test_to_del_3 <- test_to_del_3[order(test_to_del_3[[c('plot1')]]), ]

head(test_to_del_3)

all.equal(test_to_del$distance, test_to_del_3$distance) #TRUE

rm(test_to_del, test_to_del_2, to_drop_del, test_to_del_3, test_to_del_4)



# --------------------------------------- 2) Functions used for preparing tables for distance-decay models


#first, create an example on a small dataset outside the function

test_small_vegmat <- EVA_veg_forest$TyrAdr_smf$Period1
test_small_predmat <- Matched_datasets_forest$TyrAdr_smf$Period1

#PlotID is the first column; geo coords included in semi-last and last column
colnames(test_small_vegmat)[1] #"PlotID"
colnames(test_small_vegmat)[c((ncol(test_small_vegmat)-1), ncol(test_small_vegmat))] #"X_laea" "Y_laea"

#extract PlotID, the geo coords, and plot size from test_small_predmat - these fields will be used to create a tmp predData when using formatsitepair
test_small_predmat <- test_small_predmat[c("PlotID", "X_laea", "Y_laea", "Releve_area_m2", "PlotID_cov")]

#drop coords columns from the site x species matrix
test_small_vegmat <- test_small_vegmat[setdiff(colnames(test_small_vegmat), c("X_laea", "Y_laea"))]

#compute Bray-Curtis index using formatsitepair
#The table includes the weights column, as well as the geo coordinates for each plot pair. These columns can be excluded because
#the weights column is only relevant for fitting gdms, while the geo distance between plot pairs will be computed using sf.
test_small_bc <- gdm::formatsitepair(bioData = test_small_vegmat, bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                     XColumn = 'X_laea', YColumn = 'Y_laea', predData = test_small_predmat)


#select only relevant columns
test_small_bc <- test_small_bc[c('distance', 's1.Releve_area_m2', 's1.PlotID_cov', 's2.Releve_area_m2', 's2.PlotID_cov')]

#compute the geographical (Euclidean) distance between sites

#select only the PlotID and coords columns in test_small_predmat to create an sf object
test_small_sp <- st_as_sf(x = test_small_predmat[c('PlotID', 'X_laea', 'Y_laea')], coords = c('X_laea', 'Y_laea'), crs = 3035)

#--compute Euclidean distance between plot pairs
test_small_dist <- st_distance(x = test_small_sp, y = test_small_sp, by_element = F, which = 'Euclidean')

class(test_small_dist) #units

#check if a matrix is returned when the sf object has no crs assigned (see help of st_distance())
test_small_sp_v2 <- st_as_sf(x = test_small_predmat[c('PlotID', 'X_laea', 'Y_laea')], coords = c('X_laea', 'Y_laea'))

test_small_dist_v2 <- st_distance(x = test_small_sp_v2, y = test_small_sp_v2, by_element = F, which = 'Euclidean')

class(test_small_dist_v2) #Yes, a matrix array is returned (easier to work with it)

#check order of plot-pairs in distance matrix - this order must match that in table formatted for gdms

head(test_small_bc)
head(test_small_sp_v2)
head(test_small_dist_v2[1:10, 1:10])

#the lower triangle of the distance matrix is returned
st_distance(x = test_small_sp_v2[1, ], y = test_small_sp_v2[2, ], which = 'Euclidean') #403597-403706
st_distance(x = test_small_sp_v2[1, ], y = test_small_sp_v2[3, ], which = 'Euclidean') #403597-403707
st_distance(x = test_small_sp_v2[1, ], y = test_small_sp_v2[4, ], which = 'Euclidean') #403597-403710
st_distance(x = test_small_sp_v2[1, ], y = test_small_sp_v2[5, ], which = 'Euclidean') #403597-403711

head(test_small_bc[1237:1239, ])

st_distance(x = test_small_sp_v2[2, ], y = test_small_sp_v2[3, ], which = 'Euclidean') #403706-403707

#The lower.tri of the matrix matches the order produced by formatsitepair
class(test_small_dist_v2[lower.tri(test_small_dist_v2)])

#test equality between distance matrices returned by st_distance when setting (or non setting) crs in st_as_sf
all.equal(as.numeric(test_small_dist[lower.tri(test_small_dist)]), test_small_dist_v2[lower.tri(test_small_dist_v2)]) #T

#--use decostand to transform abundance matrix into 1/0 matrix

#check class of obj returned by decostand
class(decostand(x = test_small_vegmat[-1], method = 'pa')) #data.frame

#give it a look
#decostand(x = test_small_vegmat[-1], method = 'pa')[1:10, 1:10]

#create p/a version of test_small_vegmat
test_small_vegmat_pa <- cbind(test_small_vegmat[1], decostand(test_small_vegmat[-1], method = 'pa'))

identical(test_small_vegmat_pa[1], test_small_vegmat[1]) #T

#compute Jaccard and Horn-Morisita and create final output

#horn
test_small_horn <- gdm::formatsitepair(bioData = test_small_vegmat, bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                       XColumn = 'X_laea', YColumn = 'Y_laea', predData = test_small_predmat, dist = 'horn')

#jaccard
test_small_jac <- gdm::formatsitepair(bioData = test_small_vegmat_pa, bioFormat = 1, abundance = FALSE, siteColumn = 'PlotID',
                                      XColumn = 'X_laea', YColumn = 'Y_laea', predData = test_small_predmat, dist = 'jaccard')


#check the order of PlotID is maintained across datasets
identical(test_small_horn[c('s1.PlotID_cov', 's2.PlotID_cov')], test_small_jac[c('s1.PlotID_cov', 's2.PlotID_cov')]) #T

#create final output
test_final_out <- cbind(test_small_bc, horn = test_small_horn[[1]], jaccard = test_small_jac[[1]],
                        euc_dist = test_small_dist_v2[lower.tri(test_small_dist_v2)])

colnames(test_final_out)[1] <- 'bray'

#create a column including difference in plot size
test_final_out$Abs_diff_plot_size <- with(test_final_out, abs(s1.Releve_area_m2 - s2.Releve_area_m2))

#drop columns with individual plots' size
test_final_out$s1.Releve_area_m2 <- NULL
test_final_out$s2.Releve_area_m2 <- NULL

#re-order columns
test_final_out <- test_final_out[c(1, 4, 5, 6, 2, 3, 7)]

hist(test_final_out$Abs_diff_plot_size)
hist(sqrt(test_final_out$Abs_diff_plot_size))
plot(test_final_out$Abs_diff_plot_size, test_final_out$bray)


#function for computing several dissimilarity indices and geographic distance among plot pairs


get_distances_ind <- function(veg_mat, pred_mat) {
  
  #load required packages
  require(gdm)
  require(vegan)
  require(sf)
  
  #check that PlotID, X_laea, Y_laea, Releve_area_m2, and PlotID_cov are included in veg_mat and pred_mat
  check_cols <- c("PlotID", "X_laea", "Y_laea", "Releve_area_m2", "PlotID_cov")
  
  if(!isTRUE(all(check_cols[c(1, 2, 3)] %in% colnames(veg_mat)) & all(check_cols %in% colnames(pred_mat)))) stop('veg_mat must include PlotID, coords, while pred_mat also plot size and PlotID_cov')
  
  #select only check_cols in pred_mat
  pred_mat <- pred_mat[check_cols]
  
  #create spatial object for computing Euclidean distance between plot pairs
  predmat_sp <- st_as_sf(x = pred_mat[c("PlotID", "X_laea", "Y_laea")], coords = c("X_laea", "Y_laea"))
  
  #compute Euclidean distance matrix and return its upper triangle
  euc_dist <- st_distance(x = predmat_sp, y = predmat_sp, by_element = F, which = 'Euclidean')
  
  #extract lower triangle
  euc_dist <- euc_dist[lower.tri(euc_dist)]
  
  #rm predmat_sp
  rm(predmat_sp)
  
  #drop coordinates columns from veg_mat
  veg_mat <- veg_mat[setdiff(colnames(veg_mat), c("X_laea", "Y_laea"))]
  
  #compute indices
  res <- do.call(cbind, lapply(c('bray', 'horn', 'jaccard'), function(ind) {
    
    if(ind == 'jaccard') {
      
      #extract PlotID col and drop it from veg_mat
      plot_id_col <- veg_mat['PlotID']
      
      #drop PlotID col from veg_mat
      veg_mat <- veg_mat[setdiff(colnames(veg_mat), 'PlotID')]
      
      #use decostand to transform abund data to p/a
      veg_mat <- vegan::decostand(x = veg_mat, method = 'pa')
      
      #cbind plot_id_col
      veg_mat <- cbind(plot_id_col, veg_mat)
      
      #compute jaccard
      jac_mat <- gdm::formatsitepair(bioData = veg_mat, bioFormat = 1, abundance = FALSE, siteColumn = 'PlotID',
                                     XColumn = 'X_laea', YColumn = 'Y_laea', predData = pred_mat, dist = ind)
      
      #return res - only distance column, i.e. the value of the Jaccard index
      jac_mat <- jac_mat[c('distance')]
      
      #rename col
      colnames(jac_mat) <- ind
      
      return(jac_mat)
      
      } else {
        
        #compute ind
        ind_mat <- gdm::formatsitepair(bioData = veg_mat, bioFormat = 1, abundance = TRUE, siteColumn = 'PlotID',
                                     XColumn = 'X_laea', YColumn = 'Y_laea', predData = pred_mat, dist = ind)
        
        if(ind == 'horn') {
          
          #extract column including values of horn index
          ind_mat <- ind_mat[c('distance')]
          
          colnames(ind_mat) <- ind
          
          return(ind_mat)
          
          } else {
            
            #extract column including values of bray index
            colnames(ind_mat)[1] <- ind
            
            ind_mat <- ind_mat[c(ind, 's1.Releve_area_m2', 's1.PlotID_cov', 's2.Releve_area_m2', 's2.PlotID_cov')]
            
            return(ind_mat)
            
          }
        
        }
    
    }))
  
  #create col with abs difference in plot size and drop individual plot size columns
  res$Abs_diff_plot_size <- abs(res[['s1.Releve_area_m2']] - res[['s2.Releve_area_m2']])
  
  #also attach geo dist col
  res$euc_dist <- euc_dist
  
  #only keep id and plot pairs cols
  res <- res[c('bray', 'horn', 'jaccard', 'euc_dist', 's1.PlotID_cov', 's2.PlotID_cov', 'Abs_diff_plot_size')]
  
  #return result
  return(res)
  
  }


identical(x = test_final_out, y = get_distances_ind(veg_mat = EVA_veg_forest$TyrAdr_smf$Period1, pred_mat = Matched_datasets_forest$TyrAdr_smf$Period1)) #TRUE


rm(test_small_vegmat, test_small_predmat, test_small_bc, test_small_sp, test_small_dist, test_small_sp_v2,
   test_small_dist_v2, test_small_vegmat_pa, test_small_horn, test_small_jac, test_final_out)





