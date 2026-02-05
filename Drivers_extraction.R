#in this code, values of environmental drivers are obtained for all ecoregions

library(sf)
library(terra)
library(mapview)

#----------------load layers

#--climate

#import climatic stack already projected to epsg:3035
clim_stack_proj <- rast('EasyClimateData/clim_stack/clim_stack_proj3035.tif')

#rename layers in stack - note that clim_lyr_nm was created in the TopoClim_for_ecor script
nlyr(clim_stack_proj) == length(clim_lyr_nm)

names(clim_stack_proj) <- clim_lyr_nm

#--topography

#import proj stack of topography
topo_stack_proj <- rast('TopographyData/topo_stack_proj3035.tif')

#modify names of the layers
names(topo_stack_proj) <- c('Elevation', 'Roughness', 'Slope')

#--human modification index

#import hmi stack
hmi_stack_proj <- rast('HumanModificationIndex/hmi_stack_proj3035.tif')

#extract index year - strings all have the same pattern - note that hmi_nm was created in the TopoClim_for_ecor script
names(hmi_stack_proj) <- paste0('Hum_mod_ind_', substr(hmi_nm, start = 13, stop = 16))

#----------------------------------extract data for Grasslands and Forests

#for the moment, I'm extracting data for grasslands (EUNIS 1st lev 'R') and forests (EUNIS 1st lev 'T')
#because I'm conducting analyses on beta diversity patterns on these two habitat types

#----------------------------------Grasslands

#vector including names of ecoregions with minimum sample size for analyses
length(Grass_eco_minN) #14 ecoregions

Grass_sel_meta <- Grass_meta[Grass_meta$ECO_NAME %in% Grass_eco_minN, ]

#check range of Sampl_year within the ecoregions
lapply(Grass_eco_minN, function(i) range(as.integer(Grass_sel_meta[Grass_sel_meta$ECO_NAME == i, 'Sampl_year'])))


#----------------extract ext of the different ecoregions


#not needed if I'm not going to crop climatic stack
Grass_ext <- lapply(Grass_eco_minN, function(i) {
  
  plot_id <- Grass_sel_meta$PlotID[Grass_sel_meta$ECO_NAME == i]
  
  res <- ext(st_bbox(EVA_meta.sp[EVA_meta.sp$PlotID %in% plot_id, ]))
  
  res <- res + 1000
  
  return(res)
  
  })

names(Grass_ext) <- Grass_eco_minN


#----------------extract E-OBS climatic data


#probably faster if I directly use the Europe-wide layer instead of creating ecoregion-specific stacks
#crop climatic stack
#Grass_clim_stacks <- lapply(Grass_eco_minN, function(i) crop(x = clim_stack_proj, y = Grass_ext[[i]]))


#set width of temporal window for extracting climate data
temp_width <- 5

#create Cl_time_start column - start of the temporal window
#I'm doing this on the whole dataset for convenience

Grass_sel_meta$Cl_time_start <- with(Grass_sel_meta, (as.integer(Sampl_year) - (temp_width - 1)))

#create Climate_cellID column
Grass_cl_cellID <- terra::cellFromXY(object = clim_stack_proj[[1]], xy = as.matrix(Grass_sel_meta[c('X_laea', 'Y_laea')]))

#check NA
sum(is.na(Grass_cl_cellID)) #0

Grass_sel_meta$Climate_cellID <- Grass_cl_cellID

#extract climate data form cells
test_start <- proc.time()

Grass_sel_clmd <- extr_climate_eva_cells(x = Grass_sel_meta, cl_stack = clim_stack_proj,
                                       cellID_col = 'Climate_cellID', from_col = 'Cl_time_start', to_col = 'Sampl_year')

test_end <- proc.time() - test_start #approx 30 mins

#check number and position of missing values
sum(is.na(Grass_sel_clmd$Prcp)) #1226 (out of 72553)
sum(is.na(Grass_sel_clmd$Tavg)) #1219

setdiff(which(is.na(Grass_sel_clmd$Prcp)), which(is.na(Grass_sel_clmd$Tavg))) #9 cells
setdiff(which(is.na(Grass_sel_clmd$Tavg)), which(is.na(Grass_sel_clmd$Prcp))) #2 cells

#join climate data
identical(class(Grass_sel_meta$Climate_cellID), class(Grass_sel_clmd$CellID)) #T
identical(class(Grass_sel_meta$Cl_time_start), class(Grass_sel_clmd$Cl_time_start)) #T
identical(class(Grass_sel_meta$Sampl_year), class(Grass_sel_clmd$Cl_time_end)) #T

#prova <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_clmd, by = c('Climate_cellID' = 'CellID',
#                                                                         'Cl_time_start' = 'Cl_time_start',
#                                                                         'Sampl_year' = 'Cl_time_end'))
#all(mapply(identical, Grass_sel_meta, prova[colnames(Grass_sel_meta)])) #T
#rm(prova)

Grass_sel_meta <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_clmd, by = c('Climate_cellID' = 'CellID',
                                                                         'Cl_time_start' = 'Cl_time_start',
                                                                         'Sampl_year' = 'Cl_time_end'))

#check number of locations with NA for Prcp and/or Tavg
sum(is.na(Grass_sel_meta$Prcp)) #3079
sum(is.na(Grass_sel_meta$Tavg)) #3076

#plot ids for which both Prcp and Tavg are NA
Cl_missing_loc.gr <- intersect(Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Prcp))], Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Tavg))])

#unique for Prcp
Prcp_missing_loc.gr <- setdiff(Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Prcp))], Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Tavg))])

#unique for Tavg
Tavg_missing_loc.gr <- setdiff(Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Tavg))], Grass_sel_meta$PlotID[which(is.na(Grass_sel_meta$Prcp))])

#check
sum(duplicated(Cl_missing_loc.gr)) #0
sum(duplicated(Prcp_missing_loc.gr)) #0
sum(duplicated(Tavg_missing_loc.gr)) #0

#get info on locations missing cl data
Cl_missing_loc.gr <- Grass_sel_meta[Grass_sel_meta$PlotID %in% Cl_missing_loc.gr, c('PlotID', 'Prcp', 'Tavg', 'Cl_time_start', 'Sampl_year', 'X_laea', 'Y_laea')]

#check
#all(is.na(Cl_missing_loc.gr[[c('Prcp')]])) #T
#all(is.na(Cl_missing_loc.gr[[c('Tavg')]])) #T

Cl_missing_loc.gr <- fill_climate_gap(x = Cl_missing_loc.gr, cl_stack = clim_stack_proj, loc_id_col = 'PlotID',
                                   from_col = 'Cl_time_start', to_col = 'Sampl_year', rad_meters = ((res(clim_stack_proj)[1])+1))

#check remaining NAs
sum(is.na(Cl_missing_loc.gr$Prcp)) #1197
sum(is.na(Cl_missing_loc.gr$Tavg)) #1197
identical(which(is.na(Cl_missing_loc.gr$Prcp)), which(is.na(Cl_missing_loc.gr$Tavg))) #T

#fill in data
class(Grass_sel_meta$PlotID); class(Cl_missing_loc.gr$PlotID)

for(i in Cl_missing_loc.gr$PlotID) {
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Prcp'] <- Cl_missing_loc.gr[Cl_missing_loc.gr$PlotID == i, 'Prcp']
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Tavg'] <- Cl_missing_loc.gr[Cl_missing_loc.gr$PlotID == i, 'Tavg']
}

rm(i)

#check
sum(is.na(Grass_sel_meta$Prcp)) #1206
sum(is.na(Grass_sel_meta$Tavg)) #1203


#Plots in Prcp_missing_loc.gr (9 plots) and Tavg_missing_loc.gr (6 plots) were sampled in the same Sampl_year (2020 and 2014, respectively)
#Climatic data are missing for these specific years, although they are available for most of or all other years of the temporal window
#Grass_sel_meta[Grass_sel_meta$PlotID %in% Prcp_missing_loc.gr, c('Sampl_year')] #all 2020
#Grass_sel_meta[Grass_sel_meta$PlotID %in% Tavg_missing_loc.gr, c('Sampl_year')] #all 2014

#get coordinates of locations for which either Prcp or Tavg data are missing
Prcp_missing_loc.gr <- Grass_sel_meta[Grass_sel_meta$PlotID %in% Prcp_missing_loc.gr, c('X_laea', 'Y_laea', 'PlotID', 'Sampl_year')]
Tavg_missing_loc.gr <- Grass_sel_meta[Grass_sel_meta$PlotID %in% Tavg_missing_loc.gr, c('X_laea', 'Y_laea', 'PlotID', 'Sampl_year')]

prcp_tmp <- extract(clim_stack_proj[[grep(pattern = 'Prcp', x = names(clim_stack_proj))]],
                    y = Prcp_missing_loc.gr[, c('X_laea', 'Y_laea')])

#selecting years for which data are actually available between 2016 and 2020 - dropping 2020 due to NAs for all locations
prcp_tmp <- prcp_tmp[, paste0('Prcp_', 2016:2019)]

#summarise prcp values
prcp_tmp <- rowMeans(prcp_tmp)

#attach values to Prcp_missing_loc.gr
Prcp_missing_loc.gr$Prcp <- prcp_tmp

#same for Tavg
tavg_tmp <- extract(clim_stack_proj[[grep(pattern = 'Tavg', x = names(clim_stack_proj))]],
                    y = Tavg_missing_loc.gr[, c('X_laea', 'Y_laea')])

tavg_tmp <- tavg_tmp[, paste0('Tavg_', 2010:2013)]

tavg_tmp <- rowMeans(tavg_tmp)

Tavg_missing_loc.gr$Tavg <- tavg_tmp


#fill gaps in Grass_sel_meta
sum(is.na(Grass_sel_meta$Prcp)) #1206
sum(is.na(Grass_sel_meta$Tavg)) #1203

for(i in Prcp_missing_loc.gr$PlotID) {
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Prcp'] <- Prcp_missing_loc.gr[Prcp_missing_loc.gr$PlotID == i, 'Prcp']
  }

rm(i)

for(i in Tavg_missing_loc.gr$PlotID) {
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Tavg'] <- Tavg_missing_loc.gr[Tavg_missing_loc.gr$PlotID == i, 'Tavg']
}

rm(i)

#remaining NAs
sum(is.na(Grass_sel_meta$Prcp)) #1197 (out of 202595)
sum(is.na(Grass_sel_meta$Tavg)) #1197

#remove objects created to get climate data
rm(Grass_cl_cellID, test_start, Grass_sel_clmd, test_end, Cl_missing_loc.gr, Prcp_missing_loc.gr, Tavg_missing_loc.gr, prcp_tmp, tavg_tmp)


#----------------extract Geomorpho90m topography data

#create the Topo_cellID column
Grass_topo_cellID <- cellFromXY(object = topo_stack_proj[[1]], xy = as.matrix(Grass_sel_meta[c('X_laea', 'Y_laea')]))

#check NAs
sum(is.na(Grass_topo_cellID)) #0

#add Topo_cellID column to Grass_sel_meta
Grass_sel_meta$Topo_cellID <- Grass_topo_cellID

#extract topography values from xy
Grass_sel_topo <- extract(topo_stack_proj, y = unique(Grass_sel_meta$Topo_cellID))

#attach cell ID column
Grass_sel_topo$CellID <- unique(Grass_sel_meta$Topo_cellID)

#check NAs
sapply(Grass_sel_topo[c(1, 2, 3)], function(cl) sum(is.na(cl))) #Elev: 177; Rough: 202; Slope: 202

#Roughness and slope have NAs at the same position
identical(which(is.na(Grass_sel_topo$Roughness)), which(is.na(Grass_sel_topo$Slope))) #TRUE

#join topography values to Grass_sel_meta
#prova <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_topo, by = c('Topo_cellID' = 'CellID'))
#identical(Grass_sel_meta, prova[colnames(Grass_sel_meta)]) #T

Grass_sel_meta <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_topo, by = c('Topo_cellID' = 'CellID'))

#check number of locations with NAs
sapply(Grass_sel_meta[c('Elevation', 'Roughness', 'Slope')], function(cl) sum(is.na(cl))) #Elev: 431; Rough: 426; Slope: 426

#check overlap of locations missing topo values
length(intersect(which(is.na(Grass_sel_meta$Elevation)), which(is.na(Grass_sel_meta$Roughness)))) #246

#get coordinates of locations with NAs
Grass_elev_missing <- Grass_sel_meta[which(is.na(Grass_sel_meta$Elevation)), c('PlotID', 'X_laea', 'Y_laea')]

Grass_elev_missing <- data.frame(PlotID = Grass_elev_missing$PlotID,
                                 extract(x = topo_stack_proj[[c('Elevation')]], y = Grass_elev_missing[c('X_laea', 'Y_laea')],
                                         search_radius = (res(topo_stack_proj)[1] + 1)))

#select only relevant columns
Grass_elev_missing <- Grass_elev_missing[c('PlotID', 'Elevation')]

#exclude NAs
Grass_elev_missing <- Grass_elev_missing[!is.na(Grass_elev_missing$Elevation), ] #245 values

#same for slope and roughness
Grass_rghslo_missing <- Grass_sel_meta[which(is.na(Grass_sel_meta$Roughness)), c('PlotID', 'X_laea', 'Y_laea')]

Grass_rghslo_missing <- data.frame(PlotID = Grass_rghslo_missing$PlotID,
                                   extract(x = topo_stack_proj[[c('Roughness')]], y = Grass_rghslo_missing[c('X_laea', 'Y_laea')],
                                           search_radius = (res(topo_stack_proj)[1] + 1)),
                                   extract(x = topo_stack_proj[[c('Slope')]], y = Grass_rghslo_missing[c('X_laea', 'Y_laea')],
                                           search_radius = (res(topo_stack_proj)[1] + 1)))


Grass_rghslo_missing <- Grass_rghslo_missing[c('PlotID', 'Roughness', 'Slope')]

#exclude NAs
Grass_rghslo_missing <- na.omit(Grass_rghslo_missing) #248

#fill in gaps in Grass_sel_meta

#Elevation

for(i in Grass_elev_missing$PlotID) {
  
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Elevation'] <- Grass_elev_missing[Grass_elev_missing$PlotID == i, 'Elevation']
  
}

rm(i)

#check
sum(is.na(Grass_sel_meta$Elevation)) #186 (431 - 245)

#Roughness and slope

for(i in Grass_rghslo_missing$PlotID) {
  
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Roughness'] <- Grass_rghslo_missing[Grass_rghslo_missing$PlotID == i, 'Roughness']
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Slope'] <- Grass_rghslo_missing[Grass_rghslo_missing$PlotID == i, 'Slope']
  
}

rm(i)

#check
sum(is.na(Grass_sel_meta$Roughness)) #178 (426 - 248)
sum(is.na(Grass_sel_meta$Slope)) #178 (426 - 248)

rm(Grass_topo_cellID, Grass_sel_topo, Grass_elev_missing, Grass_rghslo_missing)


#----------------extract Human Modification Index values

#create the Hmi_cellID column
Grass_hmi_cellID <- cellFromXY(object = hmi_stack_proj[[1]], xy = as.matrix(Grass_sel_meta[c('X_laea', 'Y_laea')]))

#check NA
sum(is.na(Grass_hmi_cellID)) #0

#add Hmi_cellID to Grass_sel_meta
Grass_sel_meta$Hmi_cellID <- Grass_hmi_cellID

#add the Hmi_yr column, which maps Sampl_year to the closest available year of HM data
Grass_sel_meta$Hmi_yr <- get_hmi_year(x = Grass_sel_meta$Sampl_year) 

#extract hmi values from cells
test_start <- proc.time()

Grass_sel_hmi <- extr_hmi_eva_cells(x = Grass_sel_meta, hmi_stack = hmi_stack_proj, cellID_col = 'Hmi_cellID', hmi_yr_col = 'Hmi_yr') 

test_end <- proc.time() - test_start #approx 13 mins

#join hmi data to Grass_sel_meta
#prova <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_hmi, by = c('Hmi_cellID' = 'CellID', 'Hmi_yr' = 'Hmi_yr'))
#identical(prova[colnames(Grass_sel_meta)], Grass_sel_meta) #T

Grass_sel_meta <- dplyr::left_join(x = Grass_sel_meta, y = Grass_sel_hmi, by = c('Hmi_cellID' = 'CellID', 'Hmi_yr' = 'Hmi_yr'))

#check number of locations with NA
sum(is.na(Grass_sel_meta$Hmi_value)) #5326

#get coordinates of locations with NA
Grass_hmi_missing <- Grass_sel_meta[which(is.na(Grass_sel_meta$Hmi_value)), c('PlotID', 'Hmi_yr', 'X_laea', 'Y_laea')]

#fill gaps
Grass_hmi_missing <- fill_hmi_gap(x = Grass_hmi_missing, hmi_stack = hmi_stack_proj, loc_id_col = 'PlotID', rad_meters = (res(hmi_stack_proj)[1] + 1))

#check NA
sum(is.na(Grass_hmi_missing$Hmi_value)) #914

#exclude NAs
Grass_hmi_missing <- Grass_hmi_missing[!is.na(Grass_hmi_missing$Hmi_value), ]

#fill in hmi values in Grass_sel_meta

for(i in Grass_hmi_missing$PlotID) {
  
  Grass_sel_meta[Grass_sel_meta$PlotID == i, 'Hmi_value'] <- Grass_hmi_missing[Grass_hmi_missing$PlotID == i, 'Hmi_value']
  
}

rm(i)

#check NAs
sum(is.na(Grass_sel_meta$Hmi_value)) #914

rm(Grass_hmi_cellID, Grass_sel_hmi, Grass_hmi_missing)

#-----save Grass_sel_meta

#save Grass_sel_meta - note that this object still contains NAs for the environmental drivers
save(Grass_sel_meta, file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Grass_selection_meta.RData')


#----------------------------------Forests

#vector including names of ecoregions with minimum sample size for analyses
length(Forest_eco_minN) #11 ecoregions

Forest_sel_meta <- Forest_meta[Forest_meta$ECO_NAME %in% Forest_eco_minN, ]

#check range of Sampl_year within the ecoregions
lapply(Forest_eco_minN, function(i) range(as.integer(Forest_sel_meta[Forest_sel_meta$ECO_NAME == i, 'Sampl_year'])))


#----------------extract E-OBS climatic data

#set width of temporal window for extracting climate data - temp_width kept at 5 as for grasslands
#temp_width <- 5

#create Cl_time_start column - start of the temporal window
Forest_sel_meta$Cl_time_start <- with(Forest_sel_meta, (as.integer(Sampl_year) - (temp_width - 1)))

#create Climate_cellID column
Forest_cl_cellID <- terra::cellFromXY(object = clim_stack_proj[[1]], xy = as.matrix(Forest_sel_meta[c('X_laea', 'Y_laea')]))

#check NA
sum(is.na(Forest_cl_cellID)) #0

Forest_sel_meta$Climate_cellID <- Forest_cl_cellID

#extract climate data form cells
test_start <- proc.time()

Forest_sel_clmd <- extr_climate_eva_cells(x = Forest_sel_meta, cl_stack = clim_stack_proj,
                                         cellID_col = 'Climate_cellID', from_col = 'Cl_time_start', to_col = 'Sampl_year')

test_end <- proc.time() - test_start #approx 17 mins

#check number and position of missing values
sum(is.na(Forest_sel_clmd$Prcp)) #917 (out of 54674)
sum(is.na(Forest_sel_clmd$Tavg)) #841

setdiff(which(is.na(Forest_sel_clmd$Prcp)), which(is.na(Forest_sel_clmd$Tavg))) #77 cells
setdiff(which(is.na(Forest_sel_clmd$Tavg)), which(is.na(Forest_sel_clmd$Prcp))) #1 cell

#join climate data
identical(class(Forest_sel_meta$Climate_cellID), class(Forest_sel_clmd$CellID)) #T
identical(class(Forest_sel_meta$Cl_time_start), class(Forest_sel_clmd$Cl_time_start)) #T
identical(class(Forest_sel_meta$Sampl_year), class(Forest_sel_clmd$Cl_time_end)) #T


Forest_sel_meta <- dplyr::left_join(x = Forest_sel_meta, y = Forest_sel_clmd, by = c('Climate_cellID' = 'CellID',
                                                                                  'Cl_time_start' = 'Cl_time_start',
                                                                                  'Sampl_year' = 'Cl_time_end'))

#check number of locations with NA for Prcp and/or Tavg
sum(is.na(Forest_sel_meta$Prcp)) #1715
sum(is.na(Forest_sel_meta$Tavg)) #1621


#plot ids for which both Prcp and Tavg are NA
Cl_missing_loc.for <- intersect(Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Prcp))], Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Tavg))])

#unique for Prcp
Prcp_missing_loc.for <- setdiff(Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Prcp))], Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Tavg))])

#unique for Tavg
Tavg_missing_loc.for <- setdiff(Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Tavg))], Forest_sel_meta$PlotID[which(is.na(Forest_sel_meta$Prcp))])

#check
sum(duplicated(Cl_missing_loc.for)) #0
sum(duplicated(Prcp_missing_loc.for)) #0
sum(duplicated(Tavg_missing_loc.for)) #0


#get info on locations missing cl data
Cl_missing_loc.for <- Forest_sel_meta[Forest_sel_meta$PlotID %in% Cl_missing_loc.for, c('PlotID', 'Prcp', 'Tavg', 'Cl_time_start', 'Sampl_year', 'X_laea', 'Y_laea')]

#check
#all(is.na(Cl_missing_loc.for[[c('Prcp')]])) #T
#all(is.na(Cl_missing_loc.for[[c('Tavg')]])) #T

Cl_missing_loc.for <- fill_climate_gap(x = Cl_missing_loc.for, cl_stack = clim_stack_proj, loc_id_col = 'PlotID',
                                      from_col = 'Cl_time_start', to_col = 'Sampl_year', rad_meters = ((res(clim_stack_proj)[1])+1))


#check remaining NAs
sum(is.na(Cl_missing_loc.for$Prcp)) #1108
sum(is.na(Cl_missing_loc.for$Tavg)) #1050
identical(which(is.na(Cl_missing_loc.for$Prcp)), which(is.na(Cl_missing_loc.for$Tavg))) #FALSE

#fill in data
class(Forest_sel_meta$PlotID); class(Cl_missing_loc.for$PlotID)

for(i in Cl_missing_loc.for$PlotID) {
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Prcp'] <- Cl_missing_loc.for[Cl_missing_loc.for$PlotID == i, 'Prcp']
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Tavg'] <- Cl_missing_loc.for[Cl_missing_loc.for$PlotID == i, 'Tavg']
}

rm(i)

#check
sum(is.na(Forest_sel_meta$Prcp)) #1211 (1108 + 103 unique to Prcp)
sum(is.na(Forest_sel_meta$Tavg)) #1059 (1050 + 9 unique to Prcp)

#check Sampl_year of obs for which Prcp and Tavg are missing
Forest_sel_meta[Forest_sel_meta$PlotID %in% Prcp_missing_loc.for, 'Sampl_year'] #these are from different Sampl_year
Forest_sel_meta[Forest_sel_meta$PlotID %in% Tavg_missing_loc.for, 'Sampl_year'] #these are all from the same Sampl_year (2015)

#get coordinates of locations for which either Prcp or Tavg data are missing - then extract climate values for these locations

#Prcp
Prcp_missing_loc.for <- Forest_sel_meta[Forest_sel_meta$PlotID %in% Prcp_missing_loc.for, c('X_laea', 'Y_laea', 'PlotID', 'Sampl_year')]

#the prcp_tmp object was removed from the env
prcp_tmp <- extract(clim_stack_proj[[grep(pattern = 'Prcp', x = names(clim_stack_proj))]],
                    y = Prcp_missing_loc.for[, c('X_laea', 'Y_laea')])

#most of the plots have NAs for the whole temporal window (5-yrs) -> these will be assigned NA for Prcp
#the other plots have missing values for 1 or max 2 years -> these will be assigned the mean of available data over the time series

#add starting year of temp window to Prcp_missing_loc.for
Prcp_missing_loc.for$Start_yr <- (as.integer(Prcp_missing_loc.for$Sampl_year) - (t_win - 1))

missing_over_time_prcp <- do.call(rbind, lapply(seq_len(nrow(Prcp_missing_loc.for)), function(i) {
  
  #extract temporal window
  temp_wind <- Prcp_missing_loc.for[i, c('Start_yr', 'Sampl_year')]
  
  temp_wind <- seq.int(from = temp_wind[['Start_yr']], to = as.integer(temp_wind[['Sampl_year']]), by = 1)
  
  #extract data for the plot
  prcp_data <- unlist(prcp_tmp[i, paste0('Prcp_', temp_wind)])

  #count NAs
  prcp_nas <- sum(is.na(prcp_data))
  
  if(prcp_nas == t_win) prcp_val <- NA else prcp_val <- mean(prcp_data, na.rm = T)
  
  #return res
  res <- data.frame(Prcp_val = prcp_val, Prcp_nas = prcp_nas)
  
  return(res)
  
  }))

#add Prcp values
Prcp_missing_loc.for$Prcp <- missing_over_time_prcp$Prcp_val

#Tavg
Tavg_missing_loc.for <- Forest_sel_meta[Forest_sel_meta$PlotID %in% Tavg_missing_loc.for, c('X_laea', 'Y_laea', 'PlotID', 'Sampl_year')]

#tavg_tmp was deleted
tavg_tmp <- extract(clim_stack_proj[[grep(pattern = 'Tavg', x = names(clim_stack_proj))]],
                    y = Tavg_missing_loc.for[, c('X_laea', 'Y_laea')])

tavg_tmp <- tavg_tmp[, paste0('Tavg_', 2011:2015)]

#data are missing only in 2014
tavg_tmp <- rowMeans(tavg_tmp, na.rm = T)

Tavg_missing_loc.for$Tavg <- tavg_tmp

#fill gaps in Forest_sel_meta
sum(is.na(Forest_sel_meta$Prcp)) #1211
sum(is.na(Forest_sel_meta$Tavg)) #1059

for(i in Prcp_missing_loc.for$PlotID) {
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Prcp'] <- Prcp_missing_loc.for[Prcp_missing_loc.for$PlotID == i, 'Prcp']
}

rm(i)

for(i in Tavg_missing_loc.for$PlotID) {
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Tavg'] <- Tavg_missing_loc.for[Tavg_missing_loc.for$PlotID == i, 'Tavg']
}

rm(i)

#remaining NAs
sum(is.na(Forest_sel_meta$Prcp)) #1184 (out of 111815)
sum(is.na(Forest_sel_meta$Tavg)) #1050

#remove objects created to get climate data
rm(Forest_cl_cellID, test_start, Forest_sel_clmd, test_end, Cl_missing_loc.for, Prcp_missing_loc.for, Tavg_missing_loc.for, prcp_tmp, tavg_tmp, missing_over_time_prcp)


#----------------extract Geomorpho90m topography data

#create the Topo_cellID column
Forest_topo_cellID <- cellFromXY(object = topo_stack_proj[[1]], xy = as.matrix(Forest_sel_meta[c('X_laea', 'Y_laea')]))

#check NAs
sum(is.na(Forest_topo_cellID)) #0

#add Topo_cellID column to Forest_sel_meta
Forest_sel_meta$Topo_cellID <- Forest_topo_cellID

#extract topography values from xy
Forest_sel_topo <- extract(topo_stack_proj, y = unique(Forest_sel_meta$Topo_cellID))

#attach cell ID column
Forest_sel_topo$CellID <- unique(Forest_sel_meta$Topo_cellID)

#check NAs
sapply(Forest_sel_topo[c(1, 2, 3)], function(cl) sum(is.na(cl))) #Elev: 45; Rough: 70; Slope: 70

#Roughness and slope have NAs at the same position
identical(which(is.na(Forest_sel_topo$Roughness)), which(is.na(Forest_sel_topo$Slope))) #TRUE

#join topography data to Forest_sel_meta
Forest_sel_meta <- dplyr::left_join(x = Forest_sel_meta, y = Forest_sel_topo, by = c('Topo_cellID' = 'CellID'))

#check number of locations with NAs
sapply(Forest_sel_meta[c('Elevation', 'Roughness', 'Slope')], function(cl) sum(is.na(cl))) #Elev: 104; Rough: 132; Slope: 132

#check overlap of locations missing topo values
length(intersect(which(is.na(Forest_sel_meta$Elevation)), which(is.na(Forest_sel_meta$Roughness)))) #66

#get coordinates of locations with NAs
Forest_elev_missing <- Forest_sel_meta[which(is.na(Forest_sel_meta$Elevation)), c('PlotID', 'X_laea', 'Y_laea')]

Forest_elev_missing <- data.frame(PlotID = Forest_elev_missing$PlotID,
                                 extract(x = topo_stack_proj[[c('Elevation')]], y = Forest_elev_missing[c('X_laea', 'Y_laea')],
                                         search_radius = (res(topo_stack_proj)[1] + 1)))

#select only relevant columns
Forest_elev_missing <- Forest_elev_missing[c('PlotID', 'Elevation')]

#exclude NAs
Forest_elev_missing <- Forest_elev_missing[!is.na(Forest_elev_missing$Elevation), ] #47 values remaining (out 0f 104)


#same for slope and roughness
Forest_rghslo_missing <- Forest_sel_meta[which(is.na(Forest_sel_meta$Roughness)), c('PlotID', 'X_laea', 'Y_laea')]

Forest_rghslo_missing <- data.frame(PlotID = Forest_rghslo_missing$PlotID,
                                   extract(x = topo_stack_proj[[c('Roughness')]], y = Forest_rghslo_missing[c('X_laea', 'Y_laea')],
                                           search_radius = (res(topo_stack_proj)[1] + 1)),
                                   extract(x = topo_stack_proj[[c('Slope')]], y = Forest_rghslo_missing[c('X_laea', 'Y_laea')],
                                           search_radius = (res(topo_stack_proj)[1] + 1)))


Forest_rghslo_missing <- Forest_rghslo_missing[c('PlotID', 'Roughness', 'Slope')]

#exclude NAs
Forest_rghslo_missing <- na.omit(Forest_rghslo_missing) #71 values remaining out of 132

#fill in gaps in Forest_sel_meta

#Elevation

for(i in Forest_elev_missing$PlotID) {
  
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Elevation'] <- Forest_elev_missing[Forest_elev_missing$PlotID == i, 'Elevation']
  
}

rm(i)

#check
sum(is.na(Forest_sel_meta$Elevation)) #57 (104 - 47) 

#Roughness and slope

for(i in Forest_rghslo_missing$PlotID) {
  
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Roughness'] <- Forest_rghslo_missing[Forest_rghslo_missing$PlotID == i, 'Roughness']
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Slope'] <- Forest_rghslo_missing[Forest_rghslo_missing$PlotID == i, 'Slope']
  
}

rm(i)

#check
sum(is.na(Forest_sel_meta$Roughness)) #61 (132 - 71)
sum(is.na(Forest_sel_meta$Slope)) #61

rm(Forest_topo_cellID, Forest_sel_topo, Forest_elev_missing, Forest_rghslo_missing)


#----------------extract Human Modification Index values

#create the Hmi_cellID column
Forest_hmi_cellID <- cellFromXY(object = hmi_stack_proj[[1]], xy = as.matrix(Forest_sel_meta[c('X_laea', 'Y_laea')]))

#check NA
sum(is.na(Forest_hmi_cellID)) #0

#add Hmi_cellID to Forest_sel_meta
Forest_sel_meta$Hmi_cellID <- Forest_hmi_cellID

#add the Hmi_yr column, which maps Sampl_year to the closest available year of HM data
Forest_sel_meta$Hmi_yr <- get_hmi_year(x = Forest_sel_meta$Sampl_year) 

#extract hmi values from cells
test_start <- proc.time()

Forest_sel_hmi <- extr_hmi_eva_cells(x = Forest_sel_meta, hmi_stack = hmi_stack_proj, cellID_col = 'Hmi_cellID', hmi_yr_col = 'Hmi_yr') 

test_end <- proc.time() - test_start #approx 8 mins

#join hmi data to Forest_sel_meta

Forest_sel_meta <- dplyr::left_join(x = Forest_sel_meta, y = Forest_sel_hmi, by = c('Hmi_cellID' = 'CellID', 'Hmi_yr' = 'Hmi_yr'))

#check number of locations with NA
sum(is.na(Forest_sel_meta$Hmi_value)) #1744

#get coordinates of locations with NA
Forest_hmi_missing <- Forest_sel_meta[which(is.na(Forest_sel_meta$Hmi_value)), c('PlotID', 'Hmi_yr', 'X_laea', 'Y_laea')]

#fill gaps
Forest_hmi_missing <- fill_hmi_gap(x = Forest_hmi_missing, hmi_stack = hmi_stack_proj, loc_id_col = 'PlotID', rad_meters = (res(hmi_stack_proj)[1] + 1))

#check NA
sum(is.na(Forest_hmi_missing$Hmi_value)) #198

#exclude NAs
Forest_hmi_missing <- Forest_hmi_missing[!is.na(Forest_hmi_missing$Hmi_value), ]

#fill in hmi values in Forest_sel_meta

for(i in Forest_hmi_missing$PlotID) {
  
  Forest_sel_meta[Forest_sel_meta$PlotID == i, 'Hmi_value'] <- Forest_hmi_missing[Forest_hmi_missing$PlotID == i, 'Hmi_value']
  
}

rm(i)

#check NAs
sum(is.na(Forest_sel_meta$Hmi_value)) #198

rm(Forest_hmi_cellID, Forest_sel_hmi, test_start, test_end, Forest_hmi_missing)

#-----save Forest_sel_meta

#save Forest_sel_meta - note that this object still contains NAs for the environmental drivers
save(Forest_sel_meta, file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/Forest_selection_meta.RData')


