#this code provides an example on how to extract climatic and topographic data for the ecoregions
#here, I'm using a single ecoregion to test the code

#spatial data
library(sf)
library(terra)

#viz
library(ggplot2)
library(ggpubr)

#Notes

##It's probably better to project the climatic layers for whole Europe to CRS 3035 so that their resolution is the same 
##in each ecoregion - otherwise if projection is carried out separately,
##layers of different ecoregions will have different spatial resolutions, which will mostly depend on longitude (and less on latitude)

##There are 2 alternative options to extract environmental data from each EVA location:
##1) point based (not used anymore because it's inefficient): the values of the env variables are simply extracted at the coords of each location; [REMOVED FROM THE CODE]
##2) cell based: many EVA locations are included in the same cells (because points are clustered in space).
##To speed up data extraction, data are extracted for each unique combo of cell id + period (or time window).
##data are then joined to the EVA metadata data frame for the specific ecoregion

##Issue below has been already solved at this stage - using the GSHHG dataset
##there are points that fall in the sea (obviously only if the ecoregion has a sea-land border)
##one option could be to use the shp of the EU EMODnet coastline to remove points that fall too far from the coast
##however, doing this is problematic because the EMODnet shp is not closed, it has many holes, many small islands,
##and polygons (probably due to raster converted to lines) along the coast.
##So, it's probably easier to extract data using the search_radius argument for points falling at sea
##if points are too far from the coast, then it means they are wrongly located. These points will be assigned NA
##for the env variables, and will eventually be dropped out as they won't have info for predictors.

#Perhaps better to create a buffer around the boundary of the ecoregion polygon to be sure to extract data
#from areas not covered by the polygon due to its coarse resolution.
#using a buffered polygon may anyway lead to including points from neighbor ecoregions
#for the moment, I'm using the search_radius option (see notes above)

#----------------------------------extract polygon for a test ecoregion

#I'm using the "Cantabrian_mixed_forests" ecoregion as a test area
Cant_mf <- Grass_meta[Grass_meta$ECO_NAME == "Cantabrian_mixed_forests", ]

#get polygon of the Cant_mf ecoregion
Cant_mf_poly <- eu_ecoregions.proj[eu_ecoregions.proj$ECO_NAME == "Cantabrian mixed forests", ]

mapview(Cant_mf_poly) + mapview(EVA_meta.sp[EVA_meta.sp$ECO_NAME == "Cantabrian_mixed_forests", ])

#create ext/bbox of the ecoregion based on locations - this will be used to crop stacks at the extent of the ecoregions
#notice that some locations may fall outside the bbox of the ecoregion polygon and yet be assigned to that ecoregion by proximity
#therefore, using the ext/bbox of the points to crop the stacks avoids leaving these points out
#also, note that here I'm considering the extent defined by all locations associated with the Cantabrian_mixed_forests ecoregion,
#even though the code below extract data only for grasslands belonging to this ecoregion
#an alternative would be to get the extent of the locations associated with a single habitat type
Cant_mf_ext <- ext(st_bbox(EVA_meta.sp[EVA_meta.sp$ECO_NAME == 'Cantabrian_mixed_forests', ]))

#adding 1 km around the ext to be sure all points are included
Cant_mf_ext <- Cant_mf_ext + 1000

#----------------------------------extract E-OBS climatic data (downscaled at 1 km spatial resolution)

#notice that E-OBS data (same as those used by easyclimate v.4) are not available for 2023

#note that data can't be directly downloaded via the easyclimate R package because the spatial extent of the ecoregion and # of plots per ecoregion are too large
#data must be downloaded as yearly climatic layers from the Boku University's server. Climatic values from each location are then extracted from the layers
#climate data in Boku's server are currently (i.e., June 2025) available until 2022

#check range of Sampl_years
range(as.integer(Cant_mf$Sampl_year)) #1980-2021

#check number of years for which EVA data are available
length(unique(as.integer(Cant_mf$Sampl_year))) #40
length(seq(1980, 2022, 1)) #43 - this means EVA data are not available for all years

#import climatic stack already projected to epsg:3035
clim_stack_proj <- rast('EasyClimateData/clim_stack/clim_stack_proj3035.tif')

#extract name of climate var and year of the data from layer
clim_lyr_nm <- names(clim_stack_proj)

#extract number(s) (year in this case) from a string: https://www.r-bloggers.com/2024/06/extracting-numbers-from-strings-in-r/
#note that this will work also when there's more than a number in the string, e.g. 1980AA2000 -> it returns 1980, 2000
yr_from_str <- gregexpr("[0-9]+", clim_lyr_nm)
yr_from_str <- regmatches(clim_lyr_nm, yr_from_str)
yr_from_str <- as.integer(unlist(yr_from_str))

#extract name of the variable
#note that all climatic variables have 4-letter acronyms and the pattern used for file names is consistent
unique(substr(clim_lyr_nm, start = 11, stop = 14)) #"Prcp" "Tavg"

clim_var_nm <- substr(clim_lyr_nm, start = 11, stop = 14)

#create vector of names for layers
clim_lyr_nm <- if(isTRUE(length(yr_from_str) == length(clim_var_nm))) paste(clim_var_nm, yr_from_str, sep = "_")

#rename layers in stack
nlyr(clim_stack_proj) == length(clim_lyr_nm)

names(clim_stack_proj) <- clim_lyr_nm

#crop clim_stack_proj to extent of test ecoregion
#!Check snap options! Now using default option 'near'
#!!!!check misalignment between shp of ecoregions and clim_bmf!!!!
#the misalignment should not affect climatic data for analyses because data are extracted at the cell/point level
#moreover, I'm cropping at the extent of the ecoregion points (i.e., I'm not masking)
Cant_mf_clim <- crop(x = clim_stack_proj, y = Cant_mf_ext)

#--------------------------set temporal window for aggregating climate data

#final length of the temporal window
t_win <- 5

#check how many data points occur before and in 1980, and after 2022
anyNA(Cant_mf$Sampl_year) #F

sum(as.integer(Cant_mf$Sampl_year) > 2022) #0
sum(as.integer(Cant_mf$Sampl_year) < 1980) #0
sum(as.integer(Cant_mf$Sampl_year) == 1980) #37

#compute the temporal window to gather climate data for each data point outside the function for extracting the data
#should speed up the process
#the temporal window has length (t_win - 1):Sampl_year
Cant_mf$Cl_time_start <- with(Cant_mf, (as.integer(Sampl_year) - (t_win - 1)))

#check whether the first year in Cl_time_start is covered by E-OBS climatic data (it should)
min(Cant_mf$Cl_time_start) #1976


#------------function for filling gaps in climatic data due to locations falling in cells with no data (NA)

#First, some info on the search_radius argument of terra::extract:

#search_radius option info
#search_radius: positive number. To be used if y has point geometry. The positive number represents
#the max distance to seek the closest non-NA cell starting from the NA cell where the point falls in.
#in other words, the distance is a distance between cells (and not a point-to-cell distance)
#the radius should be expressed in meters in case of longlat data (n.b. imprecise for large lat span - see help)
#or in distance unit otherwise (typically also meter)

#check functioning of search_radius option
ex_lyr <- rast()
#nrow(ex_lyr) <- ncol(ex_lyr) <- 2 
ex_lyr[] <- c(NA_integer_, 10L, 20L, 30L)
terra::extract(x = ex_lyr, y = vect(matrix(c(0, 0), nrow = 1, ncol = 2), type = 'points'), search_radius = 200000)

anyNA(extract(x = Cant_mf_clim[[c(1)]], y = Cant_mf[c('X_laea', 'Y_laea')])) #T
which(is.na(extract(x = Cant_mf_clim[[c(1)]], y = Cant_mf[c('X_laea', 'Y_laea')])[[2]])) #13
terra::extract(x = Cant_mf_clim[[c(1)]], y = Cant_mf[13, c('X_laea', 'Y_laea')], search_radius = 1000)

#The function for filling gaps in climatic data can be used for both the point- and the cell-based approach

#rad_meters is the argument internally passed to terra::extract(..., search_radius)
#it should be set based on the resolution of the grid.
#clim_stack_proj has res 654.7999 X 654.7999
#as an example, it could be set so that data are searched within max 2 cells distance from the cell where the point falls
#this means a max of round(654.7999*2, 0) = 1310 meters (approx 1300)

#the function sets the search_radius option only when NA is returned for a year-specific climatic layer
#this means that if the gap occurs only for a subset of the whole period, the function
#will extract the missing info only for the subset
#This avoids assuming that NAs are spread over the whole period.

fill_climate_gap <- function(x, cl_stack, loc_id_col, from_col, to_col, rad_meters) {
  require(terra)
  #get plot ids
  loc_id <- x[[loc_id_col]]
  #get name of clim vars in cl_stack
  cl_vars <- unique(substr(x = names(cl_stack), start = 1, stop = 4))
  #fill gaps in climatic data
  #loop across plot IDs
  res <- do.call(rbind, lapply(loc_id, function(id) {
    #get id position
    loc_pos <- which(loc_id == id)
    #get coordinates of the point
    loc_xy <- x[loc_pos, c("X_laea", "Y_laea")]
    #set the time window
    loc_win <- seq.int(from = as.integer(x[loc_pos, from_col]), to = as.integer(x[loc_pos, to_col]), by = 1)
    #retrieve name of layers for which data are needed
    loc_lyrs <- paste(rep(cl_vars, each = length(loc_win)), loc_win, sep = '_')
    #get missing data - loop across the names of the layers
    res_cl <- unlist(lapply(loc_lyrs, function(nm) {
      #extract layer
      lyr_touse <- cl_stack[[nm]]
      #extract value for the point
      lyr_val <- extract(x = lyr_touse, y = loc_xy, ID = F)
      #if the value is NA use search_radius - notice that in case of a single point (ID = F), extract returns a 1 X 1 df
      if(anyNA(lyr_val)) {
        #when a search_radius is set and extraction is done for a single point, terra::extract
        #returns a 2-row df with the value of the cl var in the 2nd col (in case of a single layer extraction)
        #these lines below retrieve the value of the climatic layer, dropping other info
        lyr_val <- extract(x = lyr_touse, y = loc_xy, ID = F, search_radius = rad_meters)
        lyr_val <- data.frame(unique(lyr_val[[nm]]))
        colnames(lyr_val) <- nm
        return(lyr_val)
      } else {
        return(lyr_val)
      }
    }))
    #res_cl is a vector with names being the layer names
    #loop across variable names (Prcp and Tavg), extract their value for the different layers and compute average
    res_val <- lapply(cl_vars, function(var_nm) {
      var_val <- res_cl[grep(pattern = var_nm, x = names(res_cl))]
      #if the whole vector of values for the cl var is NA return NA, otherwise the mean
      if(all(is.na(var_val))) return(NA) else return(mean(var_val, na.rm = T))
    })
    #format output for i-esim id as a data frame
    res_out <- as.data.frame(res_val)
    colnames(res_out) <- cl_vars
    #add PlotID col, temporal window and coordinates
    res_out <- data.frame(PlotID = id, res_out, Cl_time_start = min(loc_win), Cl_time_end = max(loc_win), loc_xy)
  }))
  #call the garbage collector
  gc()
  #return output
  return(res)
  }

#debugonce(fill_climate_gap)

#--------------------------get climatic data version 2 (cell based)

#the idea is that if many points fall within the same cell, and values of Sampl_year are also the same among these plots,
#then it's likely faster to extract climate data by cell than by points

#get cell IDs for the ecoregion
Cant_mf_cellID <- terra::cellFromXY(object = Cant_mf_clim[[1]], xy = as.matrix(Cant_mf[c('X_laea', 'Y_laea')]))

#check if there are NAs in cellID
anyNA(Cant_mf_cellID) #F
which(is.na(Cant_mf_cellID)) #integer(0)

#add Climate_cellID column to EVA_data
Cant_mf$Climate_cellID <- Cant_mf_cellID

#get unique combination of cell and period (time window)
#Cant_mf_cellxprd <- unique(data.frame(CellID = Cant_mf$Climate_cellID, TWinFrom = Cant_mf$Cl_time_start, TWinTo = Cant_mf$Sampl_year))

#function for extracting climatic data from cell ID(s)

extr_climate_eva_cells <- function(x, cl_stack, cellID_col, from_col, to_col) {
  require(terra)
  if(!all(sapply(list(cellID_col, from_col, to_col), is.character))) stop('One or more column names are not chr')
  if(!all((c(cellID_col, from_col, to_col) %in% colnames(x)))) stop('Column names not found in x')
  #get unique combination of cell ids and period
  cell_and_tw <- data.frame(CellID = x[[cellID_col]], TWfrom = x[[from_col]], TWto = x[[to_col]])
  #get unique combinations
  cell_and_tw <- unique(cell_and_tw)
  #get cl vars name
  cl_vars <- unique(substr(x = names(cl_stack), start = 1, stop = 4))
  #loop across combinations and extract data
  res <- do.call(rbind, lapply(seq_len(nrow(cell_and_tw)), function(i) {
    #subset cmb data
    cmb_data <- cell_and_tw[i, ]
    #get temporal window
    tmp_win <- seq.int(from = as.integer(cmb_data[['TWfrom']]), to = as.integer(cmb_data[['TWto']]), by = 1)
    tmp_win <- paste(rep(cl_vars, each = length(tmp_win)), tmp_win, sep = '_') #tmp_win gets recycled here
    #subset stack
    tmp_stck <- cl_stack[[tmp_win]]
    #extract climatic values
    cl_vals <- terra::extract(x = tmp_stck, y = cmb_data[['CellID']])
    #format result as a data.frame
    cl_vals <- data.frame(lapply(cl_vars, function(nm) {
      avg_val <- mean(unlist(cl_vals[grep(pattern = nm, x = colnames(cl_vals))]))
      return(avg_val)
    }))
    colnames(cl_vals) <- cl_vars
    res_out <- data.frame(CellID = cmb_data[['CellID']], cl_vals, Cl_time_start = cmb_data[['TWfrom']], Cl_time_end = cmb_data[['TWto']])
    return(res_out)
  }))
  return(res)
  }

#test_start and end overwrite those used above
test_start <- proc.time()

Cant_mf_clmd <- extr_climate_eva_cells(x = Cant_mf, cl_stack = Cant_mf_clim,
                                              cellID_col = 'Climate_cellID', from_col = 'Cl_time_start', to_col = 'Sampl_year')

test_end <- proc.time() - test_start #approx 5 minutes (!!)

#check class of time start and time end cols
#cols have the same class because in extr_climate_eva_cells they are extracted from the ecoregion dataset
#and left with their class (as.integer coerces them to integer in seq.int() to derive the temporal window)
class(Cant_mf$Cl_time_start); class(Cant_mf$Sampl_year) #'numeric', 'character'
class(Cant_mf_clmd$Cl_time_start); class(Cant_mf_clmd$Cl_time_end) #'numeric', 'character'

#left_join bmf_eva with climate data
Cant_mf <- dplyr::left_join(x = Cant_mf, y = Cant_mf_clmd,
                            by = c('Climate_cellID' = 'CellID',
                                   'Cl_time_start' = 'Cl_time_start',
                                   'Sampl_year' = 'Cl_time_end'))

#get vector of plot ids for which climate data are missing
identical(which(is.na(Cant_mf$Prcp)), which(is.na(Cant_mf$Tavg))) #T

#get data for points for which climate data are missing
Cant_mf_cl_missing <- Cant_mf[is.na(Cant_mf$Tavg), c('PlotID', 'Prcp', 'Tavg', 'Cl_time_start', 'Sampl_year', 'X_laea', 'Y_laea')]

#test fill_climate_gap on a single point with NA for climatic data
ex_plot_cl_missing <- Cant_mf_cl_missing[1, ]

#test fill_climate_gap
fill_climate_gap(x = ex_plot_cl_missing, cl_stack = Cant_mf_clim, loc_id_col = 'PlotID',
                 from_col = 'Cl_time_start', to_col = 'Sampl_year', rad_meters = ((res(Cant_mf_clim)[1])+1))

#all NAs
extract(Cant_mf_clim[[paste('Prcp', seq.int(1981, 1985, 1), sep = '_')]], y = ex_plot_cl_missing[c('X_laea', 'Y_laea')], ID = F)

#check whether one obtains the same results doing the same procedure outside the function

#Prcp
mean(sapply(paste('Prcp', seq.int(1981, 1985, 1), sep = '_'), function(nm) {
  vals <- extract(Cant_mf_clim[[nm]], y = ex_plot_cl_missing[c('X_laea', 'Y_laea')], search_radius = ((res(Cant_mf_clim)[1])+1))
  vals <- unique(vals[[2]])
  return(vals)
}))

#Tavg
mean(sapply(paste('Tavg', seq.int(1981, 1985, 1), sep = '_'), function(nm) {
  vals <- extract(Cant_mf_clim[[nm]], y = ex_plot_cl_missing[c('X_laea', 'Y_laea')], search_radius = ((res(Cant_mf_clim)[1])+1))
  vals <- unique(vals[[2]])
  return(vals)
}))

#fill gaps in climatic dataset using fill_climate_gap

Cant_mf_cl_missing <- fill_climate_gap(x = Cant_mf_cl_missing, cl_stack = Cant_mf_clim, loc_id_col = 'PlotID',
                                       from_col = 'Cl_time_start', to_col = 'Sampl_year', rad_meters = ((res(Cant_mf_clim)[1])+1))

#fill in missing data
sum(duplicated(Cant_mf$PlotID)); sum(duplicated(Cant_mf_cl_missing$PlotID))
class(Cant_mf$PlotID); class(Cant_mf_cl_missing$PlotID)

for(i in Cant_mf_cl_missing$PlotID) {
  Cant_mf[Cant_mf$PlotID == i, 'Prcp'] <- Cant_mf_cl_missing[Cant_mf_cl_missing$PlotID == i, 'Prcp']
  Cant_mf[Cant_mf$PlotID == i, 'Tavg'] <- Cant_mf_cl_missing[Cant_mf_cl_missing$PlotID == i, 'Tavg']
}

rm(i)
rm(Cant_mf_cl_missing) #rm dataset of points with NAs
rm(Cant_mf_clmd) #rm dataset with climate data
rm(test_start, test_end) #rm test_start and test_end

#final check on remaining NAs in climatic fields

sum(is.na(Cant_mf$Prcp)) #1
sum(is.na(Cant_mf$Tavg)) #1

#same position
identical(which(is.na(Cant_mf$Prcp)), which(is.na(Cant_mf$Tavg))) #TRUE

#have a look at the points for which climate data are still missing
plot(Cant_mf_clim[[1]])
plot(st_geometry(EVA_meta.sp[EVA_meta.sp$PlotID %in% Cant_mf[which(is.na(Cant_mf$Prcp)), ], ]), add = T, col = 'red')


#----------------------------------extract topography using geomorpho90m

#get tile IDs for the ecoregion

#import geomorpho90m tile scheme -> can't open it using sf::st_read
#geomorpho_tile_scheme <- st_read(dsn = "/MOTIVATE/GDM_ForesteCasentinesi/TopographicData/geomorph90_Boundary.kmz")

#import proj stack of topography
topo_stack_proj <- rast('TopographyData/topo_stack_proj3035.tif')

#crop to the ecoregion extent
Cant_topo <- crop(x = topo_stack_proj, y = Cant_mf_ext)

#attach CellID column to ecoregion dataset
Cant_mf$Topo_cellID <- cellFromXY(object = Cant_topo[[1]], xy = as.matrix(Cant_mf[c('X_laea', 'Y_laea')]))

#check NAs in cellIDs
sum(is.na(Cant_mf$Topo_cellID)) #0

#number of unique CellIDs
length(unique(Cant_mf$Topo_cellID)) #2932 (out of 6219 locations)

#extract topographic data from XY
Cant_mf_topo <- extract(x = Cant_topo, y = unique(Cant_mf$Topo_cellID))

#coerce to data.frame
Cant_mf_topo <- data.frame(CellID = unique(Cant_mf$Topo_cellID), Cant_mf_topo)

#--check on NAs

#same NAs for all variables
identical(which(is.na(Cant_mf_topo$roughness_mosaic)), which(is.na(Cant_mf_topo$slope_mosaic))) #T
identical(which(is.na(Cant_mf_topo$elevation_mosaic_aligned)), which(is.na(Cant_mf_topo$slope_mosaic))) #FALSE

setdiff(which(is.na(Cant_mf_topo$elevation_mosaic_aligned)), which(is.na(Cant_mf_topo$slope_mosaic))) #15 1076 1171
setdiff(which(is.na(Cant_mf_topo$slope_mosaic)), which(is.na(Cant_mf_topo$elevation_mosaic_aligned))) #13

#number of NAs (per Cell - and not location)
sum(is.na(Cant_mf_topo$roughness_mosaic)) #10
sum(is.na(Cant_mf_topo$elevation_mosaic_aligned)) #12

#check for which cellID there is NA
topo_cells_NA <- list(SlopeRough = as.character(Cant_mf_topo$CellID[which(is.na(Cant_mf_topo$roughness_mosaic))]),
                      Altitude = as.character(Cant_mf_topo$CellID[which(is.na(Cant_mf_topo$elevation_mosaic_aligned))]))

#count number of observations falling in NA cells
topo_cells_cnt <- table(Cant_mf$Topo_cellID)

sum(topo_cells_cnt[topo_cells_NA$SlopeRough]) #17 observations lost due to missingness
sum(topo_cells_cnt[topo_cells_NA$Altitude]) #27

#visualize cellIDs for which we lose the larger amount of plots
plot(Cant_topo[[1]])
points(x = xyFromCell(object = Cant_topo[[1]],
                      cell = as.integer(unique(unlist(topo_cells_NA)))),
       col = 'red')

#--end of check on NAs

#join topo values to ecorigion dataset

#rename columns in Cant_mf_topo
c('Elevation', 'Roughness', 'Slope') %in% colnames(Cant_mf)
colnames(Cant_mf_topo)[c(2, 3, 4)] <- c('Elevation', 'Roughness', 'Slope')

Cant_mf <- dplyr::left_join(x = Cant_mf, y = Cant_mf_topo,
                            by = c('Topo_cellID' = 'CellID'))

#get points for which topography is missing
Cant_mf_topo_missing <- list(SlopeRough = which(is.na(Cant_mf$Roughness)),
                             Altitude = which(is.na(Cant_mf$Elevation)))

#check
length(which(is.na(Cant_mf$Roughness)))
all(which(is.na(Cant_mf$Roughness)) %in% unlist(Cant_mf_topo_missing))

length(which(is.na(Cant_mf$Elevation)))
all(which(is.na(Cant_mf$Elevation)) %in% unlist(Cant_mf_topo_missing))

#create list with coords of cells for which topo data are missing
Cant_mf_topo_missing <- lapply(Cant_mf_topo_missing, function(i, topo_stack) {
  
  cell_id <- unique(Cant_mf[i, c('Topo_cellID')])
  #get coordinates of cell with missing data
  cell_coords <- terra::xyFromCell(object = topo_stack, cell = cell_id)
  #attach cell_id to coords
  res <- data.frame(CellID = cell_id, cell_coords)
  
  return(res)
  
  },
  topo_stack = Cant_topo)

#function for extracting topography values in case of missingness
#using cellIDs instead of single points extraction: https://github.com/rspatial/terra/issues/1733
#search_radius can't be used if y is a cellID, therefore I extract the coordinates of the cells for which I miss data
#and I search missing values around these coordinates

#check result of using extract + search_radius for multiple points
#it looks like it is possible to query the value for multiple points in the same call
extract(x = Cant_topo[[c(2)]], Cant_mf_topo_missing$SlopeRough[c('x', 'y')], ID = F, search_radius = (res(Cant_topo)[1] + 1))
extract(x = Cant_topo[[c(3)]], Cant_mf_topo_missing$SlopeRough[c('x', 'y')], ID = F, search_radius = (res(Cant_topo)[1] + 1))
extract(x = Cant_topo[[c(1)]], Cant_mf_topo_missing$Altitude[c('x', 'y')], ID = F, search_radius = (res(Cant_topo)[1] + 1))

fill_topo_gap <- function(x, topo_lyr, rad_meters = (res(topo_lyr)[1] + 1)) {
  #get coords for missing data
  x_coords <- x[c('x', 'y')]
  #extract topo values
  topo_val <- extract(x = topo_lyr, y = x_coords, ID = F, search_radius = rad_meters)
  #keep only relevant columns and attach (original) cell id
  topo_val <- data.frame(topo_val[c(1, 2)], x['CellID'])
  #return result
  return(topo_val)
  }


#check
fill_topo_gap(x = Cant_mf_topo_missing$SlopeRough, topo_lyr = Cant_topo[[2]])

#get missing for slope and roughness
Cant_rghslo_missing <- do.call(cbind, lapply(Cant_topo[[c("roughness_mosaic", "slope_mosaic")]], fill_topo_gap, x = Cant_mf_topo_missing$SlopeRough))
Cant_rghslo_missing <- Cant_rghslo_missing[c(1, 4, 6)]
#remove NAs
Cant_rghslo_missing <- na.omit(Cant_rghslo_missing)

#get missing for elevation
Cant_elev_missing <- fill_topo_gap(x = Cant_mf_topo_missing$Altitude, topo_lyr = Cant_topo[[1]])
Cant_elev_missing <- Cant_elev_missing[c(1, 3)]
#remove NAs
Cant_elev_missing <- na.omit(Cant_elev_missing)

#check for how many locations I can actually fill the gaps
#Roughness-slope
sum(table(Cant_mf[Cant_mf$Topo_cellID %in% Cant_rghslo_missing$CellID, 'Topo_cellID'])) #11 (out of 17)
#Elevation
sum(table(Cant_mf[Cant_mf$Topo_cellID %in% Cant_elev_missing$CellID, 'Topo_cellID'])) #19 (put of 27)

#fill topo gaps

for(i in Cant_rghslo_missing$CellID) {
  Cant_mf[Cant_mf$Topo_cellID == i, 'Roughness'] <- Cant_rghslo_missing[Cant_rghslo_missing$CellID == i, 'roughness_mosaic']
  Cant_mf[Cant_mf$Topo_cellID == i, 'Slope'] <- Cant_rghslo_missing[Cant_rghslo_missing$CellID == i, 'slope_mosaic']
}

rm(i)

for(i in Cant_elev_missing$CellID) {
  Cant_mf[Cant_mf$Topo_cellID == i, 'Elevation'] <- Cant_elev_missing[Cant_elev_missing$CellID == i, 'elevation_mosaic_aligned']
}

rm(i)

#check
sum(is.na(Cant_mf$Roughness)) #6
sum(is.na(Cant_mf$Elevation)) #8

rm(Cant_mf_topo, topo_cells_NA, topo_cells_cnt, Cant_mf_topo_missing, Cant_rghslo_missing, Cant_elev_missing)

#----------------------------------extract Human Modification Index

#human modification index is available every 5 years from 1990 to 2020
#every Sampl_year in EVA is matched with the closest year of the hmi

table(Cant_mf$Sampl_year)

#1980 - 1990: assigned to 1990
#years are otherwise assigned to closest hmi year, e.g. 1991, 1992 are assigned to 1990,
#while 1993, 1994 to 1995

#write a small function that retrieves the year from which to extract hmi based on the Sampl_year field in EVA meta data
#the function is vectorised
get_hmi_year <- function(x) {
  if(is.character(x)) x <- as.integer(x)
  if(min(x) < 1980 || max(x) > 2025) stop('Data are outside plausible temporal window')
  res <- ifelse(x <= 1992, '1990',
                ifelse(x > 1992 & x <= 1997, '1995',
                       ifelse(x > 1997 & x <= 2002, '2000',
                              ifelse(x > 2002 & x <= 2007, '2005',
                                     ifelse(x > 2007 & x <= 2012, '2010',
                                            ifelse(x > 2012 & x <= 2017, '2015', '2020'))))))
  return(res)
}

#check -> looks good
cbind(A = get_hmi_year(x = seq.int(from = 1980, to = 2022, 1)), B = seq.int(from = 1980, to = 2022, 1))


#import hmi stack
hmi_stack_proj <- rast('HumanModificationIndex/hmi_stack_proj3035.tif')

#modify names of the layers
hmi_nm <- names(hmi_stack_proj)

#extract index year - strings all have the same pattern
names(hmi_stack_proj) <- paste0('Hum_mod_ind_', substr(hmi_nm, start = 13, stop = 16))

#crop it to the extent of the ecoregion
Cant_hmi <- crop(x = hmi_stack_proj, y = Cant_mf_ext)

#add the Hmi_cellID column to the ecoregion dataset
Cant_mf$Hmi_cellID <- terra::cellFromXY(object = Cant_hmi, xy = as.matrix(Cant_mf[c('X_laea', 'Y_laea')]))

#check NAs in cellID
sum(is.na(Cant_mf$Hmi_cellID)) #0

#also add the Hmi_yr column to the ecoregion dataset (this column includes the hmi year associated with every sampling year)
Cant_mf$Hmi_yr <- get_hmi_year(Cant_mf$Sampl_year)

#function for extracting value of the hmi at ecoregion locations (cell based)
extr_hmi_eva_cells <- function(x, hmi_stack, cellID_col, hmi_yr_col) {
  require(terra)
  if(!all(sapply(list(cellID_col, hmi_yr_col), is.character))) stop('One or more column names are not chr')
  if(!all((c(cellID_col, hmi_yr_col) %in% colnames(x)))) stop('Column names not found in x')
  #get combination of cell ids and hmi year
  cell_and_yr <- data.frame(CellID = x[[cellID_col]], Hmi_yr = x[[hmi_yr_col]])
  #get unique combinations
  cell_and_yr <- unique(cell_and_yr)
  #loop across cell_ids and hmi years and extract data
  res <- do.call(rbind, lapply(seq_len(nrow(cell_and_yr)), function(i) {
    cmb_data <- cell_and_yr[i, ]
    #retrieve name of the layer to be used
    lyr_nm <- grep(pattern = cmb_data[['Hmi_yr']], x = names(hmi_stack), value = T)
    #extract value of hmi at the cell
    hmi_val <- terra::extract(hmi_stack[[lyr_nm]], y = cmb_data[['CellID']])
    colnames(hmi_val) <- 'Hmi_value'
    return(hmi_val)
  }))
  #format output as data frame
  res <- data.frame(cell_and_yr, res)
  return(res)
}

#check
extr_hmi_eva_cells(x = Cant_mf[1:10, ], hmi_stack = Cant_hmi, cellID_col = 'Hmi_cellID', hmi_yr_col = 'Hmi_yr')

#test_start and test_end overwrite those above
test_start <- proc.time()

Cant_mf_hmi <- extr_hmi_eva_cells(x = Cant_mf, hmi_stack = Cant_hmi, cellID_col = 'Hmi_cellID', hmi_yr_col = 'Hmi_yr')

test_end <- proc.time() - test_start #approx 3 mins (similar as to climate data extraction)

#join values of hmi to ecoregion dataset

#check classes of column to join
class(Cant_mf$Hmi_cellID); class(Cant_mf_hmi$CellID)
class(Cant_mf$Hmi_yr); class(Cant_mf_hmi$Hmi_yr)

#join hmi data to ecoregion dataset
Cant_mf <- dplyr::left_join(x = Cant_mf, y = Cant_mf_hmi,
                            by = c('Hmi_cellID' = 'CellID',
                                   'Hmi_yr' = 'Hmi_yr'))


#check how many locations have hmi missing
length(which(is.na(Cant_mf$Hmi_value))) #104

#have a look at these plots
plot(Cant_hmi[[1]])
plot(EVA_meta.sp[EVA_meta.sp$PlotID %in% (Cant_mf$PlotID[is.na(Cant_mf$Hmi_value)]), ], col = 'red', add = T)

#function for filling gaps in hmi data

#first, get data of locations with missing hmi value
#only keep columns needed for data extraction
Cant_mf_hmi_missing <- Cant_mf[which(is.na(Cant_mf$Hmi_value)), c('PlotID', 'Hmi_yr', 'X_laea', 'Y_laea')]


fill_hmi_gap <- function(x, hmi_stack, loc_id_col, rad_meters = (res(hmi_stack)[1] + 1)) {
  require(terra)
  #get plot ids
  loc_id <- x[[loc_id_col]]
  #loop across plot ids
  res <- do.call(rbind, lapply(seq_along(loc_id), function(i) {
    #get plot coords
    loc_xy <- x[i, c('X_laea', 'Y_laea')]
    #use hmi_yr to get layer name
    hmi_yr_touse <- x[i, 'Hmi_yr']
    lyr_nm <- grep(pattern = hmi_yr_touse, x = names(hmi_stack), value = T)
    #extract value using the search_radius option
    hmi_val <- extract(x = hmi_stack[[lyr_nm]], y = loc_xy, ID = F, search_radius = rad_meters)
    hmi_val <- data.frame(Hmi_value = unique(hmi_val[[lyr_nm]]))
    return(hmi_val)
  }))
  res <- data.frame(PlotID = loc_id, Hmi_yr = x[['Hmi_yr']], res)
  return(res)
}

#check
fill_hmi_gap(x = Cant_mf_hmi_missing[1:10, ], hmi_stack = Cant_hmi, loc_id_col = 'PlotID')

#run on all locations
Cant_mf_hmi_missing <- fill_hmi_gap(x = Cant_mf_hmi_missing, hmi_stack = Cant_hmi, loc_id_col = 'PlotID')

#check number of remaining NA
sum(is.na(Cant_mf_hmi_missing$Hmi_value)) #21

#drop NAs
Cant_mf_hmi_missing <- na.omit(Cant_mf_hmi_missing)

#fill in missing data in ecoregion dataset
for(i in Cant_mf_hmi_missing$PlotID) {
  Cant_mf[Cant_mf$PlotID == i, 'Hmi_value'] <- Cant_mf_hmi_missing[Cant_mf_hmi_missing$PlotID == i, 'Hmi_value']
}

rm(i)
rm(Cant_mf_hmi_missing)

#check on remaining NAs
sum(is.na(Cant_mf$Hmi_value)) #21!


#----------------------------------extract Corine Land Cover (optional)

#library(clc)

#check R package for clc: https://cran.r-project.org/web/packages/clc/clc.pdf (vignette: https://cran.r-project.org/web/packages/clc/vignettes/clc.html)

#download CLC map for whole Europe, project it to 3034 and work with that
#it's in principle possible to download CLC maps for individual ecoregions, but 1) auxiliary files are not downloaded or
#anyway it's hard to locate them in the downloaded files; 2) projecting the whole EU layer at the beginning avoids
#working with projections made at the ecoregion extent


#create a shapefile of the ecoregion to download the CLC map from the viewer of Copernicus
#bmf_poly.shp <- cand_ecoreg_frame[cand_ecoreg_frame$ECO_NAME == 'Baltic mixed forests',]
#transform to epsg 3035, which is admitted by Copernicus viewer
#bmf_poly.shp <- st_transform(bmf_poly.shp, crs = 3035)
#write shapefile for CLC extraction - the shp must be compressed in a zip file
#bmf_poly.shp <- st_write(bmf_poly.shp, dsn = "Shapefile_ecoregions/Baltic_mixed_forest.shp")


#import CLC for whole EU: https://stackoverflow.com/questions/71831807/how-to-assign-and-match-land-cover-type-from-corine-to-a-dataframe-with-a-set-of

#bmf_clc_2018 <- rast("/MOTIVATE/GDM_EuropeanEcoregions/CLC_layers/199234/Results/u2018_clc2018_v2020_20u1_raster100m/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif")

#plot(bmf_clc_2018)



