#this code is a pipeline for processing the raster layers for climate, topography and 
#the human modification index. The problem is that these layers are available at different
#extents and resolutions and they all need to be projected to EPSG:3035

#the idea is to crop all of them to the bbox of the climatic data (assuming they are all in 4326 - they actually are) as these are only available for
#Europe. Using the bbox of EVA data could include data points for which we would eventually not have climatic info

#then, I can re-project all layers to the EPSG:3035. In this way, I project data only once (each projection implies data transformation)

#layers will be left at their original resolution (the one reported after projection) and resampled only for mapping predictions of beta diversity

#N.B.: projecting all layers at once allows to have the same resolution for all ecoregions. Ecoregion-specific projections would inevitably
#lead to layers with different resolutions.

library(sf)
library(terra)
library(mapview)


#------------first import one climatic layer to get the extent

#ESPG:4326 - using Modzilla version of climate data as newest version (rclone) has many more NAs :()
prcp_stack <- list.files(path = 'EasyClimateData/v4data_from_ftpzilla/prec/', pattern = 'YearlySum_cogeo.tif', full.names = T) 

prcp_stack <- rast(lapply(prcp_stack, rast))

#check
nlyr(prcp_stack) #53
nlyr(prcp_stack) == length(seq(from = 1970, to = 2022, by =  1)) #T

tavg_stack <- list.files(path = 'EasyClimateData/v4data_from_ftpzilla/tavg/', pattern = 'YearlyAvg_cogeo.tif', full.names = T) 

tavg_stack <- rast(lapply(tavg_stack, rast))

#check
nlyr(tavg_stack) #53

#compare prcp and tavg
compareGeom(prcp_stack, tavg_stack) #T

res(prcp_stack) #0.008333333 0.008333333

clim_ext <- ext(prcp_stack)

#project the climatic stacks to EPSG:3035

#first combine the prcp and tavg stacks
clim_stack <- c(prcp_stack, tavg_stack)

#divide the whole stack by 100 - this scaling parameter is used for storage reason
#better to do it now to avoid running the projection on unscaled values
clim_stack <- clim_stack/100

#first project a layer to be used as a template for the whole stack
tmp_cl <- terra::subset(clim_stack, 1)

tmp_cl <- terra::project(tmp_cl, y = 'epsg:3035', method = 'bilinear')

#then project the whole stack using the newly created template
clim_stack_proj <- project(x = clim_stack, y = tmp_cl, method = 'bilinear', filename = 'EasyClimateData/clim_stack/clim_stack_proj3035.tif')

#get min and max of prcp layers to check I have no negative sums of precipitation
sapply(as.list(subset(clim_stack_proj, 1:53)), minmax) #ok

rm(tmp_cl)

#------------import human modification index layers

#EPSG:4326
hmi_stack <- list.files(path = 'HumanModificationIndex/', pattern = '_AA.tif', full.names = T)

hmi_stack <- rast(lapply(hmi_stack, rast))

#check
nlyr(hmi_stack) #T
nlyr(hmi_stack) == length(seq(1990, 2020, by = 5)) #T

res(hmi_stack) #0.002694946 0.002694946

#crop to clim_ext
hmi_stack <- crop(x = hmi_stack, y = clim_ext, filename = 'HumanModificationIndex/geo_crop_HMI_AA.tif')

#projet to epsg:3035

#!!The year-specific hmi layers we deleted from disk (locally) as they were 9 GB each
#to (re)project the longlat hmi data to epsg:3035 (instead of previously used 3034)
#I'm directly loading the geo_crop_HMI_AA stack created above

hmi_stack <- rast('HumanModificationIndex/geo_crop_HMI_AA.tif')

#first of all, create a tmp layer to be used as a template for projection
tmp_hmi <- terra::subset(hmi_stack, 1)

#project tmp_hmi
tmp_hmi <- terra::project(x = tmp_hmi, y = 'epsg:3035', method = 'bilinear')

#check min and max value of tmp_hmi
terra::minmax(tmp_hmi) #ok

#then project the whole stack using tmp_hmi
hmi_stack_proj <- terra::project(x = hmi_stack, y = tmp_hmi, method = 'bilinear', filename = "HumanModificationIndex/hmi_stack_proj3035.tif")

#check min and max
sapply(as.list(hmi_stack_proj), minmax) #ok

rm(tmp_hmi)

#------------import topographic layers

#these layers were derived by creating a mosaic of the macro tiles of geomorpho90m
#the mosaic was built by Gabriel (CZU) - for more info, see README in the data folder
#the mosaic was already cropped by Gabriel, who used the same extent of clim_ext
#on Dec 24, I re-created (and re-wrote to disk) the topographic stack to which I added the altitude layer
#the topo stack already included slope and roughness

#topo stack
topo_stack <- list.files(path = 'TopographyData/', pattern = 'tif', full.names = T)

#EPSG:4326
topo_stack <- rast(lapply(topo_stack, rast))

nlyr(topo_stack) #3

res(topo_stack) #0.0008333333 0.0008333333

#project to epsg:3035

#even though I'm currently using only 3 layers (altitude, slope and roughness), I replicate
#the procedure followed above. First, I create a template to be used for projecting the
#whole stack. I'm doing this as it is likely that I'll have to add other topo variables

tmp_topo <- terra::subset(topo_stack, 2)

#project tmp_topo (roughness)
tmp_topo <- terra::project(x = tmp_topo, y = 'epsg:3035', method = 'bilinear')

#check min and max value of tmp_topo
terra::minmax(tmp_topo) #ok

#project the whole stack using tmp_topo
topo_stack_proj <- terra::project(x = topo_stack, y = tmp_topo, method = 'bilinear', filename = 'TopographyData/topo_stack_proj3035.tif')

rm(tmp_topo)


