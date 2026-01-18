#This code processes the shapefiles of ecoregions and European land borders, and EVA data to get a preliminary list of candidate ecoregions
#to be analyzed for MOTIVATE

#spatial analysis
library(sf)
library(terra)
library(mapview)
library(data.table)

#viz
library(ggplot2)
library(ggpubr)

#parallelization
library(future.apply)

#notes

#IDEA FOR ANALYSIS: an idea could be to follow a latitudinal gradient of forests or grasslands (or both)
#the ecoregions remaining after dropping plots with missing covariates or inadequate temporal resolution
#span a quite wide longitudinal and latitudinal (if including something from Baltic) gradient
#analyses focusing on the two decades and across latitude could evidence change in the effect of drivers
#on taxonomic dissimilarity over time (within ecoregions)
#and, at the same time, change in the effect of drivers across latitude.

#CRS used for spatial analysis:
#ETRS89-extended / LAEA Europe (https://spatialreference.org/ref/epsg/3035/)
epsg_proj <- 3035

#---------------------download and import the map of EU ecoregions

#for more details about the map of ecoregions, check README.txt in the map folder
eu_ecoregions <- st_read('MapOfEcoregions/Ecoregions2017/Ecoregions2017.shp')

#summary(eu_ecoregions)
#mapview(eu_ecoregions)

#we're mostly interested in 'Palearctic'
mapview(eu_ecoregions[eu_ecoregions$REALM == "Palearctic", ])

#some geometries are not valid - this can create problems when cropping
sum(!st_is_valid(eu_ecoregions)) #69

#check invalid geometries
inv_polys <- eu_ecoregions[!st_is_valid(eu_ecoregions), ]

mapview(inv_polys) #no problematic polygon overlaps with EU

#try to fix invalid polygons using st_make_valid() - otherwise I have problems when cropping (see below)
eu_ecoregions <- st_make_valid(eu_ecoregions)

sum(!st_is_valid(eu_ecoregions)) #0

#check which Palearctic region overlaps EU
grep(pattern = "PA", x = unique(eu_ecoregions$ECO_BIOME_), value = T)

mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA12", ]) #yes
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA13", ]) #nope
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA08", ]) #nope
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA05", ]) #part of it
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA10", ]) #nope
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA09", ]) #nope
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA04", ]) #yes
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA11", ]) #part of it
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA06", ]) #part of it
mapview(eu_ecoregions[eu_ecoregions$ECO_BIOME_ == "PA01", ]) #nope

#transform eu_ecoregions to EPSG:3035 (planar coordinates)

eu_ecoregions.proj <- st_transform(x = eu_ecoregions, crs = epsg_proj)

mapview(eu_ecoregions.proj)

#remove eu_ecoregions to save memory
rm(eu_ecoregions)

#---------------------download and import the GSHHG map of land borders

#for more details about the GSHHG map, check README.txt in the map folder

#Using shapefile (polygons)
#Full-resolution (f): Full resolution.  These contain the maximum resolution of this data and has not been decimated.
#Level-1: (Continental land masses and ocean islands, except Antarctica)

eu_land <- st_read(dsn = '/MOTIVATE/GDM_EuropeanEcoregions/GSHHG/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp')

summary(eu_land)
mapview(eu_land)

#check whether geometries are valid
sum(!st_is_valid(eu_land)) #1

which(!st_is_valid(eu_land)) #2245

#check invalid geometry
mapview(eu_land[2245, ]) #somewhere in Maine

#project to EPSG:3035
eu_land.proj <- st_transform(x = eu_land, crs = epsg_proj)

rm(eu_land)

#---------------------import EVA metadata

EVA_meta <- fread(file = 'C://MOTIVATE/EVA_dataset/200_Motivate20250415_notJUICE_header_with_precise_coordinates.csv', data.table = FALSE)

#remove columns that won't be used

EVA_meta <- EVA_meta[c(1:10, 14:18, 21, 23:42)]

#use only PlotID as the plot identifier column -> I'm keeping PlotID as integer as it is much faster to use integer in logical comparisons (than, eg, character)
identical(EVA_meta$PlotObservationID, EVA_meta$PlotID) #F
sapply(EVA_meta[c(1, 2)], class)
identical(as.integer(EVA_meta$PlotObservationID), EVA_meta$PlotID) #T

#coerce PlotObservationID to character, just in case I need a chr version of PlotID
EVA_meta$PlotObservationID <- as.character(EVA_meta$PlotObservationID)

#explo on column values

#check NAs
sapply(EVA_meta, function(i) sum(is.na(i)))

#in many plots, vegetation was recorded on a 1/0 scale
table(EVA_meta$`Cover abundance scale`)

#get rid of observations with missing coordinates
length(sort(unique(which(is.na(EVA_meta[c('Longitude', 'Latitude')]), arr.ind = T)[, 1]))) #182558

#same number as
nrow(EVA_meta) - sum(complete.cases(EVA_meta[c('Longitude', 'Latitude')])) #182558

#check if missing long and lat occur in the same rows
identical(which(is.na(EVA_meta$Longitude)), which(is.na(EVA_meta$Latitude))) #TRUE

EVA_meta <- EVA_meta[!is.na(EVA_meta$Longitude), ]

#check
anyNA(EVA_meta[c('Longitude', 'Latitude')]) #F

#--modify column names not to have white spaces or 'weird symbols'

#remove period '.' - using fixed = T to match period 'as is'
colnames(EVA_meta) <- gsub(pattern = '.', replacement = '', x = colnames(EVA_meta), fixed = T)

#substitute 'é' with e
colnames(EVA_meta) <-  gsub(pattern = 'é', replacement = 'e', x = colnames(EVA_meta), fixed = T)

#replace (m^2) by (m2)
colnames(EVA_meta)[10] <- 'Releve area m2'

#substitute white spaces with underscores
colnames(EVA_meta) <- gsub(pattern = ' ', replacement = '_', x = colnames(EVA_meta))

#remove parentheses
colnames(EVA_meta) <- gsub(pattern = '(', replacement = '', x = colnames(EVA_meta), fixed = T)
colnames(EVA_meta) <- gsub(pattern = ')', replacement = '', x = colnames(EVA_meta), fixed = T)

#replace slash with underscore
colnames(EVA_meta) <- gsub(pattern = '/', replacement = '_', x = colnames(EVA_meta), fixed = T)

#--modify habitat labels

#check frequency of multiple habitat assignment
as.data.frame(table(grep(pattern = ',', x = EVA_meta$Expert_system, value = T))) #4030 cases

#remove exclamation marks
#identical(grep(pattern = '!', x = EVA_meta$Expert_system), grep(pattern = '!', x = EVA_meta$Expert_system, fixed = T)) #T
EVA_meta$Expert_system <- gsub(pattern = '!', replacement = '', x = EVA_meta$Expert_system)

#empty assignments become 'X'
length(which(EVA_meta$Expert_system == '')) #749821
EVA_meta$Expert_system[which(EVA_meta$Expert_system == '')] <- 'X'
length(which(EVA_meta$Expert_system == 'X'))

#Unclassified obs become 'Y'
length(which(EVA_meta$Expert_system == '~')) #75287
unique(grep('~', EVA_meta$Expert_system, value = T))

EVA_meta$Expert_system[which(EVA_meta$Expert_system == '~')] <- 'Y'
length(which(EVA_meta$Expert_system == 'Y'))

#retrieve all cases of multiple habitat assignment to one observation
multy_ass <- grep(pattern = ',', x = EVA_meta$Expert_system, value = T)

#remove duplicated strings
multy_ass <- unique(multy_ass)

#function to compare first letters and check if they are identical or not
check_multy_hab <- function(x) {
  #split the string by comma
  str_parts <- strsplit(x = x, split = ',')
  #retrieve first letter(s)
  first_l <- sapply(str_parts, substr, start = 1, stop = 1)
  #check if all letters are the same
  if(length(unique(first_l)) == 1) return(TRUE) else return(FALSE)
  }

#match cases and corresponding logical (TRUE if the assignment is to be kept or FALSE otherwise)
multy_ass <- setNames(object = sapply(multy_ass, check_multy_hab), nm = multy_ass)

#only keep multiple habitat assignment to exclude
multy_ass <- multy_ass[!multy_ass]

#only keep the 'name' of the cases
multy_ass <- names(multy_ass)

#check how many plots would be lost due to excluding these cases
sum(EVA_meta$Expert_system %in% multy_ass) #501 -> these can be safely removed

#remove problematic multiple assignments
anyNA(EVA_meta$Expert_system) #F

EVA_meta <- EVA_meta[!EVA_meta$Expert_system %in% multy_ass, ]

#create column with first letter of habitat string for rough classification
EVA_meta$Eunis_lev1 <- substr(x = EVA_meta$Expert_system, start = 1, stop = 1)

#create column with description of rough EUNIS classification
hab_key <- setNames(c('empty', 'surface_waters', 'unclassified', 'wetlands_mires', 'coast_saltmarshes', 'inland_sparse_veg', 'grasslands', 'heathlands_scrub_tundra',
                      'man_made_hab', 'forests', 'coast_sand_cliff'),
                    nm = c('X', 'P', 'Y', 'Q', 'M', 'U', 'R', 'S', 'V', 'T', 'N'))

#create new column
all(unique(EVA_meta$Eunis_lev1) %in% names(hab_key)) #T
all(names(hab_key) %in% unique(EVA_meta$Eunis_lev1)) #T

EVA_meta$EunisVerbose_lev1 <- unname(hab_key[EVA_meta$Eunis_lev1])

#--exclude observations missing info on sampling year

class(EVA_meta$Date_of_recording) #chr

anyNA(EVA_meta$Date_of_recording) #F

unique(substr(x = EVA_meta$Date_of_recording, start = 7, stop = 10)) #contains ''

sum(EVA_meta$Date_of_recording == '') #180108

unique(nchar(EVA_meta$Date_of_recording)) #10 0 7

table(nchar(EVA_meta$Date_of_recording)) #96 cases of dates with 7 characters

unique(EVA_meta[which(nchar(EVA_meta$Date_of_recording) == 7), 'Date_of_recording']) #these all have the following format: 0:00:00

#assign NA to cases that have to be removed

EVA_meta$Date_of_recording[EVA_meta$Date_of_recording == ''] <- NA_character_

EVA_meta$Date_of_recording[which(EVA_meta$Date_of_recording == '0:00:00')] <- NA_character_

sum(is.na(EVA_meta$Date_of_recording)) #180204 (180108 + 96) [OK]

#remove NAs
EVA_meta <- EVA_meta[!is.na(EVA_meta$Date_of_recording), ]

#coerce date to a Date format and extract year
EVA_meta$Sampl_year <- EVA_meta$Date_of_recording

#replace '.' by '-'
EVA_meta$Sampl_year <- gsub(pattern = '.', replacement = '-', x = EVA_meta$Sampl_year, fixed = T)

#coerce to Date format
EVA_meta$Sampl_year <- as.Date(x = EVA_meta$Sampl_year, format = "%d-%m-%Y")

#extract year - note that format returns a character
EVA_meta$Sampl_year <- format(EVA_meta$Sampl_year, format = '%Y')

sort(unique(EVA_meta$Sampl_year)) #1890 -> 2024

#----exclude manipulated plots (mostly from ReSurvey) and duplicates

#here, I'm following the procedure suggested by Ilona - see email exchange 27-28/11/2025

#1) First exclude duplicates (and keep only 1 version of Danish dataset)
#2) Identify (truly) resampled plots (ReSurvey)
#3) Clean up wrongly assigned values in Manipulate_y_n field
#4) Exclude manipulated plots

#-- 1 - exclude duplicates of whole datasets

#use the RS_DUPL col because this field is filled, i.e. it is not empty, if the dataset is considered a duplicate
unique(EVA_meta$RS_DUPL) #empty string or a dataset name/code
anyNA(EVA_meta$RS_DUPL) #FALSE
table(EVA_meta$RS_DUPL)
length(which(EVA_meta$RS_DUPL != '')) #5081 plots

#exclude duplicates, which means I keep only empty string entries
EVA_meta <- EVA_meta[EVA_meta$RS_DUPL == '', ]

#keep only 1 version of the Danish dataset - I'm keeping the Danish dataset that is part of core EVA

#!!Note that I'm keeping both pres/abs and abundance data in core EVA for the moment

#check number of observations in Denmark Naturdata and Denmark Naturdata (abund.)
length(grep('^Denmark Naturdata$', x = EVA_meta$Dataset, value = T)) #278726
length(grep('Denmark Naturdata (abundance)', x = EVA_meta$Dataset, value = T, fixed = T)) #113091

table(EVA_meta[grep('^Denmark Naturdata$', x = EVA_meta$Dataset), 'Cover_abundance_scale']) #only presence/absence

#exclude DK_Naturdata_Res and DK_Naturdata_ResAbund
length(grep('^DK_Naturdata_Res$', x = EVA_meta$Dataset, value = T)) #158800
length(grep('DK_Naturdata_ResAbund', x = EVA_meta$Dataset, value = T)) #90171

#nrow(EVA_meta) - nrow(EVA_meta[EVA_meta$Dataset != 'DK_Naturdata_Res', ]) #158800
#nrow(EVA_meta) - nrow(EVA_meta[EVA_meta$Dataset != 'DK_Naturdata_ResAbund', ]) #90171

#remove presence/absence data
EVA_meta <- EVA_meta[EVA_meta$Dataset != 'DK_Naturdata_Res', ]
#remove abundance data
EVA_meta <- EVA_meta[EVA_meta$Dataset != 'DK_Naturdata_ResAbund', ]

#- remove all datasets associated with presence/absence data or pinpoint method

#check unique values of Cover_abundance_scale
sum(is.na(EVA_meta$Cover_abundance_scale)) #0
unique(EVA_meta$Cover_abundance_scale) #104

#drop 'Presence/Absence', 'Presentie/Absentie', 'frequency %', if still present

#Presence/Absence
sum(EVA_meta$Cover_abundance_scale == "Presence/Absence") #319197
unique(EVA_meta$Dataset[EVA_meta$Cover_abundance_scale == "Presence/Absence"]) #43 datasets, including Denmark Naturdata
table(EVA_meta$Dataset[EVA_meta$Cover_abundance_scale == "Presence/Absence"]) #majority of plots are from Denmark Naturdata (278726)

#Presentie/Absentie
sum(EVA_meta$Cover_abundance_scale == 'Presentie/Absentie') #0 (no longer there)

#frequency %
sum(EVA_meta$Cover_abundance_scale == 'frequency %') #0 (no longer there)

#check
c('Presentie/Absentie', 'frequency %') %in% unique(EVA_meta$Cover_abundance_scale) #FALSE FALSE

#drop plots associated with Presence/Absence category
EVA_meta <- EVA_meta[EVA_meta$Cover_abundance_scale != "Presence/Absence", ]

#check number of plots at this stage
nrow(EVA_meta) #1472627

#-- 2 - identify truly resampled plots (ReSurvey)

#to do so, combine the RS_CODE and ReSurvey_plot_Y_N fields
#truly resampled plots are those with a NON-EMPTY value in RS_code and 'Y' in ReSurvey_plot_Y_N
#because Ilona marks as 'N' plots that are in principle part of ReSurvey, although they have not been sampled yet

#check values of the RS_CODE field
unique(EVA_meta$RS_CODE) #empty or a string with dataset code
length(which(EVA_meta$RS_CODE == '')) #1305098

#assign 'NOT_RS' to all observations with an empty string
#grep(pattern = '^NOT_RS$', x = EVA_meta$RS_CODE, value = T) #character(0)
#grep(pattern = 'NOT_RS', x = EVA_meta$RS_CODE, value = T) #character(0)

EVA_meta$RS_CODE[EVA_meta$RS_CODE == ''] <- 'NOT_RS'

#combine RS_CODE with ReSurvey_plot_Y_N to create the 'Truly_RS' field
#this field will be T for truly RS plots and FALSE otherwise

#first check values of ReSurvey_plot_Y_N
unique(EVA_meta$ReSurvey_plot_Y_N) #'Y' 'N' ''
length(EVA_meta$ReSurvey_plot_Y_N[EVA_meta$ReSurvey_plot_Y_N == '']) #1049969
length(EVA_meta$ReSurvey_plot_Y_N[EVA_meta$ReSurvey_plot_Y_N == 'N']) #266443

#empty values in ReSurvey_plot_Y_N means 'N' (following Ilona here), as these observations are core EVA data
EVA_meta$ReSurvey_plot_Y_N[EVA_meta$ReSurvey_plot_Y_N == ''] <- 'N'

#check
#length(EVA_meta$ReSurvey_plot_Y_N[EVA_meta$ReSurvey_plot_Y_N == 'N']) #1316412

table(with(EVA_meta, ifelse(RS_CODE != '' & ReSurvey_plot_Y_N == 'Y', TRUE, FALSE)))

#create the Truly_RS field
EVA_meta$Truly_RS <- with(EVA_meta, ifelse(RS_CODE != '' & ReSurvey_plot_Y_N == 'Y', TRUE, FALSE))

#compare with ReSurvey_plot_Y_N
identical(which(EVA_meta$ReSurvey_plot_Y_N == 'Y'), which(EVA_meta$Truly_RS)) #TRUE

#-- 3 - Clean up wrongly assigned values in Manipulate_y_n field

#fix typo in Manipulate_y_n; 2 plots (with RS_PROJTYP == 'permanent') have Manipulate_y_n == 'M'
EVA_meta$Manipulate_y_n[EVA_meta$Manipulate_y_n == 'M'] <- 'N'

#DE_00039 -> DE_0039 (typo in Ilona's email)
#EVA_meta$RS_CODE[EVA_meta$RS_CODE == 'DE_0039'] #204
#EVA_meta$Dataset[EVA_meta$Dataset == 'DE_0039'] #(character(0))

#DE_0039; empty field manipulate -> 'N'
unique(EVA_meta[EVA_meta$RS_CODE == 'DE_0039', 'RS_PROJTYP']) #permanent [RIGHT]
EVA_meta[EVA_meta$Manipulate_y_n == '' & EVA_meta$RS_CODE == 'DE_0039', 'Manipulate_y_n'] <- 'N'

#Czechia_nvd (RS_CODE: 'CZ_0019_054'): Manipulate_y_n == 'M' already changed to 'N' (see above)
unique(EVA_meta[EVA_meta$Dataset == 'Czechia_nvd' & EVA_meta$RS_CODE == 'CZ_0019_054', 'RS_PROJTYP']) #permanent [RIGHT] 
unique(EVA_meta[EVA_meta$Dataset == 'Czechia_nvd' & EVA_meta$RS_CODE == 'CZ_0019_054', 'Manipulate_y_n']) #'N'

#DE_0037_29: empty field manipulate -> 'N' (keep 'Y' as they were)
unique(EVA_meta[EVA_meta$RS_CODE == 'DE_0037_29', 'RS_PROJTYP']) #permanent
unique(EVA_meta[EVA_meta$RS_CODE == 'DE_0037_29', 'Manipulate_y_n']) #this includes only 'Y' - as also pointed out by Michael

#AT_0004b: empty field manipulate -> 'N'
unique(EVA_meta[EVA_meta$RS_CODE == 'AT_0004b', 'RS_PROJTYP']) #resampling [RIGHT]
EVA_meta[EVA_meta$Manipulate_y_n == '' & EVA_meta$RS_CODE == 'AT_0004b', 'Manipulate_y_n'] <- 'N'

#IT_0001b and IT_0001e: empty field manipulate -> 'N'
unique(EVA_meta[EVA_meta$RS_CODE %in% c('IT_0001b', 'IT_0001e'), 'RS_PROJTYP']) #resampling [RIGHT]
EVA_meta[EVA_meta$RS_CODE %in% c('IT_0001b', 'IT_0001e') & EVA_meta$Manipulate_y_n == '', 'Manipulate_y_n'] <- 'N'

#RU_0003c: dataset not present in current selection
unique(EVA_meta[EVA_meta$RS_CODE == 'RU_0003c', 'RS_PROJTYP']) #character(0)
sum(EVA_meta$RS_CODE == 'RU_0003c') #0

#FRCH_0001: manipulate 'Y' -> 'N'
unique(EVA_meta[EVA_meta$RS_CODE == 'FRCH_0001', 'RS_PROJTYP']) #resampling [RIGHT]
unique(EVA_meta[EVA_meta$RS_CODE == 'FRCH_0001', 'Manipulate_y_n']) #'Y'
EVA_meta[EVA_meta$RS_CODE == 'FRCH_0001' & EVA_meta$Manipulate_y_n == 'Y', 'Manipulate_y_n'] <- 'N'

#HU_0002: resampling -> site comparison; manipulate remains 'Y'
unique(EVA_meta[EVA_meta$RS_CODE == 'HU_0002', 'RS_PROJTYP']) #resampling (TBC)
unique(EVA_meta[EVA_meta$RS_CODE == 'HU_0002', 'Manipulate_y_n']) #'N' 'Y'
EVA_meta[EVA_meta$RS_CODE == 'HU_0002' & EVA_meta$RS_PROJTYP == 'resampling', 'RS_PROJTYP'] <- 'Site comparison'

#IT_0001a to IT_0001e: manipulate 'Y' to 'N' (also change remaining empty fields to 'N')
unique(EVA_meta[EVA_meta$RS_CODE %in% paste0('IT_0001', c('a', 'b', 'c', 'd', 'e')), 'RS_PROJTYP']) #resampling [RIGHT]
EVA_meta[EVA_meta$RS_CODE %in% paste0('IT_0001', c('a', 'b', 'c', 'd', 'e')) & EVA_meta$Manipulate_y_n == '', 'Manipulate_y_n'] #character(0)
EVA_meta[EVA_meta$RS_CODE %in% paste0('IT_0001', c('a', 'b', 'c', 'd', 'e')) & EVA_meta$Manipulate_y_n == 'Y', 'Manipulate_y_n'] <- 'N'

#DE_0037_88: resampling to permanent
unique(EVA_meta[EVA_meta$RS_CODE == 'DE_0037_88', 'RS_PROJTYP']) #resampling
unique(EVA_meta[EVA_meta$RS_CODE == 'DE_0037_88', 'Manipulate_y_n']) #'Y'
EVA_meta[EVA_meta$RS_CODE == 'DE_0037_88', 'RS_PROJTYP'] <- 'permanent'


#-- 4 - exclude manipulated plots

#to do so, use the RS_PROJTYP column. In short, drop observations associated with 'permanent (man)'

#check unique values in RS_PROJTYP
unique(EVA_meta$RS_PROJTYP) #some typo to correct

#208 cases of Resampling -> resampling
EVA_meta$RS_PROJTYP[EVA_meta$RS_PROJTYP == 'Resampling'] <- 'resampling'

#4456 cases of Permanent (man) -> permanent (man)
EVA_meta$RS_PROJTYP[EVA_meta$RS_PROJTYP == 'Permanent (man)'] <- 'permanent (man)'

#check frequency of cases
table(EVA_meta$RS_PROJTYP) #44433 cases of 'permanent (man)'

#to keep control observations, combine RS_PROJTYP and Manipulate_y_n
#permanent (man) with 'N' in Manipulate_y_n are control plots

#count cases of permanent (man) with 'N' in Manipulate_y_n; these are control plots
sum(EVA_meta$RS_PROJTYP == 'permanent (man)' & EVA_meta$Manipulate_y_n == 'N') #12514

#count number of 'permanent (man)' that are in core EVA, i.e. outside ReSurvey
sum(EVA_meta$RS_PROJTYP == 'permanent (man)' & !EVA_meta$Truly_RS) #82 (these plots have not been resampled yet)
EVA_meta[which(EVA_meta$RS_PROJTYP == 'permanent (man)' & !EVA_meta$Truly_RS), ]

#for simplicity, I remove all 'permanent (man)' that are not simultaneously 'N' in Manipulate_y_n

#plot id of manipulated plots
man_plot_id <- EVA_meta[EVA_meta$RS_PROJTYP == 'permanent (man)', 'PlotID']

#plot id of control plots
con_plot_id <- EVA_meta[EVA_meta$RS_PROJTYP == 'permanent (man)' & EVA_meta$Manipulate_y_n == 'N', 'PlotID']

#keep plot id only of truly manipulated plots (exclude control plot ids)
man_plot_id <- man_plot_id[!man_plot_id %in% con_plot_id] #31919 (44433 - 12514)

#exclude manipulated plots
EVA_meta <- EVA_meta[!EVA_meta$PlotID %in% man_plot_id, ]

#check
#all(!(EVA_meta$PlotID %in% man_plot_id)) #TRUE

#nrow(EVA_meta) #1440708

#--account for location uncertainty (?)

#this column is very hard to use because precision of coordinates and relocation accuracy have been mixed up

anyNA(EVA_meta$Location_uncertainty_m) #T
sum(is.na(EVA_meta$Location_uncertainty_m)) #146459

#check if this is available only for RS data -> nope, this is available for EVA plots too
unique(EVA_meta[EVA_meta$ReSurvey_plot_Y_N == 'N', 'Location_uncertainty_m'])


#--extract plot coordinates and transform the object to an sf object

#I'll use this object to create a bbox ('area of interest') and remove plots falling outside the bbox

EVA_meta.sp <- st_as_sf(EVA_meta[c('PlotID', 'Longitude', 'Latitude')], coords = c('Longitude', 'Latitude'), crs = 4326)

#mapview(EVA_meta.sp)


#---------------create bbox to crop ecoregions including EVA data

#no points occur at weird locations (extreme longlat), but some plots seem to be anyway quite far from Continental Europe
st_bbox(EVA_meta.sp)

mapview(st_as_sf(st_as_sfc(st_bbox(EVA_meta.sp))), color = 'yellow')

#check range of longlat
range(st_coordinates(EVA_meta.sp)[, 1]) #-73.60  64.84
range(st_coordinates(EVA_meta.sp)[, 2]) #0.00000 80.14912

#try setting boundaries for creating the bbox - I want to exclude a small number of points, 1/10000 per long/lat extreme
seq_for_bbox <- seq(0, 1, by = 1/1e5)

#this should exclude the 1e-5 * 2 points at the long and lat boundaries
bbox_xminmax <- unname(quantile(st_coordinates(EVA_meta.sp)[, 1], probs = seq_for_bbox)[c(2, (length(seq_for_bbox) - 1))])
bbox_yminmax <- unname(quantile(st_coordinates(EVA_meta.sp)[, 2], probs = seq_for_bbox)[c(2, (length(seq_for_bbox) - 1))])

#check number of points excluded - this only counts those that would be out along either x- or y-axis
sum(st_coordinates(EVA_meta.sp)[, 1] < bbox_xminmax[1] | st_coordinates(EVA_meta.sp)[, 1] > bbox_xminmax[2]) #25
sum(st_coordinates(EVA_meta.sp)[, 2] < bbox_yminmax[1] | st_coordinates(EVA_meta.sp)[, 2] > bbox_yminmax[2]) #29

#check number of points included
sum((st_coordinates(EVA_meta.sp)[, 1] >= bbox_xminmax[1] & st_coordinates(EVA_meta.sp)[, 1] <= bbox_xminmax[2]) &
      (st_coordinates(EVA_meta.sp)[, 2] >= bbox_yminmax[1] & st_coordinates(EVA_meta.sp)[, 1] <= bbox_yminmax[2])) #1440671

bbox_xy <- c(setNames(bbox_xminmax, c('xmin', 'xmax')), setNames(bbox_yminmax, c('ymin', 'ymax')))

#create the bbox - note that I am creating the bbox in longlat as I am better able to find extreme locations using longitude and latitude
#I am anyway going to transform the bbox to a proj CRS because st_crop does not work properly with longlat data over vast regions
eu_box <- st_bbox(bbox_xy, crs = st_crs(4326))

#transform bbox to EPSG:3035 - this avoids using st_crop with unproj data over large extents (see code at the end of the script)
eu_box <- st_transform(eu_box, crs = epsg_proj)

#save longitude and latitude
EVA_meta.sp$Longitude <- st_coordinates(EVA_meta.sp)[ , 1]
EVA_meta.sp$Latitude <- st_coordinates(EVA_meta.sp)[ , 2]

#transform to EPSG:3035
EVA_meta.sp <- st_transform(x = EVA_meta.sp, crs = epsg_proj)

#also create a bbox polygon
eu_bbox.polyg <- st_as_sfc(eu_box)

#crop maps of ecoregions and continental lands
identical(st_crs(eu_ecoregions.proj), st_crs(eu_box)) #T
identical(st_crs(eu_land.proj), st_crs(eu_box)) #T

#crop maps

#ecoregions
eu_ecoregions.proj <- st_crop(x = eu_ecoregions.proj, y = eu_box)

#continental lands
eu_land.proj <- st_crop(x = eu_land.proj, y = eu_box)

#check results
mapview(eu_ecoregions.proj)
mapview(eu_land.proj)

#--check whether I get same duplicates (plots with same coordinate pairs) when coords are longlat vs. LAEA

#here I'm just checking that coordinates are considered as duplicates regardless of their 'formatting' (planar vs. geographic)
#I'm doing this using a sample of 1e5 points, otherwise it may take hours to run

#large sample of row positions
set.seed(94858)

row_smp <- sample(1e5, replace = F)

subs_laea <- EVA_meta.sp[row_smp, ]

duply_laea <- st_equals(subs_laea)

st_geometry(subs_laea) <- NULL

subs_longlat <- st_as_sf(subs_laea, coords = c('Longitude', 'Latitude'), crs = 4326)

duply_longlat <- st_equals(subs_longlat)

identical(duply_laea, duply_longlat) #T

rm(subs_laea, duply_laea, subs_longlat, duply_longlat)

#--dissolve eu_land.proj and create a X-meter buffer (in QGIS)

#the resulting polygon will be used to crop EVA plots, excluding those falling more than X meters off the coast
#I'm using eu_land.proj because the coastline of EU countries is delineated much more precisely than in the ecoregion layer

#dissolve eu_land.proj

diss_eu_land.proj <- st_union(eu_land.proj)

mapview(st_as_sf(diss_eu_land.proj))

#write the polygon and create buffer in QGIS

st_write(obj = diss_eu_land.proj, dsn = '/MOTIVATE/GDM_EuropeanEcoregions/GSHHG/diss_eu_land/diss_eu_land_proj.shp')

#import the land polygon with a 100-m buffer area. This was computed in QGIS (it took 5 mins and 38 secs)

eu_land_proj_100m_buf <- st_read(dsn = "/MOTIVATE/GDM_EuropeanEcoregions/GSHHG/diss_eu_land/eu_land_proj_buffer100m.shp")

class(eu_land_proj_100m_buf)

#--mask EVA points. Points falling outside the buffer area will be filtered out.

#st_geometry(EVA_meta.sp) -> geometry type: POINT

#from https://r-spatial.github.io/sf/reference/geos_binary_pred.html
tmp_inters <- st_intersects(x = EVA_meta.sp, y = eu_land_proj_100m_buf, sparse = FALSE)

tmp_inters <- apply(tmp_inters, 1, any)

length(which(!tmp_inters)) #20917 points falling off the coast + 100 m buffer

#extract points falling off the coast
tmp_off_the_coast <- EVA_meta.sp[which(!tmp_inters), ]

#write the shapefile to disk and have a look at the points in QGIS
st_write(obj = tmp_off_the_coast, dsn = "/MOTIVATE/GDM_EuropeanEcoregions/GSHHG/diss_eu_land/pts_off_the_coast.shp")

#filter out points falling outside buffer area
EVA_meta.sp <- EVA_meta.sp[tmp_inters, ]

#check that there are not points falling outside the bbox
sum(!apply(st_intersects(EVA_meta.sp, y = eu_bbox.polyg, sparse = FALSE), 1, any)) #0 [OK]

#--assign each plot to an ecoregion based on a distance criterion

#I am using st_nearest_feature() as points falling along the border of the ecoregion map
#would be not assigned to any ecoregion otherwise

EVA_meta.sp <- st_join(x = EVA_meta.sp, y = eu_ecoregions.proj[c(2, 4, 5)], left = T, join = st_nearest_feature)

anyNA(EVA_meta.sp$ECO_NAME) #FALSE

unique(EVA_meta.sp$ECO_NAME) #73 ecoregions

#remove '-' from ECO_NAME and replace white spaces by '_'

unique(grep(pattern = '-', x = EVA_meta.sp$ECO_NAME, value = T))

EVA_meta.sp$ECO_NAME <- gsub(pattern = '-', replacement = ' ', x = EVA_meta.sp$ECO_NAME)

EVA_meta.sp$ECO_NAME <- gsub(pattern = ' ', replacement = '_', x = EVA_meta.sp$ECO_NAME)

length(unique(EVA_meta.sp$ECO_NAME)) #still 73

#--filter out from EVA_meta plots not included in EVA_meta.sp

#save coordinates in LAEA
EVA_meta.sp$X_laea <- st_coordinates(EVA_meta.sp)[, 1]
EVA_meta.sp$Y_laea <- st_coordinates(EVA_meta.sp)[, 2]

#drop geometry
st_geometry(EVA_meta.sp) <- NULL

#drop out plots in EVA_meta not included in EVA_meta.sp
EVA_meta <- EVA_meta[EVA_meta$PlotID %in% EVA_meta.sp$PlotID, ]

#check
setdiff(EVA_meta$PlotID, EVA_meta.sp$PlotID) #integer(0)
setdiff(EVA_meta.sp$PlotID, EVA_meta$PlotID) #integer(0)

#left_join EVA_meta and EVA_meta.sp using PlotID
class(EVA_meta$PlotID); class(EVA_meta.sp$PlotID)

#check
sum(duplicated(EVA_meta$PlotID)) #0

EVA_meta <- dplyr::left_join(x = EVA_meta, y = EVA_meta.sp[c(1, 4, 5, 6, 7, 8)], by = 'PlotID')

#--split EVA and RS (ReSurvey) and check duplicates and true duplicates, so that I have clear numbers for EVA

#compare ReSurvey_plot_Y_N and Truly_RS, see if they provide same information
#unique(EVA_meta$ReSurvey_plot_Y_N) #only 'Y' or 'N'
identical(which(EVA_meta$ReSurvey_plot_Y_N == 'Y'), which(EVA_meta$Truly_RS)) #TRUE, which means it doesn't matter which field I use to isolate RS plots

#PlotID of RS time-series
RS_PlotID <- EVA_meta$PlotID[EVA_meta$Truly_RS]

length(RS_PlotID) #124308

#Isolate RS from EVA
RS_meta <- EVA_meta[EVA_meta$PlotID %in% RS_PlotID, ]

#check
unique(RS_meta$ReSurvey_plot_Y_N) #'Y'

#check RS_CODE
unique(RS_meta$RS_CODE)
length(unique(RS_meta$RS_CODE)) #441

#check ReSurvey_project - here, notice that the number of RS projects is lower than the RS codes
unique(RS_meta$ReSurvey_project)
length(unique(RS_meta$ReSurvey_project)) #437

#check ReSurvey_site
unique(RS_meta$ReSurvey_site)

#check ReSurvey_plot
unique(RS_meta$ReSurvey_plot)

#check classes
sapply(RS_meta[c("RS_CODE", "ReSurvey_site", "ReSurvey_plot")], class) #all chr

#create column to identify RS time-series (following header's description here - see word file)
RS_meta$ID_RS_ts <- with(RS_meta, paste(RS_CODE, ReSurvey_site, ReSurvey_plot, sep = '_'))

#check how many t-s are present
length(unique(RS_meta$ID_RS_ts)) #38659

#create column to identify RS observations (it should mean every RS sampling unit, a time-series is made of N sampling units)
#(following header's description here - see word file)

#first check ReSurvey_observation
class(RS_meta$ReSurvey_observation) #chr
sum(duplicated(RS_meta$ReSurvey_observation)) #18814
anyNA(RS_meta$ReSurvey_observation) #F
sum(RS_meta$ReSurvey_observation == '') #0

#this should not have duplicates
RS_meta$ID_RS_observ <- with(RS_meta, paste(RS_CODE, ReSurvey_site, ReSurvey_plot, ReSurvey_observation, sep = '_'))

#check duplicates
sum(duplicated(RS_meta$ID_RS_observ)) #2695 (?)
length(unique(RS_meta$ID_RS_observ)) #121613

unique(RS_meta$RS_CODE[duplicated(RS_meta$ID_RS_observ)])

#an example of duplicates for ID_RS_observ
RS_meta[RS_meta$RS_CODE == 'SI_0004' & duplicated(RS_meta$ID_RS_observ), ]
RS_meta[RS_meta$RS_CODE == 'CZ_0019_025' & duplicated(RS_meta$ID_RS_observ), ]
RS_meta[RS_meta$RS_CODE == 'HU_0002' & duplicated(RS_meta$ID_RS_observ), ]
RS_meta[RS_meta$RS_CODE == 'RS_0001' & duplicated(RS_meta$ID_RS_observ), ]

RS_meta[RS_meta$ID_RS_observ == 'SI_0004_Trans1_Trans1_1_Trans1_1_2021', ] #these have same coordinates, but were sampled on (1st of) Jan. and 12th August
RS_meta[RS_meta$ID_RS_observ == 'SI_0004_Trans1_Trans1_2_Trans1_2_2021', ] #same as above
RS_meta[RS_meta$ID_RS_observ == 'CZ_0019_025_KRAL_SNEZNIK_KRAL_SNEZNIK01_KRAL_SNEZNIK01_2019', ] #different coordinates, but same Date_of_recording
RS_meta[RS_meta$ID_RS_observ == 'HU_0002_SZAPKIS_KISMADV_RIP_KISMADV_RIP_16_5', ] #same coordinates and Date_of_recording
RS_meta[RS_meta$ID_RS_observ == 'RS_0001_Raska_Raska_2018_Raska_2018', ] #these are many

#I should now extract the PlotID of the first observation(s) of the time series, so that I can re-integrate them in EVA
#a way to do this is to keep observations associated with the first year of the time series
anyNA(RS_meta$ID_RS_ts) #FALSE

id_ts_start <- lapply(unique(RS_meta$ID_RS_ts), function(i) {
  
  #extract data for i-th time series
  dtf <- RS_meta[RS_meta$ID_RS_ts == i, ]
  #coerce Sampl_year to integer
  dtf$Sampl_year <- as.integer(dtf$Sampl_year)
  #extract first year of the time series
  frst_yr <- min(dtf$Sampl_year)
  #extract PlotID for the first year
  frst_plot_id <- dtf[dtf$Sampl_year == frst_yr, 'PlotID']
  #return result
  return(frst_plot_id)
  
  })

#mean (and median) number of 'starting' plots
mean(lengths(id_ts_start)) # 1.066479
median(lengths(id_ts_start)) #1

names(id_ts_start) <- unique(RS_meta$ID_RS_ts)

#check examples of time series
RS_meta[RS_meta$ID_RS_ts == names(id_ts_start)[1], c('PlotID', 'Sampl_year')]
RS_meta[RS_meta$ID_RS_ts == names(id_ts_start)[3], c('PlotID', 'Sampl_year')]

#check cases with more than a single starting plot
head(which(lengths(id_ts_start) > 1))
sum(lengths(id_ts_start) > 1) #430 (out of 38659)

#is it possible that RS_meta includes starting plots of time series that haven't been resurveyed yet
id_ts_start[42]
RS_meta[RS_meta$ID_RS_ts == names(id_ts_start)[42], c('PlotID', 'Sampl_year')]
#check Date_of_recording, which could be different from Sampl_year
RS_meta[RS_meta$PlotID %in% c(391, 415), ] #these plots were indeed surveyed on different days

#another example
id_ts_start[2770]
RS_meta[RS_meta$ID_RS_ts == names(id_ts_start)[2770], c('PlotID', 'Sampl_year')]
#check Date_of_recording
RS_meta[RS_meta$PlotID %in% c(6647, 6648, 758100, 758101, 758102, 758103), ] #these have different dates of recording, except for a couple of cases

#check how many plots I would recover
length(unlist(id_ts_start)) #41229 (out of 124308)

#sanity check
#sum(duplicated(unlist(id_ts_start))) #0

#there's one case with 105 plots that begin the time series
table(lengths(id_ts_start[lengths(id_ts_start) > 1]))

#have a look at this case; it could be an N-to-N resurvey strategy
RS_meta[RS_meta$PlotID %in% id_ts_start[[which(lengths(id_ts_start) == 105)]], ] #RS_PROJTYP is 'resampling' though

#keep the ids of plots that start the time series in EVA
#to do so, exclude ids of these 'starting' plots from the set of ids in RS_meta, the remaining plot ids (not starting plots) can be filtered put from EVA_meta

#unlist id_ts_start to obtain a vector of ids of 'starting plots'
id_ts_start <- unlist(id_ts_start)

id_ts_to_drop <- RS_meta$PlotID[!RS_meta$PlotID %in% id_ts_start]

#--separate EVA (this version includes starting plots from RS) and RS 

EVA_meta <- EVA_meta[!EVA_meta$PlotID %in% id_ts_to_drop, ]

nrow(EVA_meta) #1336712

#keep only ids contained in id_ts_to_drop for RS_meta
RS_meta <- RS_meta[RS_meta$PlotID %in% id_ts_to_drop, ]

nrow(RS_meta) #83079

#--save EVA_meta, RS_meta and RS_meta.sp

#save EVA_meta, RS_meta and RS_meta.sp to avoid having to re-run code above in case
#these objects have to be (re)used from this line on

save(EVA_meta, RS_meta, EVA_meta.sp, file = '/MOTIVATE/GDM_EuropeanEcoregions/tmp_obj/EVA_RS_meta.RData')

#--delete objects that are not used eventually



#--check on (spatial) duplicates

#this check will be on two types of duplicates:

#1) true duplicates: plots having exactly the same coordinates, Sampl_year, taxonomic composition and species' relative abundances
#2) duplicates: plots having exactly the same coordinates (but not necessarily Sampl_year, tax comp and species' rel abund)

#true duplicates are assessed on EVA + RS as a whole because the idea is to keep only a single observation per group of true duplicates
#duplicates are identified within each period (before / after 2000) and within each ecoregion separately because the idea is to fit
#period- and ecoregion-specific models - duplicates can still be used if they belong to different periods

#the check on duplicated coordinates is based on: https://github.com/r-spatial/sf/issues/1893

#first of all, subset EVA_meta.sp and create RS_meta.sp - these versions will only include plot ids contained in corresponding non-spatial versions

EVA_meta.sp <- EVA_meta.sp[EVA_meta.sp$PlotID %in% EVA_meta$PlotID, ]

all(colnames(EVA_meta.sp) %in% colnames(RS_meta)) #TRUE

RS_meta.sp <- RS_meta[colnames(EVA_meta.sp)] 

#make the *.sp dtfs spatial

EVA_meta.sp <- st_as_sf(x = EVA_meta.sp, coords = c('X_laea', 'Y_laea'), crs = epsg_proj)

RS_meta.sp <- st_as_sf(x = RS_meta.sp, coords = c('X_laea', 'Y_laea'), crs = epsg_proj)

#then, create a 'Period' column that splits the datasets in period1 and period2
#depending on whether plots were sampled before or after 2000

#attach Sampl_year to both datasets

identical(EVA_meta$PlotID, EVA_meta.sp$PlotID) #T
identical(RS_meta$PlotID, RS_meta.sp$PlotID) #T

EVA_meta.sp$Sampl_year <- EVA_meta$Sampl_year

RS_meta.sp$Sampl_year <- RS_meta$Sampl_year

#create Period column
EVA_meta.sp$Period <- ifelse(as.integer(EVA_meta.sp$Sampl_year) < 2000L, 'period1', 'period2')

table(EVA_meta.sp$Period) #period1: 794354; period2: 542358 

RS_meta.sp$Period <- ifelse(as.integer(RS_meta.sp$Sampl_year) < 2000L, 'period1', 'period2')

table(RS_meta.sp$Period) #period1: 11710; period2: 71369

#-- 1 - check on true duplicates

#!!this check is implemented on EVA and RS together as true duplicates should be detected across (and removed from) the whole dataset

#import EVA data on community composition
EVA_veg <- fread('C://MOTIVATE/EVA_dataset/200_Motivate20250415_taxa_concepts_unified.csv', data.table = FALSE)

#select only relevant columns
EVA_veg <- EVA_veg[c("PlotObservationID", "Matched concept corrected", "Cover %", "Layer")]

#rename columns
colnames(EVA_veg)[c(1, 2, 3)] <- c('PlotID', 'Species_name', 'Cover_perc')

#check NA
anyNA(EVA_veg) #FALSE
#check values of Layer field
sort(unique(EVA_veg$Layer)) #from 0 to 9

#check if all plots in EVA_meta.sp and RS_meta.sp are contained in EVA_veg
sum(!EVA_meta.sp$PlotID %in% unique(EVA_veg$PlotID)) #382
sum(!RS_meta.sp$PlotID %in% unique(EVA_veg$PlotID)) #729

EVA_meta.sp[which(!EVA_meta.sp$PlotID %in% unique(EVA_veg$PlotID)), ]

#an example, PlotID 1 is missing from EVA_veg
1 %in% EVA_veg$PlotID #FALSE

#and vice-versa -> although this is possible since I removed quite a few plots from EVA_/RS_meta.sp that should still be present in EVA_veg
sum(!unique(EVA_veg$PlotID) %in% EVA_meta.sp$PlotID) #1071455
sum(!unique(EVA_veg$PlotID) %in% RS_meta.sp$PlotID) #2325435

#save PlotID that are present in EVA_meta.sp and RS_meta.sp but are not present in EVA_veg
#I won't be able to use these plots since there's no data on the vegetation
EVA_plots_not_in_veg <- setdiff(EVA_meta.sp$PlotID, unique(EVA_veg$PlotID))
RS_plots_not_in_veg <- setdiff(RS_meta.sp$PlotID, unique(EVA_veg$PlotID))

length(EVA_plots_not_in_veg) + length(RS_plots_not_in_veg) #1111 (382+729)
all(!c(EVA_plots_not_in_veg, RS_plots_not_in_veg) %in% unique(EVA_veg$PlotID)) #T
sum(duplicated(EVA_plots_not_in_veg)) #0
sum(duplicated(RS_plots_not_in_veg)) #0
sum(duplicated(c(EVA_plots_not_in_veg, RS_plots_not_in_veg))) #0

#remove these plots from EVA_meta.sp and RS_meta.sp
#notice that I'll use these two objects to filter out true duplicates from EVA_meta and RS_meta
#so, it's ok that I'm removing plots_not_in_veg from the spatial objects at this stage
EVA_meta.sp <- EVA_meta.sp[!EVA_meta.sp$PlotID %in% EVA_plots_not_in_veg, ]
RS_meta.sp <- RS_meta.sp[!RS_meta.sp$PlotID %in% RS_plots_not_in_veg, ]

#subset EVA_veg to only keep records in EVA_meta + RS_meta
intersect(x = EVA_meta.sp$PlotID, y = RS_meta.sp$PlotID) #integer(0)
EVA_veg <- EVA_veg[EVA_veg$PlotID %in% c(EVA_meta.sp$PlotID, RS_meta.sp$PlotID), ]

#fine, EVA_veg has the same plots that are included in EVA_meta.sp and RS_meta.sp
length(unique(EVA_veg$PlotID)) #1418680
length(EVA_meta.sp$PlotID) + length(RS_meta.sp$PlotID) #1418680
setdiff(x = unique(EVA_veg$PlotID), y = c(EVA_meta.sp$PlotID, RS_meta.sp$PlotID)) #0
setdiff(x = c(EVA_meta.sp$PlotID, RS_meta.sp$PlotID), y = unique(EVA_veg$PlotID)) #0

#aggregate species' cover data in case there are multiple records for the same species within a plot

#set EVA_veg as a data.table for fast manipulation of the dataset
#using setDT avoids storing a copy of EVA_veg (as when using as.data.table)
setDT(EVA_veg)

#count number of plots having multiple records for the same species
nrow(EVA_veg[, if(any(duplicated(Species_name))) 1L else 0L, by = PlotID][V1 == 1L]) #219273 plots with at least one duplicated species

#get plot ids with duplicated species names (and different Layer values)
EVA_veg[, .(IDs = paste(which(duplicated(Species_name)), collapse = '-'), DiffLyr = if(uniqueN(Layer) == 1) 0L else 1L ), by = PlotID][IDs != '' & DiffLyr != 0L, 'PlotID']

#examples of these plots: 3552, 3587, 3597, 3600, 3601 
ex_plots_sp_dup <- EVA_veg[PlotID %in% c(3552, 3587, 3597, 3600, 3601)]

#order (asc) Layer field by PlotID and Species_name - this is to order cover values by increasing layer (for combine.cover function)
#however, note that layer order doesn't affect computation made by combine.cover
EVA_veg <- EVA_veg[order(Layer), .SD, by = .(PlotID, Species_name)]

#aggregate cover values in case of multiple records for the same species

#combine.cover aggregates cover data over multiple records for the same species
combine.cover <- function(x) {
  while (length(x) > 1) {
    x[2] <- x[1]+(100-x[1])*x[2]/100
    x <- x[-1]
  }
  return(x)
}

#order of values does not affect result of combine.cover
#replicate(n = 5, combine.cover(x = c(1, 30, 50, 22, 12)[sample(5)]), simplify = T)

#aggregate cover data
EVA_veg <- EVA_veg[, .(Cover_perc = combine.cover(Cover_perc)), by = .(PlotID, Species_name)]

#check - looks correct
#ex_plots_sp_dup[PlotID == 3601]
#EVA_veg[PlotID == 3601]

#identify groups of duplicated plots

#rbind EVA_meta.sp and RS_meta.sp to detect true duplicates over the whole dataset
identical(colnames(EVA_meta.sp), colnames(RS_meta.sp)) #T

#rbind EVA_meta.sp and RS_meta.sp to work on whole dataset
EVA_RS_sp <- rbind(EVA_meta.sp, RS_meta.sp)

#run st_equals to identify groups of plots with same coordinates
EVA_RS_true_duply <- st_equals(EVA_RS_sp)

#find unique groups of duplicated points

#count number of groups of plots (at least two duplicates)
sum(lengths(EVA_RS_true_duply) > 1) #831573 - note that this number counts multiple time the same combination of duplicates

#first subset EVA_RS_true_duply to only keep groups of plots (at least two duplicates)
EVA_RS_true_duply <- EVA_RS_true_duply[lengths(EVA_RS_true_duply) > 1]

#from plot position to plot id, plus collapse vector of plot ids to then find unique sets
EVA_RS_true_duply <- lapply(EVA_RS_true_duply, function(grp) {
  
  #from plot position to plot id
  pos2id <- EVA_RS_sp[grp, 'PlotID', drop = T]
  
  #collapse vector of plot id
  res_id <- paste(pos2id, collapse = '-')
  
  return(res_id)
  
})

#unlist EVA_RS_true_duply and keep only unique series
EVA_RS_true_duply <- unlist(EVA_RS_true_duply)

#check how many duplicates
sum(duplicated(EVA_RS_true_duply)) #677955

#keep only unique sets
EVA_RS_true_duply <- unique(EVA_RS_true_duply) #153618

#transform EVA_RS_true_duply to a list of groups of plot ids (coerced to integer)
EVA_RS_true_duply <- lapply(EVA_RS_true_duply, function(grp) {
  
  #split plot ids
  plot_ids <- unlist(strsplit(x = grp, split = '-'))
  
  #coerce to integer
  plot_ids <- as.integer(plot_ids)
  
  return(plot_ids)
  
  })

#check if there are list items with length < 2
sum(lengths(EVA_RS_true_duply) < 2) #0 (nope)

#number of plots
sum(lengths(EVA_RS_true_duply)) #831573 (this is the number of plots having at least a duplicate)

#now exclude from EVA_RS_true_duply the plot ids that are unique in terms of Sampl_year
#true duplicates are indeed associated with the same Sampl_year

#select only relevant columns in EVA_RS_sp and drop geometry
EVA_RS_sp_sub <- EVA_RS_sp[, c('PlotID', 'Sampl_year'), drop = T]

#drop PlotID not included in EVA_RS_true_duply
EVA_RS_sp_sub <- EVA_RS_sp_sub[EVA_RS_sp_sub$PlotID %in% unlist(EVA_RS_true_duply), ]

#create table with PlotID and corresponding group id - an id which groups plot ids that are potential true duplicates (and surely duplicates)

#name EVA_RS_true_duply
names(EVA_RS_true_duply) <- paste('Grp', seq_along(EVA_RS_true_duply), sep = '_')

#create table of PlotID and corresponding group id
EVA_RS_true_duply <- data.frame(PlotID = unlist(EVA_RS_true_duply), GrpID = rep(names(EVA_RS_true_duply), times = lengths(EVA_RS_true_duply)))

#join Sampl_year column to EVA_RS_true_duply by plot
class(EVA_RS_true_duply$PlotID) #int
class(EVA_RS_sp_sub$PlotID) #int

EVA_RS_true_duply <- dplyr::left_join(x = EVA_RS_true_duply, y = EVA_RS_sp_sub, by = 'PlotID')

#coerce EVA_RS_true_duply to a data.table to speed up data manipulation
setDT(EVA_RS_true_duply)

#get ids of groups with at least one duplicate for year (this is easier than getting plot ids per duplicated years, which would require a list)
grp_with_duply <- EVA_RS_true_duply[, .(HasDupYear = if(any(duplicated(Sampl_year))) TRUE else FALSE), by = .(GrpID)][which(HasDupYear), GrpID]

#check how many groups have at least a duplicate for year
length(grp_with_duply) #103279 (out of 153618 groups)

#check
EVA_RS_true_duply[GrpID %in% grp_with_duply[10]]

#focus only on plot ids belonging to grp_with_duply

#subset EVA_RS_true_duply to only keep groups in grp_with_duply
#sum(duplicated(grp_with_duply)) #0

EVA_RS_true_duply <- EVA_RS_true_duply[GrpID %in% grp_with_duply]

#subset EVA_veg to only keep plots in EVA_RS_true_duply

#notice that EVA_veg_sub will be used to create unique combination of Species_name and Cover_value for each plot
#these combinations will be transformed to strings to check true duplicates within each group of duplicates

EVA_veg_sub <- EVA_veg[PlotID %in% EVA_RS_true_duply[, PlotID]]

#exclude plots including only a single species
EVA_veg_more1sp <- EVA_veg_sub[, .N, by = PlotID]

EVA_veg_more1sp <- EVA_veg_more1sp[N > 1, PlotID]

length(EVA_veg_more1sp) #658869

#subset both EVA_veg_sub and EVA_RS_true_duply to keep plots in EVA_veg_more1sp
EVA_veg_sub <- EVA_veg_sub[PlotID %in% EVA_veg_more1sp]
EVA_RS_true_duply <- EVA_RS_true_duply[PlotID %in% EVA_veg_more1sp]

#create string of Species_name and Cover_perc - this will be used to compare taxonomic composition (and relative abundances) among potential true duplicates
#order is used to order rows by PlotID and, within each PlotID, by Species_name; this guarantees species are ordered the same ways within each plot
#notice that ordering takes place over the whole dataset (and not within chunks defined by PlotID)
#PlotComp is a new variable that will include the combination of Species_name and Cover_perc - the combo is coerced to a string for comparison among plots
EVA_veg_sub <- EVA_veg_sub[order(PlotID, Species_name), .(PlotComp = paste0(Species_name, ':', Cover_perc, collapse = '|')), by = PlotID]

#left join between EVA_RS_true_duply and EVA_veg_sub for attaching PlotComp to EVA_RS_true_duply (same as left_join(EVA_RS_true_duply, EVA_veg_sub, by = 'PlotID))
EVA_RS_true_duply <- EVA_veg_sub[EVA_RS_true_duply, on = 'PlotID']

#find true duplicates
#the idea is to create a new (list)column, PlotList, that includes the plot ids belonging to groups defined by unique combo of GrpID, Sampl_year, PlotComp
#any time two plots share the same GrpID, Sampl_year and have the same PlotComp, the correspdoning element in PlotList will have length > 1 
EVA_RS_true_duply <- EVA_RS_true_duply[, .(PlotList = list(PlotID)), by = .(GrpID, Sampl_year, PlotComp)]

#check how many cases of true duplicates
sum(lengths(EVA_RS_true_duply$PlotList) > 1) #8214

#subset EVA_RS_true_duply to keep only elements with true duplicates
EVA_RS_true_duply <- EVA_RS_true_duply[lengths(EVA_RS_true_duply[, PlotList]) > 1]

#number of true duplicates
length(unlist(EVA_RS_true_duply$PlotList)) #17280

#check distribution of sizes of groups of true duplicates
range(lengths(EVA_RS_true_duply$PlotList)) #2 70
table(lengths(EVA_RS_true_duply$PlotList))

#one group has 70 true duplicates (!) in 2005
EVA_RS_true_duply[lengths(EVA_RS_true_duply$PlotList) == 70] #Grp_82139

#notice that there might be multiple groups of true duplicates for the same year in Grp_82139
EVA_meta[EVA_meta$PlotID %in% unlist(EVA_RS_true_duply[lengths(EVA_RS_true_duply$PlotList) == 70, PlotList]), ] #looks like these are all in EVA

#get id of plots to remove to exclude true duplicates

#extract list of groups of plot ids of true duplicates
plot_ids_trued <- EVA_RS_true_duply[, PlotList]

#keep only first plot ids of each group
#there should remain sum(lengths(plot_ids_trued) - 1)
sum(lengths(plot_ids_trued)) #17280

plot_ids_drop <- lapply(plot_ids_trued, function(id) {
  
  to_drop <- id[-1]
  
  return(to_drop)
  
})

#check [RIGHT]
sum(lengths(plot_ids_trued) - 1) #9066
length(unlist(plot_ids_drop)) #9066

#unlist plot_ids_drop
plot_ids_drop <- unlist(plot_ids_drop)

#delete objects to free memory
#rm(EVA_RS_sp, EVA_veg_sub, EVA_RS_sp_sub, EVA_veg_more1sp, grp_with_duply, EVA_RS_true_duply)


#-- 2 - check on duplicates

#first of all drop plot_ids_drop from EVA_meta.sp and RS_meta.sp
class(EVA_meta.sp); class(RS_meta.sp) #sf data.frame

#check how many true duplicates are contained in the two data.frames
sum(EVA_meta.sp$PlotID %in% plot_ids_drop) #8859 - nearly the total number of plots
sum(RS_meta.sp$PlotID %in% plot_ids_drop) #207

#drop the plots
EVA_meta.sp <- EVA_meta.sp[!EVA_meta.sp$PlotID %in% plot_ids_drop, ]
RS_meta.sp <- RS_meta.sp[!RS_meta.sp$PlotID %in% plot_ids_drop, ]

#- write a function for obtaining groups of duplicates within each period and ecoregion

#example on a single ecoregion:

#sample(unique(RS_meta.sp$ECO_NAME), 1)

tmp_eco <- RS_meta.sp[RS_meta.sp$ECO_NAME == "Pyrenees_conifer_and_mixed_forests", ]

table(tmp_eco$Period)

#forgetting about the existence of periods for a while
tmp_eco_dup <- st_equals(tmp_eco)

#this behaves as a list obj
tmp_eco_dup[[1]] #1 67 -> these are plot positions
c(1, 67) %in% tmp_eco[['PlotID']] #FALSE
tmp_eco[c(1, 67), ]

lengths(tmp_eco_dup) #there are many cases of duplicates
length(which(lengths(tmp_eco_dup) > 1)) #1423
length(tmp_eco_dup[lengths(tmp_eco_dup) > 1]) #1423
tmp_eco_dup[lengths(tmp_eco_dup) > 1]

#check behavior of as.data.frame method
as.data.frame(tmp_eco_dup)

#first check sample size per ecoregion x period
smp_size_prd <- do.call(rbind, lapply(unique(EVA_meta.sp$ECO_NAME), function(nm) {
  
  dtf <- EVA_meta.sp[EVA_meta.sp$ECO_NAME == nm, , drop = T]
  #smp size period1
  smp_sz_prd1 <- sum(dtf$Period == 'period1')
  #smp size period2
  smp_sz_prd2 <- sum(dtf$Period == 'period2')
  #return a matrix
  res <- matrix(c(smp_sz_prd1, smp_sz_prd2), nrow = 1, dimnames = list(nm, c('Smp_prd1', 'Smp_prd2')))
  return(res)
  }))


get_duply_ids <- function(x) {
  
  #get periods (there might be ecoregions having all plots included within a period)
  #creating this vector prevents running the function in case there are no observations for a given period 
  prd_vec <- sort(unique(x[['Period']])) #sort is used make Period1 come before Period2
  
  #loop across periods and compute duply stats
  res <- lapply(prd_vec, function(prd) {
    dtf_pd <- x[x[['Period']] == prd, ]
    #check duplicates
    res_pd <- st_equals(dtf_pd) #this returns an sgbp object (sparse matrix) -> it has a data.frame method
    #subset res_pd to only keep elements with duplicates
    res_pd <- res_pd[lengths(res_pd) > 1]
    #anticipated return if res_pd has zero length
    if((length(res_pd) == 0)) return(NA)
    #from plot position (within dtf_pd!) to plot id
    res_pd <- lapply(res_pd, function(i) {
      #use vector of plot position to get corresponding PlotID
      pos_to_id <- dtf_pd[i, 'PlotID', drop = T]
      #collapse vector of plot ids to identify unique sets
      pos_to_id <- paste(pos_to_id, collapse = '-')
      return(pos_to_id)
    })
    #unlist res_pd
    res_pd <- unlist(res_pd)
    #keep only unique sets
    res_pd <- unique(res_pd)
    return(res_pd)
    })
  
  names(res) <- prd_vec
  
  return(res)
  
  }

#check - seems to work
tmp_duply <- get_duply_ids(x = tmp_eco)
tmp_duply$period1
tmp_eco[tmp_eco$PlotID %in% c(342951, 342952, 342953, 342954, 342955), ]

#check on empty dataset for a given period
get_duply_ids(x = EVA_meta.sp[EVA_meta.sp$ECO_NAME == 'Mediterranean_High_Atlas_juniper_steppe', ])

#run the function on EVA and RS

#EVA_duply and RS_duply includes the list of plot ids that belong to groups of duplicates
#these objects will be used later on to select plots for analyses

#--duplicates within EVA
EVA_duply <- lapply(unique(EVA_meta.sp$ECO_NAME), function(nm) {
  
  x <- EVA_meta.sp[EVA_meta.sp$ECO_NAME == nm, ]
  
  res <- get_duply_ids(x = x)
  
  return(res)
  
  })

names(EVA_duply) <- unique(EVA_meta.sp$ECO_NAME)

#check number of groups per period
sapply(EVA_duply, function(i) lengths(i))

#try out strsplit on an example
length(unlist(strsplit(EVA_duply[['Corsican_montane_broadleaf_and_mixed_forests']][['period2']], split = '-')))

#check NAs (if any)
sum(sapply(EVA_duply, function(i) sum(sapply(i, anyNA)))) #9 cases

#get number of duplicates plots, i.e. plots having at least a duplicate
EVA_num_duply <- lapply(EVA_duply, function(x) {
  
  #loop over periods
  res <- lapply(x, function(i) {
    #anticipated return if i is an NA
    if(length(i) == 1) {
      if(is.na(i)) return(NA)
    }
    #otherwise get length of vector containing all plot IDs
    vec_dup <- unlist(strsplit(i, split = '-'))
    #check if there are duplicates - there shouldn't be
    if(sum(duplicated(vec_dup)) > 0) stop('There should not be duplicated plot ids')
    #return length
    num_dup <- length(vec_dup)
    return(num_dup)
    
  })
  
  return(res)
  
})

#get numbers for period1 and period2 separately

#period1
EVA_duply_prd1 <- lapply(EVA_num_duply, function(i) i[['period1']])

#get rid of NULL entries, ecoregions for which there were no plots for period1
#which(sapply(EVA_duply_prd1, is.null))
EVA_duply_prd1 <- EVA_duply_prd1[!sapply(EVA_duply_prd1, is.null)]

#unlist
EVA_duply_prd1 <- unlist(EVA_duply_prd1)

#period2
EVA_duply_prd2 <- lapply(EVA_num_duply, function(i) i[['period2']])
EVA_duply_prd2 <- EVA_duply_prd2[!sapply(EVA_duply_prd2, is.null)]
EVA_duply_prd2 <- unlist(EVA_duply_prd2)

#attach number of duplicates plots to smp_size_prd
smp_size_prd <- as.data.frame(smp_size_prd)

smp_size_prd$Dup_prd1 <- unname(EVA_duply_prd1[rownames(smp_size_prd)])
smp_size_prd$Dup_prd2 <- unname(EVA_duply_prd2[rownames(smp_size_prd)])

#compute proportions
smp_size_prd$Prop_prd1 <- with(smp_size_prd, Dup_prd1/Smp_prd1)
smp_size_prd$Prop_prd2 <- with(smp_size_prd, Dup_prd2/Smp_prd2)

#sapply(EVA_num_duply, function(i) sum(unlist(i), na.rm = T))

#attach ECO_NAME col
smp_size_prd$ECO_NAME <- row.names(smp_size_prd)

row.names(smp_size_prd) <- seq_len(nrow(smp_size_prd))


#--duplicates within RS

#RS
RS_duply <- lapply(unique(RS_meta.sp$ECO_NAME), function(nm) {
  
  x <- RS_meta.sp[RS_meta.sp$ECO_NAME == nm, ]
  
  res <- get_duply_ids(x = x)
  
  return(res)
  
  })

#---------

##now keep only plot ids in EVA_meta.sp and RS_meta.sp from EVA_meta and RS_meta
class(EVA_meta); class(RS_meta) #data.frame data.frame
nrow(EVA_meta); nrow(RS_meta) #1.336.712 83.079
nrow(EVA_meta.sp); nrow(RS_meta.sp) #1.327.471 82.143

EVA_meta <- EVA_meta[EVA_meta$PlotID %in% EVA_meta.sp$PlotID, ]

RS_meta <- RS_meta[RS_meta$PlotID %in% RS_meta.sp$PlotID, ]

#--check missing data for key fields, including plot size, habitat classification (look for unclassified habitats and empty fields), period (1980-2022) 

#-plot size

#check how many observations miss data on plot size
class(EVA_meta$Releve_area_m2) #num
length(unique(EVA_meta$Releve_area_m2)) #1218 different plot sizes (!)
range(unique(EVA_meta$Releve_area_m2), na.rm = T) #0.01 9999.99
range(unique(RS_meta$Releve_area_m2), na.rm = T) #0.01 3000

#check NAs
sum(is.na(EVA_meta$Releve_area_m2)) #311.439 plots missing plot size
sum(is.na(RS_meta$Releve_area_m2)) #1879

#drop observations with missing data for plot size
EVA_meta <- EVA_meta[!is.na(EVA_meta$Releve_area_m2), ]
RS_meta <- RS_meta[!is.na(RS_meta$Releve_area_m2), ]

#-habitat classification

#check NAs, unclassified habitats (marked as Y), empty fields (marked as X)
sum(is.na(EVA_meta$Eunis_lev1)) #0
sum(is.na(RS_meta$Eunis_lev1)) #0
unique(EVA_meta$Eunis_lev1)
unique(RS_meta$Eunis_lev1)

#drop unclassified and empty fields for Eunis_lev1
length(which(EVA_meta$Eunis_lev1 %in% c("Y", "X"))) #192.881
length(which(RS_meta$Eunis_lev1 %in% c("Y", "X"))) #6091

EVA_meta <- EVA_meta[!EVA_meta$Eunis_lev1 %in% c("Y", "X"), ]
RS_meta <- RS_meta[!RS_meta$Eunis_lev1 %in% c("Y", "X"), ]

sort(table(EVA_meta$Eunis_lev1)) #R (grasslands) and T (forests) are the most abundant
sum(EVA_meta$Eunis_lev1 %in% c('R', 'T')) #516.172 relevés classified as either grassland or forest

sort(table(RS_meta$Eunis_lev1)) #R and T are still the most abundant
sum(RS_meta$Eunis_lev1 %in% c('R', 'T')) #58.485


#-period

#the min period is determined by the availability of the drivers over time
#at the moment, the biggest constraint is the human modification index, whose time series starts in 1990
#so, it's not possible to use data collected before 1980

#I'm not excluding these data at this stage, though. This can be done afterwards.

sum(as.integer(EVA_meta$Sampl_year) < 1980L | as.integer(EVA_meta$Sampl_year) > 2022L) #197.358


#--match PlotID between EVA_meta, RS_meta and their corresponding .sp versions

length(setdiff(EVA_meta$PlotID, EVA_meta.sp$PlotID)) #0
length(setdiff(EVA_meta.sp$PlotID, EVA_meta$PlotID)) #504.320

length(setdiff(RS_meta$PlotID, RS_meta.sp$PlotID)) #0
length(setdiff(RS_meta.sp$PlotID, RS_meta$PlotID)) #7970

EVA_meta.sp <- EVA_meta.sp[EVA_meta.sp$PlotID %in% EVA_meta$PlotID, ]

RS_meta.sp <- RS_meta.sp[RS_meta.sp$PlotID %in% RS_meta$PlotID, ]






#---------------subset EVA_meta to keep only grasslands


Grass_meta <- EVA_meta[which(EVA_meta$Eunis_lev1 == 'R'), ]

#check ecoregions including grasslands
unique(Grass_meta$ECO_NAME) #60

#drop obervations with Sampl_period not included between 1980 and 2022
Grass_meta <- Grass_meta[as.integer(Grass_meta$Sampl_year) <= 2022 & as.integer(Grass_meta$Sampl_year) >= 1980, ]

#check
range(as.integer(Grass_meta$Sampl_year))

#create period column
Grass_meta$Period <- ifelse(as.integer(Grass_meta$Sampl_year) < 2000L, 'period1', 'period2')

#check number of plots per ecoregion
anyNA(Grass_meta$ECO_NAME) #FALSE
table(Grass_meta$ECO_NAME)

#check number of plots per ecoregion x period
Grass_eco_prd <- as.data.frame(with(Grass_meta, table(ECO_NAME, Period)))

#check how many ecoregions have sufficient data (at least 1000 per period)
Grass_eco_nm <- unique(Grass_meta$ECO_NAME)

#coerce ECO_NAME and Period to chr in Grass_eco_prd
Grass_eco_prd$ECO_NAME <- as.character(Grass_eco_prd$ECO_NAME)
Grass_eco_prd$Period <- as.character(Grass_eco_prd$Period)

#14 ecoregions
Grass_eco_minN <- names(which(sapply(Grass_eco_nm, function(nm) all(Grass_eco_prd[Grass_eco_prd$ECO_NAME == nm, 'Freq'] >= 1000))))

Grass_eco_prd[Grass_eco_prd$ECO_NAME %in% Grass_eco_minN, ]

#check spatial distribution of ecoregions
mapview(eu_ecoregions.proj[eu_ecoregions.proj$ECO_NAME %in% gsub(pattern = '_', replacement = ' ', x = Grass_eco_minN), ])


##FROM HERE!!!!!!!


#compute Bray-Curtis in groups of duplicates - identify plots below a given threshold of Bray-Curtis and drop others
#the idea is to apply a random spatial shift to plots that are similar in terms of B-C

#this allows avoiding resampling, which may complicate things for matching

#https://stackoverflow.com/questions/10225098/understanding-exactly-when-a-data-table-is-a-reference-to-vs-a-copy-of-another

#example on a given ecoregion and period: Appenine_deciduous_montane_forests - Period 1
ex_admf <- EVA_duply$Appenine_deciduous_montane_forests$period1

#get plot ids
ex_admf <- as.integer(unlist(lapply(ex_admf, function(i) unlist(strsplit(i, split = '-')))))

#compute dissimilarity matrix

#FROM HERE!!!!!!!!!!!!!

#subset EVA_veg
#eva_veg_admf <- EVA_veg[PlotID %in% ex_admf]

#dcast to pivot_wide
#dcast(eva_veg_admf, PlotID ~ Species_name, value.var = 'Cover_perc', fill = 0)

#compute B-C for all plots and then subset the matrix to check mean dissimilarity within groups

#the idea is that building a unique dissimilarity matrix for all plots in a group is faster than creating a group specific dissimilarity matrix


#below, an example on a couple of plots

#re-order EVA_veg by PlotID (and within PlotID by Cover_perc)
EVA_veg[order(PlotID, Cover_perc)]

#filter plot ids
EVA_veg[PlotID %in% c(30014, 30015)]

#from long to wide
#tidyr::pivot_wider(data = EVA_veg[PlotID %in% c(30014, 30015)], names_from = 'Species_name', values_from = 'Cover_perc', values_fill = 0)
#dcast(DT.m1, family_id + age_mother ~ child, value.var = "dob")
dcast(EVA_veg[PlotID %in% c(30014, 30015)], PlotID ~ Species_name, value.var = 'Cover_perc', fill = 0)

#drop PlotID (?)

#compute B-C
as.matrix(vegan::vegdist(dcast(EVA_veg[PlotID %in% c(30014, 30015)], PlotID ~ Species_name, value.var = 'Cover_perc', fill = 0), 'bray'))



#-- check species with 0 cover in EVA_veg






