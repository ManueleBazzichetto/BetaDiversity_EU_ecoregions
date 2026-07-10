#This script includes R functions for extracting the values of climatic and topographic variables, and of the human modification index
#from each grassland and forest location.
#The functions are used in the Drivers_extraction script, so they must be created before running that script.


#-----------Function for filling gaps in climatic data

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



#-----------Function for extracting climatic data

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


#-----------Function for assigning hmi year to each grassland and forest location

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



#-----------Function for extracting the human modification index

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



#-----------Function for filling gaps in human modification index

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

