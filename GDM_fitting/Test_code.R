
#Test code for fitting GDMs and computing deviance partitions using the gdm R package.

library(gdm)

rm(list = ls())

gc()

#load example of tables formatted for GDMs

load('/MOTIVATE/GDM_EuropeanEcoregions/Data_for_analyses/tables_for_gdm_grassland/Alps_cmf_grass.RData')

#the object loaded in the environment is named 'tmp_list'

class(tmp_list) #list
length(tmp_list) #2
names(tmp_list) #"Period1" "Period2"

sapply(tmp_list, nrow)

head(tmp_list$Period1)

# ---------------- fit gdm

#period 1
test_run_prd1 <- gdm::gdm(data = tmp_list$Period1, geo = T)

gc()

#period 2
test_run_prd2 <- gdm::gdm(data = tmp_list$Period2, geo = T)

gc()

#check summary
summary(test_run_prd1)
summary(test_run_prd2)

#expl. deviance
test_run_prd1$explained
(1 - (test_run_prd1$gdmdeviance/test_run_prd1$nulldeviance))*100

#is geo included as a predictor?
test_run_prd1$geo

#test_run_prd1$sample
#test_run_prd1$sumCoeff


# ---------------- compute deviance partitions

#!!Notice that the total deviance explained (and unexplained) returned with the deviance partitioning considers only the contribution
#of the variables used as components/groups for the partitions.

#create list of variables' groups to be used for deviance partitioning
#note that if env = c('Roughness'), or any other fourth predictor, is included, the function will return the following error: Cannot partition more than three variables sets
#for this reason, I'm only including climatic variables, the human modification index and geographic distance (setting partSpace = T)
test_part_vars <- list(climate = c('Prcp', 'Tavg'), human = c('Hmi_value')) #notice that geo dist is included by setting the partSpace = TRUE

test_period_devpart <- gdm::gdm.partition.deviance(sitePairTable = tmp_list$Period1, varSets = test_part_vars, partSpace = TRUE)

test_period_devpart


# ---------------- extract iSplines from fitted GDMs

#extract an example of a gdm fitted object for testing

#imported as tmp_list
load('/MOTIVATE/GDM_EuropeanEcoregions/GDM_fitting/GDM_output_grassland/ItaScl_sdf_grass_gdm_out.RData')

#tmp_list (gdm output) has 2 period-specific elements
length(tmp_list)
names(tmp_list)

#each element is a list with 2 elements: "GDM_out" "Dev_out"
names(tmp_list$Period1)

#extract fitted onj for Period1
test_gdm_obj <- tmp_list$Period1$GDM_out

#extract ispline using the gdm::isplineExtract function
#test_ispl_obj is a list with 2 elements: x (matrix with 200 values of each predictor) and y (fitted value of the spline)
test_ispl_obj <- gdm::isplineExtract(model = test_gdm_obj)

length(test_ispl_obj)
names(test_ispl_obj) #"x" "y"
nrow(test_ispl_obj$x) #200
class(test_ispl_obj$x) #matrix

#test identity numeric values after stacking
identical(as.vector(test_ispl_obj$x), stack(as.data.frame(test_ispl_obj$x))[['values']]) #TRUE

#coerce each element to a data.frame and use the stack function to concatenate values
test_ispl_obj <- lapply(test_ispl_obj, function(dtf) stack(as.data.frame(dtf)))

head(test_ispl_obj$x)

#extract ind column
test_ispl_obj[[c(1, 2)]]

#coerce test_ispl_obj to a dtf
test_ispl_obj <- as.data.frame(test_ispl_obj)

#check x.ind and y.ind are identical
identical(test_ispl_obj$x.ind, test_ispl_obj$y.ind) #TRUE

#drop y.ind
test_ispl_obj$y.ind <- NULL

#coerce x.ind to a chr
test_ispl_obj$x.ind <- as.character(test_ispl_obj$x.ind)

#modify colnames
colnames(test_ispl_obj) <- c('Values_x', 'Variable_x', 'Values_y')

#add ECO_NM and Period cols
test_ispl_obj$ECO_NM <- 'ItaScl_sdf'
test_ispl_obj$Period <- 'Period1'

head(test_ispl_obj)

# --

#write a function that replicates what achieved above

#x is the fitted GDM for a given ecoregion x period combo

reformat_ispl_obj <- function(x, eco_nm, prd_nm) {
  
  require(gdm)
  
  #extract iSplines from fitted GDM
  isplines_obj <- gdm::isplineExtract(model = x)
  
  #reformat x and y elements included in isplines_obj
  isplines_obj <- lapply(names(isplines_obj), function(nm) {
    
    #coerce to a data.frame
    ispl_el <- as.data.frame(isplines_obj[[nm]])
    
    #stack ispl_el
    ispl_el <- stack(ispl_el)
    
    #coerce ind col to chr
    ispl_el$ind <- as.character(ispl_el$ind)
    
    #modify colnames
    colnames(ispl_el) <- paste(c('Values', 'Variable'), nm, sep = '_')
    
    return(ispl_el)
    
  })
  
  #if Variables_* cols match, cbind ispl_el(s)
  
  if(identical(isplines_obj[[c(1, 2)]], isplines_obj[[c(2, 2)]])) {
    
    isplines_obj <- as.data.frame(isplines_obj)
    
    #drop Variable_y, which is a duplicate
    isplines_obj$Variable_y <- NULL
    
    #add ECO_NM and Period columns
    isplines_obj$ECO_NM <- eco_nm
    
    isplines_obj$Period <- prd_nm
    
    return(isplines_obj)
    
  } else {
    
    stop("Variable names do not match between x and y objects")
    
  }
  
}


#compare results
test_ispl_obj_v2 <- reformat_ispl_obj(x = test_gdm_obj, eco_nm = 'ItaScl_sdf', prd_nm = 'Period1')

identical(test_ispl_obj, test_ispl_obj_v2) #TRUE








