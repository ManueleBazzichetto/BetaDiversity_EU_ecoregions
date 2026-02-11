#this code provides an example on how to run statistical matching on a test ecoregion,
#and on ecoregions belonging to habitat type 'R' (grasslands) and 'T' (forests)

library(data.table)
library(sf)
#library(terra)

#using statistical matching should allow obtaining samples for the two periods that
#have comparable (balanced) distribution of (observed) covariates. In turn, this should allow making inference
#on the 'effect' of time on the change in the influence of the environmental drivers
#on beta diversity spatial patterns.

#note that GDMs are correlative models, so they don't allow causal inference on the effect of the drivers
#within each period.

#the problem is that EVA data are spatially, temporally and environmentally biased.
#Comparing two samples that vary a lot in terms of environmental conditions may confound
#effect of time - differences in the influence of drivers may be driven by differences
#in samples' covariates rather than, e.g. an increase in the importance of anthropogenic
#pressures over time.

#at this stage, true duplicates, that is relevés with identical Sampl_year, coordinates and plant composition
#have been already excluded.

library(MatchIt) #installed on Oct. 24th 2025
#library(cobalt) ##installed on Oct. 24th 2025

#--CALIPER

#remember that the idea is to use a caliper to ameliorate covariates balance (both SDM and var ratio), without loosing too many observations!
#caliper increases balance at the cost of precision - therefore the best method + caliper combination is the one
#that allows reaching balance and keeping as many obs as possible

#caliper should only be used to reach balance of otherwise unbalanced covariates!

#for info on how to set a caliper,
#check https://stats.stackexchange.com/questions/621309/propensity-scores-variable-k1-matching-with-caliper-some-technical-questions
#calipers can be specified for individual covariates or for the propensity score (or both)
#when std.caliper = TRUE calipers are expressed in std dev units, otherwise in raw units
#for example, if using caliper on prop score, if std.caliper = TRUE the difference in prop score expressed in std dev units must be
#below the caliper value for the control obs to be matched with a treat obs. Otherwise the diff in raw values of the prop score
#must be lower than the caliper value.
#setting link = "linear.logit" when a caliper is specified means that the difference between logit of the prop score betw control and tr
#is actually being computed - so, the units refer to the logit-transformed prop score (linear prop score)
#note that this is the recommended way of using caliper for prop score (see link above)

#STRATEGY FOR USING CALIPER:
#For var ratio:
#Use a caliper, starting from a caliper = 2, for all covariates with var ratio > 2 - decrease it by a .1 step until var ratio is within the .5 - 2 range
#For SMD:
#If balanced is not reached for a single covariate and its functions/interactions - start with std caliper = 2 and progressively decrease it (by 0.1 std dev unit) until all covariates are below .1 abs SMD 
#If balance is not reached for a single covariate + prop score - start with caliper = 2 for both and decrease both by a 0.1 step - stop decreasing caliper when both variables are within
#the 0.10 threshold - if a good balance is achieved for one of the two, keep decreasing caliper only for the one that is still unbalanced
#If multiple or all covariates (and/or their functions and interactions) are not balanced (+ prop score) - start with caliper = 2, decrease it by a 0.1 step - stop for a covariate when is balanced - continue
#with the others until they are balanced - decrease caliper for all variables when they are involved in unbalanced interactions

#--K:1 ratio

#IMPORTANT! When using K:1 ratio, all treatment units should receive the same number of control obs, otherwise weights should be included
#in the analyses! See 'Common mistakes' section on this page: https://cran.r-project.org/web/packages/MatchIt/vignettes/estimating-effects.html

#--Do not estimate prop score in genetic matching

#to do not estimate prop score in genetic matching, set distance = 'mahalanobis' - see examples method_genetic()

#--Reproducibility

#remember to allow reproducibility in case of matching methods based on stochastic processes (genetic)

#-------------------------------------------------example on a test ecoregion

#select columns that are actually used for matching
#here, I'm using the Cant_mf dataset

Cant_mf_for_match <- Cant_mf[c('PlotID', 'Eunis_lev1', 'EunisVerbose_lev1', 'Sampl_year', 'X_laea', 'Y_laea', 'Period', 'Prcp', 'Tavg',
                               'Elevation', 'Roughness', 'Slope', 'Hmi_value')]

#check habitat types
unique(Cant_mf_for_match$Eunis_lev1)

#check range of years
range(as.integer(Cant_mf_for_match$Sampl_year))

#remove observations for which env data are missing
#here I should probably drop NAs for all env variables, including those not used for matching
#these observations will be anyway excluded from the final sample used to fit the GDM
#therefore, they should not be considered for matching
env_var_nm <- c('Prcp', 'Tavg', 'Elevation', 'Roughness', 'Slope', 'Hmi_value')

sapply(Cant_mf_for_match[env_var_nm], function(cl) sum(is.na(cl)))

Cant_mf_for_match <- Cant_mf_for_match[complete.cases(Cant_mf_for_match[env_var_nm]), ]

#check imbalance of period's sample size
length(which(Cant_mf_for_match$Period == 'period1'))
length(which(Cant_mf_for_match$Period == 'period2'))

#-------------------statistical matching

#check NAs for Period and Sampl_year
sum(is.na(Cant_mf_for_match$Period)) #0
sum(is.na(Cant_mf_for_match$Sampl_year)) #0

#Check the estimand argument, which can be used to change the labeling of tr vs. control groups
#Changing the labeling doesn't mean a different estimand is actually being targeted - it's a mere trick to set labels

#create 1/0 version of the Period column -> it is completely arbitrary which group is 1 and 0 in this case as the estimand is going to be the ATO
#setting the smallest group in terms of sample size as 1 allows easily using the k:1 ratio, because more control units (zeros) are matched to the ones, i.e. treatment unit
tr_grp <- names(which.min(table(Cant_mf_for_match$Period)))

Cant_mf_for_match$Period_bin <- ifelse(Cant_mf_for_match$Period == tr_grp, 1, 0)

#check
table(Cant_mf_for_match$Period_bin)

#check initial imbalance
Cant_mf_init_imb <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = NULL, distance = 'glm')

summary(Cant_mf_init_imb)

plot(summary(Cant_mf_init_imb))

#------nearest neighbor

#--propensity score

#logit
Cant_mf_pscore1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest', distance = 'glm')

Cant_mf_pscore1
summary(Cant_mf_pscore1)
summary(Cant_mf_pscore1, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Cant_mf_pscore1, un = FALSE, interactions = T)
plot(Cant_mf_pscore1, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1))
plot(summary(Cant_mf_pscore1, interactions = T))

#order
Cant_mf_pscore1_ord1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                          distance = 'glm', m.order = 'closest')

summary(Cant_mf_pscore1_ord1, un = F)
plot(Cant_mf_pscore1_ord1, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1_ord1, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1_ord1))
plot(summary(Cant_mf_pscore1_ord1, interactions = T))

#order + ratio -> this increases precision at the expenses of balance
Cant_mf_pscore1_ord1_ratio1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                          distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Cant_mf_pscore1_ord1_ratio1, un = F)
plot(Cant_mf_pscore1_ord1_ratio1, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1_ord1_ratio1, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1_ord1_ratio1))
plot(summary(Cant_mf_pscore1_ord1_ratio1, interactions = T))

#--Mahalanobis

Cant_mf_mahala1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest', distance = 'mahalanobis')

summary(Cant_mf_mahala1)
plot(Cant_mf_mahala1, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1))
plot(summary(Cant_mf_mahala1, interactions = T))

#order
Cant_mf_mahala1_ord1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                          distance = 'mahalanobis', m.order = 'closest')

summary(Cant_mf_mahala1_ord1, un = F)
plot(Cant_mf_mahala1_ord1, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1_ord1))
plot(summary(Cant_mf_mahala1_ord1, interactions = T))

#order + ratio -> this increases precision at the expenses of balance
Cant_mf_mahala1_ord1_ratio1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                                 distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Cant_mf_mahala1_ord1_ratio1)
plot(Cant_mf_mahala1_ord1_ratio1, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1_ord1_ratio1))
plot(summary(Cant_mf_mahala1_ord1_ratio1, interactions = T))

#--robust Mahalanobis

Cant_mf_rob_mahala1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match,
                         method = 'nearest', distance = 'robust_mahalanobis')

summary(Cant_mf_rob_mahala1)
plot(Cant_mf_rob_mahala1, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1))
plot(summary(Cant_mf_rob_mahala1, interactions = T))

#order
Cant_mf_rob_mahala1_ord1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                          distance = 'robust_mahalanobis', m.order = 'closest')

summary(Cant_mf_rob_mahala1_ord1, un = F)
plot(Cant_mf_rob_mahala1_ord1, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1_ord1))
plot(summary(Cant_mf_rob_mahala1_ord1, interactions = T))


#order + ratio -> this increases precision at the expenses of balance
Cant_mf_rob_mahala1_ord1_ratio1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'nearest',
                                 distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Cant_mf_rob_mahala1_ord1_ratio1)
plot(Cant_mf_rob_mahala1_ord1_ratio1, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1_ord1_ratio1))
plot(summary(Cant_mf_rob_mahala1_ord1_ratio1, interactions = T))

#------genetic matching

#see genoud (from Matching) to set pop.size: https://stackoverflow.com/questions/48500267/pop-size-argument-in-genmatch-respectively-genoud

#using GMD without prop.score
Cant_mf_genetic1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match,
                            method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Cant_mf_genetic1)
plot(Cant_mf_genetic1, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic1))
plot(summary(Cant_mf_genetic1, interactions = T))

#check https://kosukeimai.github.io/MatchIt/reference/method_genetic.html for more details
#try out with distance = some method for prop. score -> this computes the prop. score and adds it to the covariates used to
#compute the GMD

#using GMD with prop.score
Cant_mf_genetic2 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Cant_mf_genetic2)
plot(Cant_mf_genetic2, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic2))
plot(summary(Cant_mf_genetic2, interactions = T))

#order: it is not possible to set m.order = 'closest'

#ratio
Cant_mf_genetic1_ratio1 <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_for_match,
                             method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Cant_mf_genetic1_ratio1)
plot(Cant_mf_genetic1_ratio1, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic1_ratio1))
plot(summary(Cant_mf_genetic1_ratio1, interactions = T))


#------extract matched data

#extract matched data:
#match_data keeps the original order of the data
#a distance column is added when propensity score is used as the distance (see help)
#weights should include matching weights by default (though in my case weights should be all 1)
#subclass - not sure what this means when using nearest neighbor

Cant_mf_matched <- match_data(Cant_mf_mahala1_ord1_ratio1, data = Cant_mf_for_match)

#check if the distance field is returned when explicitly setting the distance column to a name
head(match_data(Cant_mf_mahala1_ord1_ratio1, data = Cant_mf_for_match, distance = 'Mahalanobis')) #nope
#check if a distance field is returned when prop score is actually estimated/used in the matchit call
head(match_data(Cant_mf_pscore1)) #yes

#weights are all 1
unique(Cant_mf_matched$weights) #all 1

#the subclass field is a factor - values are all 1 as if every observation was a subclass
head(Cant_mf_matched$subclass)
all(table(unique(as.integer(Cant_mf_matched$subclass))) == 1) #TRUE


#the weights and subclass fields can be safely excluded - these columns won't be used in the GDM
Cant_mf_matched <- Cant_mf_matched[!colnames(Cant_mf_matched) %in% c('weights', 'subclass')]

#subset EVA_veg to only include PlotID within Cant_mf_matched
Cant_mf_veg <- as.data.frame(EVA_veg[PlotID %in% Cant_mf_matched$PlotID])

#check number of species with 0 Cover_perc
sum(Cant_mf_veg$Cover_perc == 0) #0

#save data for fitting the GDM
save(Cant_mf_matched, Cant_mf_veg, file = '/Users/Manuele/Desktop/Cant_mf_GDM_data.RData')

#check whether matching stats can be extracted from summary(matchit obj)

#the summary object is a list, which includes, among other objects, the sum.matched obj, that is
#"a matrix of balance statistics for each covariate in the matched sample" (from help of summary() function)
sum_ex <- summary(Cant_mf_mahala1_ord1_ratio1, interactions = T)

#from this object, I can check if SMD for all covariates is < .1, check if the variance ratio is between 0.5 and 2
#and compute the average SMD as a means to compare different methods
sum_ex$sum.matched

#from this object I can compare sample size of matched units when a caliper is used
sum_ex$nn

#from the object returned by matchit, I can check if weights are either 1 or 0
all(Cant_mf_mahala1_ord1_ratio1$weights %in% c(0, 1))

#this can be used to check that matching was performed without replacement
!Cant_mf_mahala1_ord1_ratio1$info$replace[1]

#write a function to make these checks and extract summary stats to compare the performance of different matching methods

#note that the GWD stats is affected by the number of covariates and, more specifically, by the inclusion of the prop score

check_match_performance <- function(mobj, smd_thr = 0.1, var_thrs = c(0.5, 2)) {
  
  #check matching was performed w/o replacement
  if(mobj$info$replace[1]) stop('Matching was performed with replacement')
  
  #check that all weights are either 0 (unmatched obs) or 1 (matched obs)
  #this should indicate whether tr obs were all matched to the same number of control obs
  if(!all(mobj$weights %in% c(0, 1))) stop("Matching weights are not all 0 or 1")
  
  #check if a caliper was used
  cal_used <- 'caliper' %in% names(mobj)
  
  #extract matrix of summary stats
  m_sum <- summary(mobj, interactions = T)
  
  m_stats <- m_sum$sum.matched
  
  smd <- m_stats[, 'Std. Mean Diff.']
  
  #check smd is always below ths
  if(!(all(smd <= smd_thr))) stop(paste('SMD > than', smd_thr, 'detected!', sep = ' '))
  
  #compute average of abs smd
  avg_smd <- mean(abs(smd))
  
  #compute generalized weighted distance
  gwd_smd <- sum(smd)
  
  #compute max abs std diff
  max_smd <- max(abs(smd))
  
  #check if var ratio of all covariates is within thrs
  var_rat <- m_stats[, 'Var. Ratio']
  
  var_rat <- all(var_rat >= var_thrs[1] & var_rat <= var_thrs[2])
  
  #return result
  res <- data.frame(Avg_smd = avg_smd, Gwd_smd = gwd_smd, Max_smd = max_smd, Var_rat = var_rat)
  
  #if a caliper was used, check number of matched units
  if(cal_used) {
    
    message('Caliper detected!')
    
    #extract matrix with sample sizes
    nn_mat <- m_sum$nn
    
    #sum sample sizes of Control and Treated
    nn_mat <- sum(nn_mat['Matched', ])
    
    res$Smp_sz <- nn_mat
    
  }
  
  return(res)
  
  }

#test
check_match_performance(mobj = Cant_mf_mahala1_ord1_ratio1, smd_thr = .05)

#-------------------------------------------------statistical matching of grasslands

grass_col_to_use <- c('PlotID', 'Eunis_lev1', 'EunisVerbose_lev1', 'Sampl_year', 'X_laea', 'Y_laea', 'Period', 'Prcp', 'Tavg',
                      'Elevation', 'Roughness', 'Slope', 'Hmi_value')

exists('env_var_nm')
all(env_var_nm %in% colnames(Grass_sel_meta)) 

#drop NAs for environmental drivers
sapply(Grass_sel_meta[env_var_nm], function(cl) sum(is.na(cl)))

grass_loc_to_drop <- complete.cases(Grass_sel_meta[env_var_nm])

sum(!grass_loc_to_drop) #2135 (out of 202595)

Grass_sel_meta <- Grass_sel_meta[grass_loc_to_drop, ]

#check abundance of observations for the different ecoregions
table(Grass_sel_meta$ECO_NAME)
table(Grass_sel_meta$ECO_NAME, Grass_sel_meta$Period) #only Tyrrhenian_Adriatic_sclerophyllous_and_mixed_forests has less than 1k obs in period1

#check habitat type
unique(Grass_sel_meta$Eunis_lev1) #'R'

#check NAs in Sampl_year and Period
sum(is.na(Grass_sel_meta$Sampl_year)) #0
sum(is.na(Grass_sel_meta$Period)) #0


# ------  1. Alps_conifer_and_mixed_forests - GENETIC ------

Alps_cmf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Alps_conifer_and_mixed_forests', grass_col_to_use] 

names(which.min(table(Alps_cmf_grass$Period))) #period1

Alps_cmf_grass$Period_bin <- ifelse(Alps_cmf_grass$Period == "period1", 1, 0)

#check
table(Alps_cmf_grass$Period_bin)


#------check initial imbalance

Alps_cmf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = NULL, distance = 'glm')

summary(Alps_cmf_init_gr)
plot(summary(Alps_cmf_init_gr))
plot(Alps_cmf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Alps_cmf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest', distance = 'glm')

Alps_cmf_pscore1_gr
summary(Alps_cmf_pscore1_gr, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Alps_cmf_pscore1_gr, un = FALSE, interactions = T)
plot(Alps_cmf_pscore1_gr, type = 'jitter', interactive = F)
plot(Alps_cmf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_pscore1_gr))
plot(summary(Alps_cmf_pscore1_gr, interactions = T))
plot(summary(Alps_cmf_pscore1_gr, interactions = T, un = F))

#order
Alps_cmf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                distance = 'glm', m.order = 'closest')

summary(Alps_cmf_pscore1_ord1_gr, un = F)
plot(Alps_cmf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Alps_cmf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_pscore1_ord1_gr, un = F))
plot(summary(Alps_cmf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                       distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Alps_cmf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(Alps_cmf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(Alps_cmf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(Alps_cmf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on all covariates and prop score
Alps_cmf_pscore1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                           distance = 'glm', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                           caliper = c(1.6, Elevation = 1.1, Roughness = 1.1, Slope = 1.2), link = "linear.logit")

summary(Alps_cmf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T) #78 tr obs are not matched due to the caliper (besides control obs) - Matched (ESS) < Matched
plot(summary(Alps_cmf_pscore1_ord1_ratio1_cal1_gr, interactions = T, un = F))

#--Mahalanobis

Alps_cmf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Alps_cmf_mahala1_gr)
plot(Alps_cmf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_mahala1_gr, un = F))
plot(summary(Alps_cmf_mahala1_gr, interactions = T, un = F))

#order
Alps_cmf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                distance = 'mahalanobis', m.order = 'closest')

summary(Alps_cmf_mahala1_ord1_gr, un = F)
plot(Alps_cmf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_mahala1_ord1_gr, un = F))
plot(summary(Alps_cmf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                       distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Alps_cmf_mahala1_ord1_ratio1_gr, un = T, interactions = T)
plot(Alps_cmf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Alps_cmf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on all covariates
Alps_cmf_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                           distance = 'mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                           caliper = c(Elevation = .9, Roughness = .9, Slope = .9))

summary(Alps_cmf_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #101 tr obs are not matched due to the caliper (besides control obs) - Matched (ESS) < Matched
plot(summary(Alps_cmf_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = F))


#--robust Mahalanobis

Alps_cmf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass,
                               method = 'nearest', distance = 'robust_mahalanobis')

summary(Alps_cmf_rob_mahala1_gr)
plot(Alps_cmf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_rob_mahala1_gr, un = F))
plot(summary(Alps_cmf_rob_mahala1_gr, interactions = T, un = F))

#order
Alps_cmf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                    distance = 'robust_mahalanobis', m.order = 'closest')

summary(Alps_cmf_rob_mahala1_ord1_gr, un = F)
plot(Alps_cmf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_rob_mahala1_ord1_gr, un = F))
plot(summary(Alps_cmf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                           distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Alps_cmf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Alps_cmf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Alps_cmf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on all covariates
Alps_cmf_rob_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'nearest',
                                               distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                               caliper = c(Elevation = .9, Roughness = .9, Slope = .9))

summary(Alps_cmf_rob_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #110 tr obs are not matched due to the caliper (besides control obs) - Matched (ESS) < Matched
plot(summary(Alps_cmf_rob_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
Alps_cmf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Alps_cmf_genetic1_gr)
plot(Alps_cmf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_genetic1_gr, un = F))
plot(summary(Alps_cmf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
Alps_cmf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Alps_cmf_genetic2_gr)
plot(Alps_cmf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_genetic2_gr, un = F))
plot(summary(Alps_cmf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
Alps_cmf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass,
                                   method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Alps_cmf_genetic1_ratio1_gr, un = F, interactions = T)
plot(Alps_cmf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Alps_cmf_genetic1_ratio1_gr, un = F))
plot(summary(Alps_cmf_genetic1_ratio1_gr, interactions = T, un = F))

#using caliper on all covariates
Alps_cmf_genetic1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Alps_cmf_grass,
                                       method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis',
                                       std.caliper = TRUE, caliper = c(Elevation = .8, Roughness = .9, Slope = .9))

summary(Alps_cmf_genetic1_ratio1_cal1_gr, un = F, interactions = T) #258 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Alps_cmf_genetic1_ratio1_cal1_gr, interactions = T, un = F))

#--check matching performance

#no comparison, since genetic is the only method that provided good balance and effective K:1 ratio
Alps_gr_perf <- check_match_performance(mobj = Alps_cmf_genetic1_ratio1_cal1_gr)


# ------ 2. Baltic_mixed_forests - GENETIC -------

Baltic_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Baltic_mixed_forests', grass_col_to_use]

names(which.min(table(Baltic_mf_grass$Period))) #period2

Baltic_mf_grass$Period_bin <- ifelse(Baltic_mf_grass$Period == 'period2', 1, 0)

table(Baltic_mf_grass$Period_bin)

#------check initial imbalance

Baltic_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = NULL, distance = 'glm')

summary(Baltic_mf_init_gr)
plot(summary(Baltic_mf_init_gr))
plot(Baltic_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Baltic_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest', distance = 'glm')

Baltic_mf_pscore1_gr
summary(Baltic_mf_pscore1_gr, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Baltic_mf_pscore1_gr, un = FALSE, interactions = T)
plot(Baltic_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(Baltic_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_pscore1_gr))
plot(summary(Baltic_mf_pscore1_gr, un = F))
plot(summary(Baltic_mf_pscore1_gr, interactions = T))
plot(summary(Baltic_mf_pscore1_gr, interactions = T, un = F))

#order
Baltic_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                    distance = 'glm', m.order = 'closest')

summary(Baltic_mf_pscore1_ord1_gr, un = F)
plot(Baltic_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Baltic_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_pscore1_ord1_gr, un = F))
plot(summary(Baltic_mf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Baltic_mf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                           distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Baltic_mf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(Baltic_mf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(Baltic_mf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(Baltic_mf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#using a caliper for Elevation - w/o caliper (order + ratio) all covariates (and their functions) showing bad balance are related to Elevation
#(except Roughness, which is slightly beyond .10)
Baltic_mf_pscore1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                                 distance = 'glm', m.order = 'closest', ratio = 2,
                                                 std.caliper = TRUE, caliper = c(Elevation = 1.8))

summary(Baltic_mf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T) #0 tr observations are unmatched - but Matched (ESS) < Matched
plot(summary(Baltic_mf_pscore1_ord1_ratio1_cal1_gr, interactions = T, un = F))

#--Mahalanobis

Baltic_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Baltic_mf_mahala1_gr)
plot(Baltic_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_mahala1_gr, un = F))
plot(summary(Baltic_mf_mahala1_gr, interactions = T, un = F))

#order
Baltic_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest')

summary(Baltic_mf_mahala1_ord1_gr, un = F)
plot(Baltic_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_mahala1_ord1_gr, un = F))
plot(summary(Baltic_mf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Baltic_mf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                           distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Baltic_mf_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Baltic_mf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Baltic_mf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using a caliper for Elevation
Baltic_mf_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                                 distance = 'mahalanobis', m.order = 'closest', ratio = 2,
                                                 std.caliper = TRUE, caliper = c(Elevation = 1.9))

summary(Baltic_mf_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #0 tr observations are unmatched - but Matched (ESS) < Matched
plot(summary(Baltic_mf_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = F))

#--robust Mahalanobis

Baltic_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(Baltic_mf_rob_mahala1_gr)
plot(Baltic_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_rob_mahala1_gr, un = F))
plot(summary(Baltic_mf_rob_mahala1_gr, interactions = T, un = F))

#order
Baltic_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest')

summary(Baltic_mf_rob_mahala1_ord1_gr, un = F)
plot(Baltic_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(Baltic_mf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
Baltic_mf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                               distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Baltic_mf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Baltic_mf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Baltic_mf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using a caliper for Elevation
Baltic_mf_rob_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2,
                                                std.caliper = TRUE, caliper = c(Elevation = 1.6))

summary(Baltic_mf_rob_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #2 tr observations are unmatched - but Matched (ESS) < Matched
plot(summary(Baltic_mf_rob_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
Baltic_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Baltic_mf_genetic1_gr)
plot(Baltic_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_genetic1_gr, un = F))
plot(summary(Baltic_mf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
Baltic_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Baltic_mf_genetic2_gr)
plot(Baltic_mf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_genetic2_gr, un = F))
plot(summary(Baltic_mf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
Baltic_mf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass,
                                       method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Baltic_mf_genetic1_ratio1_gr, un = F, interactions = T)
plot(Baltic_mf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Baltic_mf_genetic1_ratio1_gr, un = F))
plot(summary(Baltic_mf_genetic1_ratio1_gr, interactions = T, un = F))

#using a caliper for Elevation
Baltic_mf_genetic1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Baltic_mf_grass,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis',
                                        std.caliper = TRUE, caliper = c(Elevation = 1.5))

summary(Baltic_mf_genetic1_ratio1_cal1_gr, un = F, interactions = T) #30 tr observations are unmatched - no diff between Matched (ESS) and Matched
plot(summary(Baltic_mf_genetic1_ratio1_cal1_gr, interactions = T, un = F))

#--check matching performance

#same as for Alps_cmf
Baltic_gr_perf <- check_match_performance(mobj = Baltic_mf_genetic1_ratio1_cal1_gr)

# ------ 3. Cantabrian_mixed_forests - GENETIC ------

Cant_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Cantabrian_mixed_forests', grass_col_to_use]

names(which.min(table(Cant_mf_grass$Period))) #period2

Cant_mf_grass$Period_bin <- ifelse(Cant_mf_grass$Period == 'period2', 1, 0)

table(Cant_mf_grass$Period_bin)

#------check initial imbalance

Cant_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = NULL, distance = 'glm')

summary(Cant_mf_init_gr)
plot(summary(Cant_mf_init_gr))
plot(Cant_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Cant_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest', distance = 'glm')

summary(Cant_mf_pscore1_gr)
plot(Cant_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1_gr))
plot(summary(Cant_mf_pscore1_gr, un = F))
plot(summary(Cant_mf_pscore1_gr, interactions = T))
plot(summary(Cant_mf_pscore1_gr, interactions = T, un = F))

#order
Cant_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                     distance = 'glm', m.order = 'closest')

summary(Cant_mf_pscore1_ord1_gr, un = F)
plot(Cant_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1_ord1_gr, un = F))
plot(summary(Cant_mf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Cant_mf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                            distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Cant_mf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(Cant_mf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(Cant_mf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(Cant_mf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#--Mahalanobis

Cant_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Cant_mf_mahala1_gr)
plot(Cant_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1_gr, un = F))
plot(summary(Cant_mf_mahala1_gr, interactions = T, un = F))

#order
Cant_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest')

summary(Cant_mf_mahala1_ord1_gr, un = F)
plot(Cant_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1_ord1_gr, un = F))
plot(summary(Cant_mf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Cant_mf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Cant_mf_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Cant_mf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Cant_mf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#--robust Mahalanobis

Cant_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(Cant_mf_rob_mahala1_gr)
plot(Cant_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1_gr, un = F))
plot(summary(Cant_mf_rob_mahala1_gr, interactions = T, un = F))

#order
Cant_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                         distance = 'robust_mahalanobis', m.order = 'closest')

summary(Cant_mf_rob_mahala1_ord1_gr, un = F)
plot(Cant_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(Cant_mf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
Cant_mf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Cant_mf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Cant_mf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Cant_mf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
Cant_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass,
                               method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Cant_mf_genetic1_gr)
plot(Cant_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic1_gr, un = F))
plot(summary(Cant_mf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
Cant_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Cant_mf_genetic2_gr)
plot(Cant_mf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic2_gr, un = F))
plot(summary(Cant_mf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
Cant_mf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Cant_mf_grass,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Cant_mf_genetic1_ratio1_gr, un = F, interactions = T)
plot(Cant_mf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Cant_mf_genetic1_ratio1_gr, un = F))
plot(summary(Cant_mf_genetic1_ratio1_gr, interactions = T, un = F))

#--check matching performance

#Genetic provided good balance, without the need of caliper
Cant_gr_perf <- check_match_performance(Cant_mf_genetic1_ratio1_gr)

# ------ 4. Carpathian_montane_forests - PScore ------

#balanced sample size between periods - don't test ratio = 2

Carp_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Carpathian_montane_forests', grass_col_to_use]

names(which.min(table(Carp_mf_grass$Period))) #period1

Carp_mf_grass$Period_bin <- ifelse(Carp_mf_grass$Period == 'period1', 1, 0)

table(Carp_mf_grass$Period_bin)

#------check initial imbalance

Carp_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = NULL, distance = 'glm')

summary(Carp_mf_init_gr)
plot(summary(Carp_mf_init_gr))
plot(Carp_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Carp_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest', distance = 'glm')

summary(Carp_mf_pscore1_gr)
plot(Carp_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(Carp_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_pscore1_gr))
plot(summary(Carp_mf_pscore1_gr, un = F))
plot(summary(Carp_mf_pscore1_gr, interactions = T))
plot(summary(Carp_mf_pscore1_gr, interactions = T, un = F))

#order
Carp_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                   distance = 'glm', m.order = 'closest')

summary(Carp_mf_pscore1_ord1_gr, un = F, interactions = T)
plot(Carp_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Carp_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_pscore1_ord1_gr, un = F))
plot(summary(Carp_mf_pscore1_ord1_gr, interactions = T, un = F))

#using caliper on all covariates (+ prop score)
Carp_mf_pscore1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                   distance = 'glm', m.order = 'closest',
                                   std.caliper = TRUE, caliper = c(1.8, Elevation = 1.9, Roughness = 1.8, Slope = 1.4),
                                   link = "linear.logit")

summary(Carp_mf_pscore1_ord1_cal1_gr, un = F, interactions = T) #624 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Carp_mf_pscore1_ord1_cal1_gr, interactions = T, un = F))

#--Mahalanobis

Carp_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Carp_mf_mahala1_gr)
plot(Carp_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_mahala1_gr, un = F))
plot(summary(Carp_mf_mahala1_gr, interactions = T, un = F))

#order
Carp_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                   distance = 'mahalanobis', m.order = 'closest')

summary(Carp_mf_mahala1_ord1_gr, un = F, interactions = T)
plot(Carp_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_mahala1_ord1_gr, un = F))
plot(summary(Carp_mf_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on all covariates
Carp_mf_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                   distance = 'mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                   caliper = c(Elevation = 1.6, Roughness = 1.7, Slope = 1.6))

summary(Carp_mf_mahala1_ord1_cal1_gr, un = F, interactions = T) #630 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Carp_mf_mahala1_ord1_cal1_gr, interactions = T, un = F))

#--robust Mahalanobis

Carp_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass,
                                  method = 'nearest', distance = 'robust_mahalanobis')

summary(Carp_mf_rob_mahala1_gr)
plot(Carp_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_rob_mahala1_gr, un = F))
plot(summary(Carp_mf_rob_mahala1_gr, interactions = T, un = F))

#order
Carp_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                       distance = 'robust_mahalanobis', m.order = 'closest')

summary(Carp_mf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(Carp_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(Carp_mf_rob_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on all covariates
Carp_mf_rob_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'nearest',
                                       distance = 'robust_mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                       caliper = c(Elevation = 1.5, Roughness = 1.5, Slope = 1.5))

summary(Carp_mf_rob_mahala1_ord1_cal1_gr, un = F, interactions = T) #605 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Carp_mf_rob_mahala1_ord1_cal1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
Carp_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass,
                               method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Carp_mf_genetic1_gr, un = F, interactions = T)
plot(Carp_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_genetic1_gr, un = F))
plot(summary(Carp_mf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
Carp_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Carp_mf_genetic2_gr)
plot(Carp_mf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Carp_mf_genetic2_gr, un = F))
plot(summary(Carp_mf_genetic2_gr, interactions = T, un = F))

#using caliper on all covariates
Carp_mf_genetic1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Carp_mf_grass,
                               method = 'genetic', pop.size = 100, distance = 'mahalanobis', std.caliper = TRUE,
                               caliper = c(Elevation = 1.1, Roughness = 1.2, Slope = 1.1))

summary(Carp_mf_genetic1_cal1_gr, un = F, interactions = T) #700 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Carp_mf_genetic1_cal1_gr, interactions = T, un = F))

#--check matching performance

#compare all methods (with caliper)
Carp_gr_perf <- do.call(rbind, lapply(list(PScore = Carp_mf_pscore1_ord1_cal1_gr, Mah = Carp_mf_mahala1_ord1_cal1_gr,
                           RMah = Carp_mf_rob_mahala1_ord1_cal1_gr, Genetic = Carp_mf_genetic1_cal1_gr),
                      check_match_performance))

Carp_gr_perf #Pscore provides a good balance and number of matched units

# ------ 5. Celtic_broadleaf_forests - Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

Celtic_bf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Celtic_broadleaf_forests', grass_col_to_use]

names(which.min(table(Celtic_bf_grass$Period))) #period1

Celtic_bf_grass$Period_bin <- ifelse(Celtic_bf_grass$Period == 'period1', 1, 0)

table(Celtic_bf_grass$Period_bin)

#------check initial imbalance

Celtic_bf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = NULL, distance = 'glm')

summary(Celtic_bf_init_gr)
plot(summary(Celtic_bf_init_gr))
plot(summary(Celtic_bf_init_gr, interactions = T))
plot(Celtic_bf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Celtic_bf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest', distance = 'glm')

summary(Celtic_bf_pscore1_gr)
plot(Celtic_bf_pscore1_gr, type = 'jitter', interactive = F)
plot(Celtic_bf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_pscore1_gr))
plot(summary(Celtic_bf_pscore1_gr, un = F))
plot(summary(Celtic_bf_pscore1_gr, interactions = T))
plot(summary(Celtic_bf_pscore1_gr, interactions = T, un = F))

#order
Celtic_bf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                   distance = 'glm', m.order = 'closest')

summary(Celtic_bf_pscore1_ord1_gr, un = F, interactions = T)
plot(Celtic_bf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Celtic_bf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_pscore1_ord1_gr, un = F))
plot(summary(Celtic_bf_pscore1_ord1_gr, interactions = T, un = F))

#using a caliper on Elevation for SMD and on Slope for var ratio
Celtic_bf_pscore1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                     distance = 'glm', m.order = 'closest', std.caliper = T,
                                     caliper = c(Elevation = 1.5, Slope = 2))

summary(Celtic_bf_pscore1_ord1_cal1_gr, un = F, interactions = T) #157 tr obs are not matched due to the caliper (besides control obs)
plot(summary(Celtic_bf_pscore1_ord1_cal1_gr, interactions = T, un = F))


#--Mahalanobis

Celtic_bf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Celtic_bf_mahala1_gr)
plot(Celtic_bf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_mahala1_gr, un = F))
plot(summary(Celtic_bf_mahala1_gr, interactions = T, un = F))

#order
Celtic_bf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                   distance = 'mahalanobis', m.order = 'closest')

summary(Celtic_bf_mahala1_ord1_gr, un = F, interactions = T)
plot(Celtic_bf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_mahala1_ord1_gr, un = F))
plot(summary(Celtic_bf_mahala1_ord1_gr, interactions = T, un = F))

#using a caliper for Elevation - same as for prop score with caliper; 1.6 is not enough to ameliorate balance for Elevation, while 1.5 does the job
Celtic_bf_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest', std.caliper = TRUE, caliper = c(Elevation = 1.5))

summary(Celtic_bf_mahala1_ord1_cal1_gr, un = F, interactions = T) #84 tr observations are unmatched
plot(summary(Celtic_bf_mahala1_ord1_cal1_gr, interactions = T, un = F))


#--robust Mahalanobis

Celtic_bf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass,
                                  method = 'nearest', distance = 'robust_mahalanobis')

summary(Celtic_bf_rob_mahala1_gr)
plot(Celtic_bf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_rob_mahala1_gr, un = F))
plot(summary(Celtic_bf_rob_mahala1_gr, interactions = T, un = F))

#order
Celtic_bf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                       distance = 'robust_mahalanobis', m.order = 'closest')

summary(Celtic_bf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(Celtic_bf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_rob_mahala1_ord1_gr, un = F))
plot(summary(Celtic_bf_rob_mahala1_ord1_gr, interactions = T, un = F))

#using a caliper for Elevation - interestingly, caliper = 1.5 (not even 1.4, actually) is not enough to ameliorate balance of Elevation^2
Celtic_bf_rob_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'nearest',
                                              distance = 'robust_mahalanobis', m.order = 'closest',
                                              std.caliper = TRUE, caliper = c(Elevation = 1.3))

summary(Celtic_bf_rob_mahala1_ord1_cal1_gr, un = F, interactions = T) #112 tr observations are unmatched
plot(summary(Celtic_bf_rob_mahala1_ord1_cal1_gr, interactions = T, un = F))


#------genetic matching

#using GMD without prop.score
Celtic_bf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Celtic_bf_genetic1_gr, un = F, interactions = T)
plot(Celtic_bf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_genetic1_gr, un = F))
plot(summary(Celtic_bf_genetic1_gr, interactions = T, un = F))

#using a caliper for Elevation
Celtic_bf_genetic1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis',
                                 std.caliper = TRUE, caliper = c(Elevation = 1.1))

summary(Celtic_bf_genetic1_cal1_gr, un = FALSE, interactions = T) #243 unmatched tr obs
plot(summary(Celtic_bf_genetic1_cal1_gr, interactions = T, un = F))


#using GMD with prop.score
Celtic_bf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Celtic_bf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Celtic_bf_genetic2_gr)
plot(Celtic_bf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Celtic_bf_genetic2_gr, un = F))
plot(summary(Celtic_bf_genetic2_gr, interactions = T, un = F))

#--check matching performance

#compare all methods (with caliper)
Celtic_gr_perf <- do.call(rbind, lapply(list(PScore = Celtic_bf_pscore1_ord1_cal1_gr,
                                             Mah = Celtic_bf_mahala1_ord1_cal1_gr,
                                             RMah = Celtic_bf_rob_mahala1_ord1_cal1_gr,
                                             Genetic = Celtic_bf_genetic1_cal1_gr), check_match_performance))

Celtic_gr_perf #Mah and RMah are rather equivalent if not considering Gwd - Mah retains more observations


# ------ 6. Central_European_mixed_forests - GENETIC ------

#balanced sample size between periods - don't test ratio = 2

CenEu_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Central_European_mixed_forests', grass_col_to_use]

names(which.min(table(CenEu_mf_grass$Period))) #period2

CenEu_mf_grass$Period_bin <- ifelse(CenEu_mf_grass$Period == 'period2', 1, 0)

table(CenEu_mf_grass$Period_bin)

#------check initial imbalance

CenEu_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = NULL, distance = 'glm')

summary(CenEu_mf_init_gr)
plot(summary(CenEu_mf_init_gr))
plot(summary(CenEu_mf_init_gr, interactions = T))
plot(CenEu_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
CenEu_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'nearest', distance = 'glm')

summary(CenEu_mf_pscore1_gr)
plot(CenEu_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(CenEu_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_pscore1_gr))
plot(summary(CenEu_mf_pscore1_gr, un = F))
plot(summary(CenEu_mf_pscore1_gr, interactions = T))
plot(summary(CenEu_mf_pscore1_gr, interactions = T, un = F))

#order
CenEu_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'nearest',
                                     distance = 'glm', m.order = 'closest')

summary(CenEu_mf_pscore1_ord1_gr, un = F, interactions = T)
plot(CenEu_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(CenEu_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_pscore1_ord1_gr, un = F))
plot(summary(CenEu_mf_pscore1_ord1_gr, interactions = T, un = F))

#--Mahalanobis

CenEu_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(CenEu_mf_mahala1_gr, un = F, interactions = T)
plot(CenEu_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_mahala1_gr, un = F))
plot(summary(CenEu_mf_mahala1_gr, interactions = T, un = F))

#order
CenEu_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest')

summary(CenEu_mf_mahala1_ord1_gr, un = F, interactions = T)
plot(CenEu_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_mahala1_ord1_gr, un = F))
plot(summary(CenEu_mf_mahala1_ord1_gr, interactions = T, un = F))

#--robust Mahalanobis

CenEu_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(CenEu_mf_rob_mahala1_gr)
plot(CenEu_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_rob_mahala1_gr, un = F))
plot(summary(CenEu_mf_rob_mahala1_gr, interactions = T, un = F))

#order
CenEu_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'nearest',
                                         distance = 'robust_mahalanobis', m.order = 'closest')

summary(CenEu_mf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(CenEu_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(CenEu_mf_rob_mahala1_ord1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
CenEu_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(CenEu_mf_genetic1_gr, un = F, interactions = T)
plot(CenEu_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(CenEu_mf_genetic1_gr, un = F))
plot(summary(CenEu_mf_genetic1_gr, interactions = T, un = F))

#I'm not running the one below as 1) it takes quite a lot to run; 2) matching of previous methods is already very good
#using GMD with prop.score
#CenEu_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = CenEu_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

#summary(CenEu_mf_genetic2_gr)
#plot(CenEu_mf_genetic2_gr, type = 'density', interactive = F)
#plot(summary(CenEu_mf_genetic2_gr, un = F))
#plot(summary(CenEu_mf_genetic2_gr, interactions = T, un = F))

#--check matching performance

#compare all methods
CenEu_gr_perf <- do.call(rbind, lapply(list(PScore = CenEu_mf_pscore1_ord1_gr,
                                            Mah = CenEu_mf_mahala1_ord1_gr,
                                            RMah = CenEu_mf_rob_mahala1_ord1_gr,
                                            Genetic = CenEu_mf_genetic1_gr),
                                       check_match_performance))
CenEu_gr_perf #Genetic


# ------ 7. English_Lowlands_beech_forests - GENETIC ------

EngLow_bf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'English_Lowlands_beech_forests', grass_col_to_use]

names(which.min(table(EngLow_bf_grass$Period))) #period2

EngLow_bf_grass$Period_bin <- ifelse(EngLow_bf_grass$Period == 'period2', 1, 0)

table(EngLow_bf_grass$Period_bin)

#------check initial imbalance

EngLow_bf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = NULL, distance = 'glm')

summary(EngLow_bf_init_gr)
plot(summary(EngLow_bf_init_gr))
plot(EngLow_bf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
EngLow_bf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest', distance = 'glm')

summary(EngLow_bf_pscore1_gr) 
summary(EngLow_bf_pscore1_gr, un = FALSE, interactions = T)
plot(EngLow_bf_pscore1_gr, type = 'jitter', interactive = F)
plot(EngLow_bf_pscore1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_pscore1_gr))
plot(summary(EngLow_bf_pscore1_gr, un = F))
plot(summary(EngLow_bf_pscore1_gr, interactions = T))
plot(summary(EngLow_bf_pscore1_gr, interactions = T, un = F))

#order
EngLow_bf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                     distance = 'glm', m.order = 'closest')

summary(EngLow_bf_pscore1_ord1_gr, un = F)
plot(EngLow_bf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(EngLow_bf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_pscore1_ord1_gr, un = F))
plot(summary(EngLow_bf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EngLow_bf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                            distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(EngLow_bf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(EngLow_bf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(EngLow_bf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(EngLow_bf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on prop score (for SMD and var ratio), and on Roughness and Slope (for var ratio)
EngLow_bf_pscore1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                            distance = 'glm', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                            caliper = c(2, Roughness = 2, Slope = 2),
                                            link = 'linear.logit')

summary(EngLow_bf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T) #65 tr obs unmatched - ESS < SS
plot(summary(EngLow_bf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T))

#--Mahalanobis

EngLow_bf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest', distance = 'mahalanobis')

summary(EngLow_bf_mahala1_gr)
plot(EngLow_bf_mahala1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_mahala1_gr, un = F))
plot(summary(EngLow_bf_mahala1_gr, interactions = T, un = F))

#order
EngLow_bf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest')

summary(EngLow_bf_mahala1_ord1_gr, un = F)
plot(EngLow_bf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_mahala1_ord1_gr, un = F))
plot(summary(EngLow_bf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EngLow_bf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EngLow_bf_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(EngLow_bf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(EngLow_bf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on Roughness (for SMD) and on Slope (for var ratio)
EngLow_bf_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                            caliper = c(Roughness = 2, Slope = 2))

summary(EngLow_bf_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = FALSE) #4 tr obs unmatched - ESS < SS
plot(summary(EngLow_bf_mahala1_ord1_ratio1_cal1_gr, interactions = T, un = FALSE))

#--robust Mahalanobis

EngLow_bf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(EngLow_bf_rob_mahala1_gr)
plot(EngLow_bf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_rob_mahala1_gr, un = F))
plot(summary(EngLow_bf_rob_mahala1_gr, interactions = T, un = F))

#order
EngLow_bf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                         distance = 'robust_mahalanobis', m.order = 'closest')

summary(EngLow_bf_rob_mahala1_ord1_gr, un = F)
plot(EngLow_bf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_rob_mahala1_ord1_gr, un = F))
plot(summary(EngLow_bf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
EngLow_bf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EngLow_bf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(EngLow_bf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(EngLow_bf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on Elevation (for SMD) and Roughness and Slope (for var ratio)
EngLow_bf_rob_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                                caliper = c(Elevation = 1.5, Roughness = 2, Slope = 2))

summary(EngLow_bf_rob_mahala1_ord1_ratio1_cal1_gr, interactions = TRUE, un = FALSE) #6 tr obs unmatched - ESS < SS
plot(summary(EngLow_bf_rob_mahala1_ord1_ratio1_cal1_gr, interactions = TRUE, un = FALSE))

#------genetic matching

#using GMD without prop.score
EngLow_bf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(EngLow_bf_genetic1_gr)
plot(EngLow_bf_genetic1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_genetic1_gr, un = F))
plot(summary(EngLow_bf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
EngLow_bf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(EngLow_bf_genetic2_gr)
plot(EngLow_bf_genetic2_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_genetic2_gr, un = F))
plot(summary(EngLow_bf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
EngLow_bf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(EngLow_bf_genetic1_ratio1_gr, un = F, interactions = TRUE)
plot(EngLow_bf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EngLow_bf_genetic1_ratio1_gr, un = F))
plot(summary(EngLow_bf_genetic1_ratio1_gr, interactions = T, un = F))

#using caliper on Roughness and Slope for var ratio
EngLow_bf_genetic1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EngLow_bf_grass,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis',
                                        std.caliper = TRUE, caliper = c(Roughness = 2, Slope = 2))

summary(EngLow_bf_genetic1_ratio1_cal1_gr, interactions = T, un = F) #17 tr obs unmatched
plot(summary(EngLow_bf_genetic1_ratio1_cal1_gr, interactions = T, un = F))


#--check matching performance

#Genetic is the only method that provided good balance and effective K:1 ratio
EngLow_gr_perf <- check_match_performance(EngLow_bf_genetic1_ratio1_cal1_gr) #Genetic

# ------ 8. European_Atlantic_mixed_forests - GENETIC ------

EuAtl_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'European_Atlantic_mixed_forests', grass_col_to_use]

names(which.min(table(EuAtl_mf_grass$Period))) #period1

EuAtl_mf_grass$Period_bin <- ifelse(EuAtl_mf_grass$Period == 'period1', 1, 0)

table(EuAtl_mf_grass$Period_bin)

#------check initial imbalance

EuAtl_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = NULL, distance = 'glm')

summary(EuAtl_mf_init_gr)
plot(summary(EuAtl_mf_init_gr))
plot(EuAtl_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
EuAtl_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest', distance = 'glm')

summary(EuAtl_mf_pscore1_gr) 
summary(EuAtl_mf_pscore1_gr, un = FALSE, interactions = T)
plot(EuAtl_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(EuAtl_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_pscore1_gr))
plot(summary(EuAtl_mf_pscore1_gr, un = F))
plot(summary(EuAtl_mf_pscore1_gr, interactions = T))
plot(summary(EuAtl_mf_pscore1_gr, interactions = T, un = F))

#order
EuAtl_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                     distance = 'glm', m.order = 'closest')

summary(EuAtl_mf_pscore1_ord1_gr, un = F)
plot(EuAtl_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(EuAtl_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_pscore1_ord1_gr, un = F))
plot(summary(EuAtl_mf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                            distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(EuAtl_mf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(EuAtl_mf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(EuAtl_mf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
EuAtl_mf_pscore1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                           distance = 'glm', m.order = 'closest', ratio = 2,
                                           std.caliper = T, caliper = c(Elevation = 2))

summary(EuAtl_mf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T) #24 tr obs unmatched
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_cal1_gr, un = F, interactions = T))


#--Mahalanobis

EuAtl_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(EuAtl_mf_mahala1_gr)
plot(EuAtl_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_mahala1_gr, un = F))
plot(summary(EuAtl_mf_mahala1_gr, interactions = T, un = F))

#order
EuAtl_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest')

summary(EuAtl_mf_mahala1_ord1_gr, un = F)
plot(EuAtl_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_mahala1_ord1_gr, un = F))
plot(summary(EuAtl_mf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EuAtl_mf_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(EuAtl_mf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
EuAtl_mf_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                           distance = 'mahalanobis', m.order = 'closest', ratio = 2,
                                           std.caliper = T, caliper = c(Elevation = 2))

summary(EuAtl_mf_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #24 tr obs unmatched
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T))


#--robust Mahalanobis

EuAtl_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(EuAtl_mf_rob_mahala1_gr)
plot(EuAtl_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_rob_mahala1_gr, un = F))
plot(summary(EuAtl_mf_rob_mahala1_gr, interactions = T, un = F))

#order
EuAtl_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                         distance = 'robust_mahalanobis', m.order = 'closest')

summary(EuAtl_mf_rob_mahala1_ord1_gr, un = F)
plot(EuAtl_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(EuAtl_mf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EuAtl_mf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(EuAtl_mf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'nearest',
                                               distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2,
                                               std.caliper = T, caliper = c(Elevation = 2))

summary(EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T) #24 tr obs unmatched
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_gr, un = F, interactions = T))


#------genetic matching

#using GMD without prop.score
EuAtl_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(EuAtl_mf_genetic1_gr)
plot(EuAtl_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_genetic1_gr, un = F))
plot(summary(EuAtl_mf_genetic1_gr, interactions = T, un = F))

#no need to run this below, given quality of matching achieved without prop score
#using GMD with prop.score
#EuAtl_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

#summary(EuAtl_mf_genetic2_gr)
#plot(EuAtl_mf_genetic2_gr, type = 'density', interactive = F)
#plot(summary(EuAtl_mf_genetic2_gr, un = F))
#plot(summary(EuAtl_mf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
EuAtl_mf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(EuAtl_mf_genetic1_ratio1_gr, un = F, interactions = T)
plot(EuAtl_mf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(EuAtl_mf_genetic1_ratio1_gr, un = F))
plot(summary(EuAtl_mf_genetic1_ratio1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
EuAtl_mf_genetic1_ratio1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = EuAtl_mf_grass,
                                       method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis',
                                       std.caliper = TRUE, caliper = c(Elevation = 2))


summary(EuAtl_mf_genetic1_ratio1_cal1_gr, un = F, interactions = T) #24 tr obs unmatched
plot(summary(EuAtl_mf_genetic1_ratio1_cal1_gr, un = F, interactions = T))

#--check matching performance

#compare all methods
EuAtl_gr_perf <- do.call(rbind, lapply(list(PScore = EuAtl_mf_pscore1_ord1_ratio1_cal1_gr,
                                            Mah = EuAtl_mf_mahala1_ord1_ratio1_cal1_gr,
                                            RMah = EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_gr,
                                            Genetic = EuAtl_mf_genetic1_ratio1_cal1_gr),
                                       check_match_performance))

EuAtl_gr_perf #Genetic


# ------ 9. Italian_sclerophyllous_and_semi_deciduous_forests - Robust Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

ItaScl_sdf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Italian_sclerophyllous_and_semi_deciduous_forests', grass_col_to_use]

names(which.min(table(ItaScl_sdf_grass$Period))) #period2

ItaScl_sdf_grass$Period_bin <- ifelse(ItaScl_sdf_grass$Period == 'period2', 1, 0)

table(ItaScl_sdf_grass$Period_bin)

#------check initial imbalance

ItaScl_sdf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = NULL, distance = 'glm')

summary(ItaScl_sdf_init_gr)
plot(summary(ItaScl_sdf_init_gr))
plot(summary(ItaScl_sdf_init_gr, interactions = T))
plot(ItaScl_sdf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
ItaScl_sdf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'nearest', distance = 'glm')

summary(ItaScl_sdf_pscore1_gr)
plot(ItaScl_sdf_pscore1_gr, type = 'jitter', interactive = F)
plot(ItaScl_sdf_pscore1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_pscore1_gr))
plot(summary(ItaScl_sdf_pscore1_gr, un = F))
plot(summary(ItaScl_sdf_pscore1_gr, interactions = T))
plot(summary(ItaScl_sdf_pscore1_gr, interactions = T, un = F))

#order
ItaScl_sdf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'nearest',
                                    distance = 'glm', m.order = 'closest')

summary(ItaScl_sdf_pscore1_ord1_gr, un = F, interactions = T)
plot(ItaScl_sdf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(ItaScl_sdf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_pscore1_ord1_gr, un = F))
plot(summary(ItaScl_sdf_pscore1_ord1_gr, interactions = T, un = F))

#--Mahalanobis

ItaScl_sdf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'nearest', distance = 'mahalanobis')

summary(ItaScl_sdf_mahala1_gr)
plot(ItaScl_sdf_mahala1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_mahala1_gr, un = F))
plot(summary(ItaScl_sdf_mahala1_gr, interactions = T, un = F))

#order
ItaScl_sdf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest')

summary(ItaScl_sdf_mahala1_ord1_gr, un = F, interactions = T)
plot(ItaScl_sdf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_mahala1_ord1_gr, un = F))
plot(summary(ItaScl_sdf_mahala1_ord1_gr, interactions = T, un = F))

#--robust Mahalanobis

ItaScl_sdf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(ItaScl_sdf_rob_mahala1_gr)
plot(ItaScl_sdf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_rob_mahala1_gr, un = F))
plot(summary(ItaScl_sdf_rob_mahala1_gr, interactions = T, un = F))

#order
ItaScl_sdf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest')

summary(ItaScl_sdf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(ItaScl_sdf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_rob_mahala1_ord1_gr, un = F))
plot(summary(ItaScl_sdf_rob_mahala1_ord1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
ItaScl_sdf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(ItaScl_sdf_genetic1_gr, un = F, interactions = T)
plot(ItaScl_sdf_genetic1_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_genetic1_gr, un = F))
plot(summary(ItaScl_sdf_genetic1_gr, interactions = T, un = F))

#using GMD with prop.score
ItaScl_sdf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = ItaScl_sdf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(ItaScl_sdf_genetic2_gr)
plot(ItaScl_sdf_genetic2_gr, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_genetic2_gr, un = F))
plot(summary(ItaScl_sdf_genetic2_gr, interactions = T, un = F))

#--check matching performance

#comparison is between Mahalanobis and Robust Mahalanobis
ItaScl_gr_perf <- do.call(rbind, lapply(list(Mah = ItaScl_sdf_mahala1_ord1_gr, RMah = ItaScl_sdf_rob_mahala1_ord1_gr),
                                        check_match_performance))

ItaScl_gr_perf #RMah


# ------ 10. North_Atlantic_moist_mixed_forests - GENETIC ------

#balanced sample size between periods - don't test ratio = 2

NAtl_mmf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'North_Atlantic_moist_mixed_forests', grass_col_to_use]

names(which.min(table(NAtl_mmf_grass$Period))) #period1

NAtl_mmf_grass$Period_bin <- ifelse(NAtl_mmf_grass$Period == 'period1', 1, 0)

table(NAtl_mmf_grass$Period_bin)

#------check initial imbalance

NAtl_mmf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = NULL, distance = 'glm')

summary(NAtl_mmf_init_gr)
plot(summary(NAtl_mmf_init_gr))
plot(summary(NAtl_mmf_init_gr, interactions = T))
plot(NAtl_mmf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
NAtl_mmf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest', distance = 'glm')

summary(NAtl_mmf_pscore1_gr)
plot(NAtl_mmf_pscore1_gr, type = 'jitter', interactive = F)
plot(NAtl_mmf_pscore1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_pscore1_gr))
plot(summary(NAtl_mmf_pscore1_gr, un = F))
plot(summary(NAtl_mmf_pscore1_gr, interactions = T))
plot(summary(NAtl_mmf_pscore1_gr, interactions = T, un = F))

#order
NAtl_mmf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                      distance = 'glm', m.order = 'closest')

summary(NAtl_mmf_pscore1_ord1_gr, un = F, interactions = T)
plot(NAtl_mmf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(NAtl_mmf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_pscore1_ord1_gr, un = F))
plot(summary(NAtl_mmf_pscore1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation (for var ratio) and on prop score (for SMD) - Roughness and Slope were included as a caliper of .9 made the unbalanced
NAtl_mmf_pscore1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                    distance = 'glm', m.order = 'closest', std.caliper = TRUE,
                                    caliper = c(1, Elevation = .5, Roughness = 1.9, Slope = 1.9),
                                    link = 'linear.logit')

summary(NAtl_mmf_pscore1_ord1_cal1_gr, un = F, interactions = T) #127 tr obs unmatched
plot(summary(NAtl_mmf_pscore1_ord1_cal1_gr, un = F, interactions = T))


#--Mahalanobis

NAtl_mmf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest', distance = 'mahalanobis')

summary(NAtl_mmf_mahala1_gr)
plot(NAtl_mmf_mahala1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_mahala1_gr, un = F))
plot(summary(NAtl_mmf_mahala1_gr, interactions = T, un = F))

#order
NAtl_mmf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                      distance = 'mahalanobis', m.order = 'closest')

summary(NAtl_mmf_mahala1_ord1_gr, un = F, interactions = T)
plot(NAtl_mmf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_mahala1_ord1_gr, un = F))
plot(summary(NAtl_mmf_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
NAtl_mmf_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest',
                                    std.caliper = TRUE, caliper = c(Elevation = 1.4))

summary(NAtl_mmf_mahala1_ord1_cal1_gr, un = F, interactions = T) #5 tr obs unmatched
plot(summary(NAtl_mmf_mahala1_ord1_cal1_gr, un = F, interactions = T))

#--robust Mahalanobis

NAtl_mmf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass,
                                     method = 'nearest', distance = 'robust_mahalanobis')

summary(NAtl_mmf_rob_mahala1_gr)
plot(NAtl_mmf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_rob_mahala1_gr, un = F))
plot(summary(NAtl_mmf_rob_mahala1_gr, interactions = T, un = F))

#order
NAtl_mmf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                          distance = 'robust_mahalanobis', m.order = 'closest')

summary(NAtl_mmf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(NAtl_mmf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_rob_mahala1_ord1_gr, un = F))
plot(summary(NAtl_mmf_rob_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
NAtl_mmf_rob_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest',
                                        std.caliper = TRUE, caliper = c(Elevation = 2))

summary(NAtl_mmf_rob_mahala1_ord1_cal1_gr, un = F, interactions = T) #4 tr obs unmatched
plot(summary(NAtl_mmf_rob_mahala1_ord1_cal1_gr, un = F, interactions = T))


#------genetic matching

#using GMD without prop.score
NAtl_mmf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass,
                                  method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(NAtl_mmf_genetic1_gr, un = F, interactions = T)
plot(NAtl_mmf_genetic1_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_genetic1_gr, un = F))
plot(summary(NAtl_mmf_genetic1_gr, interactions = T, un = F))

#using caliper on Elevation for var ratio
NAtl_mmf_genetic1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis',
                                std.caliper = TRUE, caliper = c(Elevation = 1.7))

summary(NAtl_mmf_genetic1_cal1_gr, un = F, interactions = T) #4 tr obs unmatched
plot(summary(NAtl_mmf_genetic1_cal1_gr, un = F, interactions = T))


#using GMD with prop.score
NAtl_mmf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = NAtl_mmf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(NAtl_mmf_genetic2_gr)
plot(NAtl_mmf_genetic2_gr, type = 'density', interactive = F)
plot(summary(NAtl_mmf_genetic2_gr, un = F))
plot(summary(NAtl_mmf_genetic2_gr, interactions = T, un = F))

#--check matching performance

#compare performance of Mahalanobis, Robust Mahalanobis and genetic
NAtl_gr_perf <- do.call(rbind, lapply(list(PScore = NAtl_mmf_pscore1_ord1_cal1_gr,
                                           Mah = NAtl_mmf_mahala1_ord1_cal1_gr,
                                           RMah = NAtl_mmf_rob_mahala1_ord1_cal1_gr,
                                           Genetic = NAtl_mmf_genetic1_cal1_gr),
                                      check_match_performance))

NAtl_gr_perf #Genetic (RMah provided similar performance)

# ------ 11. Pannonian_mixed_forests - Robust Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

Pan_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Pannonian_mixed_forests', grass_col_to_use]

names(which.min(table(Pan_mf_grass$Period))) #period1

Pan_mf_grass$Period_bin <- ifelse(Pan_mf_grass$Period == 'period1', 1, 0)

table(Pan_mf_grass$Period_bin)

#------check initial imbalance

Pan_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = NULL, distance = 'glm')

summary(Pan_mf_init_gr)
plot(summary(Pan_mf_init_gr))
plot(summary(Pan_mf_init_gr, interactions = T))
plot(Pan_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Pan_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest', distance = 'glm')

summary(Pan_mf_pscore1_gr)
plot(Pan_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(Pan_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_pscore1_gr))
plot(summary(Pan_mf_pscore1_gr, un = F))
plot(summary(Pan_mf_pscore1_gr, interactions = T))
plot(summary(Pan_mf_pscore1_gr, interactions = T, un = F))

#order
Pan_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                    distance = 'glm', m.order = 'closest')

summary(Pan_mf_pscore1_ord1_gr, un = F, interactions = T)
plot(Pan_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Pan_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_pscore1_ord1_gr, un = F))
plot(summary(Pan_mf_pscore1_ord1_gr, interactions = T, un = F))

#caliper on Elevation and propensity score for SMD
Pan_mf_pscore1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                  distance = 'glm', m.order = 'closest', std.caliper = TRUE,
                                  caliper = c(1.9, Elevation = 1.9), link = "linear.logit")

summary(Pan_mf_pscore1_ord1_cal1_gr, un = F, interactions = T) #336 tr obs unmatched
plot(summary(Pan_mf_pscore1_ord1_cal1_gr, un = F, interactions = TRUE))

#--Mahalanobis

Pan_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Pan_mf_mahala1_gr)
plot(Pan_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_mahala1_gr, un = F))
plot(summary(Pan_mf_mahala1_gr, interactions = T, un = F))

#order
Pan_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest')

summary(Pan_mf_mahala1_ord1_gr, un = F, interactions = T)
plot(Pan_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_mahala1_ord1_gr, un = F))
plot(summary(Pan_mf_mahala1_ord1_gr, interactions = T, un = F))

#using a caliper for Elevation - 2 was too large to ameliorate balance for Elevation, while 1.9 was sufficient
Pan_mf_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                       distance = 'mahalanobis', m.order = 'closest',
                                       std.caliper = TRUE, caliper = c(Elevation = 1.9))

summary(Pan_mf_mahala1_ord1_cal1_gr, un = F, interactions = T) #269 tr observations are unmatched
plot(summary(Pan_mf_mahala1_ord1_cal1_gr, interactions = T, un = F))


#--robust Mahalanobis

Pan_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(Pan_mf_rob_mahala1_gr)
plot(Pan_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_rob_mahala1_gr, un = F))
plot(summary(Pan_mf_rob_mahala1_gr, interactions = T, un = F))

#order
Pan_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest')

summary(Pan_mf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(Pan_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(Pan_mf_rob_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation, Roughness and Slope
Pan_mf_rob_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'nearest',
                                      distance = 'robust_mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                      caliper = c(Elevation = 1.6, Roughness = 2, Slope = 2))

summary(Pan_mf_rob_mahala1_ord1_cal1_gr, un = FALSE, interactions = T)
plot(summary(Pan_mf_rob_mahala1_ord1_cal1_gr, un = FALSE, interactions = TRUE))


#------genetic matching

#using GMD without prop.score
Pan_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Pan_mf_genetic1_gr, un = F, interactions = T)
plot(Pan_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_genetic1_gr, un = F))
plot(summary(Pan_mf_genetic1_gr, interactions = T, un = F))

#using a caliper for Elevation
Pan_mf_genetic1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass,
                              method = 'genetic', pop.size = 100, distance = 'mahalanobis',
                              std.caliper = TRUE, caliper = c(Elevation = 1.1))

summary(Pan_mf_genetic1_cal1_gr, un = F, interactions = T) #376 tr observations are unmatched
plot(summary(Pan_mf_genetic1_cal1_gr, interactions = T, un = F))

#using GMD with prop.score
Pan_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Pan_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

summary(Pan_mf_genetic2_gr)
plot(Pan_mf_genetic2_gr, type = 'density', interactive = F)
plot(summary(Pan_mf_genetic2_gr, un = F))
plot(summary(Pan_mf_genetic2_gr, interactions = T, un = F))

#--check matching performance

#compare all methods (with caliper)
Pan_gr_perf <- do.call(rbind, lapply(list(PScore = Pan_mf_pscore1_ord1_cal1_gr,
                                          Mah = Pan_mf_mahala1_ord1_cal1_gr,
                                          RMah = Pan_mf_rob_mahala1_ord1_cal1_gr,
                                          Genetic = Pan_mf_genetic1_cal1_gr),
                                     check_match_performance))

Pan_gr_perf #RMah provides good balance and a sample size in line with other matching approaches

# ------ 12. Sarmatic_mixed_forests - GENETIC ------

Sar_mf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Sarmatic_mixed_forests', grass_col_to_use]

names(which.min(table(Sar_mf_grass$Period))) #period1

Sar_mf_grass$Period_bin <- ifelse(Sar_mf_grass$Period == 'period1', 1, 0)

table(Sar_mf_grass$Period_bin)

#------check initial imbalance

Sar_mf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = NULL, distance = 'glm')

summary(Sar_mf_init_gr)
plot(summary(Sar_mf_init_gr))
plot(Sar_mf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Sar_mf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest', distance = 'glm')

summary(Sar_mf_pscore1_gr) 
summary(Sar_mf_pscore1_gr, un = FALSE, interactions = T)
plot(Sar_mf_pscore1_gr, type = 'jitter', interactive = F)
plot(Sar_mf_pscore1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_pscore1_gr))
plot(summary(Sar_mf_pscore1_gr, un = F))
plot(summary(Sar_mf_pscore1_gr, interactions = T))
plot(summary(Sar_mf_pscore1_gr, interactions = T, un = F))

#order
Sar_mf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                    distance = 'glm', m.order = 'closest')

summary(Sar_mf_pscore1_ord1_gr, un = F)
plot(Sar_mf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(Sar_mf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_pscore1_ord1_gr, un = F))
plot(summary(Sar_mf_pscore1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Sar_mf_pscore1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                           distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Sar_mf_pscore1_ord1_ratio1_gr, un = F, interactions = T)
plot(Sar_mf_pscore1_ord1_ratio1_gr, type = 'jitter', interactive = F)
plot(Sar_mf_pscore1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_pscore1_ord1_ratio1_gr, un = F))
plot(summary(Sar_mf_pscore1_ord1_ratio1_gr, interactions = T, un = F))

#--Mahalanobis

Sar_mf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest', distance = 'mahalanobis')

summary(Sar_mf_mahala1_gr)
plot(Sar_mf_mahala1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_mahala1_gr, un = F))
plot(summary(Sar_mf_mahala1_gr, interactions = T, un = F))

#order
Sar_mf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest')

summary(Sar_mf_mahala1_ord1_gr, un = F)
plot(Sar_mf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_mahala1_ord1_gr, un = F))
plot(summary(Sar_mf_mahala1_ord1_gr, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Sar_mf_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                           distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Sar_mf_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Sar_mf_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Sar_mf_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#--robust Mahalanobis

Sar_mf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(Sar_mf_rob_mahala1_gr)
plot(Sar_mf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_rob_mahala1_gr, un = F))
plot(summary(Sar_mf_rob_mahala1_gr, interactions = T, un = F))

#order
Sar_mf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest')

summary(Sar_mf_rob_mahala1_ord1_gr, un = F)
plot(Sar_mf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_rob_mahala1_ord1_gr, un = F))
plot(summary(Sar_mf_rob_mahala1_ord1_gr, interactions = T, un = F))


#order + ratio -> this increases precision at the expenses of balance
Sar_mf_rob_mahala1_ord1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'nearest',
                                               distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Sar_mf_rob_mahala1_ord1_ratio1_gr, un = F, interactions = T)
plot(Sar_mf_rob_mahala1_ord1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_rob_mahala1_ord1_ratio1_gr, un = F))
plot(summary(Sar_mf_rob_mahala1_ord1_ratio1_gr, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score
Sar_mf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Sar_mf_genetic1_gr)
plot(Sar_mf_genetic1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_genetic1_gr, un = F))
plot(summary(Sar_mf_genetic1_gr, interactions = T, un = F))

#no need to run this below, given quality of matching achieved without prop score
#using GMD with prop.score
#Sar_mf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

#summary(Sar_mf_genetic2_gr)
#plot(Sar_mf_genetic2_gr, type = 'density', interactive = F)
#plot(summary(Sar_mf_genetic2_gr, un = F))
#plot(summary(Sar_mf_genetic2_gr, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'

#ratio
Sar_mf_genetic1_ratio1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = Sar_mf_grass,
                                       method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Sar_mf_genetic1_ratio1_gr, un = F, interactions = T)
plot(Sar_mf_genetic1_ratio1_gr, type = 'density', interactive = F)
plot(summary(Sar_mf_genetic1_ratio1_gr, un = F))
plot(summary(Sar_mf_genetic1_ratio1_gr, interactions = T, un = F))

#--check matching performance

#compare all methods
Sar_gr_perf <- do.call(rbind, lapply(list(PScore = Sar_mf_pscore1_ord1_ratio1_gr,
                                          Mah = Sar_mf_mahala1_ord1_ratio1_gr,
                                          RMah = Sar_mf_rob_mahala1_ord1_ratio1_gr,
                                          Genetic = Sar_mf_genetic1_ratio1_gr), check_match_performance))

Sar_gr_perf #Genetic


# ------ 13. Western_European_broadleaf_forests - Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

WesEu_bf_grass <- Grass_sel_meta[Grass_sel_meta$ECO_NAME == 'Western_European_broadleaf_forests', grass_col_to_use]

names(which.min(table(WesEu_bf_grass$Period))) #period2

WesEu_bf_grass$Period_bin <- ifelse(WesEu_bf_grass$Period == 'period2', 1, 0)

table(WesEu_bf_grass$Period_bin)

#------check initial imbalance

WesEu_bf_init_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = NULL, distance = 'glm')

summary(WesEu_bf_init_gr)
plot(summary(WesEu_bf_init_gr))
plot(summary(WesEu_bf_init_gr, interactions = T))
plot(WesEu_bf_init_gr, type = 'density')

#------nearest neighbor

#--propensity score

#logit
WesEu_bf_pscore1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest', distance = 'glm')

summary(WesEu_bf_pscore1_gr)
plot(WesEu_bf_pscore1_gr, type = 'jitter', interactive = F)
plot(WesEu_bf_pscore1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_pscore1_gr))
plot(summary(WesEu_bf_pscore1_gr, un = F))
plot(summary(WesEu_bf_pscore1_gr, interactions = T))
plot(summary(WesEu_bf_pscore1_gr, interactions = T, un = F))

#order
WesEu_bf_pscore1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                  distance = 'glm', m.order = 'closest')

summary(WesEu_bf_pscore1_ord1_gr, un = F, interactions = T)
plot(WesEu_bf_pscore1_ord1_gr, type = 'jitter', interactive = F)
plot(WesEu_bf_pscore1_ord1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_pscore1_ord1_gr, un = F))
plot(summary(WesEu_bf_pscore1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation and prop score (for SMD) and on Roughness (for var ratio)
WesEu_bf_pscore1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                         distance = 'glm', m.order = 'closest',
                                         std.caliper = TRUE, caliper = c(1.6, Elevation = 2, Roughness = 2),
                                         link = "linear.logit")

summary(WesEu_bf_pscore1_ord1_cal1_gr, un = F, interactions = T) #1349 tr observations are unmatched
plot(summary(WesEu_bf_pscore1_ord1_cal1_gr, interactions = T, un = F))

#--Mahalanobis

WesEu_bf_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest', distance = 'mahalanobis')

summary(WesEu_bf_mahala1_gr)
plot(WesEu_bf_mahala1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_mahala1_gr, un = F))
plot(summary(WesEu_bf_mahala1_gr, interactions = T, un = F))

#order
WesEu_bf_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                  distance = 'mahalanobis', m.order = 'closest')

summary(WesEu_bf_mahala1_ord1_gr, un = F, interactions = T)
plot(WesEu_bf_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_mahala1_ord1_gr, un = F))
plot(summary(WesEu_bf_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation (for SMD) and on Roughness (for var ratio)
WesEu_bf_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest',
                                    std.caliper = TRUE, caliper = c(Elevation = 1.2, Roughness = 2))

summary(WesEu_bf_mahala1_ord1_cal1_gr, un = F, interactions = T) #1018 tr observations are unmatched
plot(summary(WesEu_bf_mahala1_ord1_cal1_gr, interactions = T, un = F))


#--robust Mahalanobis

WesEu_bf_rob_mahala1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass,
                                 method = 'nearest', distance = 'robust_mahalanobis')

summary(WesEu_bf_rob_mahala1_gr)
plot(WesEu_bf_rob_mahala1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_rob_mahala1_gr, un = F))
plot(summary(WesEu_bf_rob_mahala1_gr, interactions = T, un = F))

#order
WesEu_bf_rob_mahala1_ord1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                      distance = 'robust_mahalanobis', m.order = 'closest')

summary(WesEu_bf_rob_mahala1_ord1_gr, un = F, interactions = T)
plot(WesEu_bf_rob_mahala1_ord1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_rob_mahala1_ord1_gr, un = F))
plot(summary(WesEu_bf_rob_mahala1_ord1_gr, interactions = T, un = F))

#using caliper on Elevation (for SMD) and on Roughness (for var ratio)
WesEu_bf_rob_mahala1_ord1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest',
                                        std.caliper = TRUE, caliper = c(Elevation = 1.2, Roughness = 2))

summary(WesEu_bf_rob_mahala1_ord1_cal1_gr, un = F, interactions = T) #1089 tr observations are unmatched
plot(summary(WesEu_bf_rob_mahala1_ord1_cal1_gr, interactions = T, un = F))


#------genetic matching

#using GMD without prop.score
WesEu_bf_genetic1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass,
                              method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(WesEu_bf_genetic1_gr, un = F, interactions = T)
plot(WesEu_bf_genetic1_gr, type = 'density', interactive = F)
plot(summary(WesEu_bf_genetic1_gr, un = F))
plot(summary(WesEu_bf_genetic1_gr, interactions = T, un = F))

#not running it as it takes super long and results are not going to change much from genetic1
#using GMD with prop.score
#WesEu_bf_genetic2_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass, method = 'genetic', distance = 'glm', pop.size = 100)

#summary(WesEu_bf_genetic2_gr)
#plot(WesEu_bf_genetic2_gr, type = 'density', interactive = F)
#plot(summary(WesEu_bf_genetic2_gr, un = F))
#plot(summary(WesEu_bf_genetic2_gr, interactions = T, un = F))

#using caliper on Elevation (for SMD) and on Roughness (for var ratio)
#Slope and Roughness exhibited SMD > .1 after setting caliper for Elevation and Roughness (Roughness originally for var ratio) at 2
#that is why a caliper for Slope was included
WesEu_bf_genetic1_cal1_gr <- matchit(Period_bin ~ Elevation + Roughness + Slope, data = WesEu_bf_grass,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis',
                                std.caliper = TRUE, caliper = c(Elevation = 1.6, Roughness = .6, Slope = .7))

summary(WesEu_bf_genetic1_cal1_gr, un = F, interactions = T)
plot(summary(WesEu_bf_genetic1_cal1_gr, un = F))
plot(summary(WesEu_bf_genetic1_cal1_gr, un = F, interactions = T))

#--check matching performance

#compare all methods
Wes_gr_perf <- do.call(rbind, lapply(list(PScore = WesEu_bf_pscore1_ord1_cal1_gr,
                                          Mah = WesEu_bf_mahala1_ord1_cal1_gr,
                                          RMah = WesEu_bf_rob_mahala1_ord1_cal1_gr,
                                          Genetic = WesEu_bf_genetic1_cal1_gr), check_match_performance))

Wes_gr_perf #Mah provides good balance performance and large sample size (very similar performance than RMah)




# ----
# ----

# ------ End of stat matching for grasslands ------

#----extract matched datasets and prepare them for GDMs

#Alps_conifer_and_mixed_forests: Alps_cmf_genetic1_ratio1_cal1_gr
#Baltic_mixed_forests: Baltic_mf_genetic1_ratio1_cal1_gr
#Cantabrian_mixed_forests: Cant_mf_genetic1_ratio1_gr
#Carpathian_montane_forests: Carp_mf_pscore1_ord1_cal1_gr
#Celtic_broadleaf_forests: Celtic_bf_mahala1_ord1_cal1_gr
#Central_European_mixed_forests: CenEu_mf_genetic1_gr
#English_Lowlands_beech_forests: EngLow_bf_genetic1_ratio1_cal1_gr
#European_Atlantic_mixed_forests: EuAtl_mf_genetic1_ratio1_cal1_gr
#Italian_sclerophyllous_and_semi_deciduous_forests: ItaScl_sdf_rob_mahala1_ord1_gr
#North_Atlantic_moist_mixed_forests: NAtl_mmf_genetic1_cal1_gr
#Pannonian_mixed_forests: Pan_mf_rob_mahala1_ord1_cal1_gr
#Sarmatic_mixed_forests: Sar_mf_genetic1_ratio1_gr
#Western_European_broadleaf_forests: WesEu_bf_mahala1_ord1_cal1_gr

#FROM HERE! Create Matched_datasets_grass including WesEu_bf_mahala1_ord1_cal1_gr

Matched_datasets_grass <- list(Alps_cmf = Alps_cmf_genetic1_ratio1_cal1_gr,
                               Baltic_mf = Baltic_mf_genetic1_ratio1_cal1_gr,
                               Cant_mf = Cant_mf_genetic1_ratio1_gr,
                               Carp_mf = Carp_mf_pscore1_ord1_cal1_gr,
                               Celtic_bf = Celtic_bf_mahala1_ord1_cal1_gr,
                               CenEu_mf = CenEu_mf_genetic1_gr,
                               EngLow_bf = EngLow_bf_genetic1_ratio1_cal1_gr,
                               EuAtl_mf = EuAtl_mf_genetic1_ratio1_cal1_gr,
                               ItaScl_sdf = ItaScl_sdf_rob_mahala1_ord1_gr,
                               NAtl_mmf = NAtl_mmf_genetic1_cal1_gr,
                               Pan_mf = Pan_mf_rob_mahala1_ord1_cal1_gr,
                               Sar_mf = Sar_mf_genetic1_ratio1_gr)


Matched_datasets_grass <- lapply(Matched_datasets_grass, function(mobj) {
  
  #extract matched data - drop 'matchdata' class because this obj type is not accepted by gdm package
  dtf <- as.data.frame(MatchIt::match_data(mobj))
  
  #check weights are all 1
  if(!all(dtf[['weights']] == 1)) stop('Not all weights are 1')
  
  #drop weights, subclass and distance columns
  dtf <- dtf[setdiff(colnames(dtf), c('distance', 'weights', 'subclass'))]
  
  #split data frame in period 1 and 2
  dtf_pr1 <- dtf[which(dtf[['Period']] == 'period1'), ]
  
  dtf_pr2 <- dtf[which(dtf[['Period']] == 'period2'), ]
  
  #res
  res <- list(Period1 = dtf_pr1, Period2 =  dtf_pr2)
  
  return(res)
  
  })


#check all datasets have min sample size (1k observations per period)
all(sapply(Matched_datasets_grass, function(eco) all(sapply(eco, function(prd) nrow(prd) >= 1000L))))


#save data to run GDM in another R project
save(Matched_datasets_grass, EVA_veg, file = '/Temporary_proj_run_GDM/Tmp_data_for_GDM.RData')

#save EVA_duply - these plot id will be used to filter out observations to fit GDMs
save(EVA_duply, file = '/Temporary_proj_run_GDM/EVA_duplicates.RData')


#-------------------------------------------------statistical matching of forests

#use grass_col_to_use because columns to select in Forest_meta are the same

all(env_var_nm %in% colnames(Forest_sel_meta)) #T

#drop NAs for environmental drivers
sapply(Forest_sel_meta[env_var_nm], function(cl) sum(is.na(cl)))

forest_loc_to_drop <- complete.cases(Forest_sel_meta[env_var_nm])

sum(!forest_loc_to_drop) #1375 (out of 111815)

Forest_sel_meta <- Forest_sel_meta[forest_loc_to_drop, ]

#check abundance of observations for the different ecoregions
table(Forest_sel_meta$ECO_NAME)
table(Forest_sel_meta$ECO_NAME, Forest_sel_meta$Period) #only Aegean_and_Western_Turkey_sclerophyllous_and_mixed_forests has less than 1k obs

#check habitat type
unique(Forest_sel_meta$Eunis_lev1) #'T'

#check NAs in Sampl_year and Period
sum(is.na(Forest_sel_meta$Sampl_year)) #0
sum(is.na(Forest_sel_meta$Period)) #0

#check plot size
sum(is.na(Forest_sel_meta$Releve_area_m2)) #0
range(Forest_sel_meta$Releve_area_m2) #1 1000
table(Forest_sel_meta$Releve_area_m2)

#drop observations associated with Releve_area_m2 < 50 m2
sum(Forest_sel_meta$Releve_area_m2 < 50) #10201

Forest_sel_meta <- Forest_sel_meta[Forest_sel_meta$Releve_area_m2 >= 50, ] #100239

#check abundance of observations for the different ecoregions (per period)
table(Forest_sel_meta$ECO_NAME, Forest_sel_meta$Period) #only Aegean_and_Western_Turkey_sclerophyllous_and_mixed_forests has less than 1k obs

#drop observations for Aegean_and_Western_Turkey_sclerophyllous_and_mixed_forests
Forest_sel_meta <- Forest_sel_meta[Forest_sel_meta$ECO_NAME != 'Aegean_and_Western_Turkey_sclerophyllous_and_mixed_forests', ] #98064

#check differences in ecoregions between forests and grasslands
setdiff(unique(Forest_sel_meta$ECO_NAME), unique(Grass_sel_meta$ECO_NAME)) #Dinaric_Mountains_mixed_forests (and TyrrAdriaticSclerophyllous which was kept in Grass_sel_meta but not matched)
setdiff(unique(Grass_sel_meta$ECO_NAME), unique(Forest_sel_meta$ECO_NAME)) #5 more

for_sel_econm <- unique(Forest_sel_meta$ECO_NAME)

for_sel_econm <- setNames(for_sel_econm, nm = c('WesEu_bf', 'Alps_cmf', 'Pan_mf', 'EuAtl_mf', 'Sar_mf', 'CenEu_mf', 'Carp_mf',
                                                'ItaScl_sdf', 'DinMon_mf', 'TyrAdr_smf'))

#check
all(for_sel_econm %in% unique(Forest_sel_meta$ECO_NAME)) #T

#order for_sel_econm alphabetically
for_sel_econm <- sort(for_sel_econm)

#select columns for matching and add Period_bin column 
For_sel_ecor <- lapply(for_sel_econm, function(eco_nm) {
  
  #select grass_col_to_use
  dtf <- Forest_sel_meta[Forest_sel_meta$ECO_NAME == eco_nm, grass_col_to_use]
  
  #add Period_bin col
  prd_tr <- table(dtf[['Period']])
  prd_tr <- names(prd_tr[which.min(prd_tr)])
  
  dtf$Period_bin <- ifelse(dtf$Period == prd_tr, 1, 0)
  
  return(dtf)
  
  })


#check
lapply(For_sel_ecor, function(i) table(i[['Period_bin']]))

#check correlation between topographic variables before running the statistical matching

#roughness and slope are super correlated - removing slope fixes the problem
lapply(For_sel_ecor, function(i) {
  
  require(car)
  
  dtf <- i[c('Period_bin', 'Elevation', 'Roughness', 'Slope')]
  
  vif_vals <- car::vif(lm(formula(dtf), data = dtf))
  
  return(vif_vals)
  
})

#the formula for matchit will always be the same, so it's worth creating an object for it
formula_for_matchit <- as.formula("Period_bin ~ Elevation + Roughness")

#create a vector of seeds for genetic algorithm
set.seed(5768)

forest_seeds <- round(runif(n = length(For_sel_ecor), min = 1, max = 1000), digits = 0)

forest_seeds <- setNames(object = forest_seeds, nm = names(For_sel_ecor))

# ------  1. Alps_conifer_and_mixed_forests - PScore ------

table(For_sel_ecor$Alps_cmf$Period)
table(For_sel_ecor$Alps_cmf$Period_bin)

#------check initial imbalance

Alps_cmf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = NULL, distance = 'glm')

summary(Alps_cmf_init_for)
plot(summary(Alps_cmf_init_for))
plot(Alps_cmf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Alps_cmf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = 'nearest', distance = 'glm')

Alps_cmf_pscore1_for
summary(Alps_cmf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Alps_cmf_pscore1_for, un = FALSE, interactions = T)
plot(Alps_cmf_pscore1_for, type = 'jitter', interactive = F)
plot(Alps_cmf_pscore1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_pscore1_for))
plot(summary(Alps_cmf_pscore1_for, interactions = T))
plot(summary(Alps_cmf_pscore1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_pscore1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = 'nearest',
                                           distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(Alps_cmf_pscore1_ord1_ratio1_for, un = F, interactions = T)
plot(Alps_cmf_pscore1_ord1_ratio1_for, type = 'jitter', interactive = F)
plot(Alps_cmf_pscore1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_pscore1_ord1_ratio1_for, un = F))
plot(summary(Alps_cmf_pscore1_ord1_ratio1_for, interactions = T, un = F))

#--Mahalanobis

Alps_cmf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = 'nearest', distance = 'mahalanobis')

summary(Alps_cmf_mahala1_for)
plot(Alps_cmf_mahala1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_mahala1_for, un = F))
plot(summary(Alps_cmf_mahala1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_mahala1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = 'nearest',
                                           distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Alps_cmf_mahala1_ord1_ratio1_for, un = T, interactions = T)
plot(Alps_cmf_mahala1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_mahala1_ord1_ratio1_for, un = F))
plot(summary(Alps_cmf_mahala1_ord1_ratio1_for, interactions = T, un = F))

#--robust Mahalanobis

Alps_cmf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(Alps_cmf_rob_mahala1_for)
plot(Alps_cmf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_rob_mahala1_for, un = F))
plot(summary(Alps_cmf_rob_mahala1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
Alps_cmf_rob_mahala1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf, method = 'nearest',
                                               distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(Alps_cmf_rob_mahala1_ord1_ratio1_for, un = F, interactions = T)
plot(Alps_cmf_rob_mahala1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_rob_mahala1_ord1_ratio1_for, un = F))
plot(summary(Alps_cmf_rob_mahala1_ord1_ratio1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['Alps_cmf']])
Alps_cmf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Alps_cmf_genetic1_for)
plot(Alps_cmf_genetic1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_genetic1_for, un = F))
plot(summary(Alps_cmf_genetic1_for, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'
#ratio
set.seed(forest_seeds[['Alps_cmf']])
Alps_cmf_genetic1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Alps_cmf,
                                       method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(Alps_cmf_genetic1_ratio1_for, un = F, interactions = T)
plot(Alps_cmf_genetic1_ratio1_for, type = 'density', interactive = F)
plot(summary(Alps_cmf_genetic1_ratio1_for, un = F))
plot(summary(Alps_cmf_genetic1_ratio1_for, interactions = T, un = F))

#--check matching performance

#excluding genetic since it does not provide good balance in case of k:1 ratio
Alps_for_perf <- do.call(rbind, lapply(list(PScore = Alps_cmf_pscore1_ord1_ratio1_for,
                                            Mah = Alps_cmf_mahala1_ord1_ratio1_for,
                                            RMah = Alps_cmf_rob_mahala1_ord1_ratio1_for), check_match_performance))

Alps_for_perf #PScore


# ------ 2. Carpathian_montane_forests - Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$Carp_mf$Period)
table(For_sel_ecor$Carp_mf$Period_bin)

#------check initial imbalance

Carp_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = NULL, distance = 'glm')

summary(Carp_mf_init_for)
plot(summary(Carp_mf_init_for))
plot(Carp_mf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Carp_mf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = 'nearest', distance = 'glm')

summary(Carp_mf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Carp_mf_pscore1_for, un = FALSE, interactions = T)
plot(Carp_mf_pscore1_for, type = 'jitter', interactive = F)
plot(Carp_mf_pscore1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_pscore1_for))
plot(summary(Carp_mf_pscore1_for, interactions = T))
plot(summary(Carp_mf_pscore1_for, interactions = T, un = F))

#order
Carp_mf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = 'nearest',
                                   distance = 'glm', m.order = 'closest')

summary(Carp_mf_pscore1_ord1_for, un = F, interactions = T)
plot(Carp_mf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(Carp_mf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_pscore1_ord1_for, un = F))
plot(summary(Carp_mf_pscore1_ord1_for, interactions = T, un = F))

#--Mahalanobis

Carp_mf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = 'nearest', distance = 'mahalanobis')

summary(Carp_mf_mahala1_for)
plot(Carp_mf_mahala1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_mahala1_for, un = F))
plot(summary(Carp_mf_mahala1_for, interactions = T, un = F))

#order
Carp_mf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = 'nearest',
                                   distance = 'mahalanobis', m.order = 'closest')

summary(Carp_mf_mahala1_ord1_for, un = F, interactions = T)
plot(Carp_mf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_mahala1_ord1_for, un = F))
plot(summary(Carp_mf_mahala1_ord1_for, interactions = T, un = F))

#--robust Mahalanobis

Carp_mf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(Carp_mf_rob_mahala1_for)
plot(Carp_mf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_rob_mahala1_for, un = F))
plot(summary(Carp_mf_rob_mahala1_for, interactions = T, un = F))

#order
Carp_mf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf, method = 'nearest',
                                       distance = 'robust_mahalanobis', m.order = 'closest')

summary(Carp_mf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(Carp_mf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_rob_mahala1_ord1_for, un = F))
plot(summary(Carp_mf_rob_mahala1_ord1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['Carp_mf']])
Carp_mf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Carp_mf,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Carp_mf_genetic1_for)
plot(Carp_mf_genetic1_for, type = 'density', interactive = F)
plot(summary(Carp_mf_genetic1_for, un = F))
plot(summary(Carp_mf_genetic1_for, interactions = T, un = F))

#--check matching performance

#
Carp_for_perf <- do.call(rbind, lapply(list(PScore = Carp_mf_pscore1_ord1_for,
                                            Mah = Carp_mf_mahala1_ord1_for,
                                            RMah = Carp_mf_rob_mahala1_for,
                                            Genetic = Carp_mf_genetic1_for), check_match_performance))

Carp_for_perf #PScore and Mah provide very close performances - choosing Mah only becaause it has slightly lower Avg_smd


# ------ 3. Central_European_mixed_forests - PScore ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$CenEu_mf$Period)
table(For_sel_ecor$CenEu_mf$Period_bin)

#------check initial imbalance

CenEu_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = NULL, distance = 'glm')

summary(CenEu_mf_init_for)
plot(summary(CenEu_mf_init_for))
plot(CenEu_mf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
CenEu_mf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest', distance = 'glm')

summary(CenEu_mf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(CenEu_mf_pscore1_for, un = FALSE, interactions = T)
plot(CenEu_mf_pscore1_for, type = 'jitter', interactive = F)
plot(CenEu_mf_pscore1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_pscore1_for))
plot(summary(CenEu_mf_pscore1_for, interactions = T))
plot(summary(CenEu_mf_pscore1_for, interactions = T, un = F))

#order
CenEu_mf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                    distance = 'glm', m.order = 'closest')

summary(CenEu_mf_pscore1_ord1_for, un = F, interactions = T)
plot(CenEu_mf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(CenEu_mf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_pscore1_ord1_for, un = F))
plot(summary(CenEu_mf_pscore1_ord1_for, interactions = T, un = F))

#using caliper on all covariates (+ prop score)
CenEu_mf_pscore1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                        distance = 'glm', m.order = 'closest',
                                        std.caliper = TRUE, caliper = c(2, Elevation = 2, Roughness = 2),
                                        link = "linear.logit")

summary(CenEu_mf_pscore1_ord1_cal1_for, un = F, interactions = T) #377 tr obs are not matched due to the caliper (besides control obs)
plot(summary(CenEu_mf_pscore1_ord1_cal1_for, interactions = T, un = F))

#--Mahalanobis

CenEu_mf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest', distance = 'mahalanobis')

summary(CenEu_mf_mahala1_for)
plot(CenEu_mf_mahala1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_mahala1_for, un = F))
plot(summary(CenEu_mf_mahala1_for, interactions = T, un = F))

#order
CenEu_mf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                    distance = 'mahalanobis', m.order = 'closest')

summary(CenEu_mf_mahala1_ord1_for, un = F, interactions = T)
plot(CenEu_mf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_mahala1_ord1_for, un = F))
plot(summary(CenEu_mf_mahala1_ord1_for, interactions = T, un = F))

#using caliper on all covariates
CenEu_mf_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                        distance = 'mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                        caliper = c(Elevation = 2, Roughness = 2))

summary(CenEu_mf_mahala1_ord1_cal1_for, un = F, interactions = T) #341 tr obs are not matched due to the caliper (besides control obs)
plot(summary(CenEu_mf_mahala1_ord1_cal1_for, interactions = T, un = F))

#--robust Mahalanobis

CenEu_mf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf,
                                   method = 'nearest', distance = 'robust_mahalanobis')

summary(CenEu_mf_rob_mahala1_for)
plot(CenEu_mf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_rob_mahala1_for, un = F))
plot(summary(CenEu_mf_rob_mahala1_for, interactions = T, un = F))

#order
CenEu_mf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                        distance = 'robust_mahalanobis', m.order = 'closest')

summary(CenEu_mf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(CenEu_mf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_rob_mahala1_ord1_for, un = F))
plot(summary(CenEu_mf_rob_mahala1_ord1_for, interactions = T, un = F))

#using caliper on all covariates
CenEu_mf_rob_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf, method = 'nearest',
                                            distance = 'robust_mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                            caliper = c(Elevation = 2, Roughness = 2))

summary(CenEu_mf_rob_mahala1_ord1_cal1_for, un = F, interactions = T) #303 tr obs are not matched due to the caliper (besides control obs)
plot(summary(CenEu_mf_rob_mahala1_ord1_cal1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['CenEu_mf']])
CenEu_mf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf,
                                method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(CenEu_mf_genetic1_for)
plot(CenEu_mf_genetic1_for, type = 'density', interactive = F)
plot(summary(CenEu_mf_genetic1_for, un = F))
plot(summary(CenEu_mf_genetic1_for, interactions = T, un = F))

#using caliper on all covariates
set.seed(forest_seeds[['CenEu_mf']])
CenEu_mf_genetic1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$CenEu_mf,
                                    method = 'genetic', pop.size = 100, distance = 'mahalanobis', std.caliper = TRUE,
                                    caliper = c(Elevation = 1.5, Roughness = 2))

summary(CenEu_mf_genetic1_cal1_for, un = F, interactions = T) #189 tr obs are not matched due to the caliper (besides control obs)
plot(summary(CenEu_mf_genetic1_cal1_for, interactions = T, un = F))


#--check matching performance

#all methods with caliper
CenEu_for_perf <- do.call(rbind, lapply(list(PScore = CenEu_mf_pscore1_ord1_cal1_for,
                                            Mah = CenEu_mf_mahala1_ord1_cal1_for,
                                            RMah = CenEu_mf_rob_mahala1_ord1_cal1_for,
                                            Genetic = CenEu_mf_genetic1_cal1_for), check_match_performance))

CenEu_for_perf #PScore and Mah provide very similar performances - selecting PScore only because it has a lower Max_smd

# ------ 4. Dinaric_Mountains_mixed_forests - PScore ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$DinMon_mf$Period)
table(For_sel_ecor$DinMon_mf$Period_bin)

#------check initial imbalance

DinMon_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = NULL, distance = 'glm')

summary(DinMon_mf_init_for)
plot(summary(DinMon_mf_init_for))
plot(DinMon_mf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
DinMon_mf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = 'nearest', distance = 'glm')

summary(DinMon_mf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(DinMon_mf_pscore1_for, un = FALSE, interactions = T)
plot(DinMon_mf_pscore1_for, type = 'jitter', interactive = F)
plot(DinMon_mf_pscore1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_pscore1_for))
plot(summary(DinMon_mf_pscore1_for, interactions = T))
plot(summary(DinMon_mf_pscore1_for, interactions = T, un = F))

#order
DinMon_mf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = 'nearest',
                                     distance = 'glm', m.order = 'closest')

summary(DinMon_mf_pscore1_ord1_for, un = F, interactions = T)
plot(DinMon_mf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(DinMon_mf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_pscore1_ord1_for, un = F))
plot(summary(DinMon_mf_pscore1_ord1_for, interactions = T, un = F))

#--Mahalanobis

DinMon_mf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = 'nearest', distance = 'mahalanobis')

summary(DinMon_mf_mahala1_for)
plot(DinMon_mf_mahala1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_mahala1_for, un = F))
plot(summary(DinMon_mf_mahala1_for, interactions = T, un = F))

#order
DinMon_mf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = 'nearest',
                                     distance = 'mahalanobis', m.order = 'closest')

summary(DinMon_mf_mahala1_ord1_for, un = F, interactions = T)
plot(DinMon_mf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_mahala1_ord1_for, un = F))
plot(summary(DinMon_mf_mahala1_ord1_for, interactions = T, un = F))

#--robust Mahalanobis

DinMon_mf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(DinMon_mf_rob_mahala1_for)
plot(DinMon_mf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_rob_mahala1_for, un = F))
plot(summary(DinMon_mf_rob_mahala1_for, interactions = T, un = F))

#order
DinMon_mf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf, method = 'nearest',
                                         distance = 'robust_mahalanobis', m.order = 'closest')

summary(DinMon_mf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(DinMon_mf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_rob_mahala1_ord1_for, un = F))
plot(summary(DinMon_mf_rob_mahala1_ord1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['DinMon_mf']])
DinMon_mf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$DinMon_mf,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(DinMon_mf_genetic1_for)
plot(DinMon_mf_genetic1_for, type = 'density', interactive = F)
plot(summary(DinMon_mf_genetic1_for, un = F))
plot(summary(DinMon_mf_genetic1_for, interactions = T, un = F))

#--check matching performance

#all methods with caliper
DinMon_for_perf <- do.call(rbind, lapply(list(PScore = DinMon_mf_pscore1_ord1_for,
                                             Mah = DinMon_mf_mahala1_ord1_for,
                                             RMah = DinMon_mf_rob_mahala1_ord1_for,
                                             Genetic = DinMon_mf_genetic1_for), check_match_performance))

DinMon_for_perf #PScore


# ------ 5. European_Atlantic_mixed_forests - GENETIC ------

table(For_sel_ecor$EuAtl_mf$Period)
table(For_sel_ecor$EuAtl_mf$Period_bin)

#------check initial imbalance

EuAtl_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = NULL, distance = 'glm')

summary(EuAtl_mf_init_for)
plot(summary(EuAtl_mf_init_for))
plot(EuAtl_mf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
EuAtl_mf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest', distance = 'glm')

EuAtl_mf_pscore1_for
summary(EuAtl_mf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(EuAtl_mf_pscore1_for, un = FALSE, interactions = T)
plot(EuAtl_mf_pscore1_for, type = 'jitter', interactive = F)
plot(EuAtl_mf_pscore1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_pscore1_for))
plot(summary(EuAtl_mf_pscore1_for, interactions = T))
plot(summary(EuAtl_mf_pscore1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_pscore1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                            distance = 'glm', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr


summary(EuAtl_mf_pscore1_ord1_ratio1_for, un = F, interactions = T)
plot(EuAtl_mf_pscore1_ord1_ratio1_for, type = 'jitter', interactive = F)
plot(EuAtl_mf_pscore1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_for, un = F))
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_for, interactions = T, un = F))

#caliper on all covariates and prop score (for both SMD and Var Ratio)
EuAtl_mf_pscore1_ord1_ratio1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                                 distance = 'glm', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                                 link = 'linear.logit', caliper = c(2, Elevation = 2, Roughness = 2))

summary(EuAtl_mf_pscore1_ord1_ratio1_cal1_for, un = F, interactions = T) #92 tr obs unmatched - Matched ESS < Matched
plot(summary(EuAtl_mf_pscore1_ord1_ratio1_cal1_for, un = F, interactions = T))


#--Mahalanobis

EuAtl_mf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest', distance = 'mahalanobis')

summary(EuAtl_mf_mahala1_for)
plot(EuAtl_mf_mahala1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_mahala1_for, un = F))
plot(summary(EuAtl_mf_mahala1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_mahala1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EuAtl_mf_mahala1_ord1_ratio1_for, un = T, interactions = T)
plot(EuAtl_mf_mahala1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_for, un = F))
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_for, interactions = T, un = F))

#caliper on all covariates (for both SMD and Var Ratio)
EuAtl_mf_mahala1_ord1_ratio1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                                 distance = 'mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                                 caliper = c(Elevation = 2, Roughness = 2))

summary(EuAtl_mf_mahala1_ord1_ratio1_cal1_for, un = FALSE, interactions = TRUE) #114 tr obs unmatched - Matched ESS < Matched
plot(summary(EuAtl_mf_mahala1_ord1_ratio1_cal1_for, un = FALSE, interactions = TRUE))

#--robust Mahalanobis

EuAtl_mf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf,
                                    method = 'nearest', distance = 'robust_mahalanobis')

summary(EuAtl_mf_rob_mahala1_for)
plot(EuAtl_mf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_rob_mahala1_for, un = F))
plot(summary(EuAtl_mf_rob_mahala1_for, interactions = T, un = F))

#order + ratio -> this increases precision at the expenses of balance
EuAtl_mf_rob_mahala1_ord1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2) #2 control units matched to 1 tr

summary(EuAtl_mf_rob_mahala1_ord1_ratio1_for, un = F, interactions = T)
plot(EuAtl_mf_rob_mahala1_ord1_ratio1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_for, un = F))
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_for, interactions = T, un = F))

#caliper on all covariates (for both SMD and Var Ratio)
EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf, method = 'nearest',
                                                 distance = 'robust_mahalanobis', m.order = 'closest', ratio = 2, std.caliper = TRUE,
                                                 caliper = c(Elevation = 2, Roughness = 2))

summary(EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_for, un = FALSE, interactions = TRUE) #101 tr obs unmatched - Matched ESS < Matched
plot(summary(EuAtl_mf_rob_mahala1_ord1_ratio1_cal1_for, un = FALSE, interactions = TRUE))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['EuAtl_mf']])
EuAtl_mf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf,
                                 method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(EuAtl_mf_genetic1_for)
plot(EuAtl_mf_genetic1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_genetic1_for, un = F))
plot(summary(EuAtl_mf_genetic1_for, interactions = T, un = F))

#order: it is not possible to set m.order = 'closest'
#ratio
set.seed(forest_seeds[['EuAtl_mf']])
EuAtl_mf_genetic1_ratio1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis')


summary(EuAtl_mf_genetic1_ratio1_for, un = F, interactions = T)
plot(EuAtl_mf_genetic1_ratio1_for, type = 'density', interactive = F)
plot(summary(EuAtl_mf_genetic1_ratio1_for, un = F))
plot(summary(EuAtl_mf_genetic1_ratio1_for, interactions = T, un = F))

#caliper on all covariates (for both SMD and Var Ratio)
set.seed(forest_seeds[['EuAtl_mf']])
EuAtl_mf_genetic1_ratio1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$EuAtl_mf,
                                        method = 'genetic', ratio = 2, pop.size = 100, distance = 'mahalanobis',
                                        std.caliper = TRUE, caliper = c(Elevation = 2, Roughness = 2))

summary(EuAtl_mf_genetic1_ratio1_cal1_for, un = FALSE, interactions = TRUE) #363 tr obs unmatched
plot(summary(EuAtl_mf_genetic1_ratio1_cal1_for, un = FALSE, interactions = TRUE))

#--check matching performance

#
EuAtl_for_perf <- check_match_performance(mobj = EuAtl_mf_genetic1_ratio1_cal1_for) 

EuAtl_for_perf #Genetic is the only method providing balance and allowing to match each treatm obs to the same number of control obs (2)


# ------ 6. Italian_sclerophyllous_and_semi_deciduous_forests - PScore ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$ItaScl_sdf$Period)
table(For_sel_ecor$ItaScl_sdf$Period_bin)

#------check initial imbalance

ItaScl_sdf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = NULL, distance = 'glm')

summary(ItaScl_sdf_init_for)
plot(summary(ItaScl_sdf_init_for))
plot(ItaScl_sdf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
ItaScl_sdf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest', distance = 'glm')

summary(ItaScl_sdf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(ItaScl_sdf_pscore1_for, un = FALSE, interactions = T)
plot(ItaScl_sdf_pscore1_for, type = 'jitter', interactive = F)
plot(ItaScl_sdf_pscore1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_pscore1_for))
plot(summary(ItaScl_sdf_pscore1_for, interactions = T))
plot(summary(ItaScl_sdf_pscore1_for, interactions = T, un = F))

#order
ItaScl_sdf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                      distance = 'glm', m.order = 'closest')

summary(ItaScl_sdf_pscore1_ord1_for, un = F, interactions = T)
plot(ItaScl_sdf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(ItaScl_sdf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_pscore1_ord1_for, un = F))
plot(summary(ItaScl_sdf_pscore1_ord1_for, interactions = T, un = F))

#using caliper on roughness and prop score
ItaScl_sdf_pscore1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                          distance = 'glm', m.order = 'closest',
                                          std.caliper = TRUE, caliper = c(2, Roughness = 2),
                                          link = "linear.logit")

summary(ItaScl_sdf_pscore1_ord1_cal1_for, un = F, interactions = T) #72 tr obs are not matched due to the caliper (besides control obs)
plot(summary(ItaScl_sdf_pscore1_ord1_cal1_for, interactions = T, un = F))

#--Mahalanobis

ItaScl_sdf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest', distance = 'mahalanobis')

summary(ItaScl_sdf_mahala1_for)
plot(ItaScl_sdf_mahala1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_mahala1_for, un = F))
plot(summary(ItaScl_sdf_mahala1_for, interactions = T, un = F))

#order
ItaScl_sdf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                      distance = 'mahalanobis', m.order = 'closest')

summary(ItaScl_sdf_mahala1_ord1_for, un = F, interactions = T)
plot(ItaScl_sdf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_mahala1_ord1_for, un = F))
plot(summary(ItaScl_sdf_mahala1_ord1_for, interactions = T, un = F))

#using caliper on roughness
ItaScl_sdf_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                          distance = 'mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                          caliper = c(Roughness = 1.9))

summary(ItaScl_sdf_mahala1_ord1_cal1_for, un = F, interactions = T) #38 tr obs are not matched due to the caliper (besides control obs)
plot(summary(ItaScl_sdf_mahala1_ord1_cal1_for, interactions = T, un = F))

#--robust Mahalanobis

ItaScl_sdf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf,
                                     method = 'nearest', distance = 'robust_mahalanobis')

summary(ItaScl_sdf_rob_mahala1_for)
plot(ItaScl_sdf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_rob_mahala1_for, un = F))
plot(summary(ItaScl_sdf_rob_mahala1_for, interactions = T, un = F))

#order
ItaScl_sdf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                          distance = 'robust_mahalanobis', m.order = 'closest')

summary(ItaScl_sdf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(ItaScl_sdf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_rob_mahala1_ord1_for, un = F))
plot(summary(ItaScl_sdf_rob_mahala1_ord1_for, interactions = T, un = F))

#using caliper on roughness
ItaScl_sdf_rob_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf, method = 'nearest',
                                              distance = 'robust_mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                              caliper = c(Roughness = 2))

summary(ItaScl_sdf_rob_mahala1_ord1_cal1_for, un = F, interactions = T) #55 tr obs are not matched due to the caliper (besides control obs)
plot(summary(ItaScl_sdf_rob_mahala1_ord1_cal1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['ItaScl_sdf']])
ItaScl_sdf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf,
                                  method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(ItaScl_sdf_genetic1_for)
plot(ItaScl_sdf_genetic1_for, type = 'density', interactive = F)
plot(summary(ItaScl_sdf_genetic1_for, un = F))
plot(summary(ItaScl_sdf_genetic1_for, interactions = T, un = F))

#using caliper on roughness
set.seed(forest_seeds[['ItaScl_sdf']])
ItaScl_sdf_genetic1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$ItaScl_sdf,
                                      method = 'genetic', pop.size = 100, distance = 'mahalanobis', std.caliper = TRUE,
                                      caliper = c(Roughness = 2))

summary(ItaScl_sdf_genetic1_cal1_for, un = F, interactions = T) #39 tr obs are not matched due to the caliper (besides control obs)
plot(summary(ItaScl_sdf_genetic1_cal1_for, interactions = T, un = F))

#--check matching performance

#all methods with caliper
ItaScl_for_perf <- do.call(rbind, lapply(list(PScore = ItaScl_sdf_pscore1_ord1_cal1_for,
                                             Mah = ItaScl_sdf_mahala1_ord1_cal1_for,
                                             RMah = ItaScl_sdf_rob_mahala1_ord1_cal1_for,
                                             Genetic = ItaScl_sdf_genetic1_cal1_for), check_match_performance))

ItaScl_for_perf #PScore

# ------ 7. Pannonian_mixed_forests - GENETIC ------

#Sample size period 2 is slightly lower than (sample size period 1) * 2 - I'm not testing ratio = 2 for consistency with what done above
#if larger sample size is not at least twice than lower sample size, ratio = 2 is not tested

table(For_sel_ecor$Pan_mf$Period)
table(For_sel_ecor$Pan_mf$Period_bin)

#------check initial imbalance

Pan_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = NULL, distance = 'glm')

summary(Pan_mf_init_for)
plot(summary(Pan_mf_init_for))
plot(Pan_mf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
Pan_mf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = 'nearest', distance = 'glm')

summary(Pan_mf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(Pan_mf_pscore1_for, un = FALSE, interactions = T)
plot(Pan_mf_pscore1_for, type = 'jitter', interactive = F)
plot(Pan_mf_pscore1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_pscore1_for))
plot(summary(Pan_mf_pscore1_for, interactions = T))
plot(summary(Pan_mf_pscore1_for, interactions = T, un = F))

#order
Pan_mf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = 'nearest',
                                       distance = 'glm', m.order = 'closest')

summary(Pan_mf_pscore1_ord1_for, un = F, interactions = T)
plot(Pan_mf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(Pan_mf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_pscore1_ord1_for, un = F))
plot(summary(Pan_mf_pscore1_ord1_for, interactions = T, un = F))

#--Mahalanobis

Pan_mf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = 'nearest', distance = 'mahalanobis')

summary(Pan_mf_mahala1_for)
plot(Pan_mf_mahala1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_mahala1_for, un = F))
plot(summary(Pan_mf_mahala1_for, interactions = T, un = F))

#order
Pan_mf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = 'nearest',
                                       distance = 'mahalanobis', m.order = 'closest')

summary(Pan_mf_mahala1_ord1_for, un = F, interactions = T)
plot(Pan_mf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_mahala1_ord1_for, un = F))
plot(summary(Pan_mf_mahala1_ord1_for, interactions = T, un = F))

#--robust Mahalanobis

Pan_mf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf,
                                      method = 'nearest', distance = 'robust_mahalanobis')

summary(Pan_mf_rob_mahala1_for)
plot(Pan_mf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_rob_mahala1_for, un = F))
plot(summary(Pan_mf_rob_mahala1_for, interactions = T, un = F))

#order
Pan_mf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf, method = 'nearest',
                                           distance = 'robust_mahalanobis', m.order = 'closest')

summary(Pan_mf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(Pan_mf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_rob_mahala1_ord1_for, un = F))
plot(summary(Pan_mf_rob_mahala1_ord1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['Pan_mf']])
Pan_mf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$Pan_mf,
                                   method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(Pan_mf_genetic1_for)
plot(Pan_mf_genetic1_for, type = 'density', interactive = F)
plot(summary(Pan_mf_genetic1_for, un = F))
plot(summary(Pan_mf_genetic1_for, interactions = T, un = F))

#--check matching performance

#
Pan_for_perf <- do.call(rbind, lapply(list(PScore = Pan_mf_pscore1_ord1_for,
                                              Mah = Pan_mf_mahala1_ord1_for,
                                              RMah = Pan_mf_rob_mahala1_ord1_for,
                                              Genetic = Pan_mf_genetic1_for), check_match_performance))

Pan_for_perf #Genetic

# ------ 8. Sarmatic_mixed_forests - NO NEED OF MATCHING ------

table(For_sel_ecor$Sar_mf$Period)
table(For_sel_ecor$Sar_mf$Period_bin)

#------check initial imbalance

Sar_mf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$Sar_mf, method = NULL, distance = 'glm')

#No need to match data for ratio 1:1
summary(Sar_mf_init_for)
plot(summary(Sar_mf_init_for))
plot(summary(Sar_mf_init_for, interactions = T))
plot(Sar_mf_init_for, type = 'density')

# ------ 9. Tyrrhenian_Adriatic_sclerophyllous_and_mixed_forests - Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$TyrAdr_smf$Period)
table(For_sel_ecor$TyrAdr_smf$Period_bin)

#------check initial imbalance

TyrAdr_smf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = NULL, distance = 'glm')

summary(TyrAdr_smf_init_for)
plot(summary(TyrAdr_smf_init_for))
plot(TyrAdr_smf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
TyrAdr_smf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest', distance = 'glm')

summary(TyrAdr_smf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(TyrAdr_smf_pscore1_for, un = FALSE, interactions = T)
plot(TyrAdr_smf_pscore1_for, type = 'jitter', interactive = F)
plot(TyrAdr_smf_pscore1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_pscore1_for))
plot(summary(TyrAdr_smf_pscore1_for, interactions = T))
plot(summary(TyrAdr_smf_pscore1_for, interactions = T, un = F))

#order
TyrAdr_smf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                      distance = 'glm', m.order = 'closest')

summary(TyrAdr_smf_pscore1_ord1_for, un = F, interactions = T)
plot(TyrAdr_smf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(TyrAdr_smf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_pscore1_ord1_for, un = F))
plot(summary(TyrAdr_smf_pscore1_ord1_for, interactions = T, un = F))

#using caliper on roughness (for Var Ratio - see var ratio of roughness^2)
TyrAdr_smf_pscore1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                            distance = 'glm', m.order = 'closest',
                                            std.caliper = TRUE, caliper = c(Roughness = 2))

summary(TyrAdr_smf_pscore1_ord1_cal1_for, un = F, interactions = T) #22 tr obs are not matched due to the caliper (besides control obs)
plot(summary(TyrAdr_smf_pscore1_ord1_cal1_for, interactions = T, un = F))

#--Mahalanobis

TyrAdr_smf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest', distance = 'mahalanobis')

summary(TyrAdr_smf_mahala1_for)
plot(TyrAdr_smf_mahala1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_mahala1_for, un = F))
plot(summary(TyrAdr_smf_mahala1_for, interactions = T, un = F))

#order
TyrAdr_smf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                      distance = 'mahalanobis', m.order = 'closest')

summary(TyrAdr_smf_mahala1_ord1_for, un = F, interactions = T)
plot(TyrAdr_smf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_mahala1_ord1_for, un = F))
plot(summary(TyrAdr_smf_mahala1_ord1_for, interactions = T, un = F))

#using caliper on roughness
TyrAdr_smf_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                            distance = 'mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                            caliper = c(Roughness = 2))

summary(TyrAdr_smf_mahala1_ord1_cal1_for, un = F, interactions = T) #18 tr obs are not matched due to the caliper (besides control obs)
plot(summary(TyrAdr_smf_mahala1_ord1_cal1_for, interactions = T, un = F))

#--robust Mahalanobis

TyrAdr_smf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf,
                                     method = 'nearest', distance = 'robust_mahalanobis')

summary(TyrAdr_smf_rob_mahala1_for)
plot(TyrAdr_smf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_rob_mahala1_for, un = F))
plot(summary(TyrAdr_smf_rob_mahala1_for, interactions = T, un = F))

#order
TyrAdr_smf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                          distance = 'robust_mahalanobis', m.order = 'closest')

summary(TyrAdr_smf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(TyrAdr_smf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_rob_mahala1_ord1_for, un = F))
plot(summary(TyrAdr_smf_rob_mahala1_ord1_for, interactions = T, un = F))

#using caliper on roughness
TyrAdr_smf_rob_mahala1_ord1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf, method = 'nearest',
                                                distance = 'robust_mahalanobis', m.order = 'closest', std.caliper = TRUE,
                                                caliper = c(Roughness = 2))

summary(TyrAdr_smf_rob_mahala1_ord1_cal1_for, un = F, interactions = T) #17 tr obs are not matched due to the caliper (besides control obs)
plot(summary(TyrAdr_smf_rob_mahala1_ord1_cal1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['TyrAdr_smf']])
TyrAdr_smf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf,
                                  method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(TyrAdr_smf_genetic1_for, interactions = T, un = F)
plot(TyrAdr_smf_genetic1_for, type = 'density', interactive = F)
plot(summary(TyrAdr_smf_genetic1_for, un = F))
plot(summary(TyrAdr_smf_genetic1_for, interactions = T, un = F))

#using caliper on roughness
set.seed(forest_seeds[['TyrAdr_smf']])
TyrAdr_smf_genetic1_cal1_for <- matchit(formula_for_matchit, data = For_sel_ecor$TyrAdr_smf,
                                        method = 'genetic', pop.size = 100, distance = 'mahalanobis', std.caliper = TRUE,
                                        caliper = c(Roughness = 2))

summary(TyrAdr_smf_genetic1_cal1_for, un = F, interactions = T) #15 tr obs are not matched due to the caliper (besides control obs)
plot(summary(TyrAdr_smf_genetic1_cal1_for, interactions = T, un = F))

#--check matching performance

#all methods with caliper
TyrAdr_for_perf <- do.call(rbind, lapply(list(PScore = TyrAdr_smf_pscore1_ord1_cal1_for,
                                              Mah = TyrAdr_smf_mahala1_ord1_cal1_for,
                                              RMah = TyrAdr_smf_rob_mahala1_ord1_cal1_for,
                                              Genetic = TyrAdr_smf_genetic1_cal1_for), check_match_performance))

TyrAdr_for_perf #Mahalanobis

# ------ 10. Western_European_broadleaf_forests - Mahalanobis ------

#balanced sample size between periods - don't test ratio = 2

table(For_sel_ecor$WesEu_bf$Period)
table(For_sel_ecor$WesEu_bf$Period_bin)

#------check initial imbalance

WesEu_bf_init_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = NULL, distance = 'glm')

summary(WesEu_bf_init_for)
plot(summary(WesEu_bf_init_for))
plot(WesEu_bf_init_for, type = 'density')

#------nearest neighbor

#--propensity score

#logit
WesEu_bf_pscore1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = 'nearest', distance = 'glm')

summary(WesEu_bf_pscore1_for, un = FALSE) #this prevents comparison pre-matching to be printed
summary(WesEu_bf_pscore1_for, un = FALSE, interactions = T)
plot(WesEu_bf_pscore1_for, type = 'jitter', interactive = F)
plot(WesEu_bf_pscore1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_pscore1_for))
plot(summary(WesEu_bf_pscore1_for, interactions = T))
plot(summary(WesEu_bf_pscore1_for, interactions = T, un = F))

#order
WesEu_bf_pscore1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = 'nearest',
                                       distance = 'glm', m.order = 'closest')

summary(WesEu_bf_pscore1_ord1_for, un = F, interactions = T)
plot(WesEu_bf_pscore1_ord1_for, type = 'jitter', interactive = F)
plot(WesEu_bf_pscore1_ord1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_pscore1_ord1_for, un = F))
plot(summary(WesEu_bf_pscore1_ord1_for, interactions = T, un = F))

#--Mahalanobis

WesEu_bf_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = 'nearest', distance = 'mahalanobis')

summary(WesEu_bf_mahala1_for)
plot(WesEu_bf_mahala1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_mahala1_for, un = F))
plot(summary(WesEu_bf_mahala1_for, interactions = T, un = F))

#order
WesEu_bf_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = 'nearest',
                                       distance = 'mahalanobis', m.order = 'closest')

summary(WesEu_bf_mahala1_ord1_for, un = F, interactions = T)
plot(WesEu_bf_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_mahala1_ord1_for, un = F))
plot(summary(WesEu_bf_mahala1_ord1_for, interactions = T, un = F))

#--robust Mahalanobis

WesEu_bf_rob_mahala1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf,
                                      method = 'nearest', distance = 'robust_mahalanobis')

summary(WesEu_bf_rob_mahala1_for)
plot(WesEu_bf_rob_mahala1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_rob_mahala1_for, un = F))
plot(summary(WesEu_bf_rob_mahala1_for, interactions = T, un = F))

#order
WesEu_bf_rob_mahala1_ord1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf, method = 'nearest',
                                           distance = 'robust_mahalanobis', m.order = 'closest')

summary(WesEu_bf_rob_mahala1_ord1_for, un = F, interactions = T)
plot(WesEu_bf_rob_mahala1_ord1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_rob_mahala1_ord1_for, un = F))
plot(summary(WesEu_bf_rob_mahala1_ord1_for, interactions = T, un = F))

#------genetic matching

#using GMD without prop.score - note that distance is set to mahalanobis so that prop score is not estimated - see examples method_genetic() 
set.seed(forest_seeds[['WesEu_bf']])
WesEu_bf_genetic1_for <- matchit(formula_for_matchit, data = For_sel_ecor$WesEu_bf,
                                   method = 'genetic', pop.size = 100, distance = 'mahalanobis')

summary(WesEu_bf_genetic1_for, interactions = T, un = F)
plot(WesEu_bf_genetic1_for, type = 'density', interactive = F)
plot(summary(WesEu_bf_genetic1_for, un = F))
plot(summary(WesEu_bf_genetic1_for, interactions = T, un = F))

#--check matching performance

#all methods with caliper
WesEu_for_perf <- do.call(rbind, lapply(list(PScore = WesEu_bf_pscore1_ord1_for,
                                              Mah = WesEu_bf_mahala1_ord1_for,
                                              RMah = WesEu_bf_rob_mahala1_ord1_for,
                                              Genetic = WesEu_bf_genetic1_for), check_match_performance))

WesEu_for_perf #Mahalanobis


# ----
# ----

# ------ End of stat matching for forests ------

#----extract matched datasets and prepare them for GDMs

#Alps_conifer_and_mixed_forests: Alps_cmf_pscore1_ord1_ratio1_for
#Carpathian_montane_forests: Carp_mf_mahala1_ord1_for
#Central_European_mixed_forests: CenEu_mf_pscore1_ord1_cal1_for
#Dinaric_Mountains_mixed_forests: DinMon_mf_pscore1_ord1_for
#European_Atlantic_mixed_forests: EuAtl_mf_genetic1_ratio1_cal1_for
#Italian_sclerophyllous_and_semi_deciduous_forests: ItaScl_sdf_pscore1_ord1_cal1_for
#Pannonian_mixed_forests: Pan_mf_genetic1_for
#Sarmatic_mixed_forests: MATCHING NOT NEEDED
#Tyrrhenian_Adriatic_sclerophyllous_and_mixed_forests: TyrAdr_smf_mahala1_ord1_cal1_for
#Western_European_broadleaf_forests: WesEu_bf_mahala1_ord1_for


Matched_datasets_forest <- list(Alps_cmf_pscore1_ord1_ratio1_for,
                                Carp_mf_mahala1_ord1_for,
                                CenEu_mf_pscore1_ord1_cal1_for,
                                DinMon_mf_pscore1_ord1_for,
                                EuAtl_mf_genetic1_ratio1_cal1_for,
                                ItaScl_sdf_pscore1_ord1_cal1_for,
                                Pan_mf_genetic1_for,
                                TyrAdr_smf_mahala1_ord1_cal1_for,
                                WesEu_bf_mahala1_ord1_for)

names(Matched_datasets_forest) <- names(for_sel_econm[-8]) 


Matched_datasets_forest <- lapply(Matched_datasets_forest, function(mobj) {
  
  #extract matched data - drop 'matchdata' class because this obj type is not accepted by gdm package
  dtf <- as.data.frame(MatchIt::match_data(mobj))
  
  #check weights are all 1
  if(!all(dtf[['weights']] == 1)) stop('Not all weights are 1')
  
  #drop weights, subclass and distance columns
  dtf <- dtf[setdiff(colnames(dtf), c('distance', 'weights', 'subclass'))]
  
  #split data frame in period 1 and 2
  dtf_pr1 <- dtf[which(dtf[['Period']] == 'period1'), ]
  
  dtf_pr2 <- dtf[which(dtf[['Period']] == 'period2'), ]
  
  #res
  res <- list(Period1 = dtf_pr1, Period2 =  dtf_pr2)
  
  return(res)
  
})

#attach dataset for Sarmatic_mixed_forests, which did not need matching

identical(colnames(For_sel_ecor$Sar_mf), colnames(Matched_datasets_forest$Alps_cmf$Period1)) #TRUE
'Sar_mf' %in% names(Matched_datasets_forest)

Matched_datasets_forest[['Sar_mf']] <- list(Period1 = For_sel_ecor$Sar_mf[which(For_sel_ecor$Sar_mf$Period == 'period1'), ],
                                                            Period2 = For_sel_ecor$Sar_mf[which(For_sel_ecor$Sar_mf$Period == 'period2'), ])

#re-order alphabetically
Matched_datasets_forest <- Matched_datasets_forest[sort(names(Matched_datasets_forest))]

#check all datasets have min sample size (1k observations per period)
all(sapply(Matched_datasets_forest, function(eco) all(sapply(eco, function(prd) nrow(prd) >= 1000L))))


#save data to run GDM in another R project
save(Matched_datasets_forest, file = '/Temporary_proj_run_GDM/Tmp_data_for_GDM_forest.RData')

#save EVA_veg to run GDM in another R project.
#Notice that EVA_veg was already saved above. However, reloading the obj saved above in the Rprj where I am currently running the GDMs
#will overwrite the list of matched datasets for grasslands. !!!!Next time, save EVA_veg separately, so that I don't have to export it twice
save(EVA_veg, file = '/Temporary_proj_run_GDM/EVA_veg_data.RData')


#-------------------------------------------------create maps showing ecoregions selected for analyses

library(paletteer)

sel_ecor_analysis <- list(Grass = unique(Grass_sel_meta$ECO_NAME),
                          Forest = unique(Forest_sel_meta$ECO_NAME))

length(sel_ecor_analysis$Grass) #14 - drop Tyrrhenian_Adriatic_sclerophyllous_and_mixed_forests that did not meet min smp size requirement (for grass)

sel_ecor_analysis$Grass <- sel_ecor_analysis$Grass[!sel_ecor_analysis$Grass %in% 'Tyrrhenian_Adriatic_sclerophyllous_and_mixed_forests'] 

#drop "_" and add '-' in ecor names to match ECO_NAME in eu_ecoregions.proj
sel_ecor_analysis$Grass <- gsub(pattern = '_', replacement = ' ', x = sel_ecor_analysis$Grass)
sel_ecor_analysis$Forest <- gsub(pattern = '_', replacement = ' ', x = sel_ecor_analysis$Forest)

sel_ecor_analysis$Grass[9] <- "Italian sclerophyllous and semi-deciduous forests"
sel_ecor_analysis$Forest[c(8, 10)] <- c("Italian sclerophyllous and semi-deciduous forests", "Tyrrhenian-Adriatic sclerophyllous and mixed forests")

#create palette for ecoregions
sel_ecor_palette <- unique(unlist(sel_ecor_analysis)) 

sel_ecor_palette <- setNames(as.character(paletteer::paletteer_d("ggsci::default_igv", n = length(sel_ecor_palette))),
                             nm = sel_ecor_palette)

all(names(sel_ecor_palette) %in% eu_ecoregions.proj$ECO_NAME) #TRUE

#create bbox including all ecoregions (both grass and for)
sel_ecor_bbox <- st_bbox(obj = eu_ecoregions.proj[eu_ecoregions.proj$ECO_NAME %in% names(sel_ecor_palette), ])

#not using the palette - I am simply coloring grass in lightgreen and forests in dark green

#grass
plot(st_geometry(st_crop(x = eu_ecoregions.proj, y = sel_ecor_bbox)))
plot(st_geometry(eu_ecoregions.proj[eu_ecoregions.proj$ECO_NAME %in% sel_ecor_analysis$Grass, ]),
     col = 'lightgreen', add = T)


#forest
plot(st_geometry(st_crop(x = eu_ecoregions.proj, y = sel_ecor_bbox)))
plot(st_geometry(eu_ecoregions.proj[eu_ecoregions.proj$ECO_NAME %in% sel_ecor_analysis$Forest, ]),
     col = 'darkgreen', add = T)



