This repository contains the R code required to reproduce the results presented in the manuscript: *Weak evidence of biotic homogenisation in European grasslands and forests despite changes in the influence of climate and human pressure over the past four decades*.



The analyses aim to assess taxonomic homogenisation and differentiation in European grasslands and forests, and to investigate temporal changes in the relative importance of climate and human pressure in shaping spatial taxonomic dissimilariy over the past four decades.



The R scripts are organized as outlines in the scheme below and should be run in the order presented:



The first set of scripts, available in the repository's root directory, includes code for pre-processing vegetation-plot data from the [**European Vegetation Archive**](https://doi.org/10.1111/avsc.12191) and [**ReSurveyEurope**](https://doi.org/10.1111/jvs.13235) for subsequent analyses. These scripts implement the following steps:



1. **EVA\_in\_ecoregions**: pre-processing of EVA and ReSurveyEurope data to generate grassland and forest datasets. These datasets include only observations that meet the criteria for inclusion in the analyses;
2. **Spat\_oper\_envlayers**: processing of geospatial raster layers containing climatic, topographic and human pressure data. These layers are subsequently used in the **Drivers\_extraction** script to extract climatic, human pressure and topographic variables from the location of each vegetation plot. The R functions used in **Drivers\_extraction** are available in the **TopoClim\_for\_ecor** script;
3. **Match\_data\_ecor**: statistical matching of period-specific plot data based on topographic variables, performed separately for each ecoregion. Results of the statistical matching are summarised in tables generated using the code in the **Match\_summary\_tabs** script.



Following this initial data pre-processing phase, additional sets of scripts, organised into task-specific folders, are used to carry out the following analyses:



1. **Data\_for\_analyses** folder: includes scripts for formatting data for temporal comparisons of spatial taxonomic dissimilarity (including distance-decay modelling) and for fitting Generalised Dissimilarity Models;
2. **DistDecay\_beta**: includes scripts for conducting the temporal comparisons and fitting distance-decay models. This folder also contains the R code used to extract, plot and summarise the results in tables;
3. **GDM\_fitting**: includes scripts for fitting Generalised Dissimilarity Models and for extracting, plotting and summarising the results in tables.