This repository includes the R code to replicate the results presented in the manuscript: *Weak evidence of biotic homogenisation in European grasslands and forests despite changes in the influence of climate and human pressure over the past four decades*.



The analyses aim to assess taxonomic homogenisation and differentiation in European grasslands and forests, and to investigate temporal changes in the relative importance of climate and human pressure in shaping spatial taxonomic dissimilariy over the last four decades.



The R scripts are organized (and should be used) as presented in the scheme below:



A first set of scripts, which can be found on the repository home, includes code for pre-processing vegetation-plot data from the [*European Vegetation Archive*](https://doi.org/10.1111/avsc.12191) and [*ReSurveyEurope*](https://doi.org/10.1111/jvs.13235) for subsequent analyses. The scripts implement the following steps:



1. EVA\_in\_ecoregions: pre-processing of EVA and ReSurveyEurope data to obtain grassland and forest datasets. These datasets include observations complying with requirements for being used in the analyses;
2. Spat\_oper\_envlayers: processing of geospatial raster layers associated with climatic, topographic and human pressure data. These layers are used in the Drivers\_extraction script to retrieve the value of the climatic, human pressure and topographic variables from each vegetation-plot location. The R functions used in Drivers\_extraction can be found in the TopoClim\_for\_ecor script;
3. Match\_data\_ecor: performing, separately for each ecoregion, the statistical matching of period-specific plot data based on topography. Results of the statistical matching are summarised in tables created using the code included the Match\_summary\_tabs script.



After this first phase of data pre-processing, other sets of scripts, contained in task-specific folders, are used to implement the following steps:



1. Data\_for\_analyses folder: includes scripts for formatting data for temporal comparisons of spatial taxonomic dissimilarity (including distance-decay modelling) and for fitting Generalised Dissimilarity Models;
2. DistDecay\_beta: includes scripts for running the temporal comparisons and fitting the distance-decay models. In addition, this folder includes the R code for extracting, plotting and summarising the results in tables;
3. GDM\_fitting: includes scripts for fitting the Generalised Dissimilarity Models and for extracting, plotting and summarising the results in tables.

