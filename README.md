This repository includes the R code to reproduce results presented in the manuscript: \_Weak evidence of biotic homogenisation in European grasslands and forests despite changes in the influence of climate and human pressure over the past four decades\_.



The manuscript aims to assess taxonomic homogenisation and differentiation in European grasslands and forests, and to investigate temporal changes in the relative importance of climate and human pressure in shaping spatial taxonomic dissimilariy possibly occurred over the last four decades.



The R scripts are organized (and should be used) as presented in the scheme below:



A first set of scripts, which can be found on the repository home, pre-process vegetation-plot observations from the \[\*\*European Vegetation Archive\*\*](https://doi.org/10.1111/avsc.12191) and \[\*\*ReSurveyEurope\*\*](https://doi.org/10.1111/jvs.13235) for subsequent analyses. Specifically, the scripts implement the following steps:



1. EVA\_in\_ecoregions: pre-process EVA and ReSurveyEurope to obtain grassland and forest datasets including observations complying with requirements needed for the analyses;
2. Spat\_oper\_envlayers: process raster layers associated with climatic, topographic and human pressure data -> these layers are then used in the Drivers\_extraction script to extract climate, human pressure and topography values at each plot location (using R functions included in the TopoClim\_for\_ecor script);
3. Match\_data\_ecor: perform, separately for each ecoregion, statistical matching of period-specific plot data based on topography. Results of the statistical matching are summarised in tables created using code included the Match\_summary\_tabs script.



After this first phase of plot data pre-processing, other scripts, contained in task-specific folders, are used to implement the following steps:



1. Data\_for\_analyses folder: includes scripts for formatting data for temporal comparisons of spatial taxonomic dissimilarity (including distance-decay modelling) and for fitting Generalised Dissimilarity Models;
2. DistDecay\_beta: includes scripts for running the temporal comparisons and fitting distance-decay models. This folder includes also the R code for extracting, plotting and summarising in tables the results;
3. GDM\_fitting: includes scripts for fitting the Generalised Dissimilarity Models and for extracting, plotting and summarising in tables the results.

