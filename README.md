This repository contains scripts to perform simulations on VIRION (VIrtual Ralstonia-tomato plant interactION).
VIRION combines a tomato metabolic model and a bacterial pathogen metabolic model.

The tomato plant model is based on VYTOP, a metabolic model of a whole plant developed with the cell genome-scale metabolic model of tomato Sl2183.
For more informations on the plant model, see:
- the associated article: Gerlin et al., 2022, Plant Physiol https://doi.org/10.1093/plphys/kiab548
- the github repository: https://github.com/lgerlin/slyc-metabolic-model/

The pathogen metabolic model is based on iRP147, the genome-scale metabolic model of R. pseudosolanacearum.
For more information on the pathogen model, see the associated article: Peyraud et al., 2016, PLoS Path https://doi.org/10.1371/journal.ppat.1005939

The scripts require Python (preferably 3.5), CPLEX Python API developed by IBM (free for academics), and the python packages lxml, matplotlib, pickle, numpy and pandas.

The main script to reproduce the simulations presented in the paper is FBA_VIRION_multiple_cfus.py.
It requires parsing the plant and pathogen sbml through parser_sbml_Ralsto.py (pathogen), parser_sbml_VYTOP.py (plant) and merging the two sbml through parser_merge_VIRION.py. To simulate the plant - pathogen interaction, a first simulation without the pathogen is done using the script FBA_VYTOP.py.

When running FBA_VIRION_multiple_cfus.py, you can switch, line 48, between 8 scenarios depending which result of the paper you want to reproduce:
1-PhotonsOnly
has no additional constraints: only photon uptake is taken from the fluxes of a healthy plant
2-PhotonsNitrogen
constrains photons + N uptake using values from the fluxes of a healthy plant
3-PhotonsIron
constrains photons + iron uptake using values from the fluxes of a healthy plant
4-PhotonsTranspiration
constrains photons + reduces xylem fluxes according to transpiration reduction
5-All
constrains photons + reduces xylem fluxes according to transpiration reduction
+ iron/nitrogen uptake using values from the fluxes of a healthy plant
6-All-Stem
authorizes nutrient exchanges from stem to xylem with the constrainst defined in 5
7-PhotonsOnly-Stem
authorizes nutrient exchanges from stem to xylem with the constraints defined in 1
8-PhotonsXylem
limits xylem elements from the amount in healthy plant, with the constraints defined for 1

The input folder contains sbml files for both the pathogen and the plant, and specific files for whole plant modeling as in VYTOP.

The output folder contains simulations results in different formats.

Three analysis scripts (analysis_bottlenecks.py, analaysis_stem_hijacking.py, analysis_putrescine.py) contain the code to reproduce the figures shown in the article.