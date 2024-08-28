![image of a smallmouth bass](https://github.com/meairey/literate-potato/blob/main/Graphics/SMB_image.jpg?raw=true)


# Smallmouth Bass Long-term Changepoint Analysis
Detmer et al., in review @ the Journal of Applied Ecology

## Study framework 
Little Moose Lake is the site of a longterm smallmouth bass removal project. The fish community within this lake has been monitored for the last two decades. The goal of this analysis is to examin long-term trends in catch per unit effort (CPUE) and size of the fish community in response to that bass removal.

This is primarily conducted through the use of a changepoint analysis through the package `ecp` and a zero-inflated regression through the package `pscl`.  


## Data availability

The data folder is included in the `.gitignore`. Please contact `ma2276@cornell.edu` with questions about data availability. 

## Repository structure

The `Analysis` folder contains the `.R` script for generating changepoints and regressions. The data that goes into this script is generated in the `Documentation_DataCleaning` folder and uses functions as created in the `Function_Source_Files` folder. Graphics for the manuscript are kept in the `Graphics` folder. The `Documentation_DataCleaning` folder also contains information that is included in the supplemental section of this manuscript regarding sampling frequency and standard protocol deviations.

