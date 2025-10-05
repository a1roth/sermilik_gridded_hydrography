# sermilik_gridded_hydrography
Code and data for Sermilik Fjord gridded hydrography using objective mapping method (version 1: June 2025).

October 2025 updates: Depth units have been changed to meters for all gridded data and the climatology (Units previously km, but described as m). Now description and values match. 

## Project Description
Full description in Roth, A. et al., (in review), A fjord dataset for multi-disciplinary applications: Thirteen years of ocean observations in Sermilik Fjord, Southeast Greenland, Earth Systems Science Data (DOI coming soon).

This code is for creating standardized, along-fjord gridded fields of temperature and salinity from discrete CTD and XCTD profiles in Sermilik Fjord. An objective mapping method (also known as optimum or optimal interpolation) is used. Gridded fields are useful for calculating fjord-wide averages (eg. climatologies), forcings or validation for modeling studies, and environmental context for other discrete measurements in the fjord (biogeochemical bottle samples).   

<img width="1590" alt="obmap_overview" src="https://github.com/user-attachments/assets/db79abc2-c906-4a37-ac99-3d29bd554b31" />

## Key Output Files Available For Use
**SermilikGriddedSections_v1.nc**: includes all available gridded along-fjord temperature and salinity fields for 2009 - 2023, along with the mapping relative error related to each gridded field (see Roth et al. (2025) for more information). This dataset is best used to create a timeseries of a specific location in the fjord or to manipulate all the gridded fields together (eg. calculating averages or a climatology). This data can be used by modeling studies for forcing and validation, intercomparison studies with other fjords, and by biogeochemical or ecosystem related studies that require the hydrographic context of Sermilik Fjord.

**SermilikGriddedClimatology_v1.nc**: includes the summer climatology gridded along-fjord temperature and salinity fields for direct use. The root mean square deviation fields for each variable are also included. These are the data displayed in Figure 8 of Roth et al. (2025). These data can be used to extract spatial mean values of the climatology from particular regions of the fjord (eg. the ice mélange region) or the deep (> 400 m depth) Atlantic Water region. 

**SermilikSection_YEAR.nc**: These files for each year contain the manually selected, discrete CTD and XCTD profiles used as input data to generate the gridded fields. These data are NOT the complete dataset for the field campaign of each year in Sermilik Fjord. The complete profile datasets for each individual field campaign are archived separately at [Arctic Data Center: Sermilik Fjord Hydrography Data Portal]([https://pages.github.com/)(https://arcticdata.io/catalog/portals/sermilik/Data)). The data in these files were chosen to best represent a synoptic along-fjord section along the thalweg line of the fjord. Each file also includes the residuals of the original profile data compared with the final gridded fields. Use these files if you are interested in the input data used to generate the gridded fields or interested in the comparison of the original profile data to the final gridded fields.

## Prerequisites
- Code is provided in MATLAB. MATLAB 2020b or later is suggested.
-  

## Table of Contents
- [Input Data Description](#InputDataDescription)
- [Workflow](#Workflow)
- [Parameters](#Parameters)
- [Output data](#Outputdata)
- [Calculating a climatology](#Calculatingaclimatology)
- [Adapting code to other fjord profile datasets](#Adaptingcodetootherfjordprofiledatasets)

## Input Data Description

### Thalweg section

### Background Grids

## Workflow

## Dependecies

## Output Data Examples

## Calculating a Climatology

## Last Updated
This README was last updated on [Date].

## License
This project is licensed under
