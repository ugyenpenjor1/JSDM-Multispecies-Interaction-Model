R scripts and data to run multispecies occupancy model for two or more interacting species of [Rota *et al.,* 2016 Methods in Ecology and Evolution] (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12587)

This repo will demonstrate to load and run the model. The codes are annotated for readability.

For convenience and to remove clutter, the data have been cleaned but for security reasons (because the study species are globally endangered) we do not provide GPS locations. 

Specific files included are:

interactionModelData_Penjor2021_dcat.RData = raw data containing species' detection history and covariates

intMod.txt = model specification in BUGS language

modelScript.R = script to load the data and run the model

Details: 

1. Total number of sites = 849
2. Total number of site covariates = 5 (settlement, forest, sambar deer, gaur, serow)
3. Total number of detection covariates = 2 (Trail and Disturb)
4. Settlement density (continuous covariate within a 4 km circular radius of camera trap location)
5. Forest cover (continuous covariate per cent cover within a 4 km radius of camera trap location)
6. Relative abundance (continuous covariate) of i) sambar deer (*Rusa unicolor*); ii) gaur (*Bos gaurus*), and iii) serow (*Capricornis thar*)

Trail = (binary) coded 1/0 for camera traps placed on- or off-trail
Disturb = daily encounter rates of humans, livestock and domestic dogs at each camera station

# The manuscript is *in review*. Please contact ugyenpenjor.bt@gmail.com to learn more about the project.
