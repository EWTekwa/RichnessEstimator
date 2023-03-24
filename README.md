# RichnessEstimator
Code and data for:
An improved species richness estimator using spatial abundance data
Eden Tekwa1,2*, Matthew A. Whalen1,2, Patrick T. Martone3, Mary I. O’Connor1

1Department of Zoology, University of British Columbia, Vancouver, BC, Canada.
2Hakai Institute, Heriot Bay, BC, Canada
3Department of Botany, University of British Columbia, Vancouver, BC, Canada.
*Corresponding author email: ewtekwa@gmail.com, ORCID: 0000-0003-2971-6128

Note: The following documentation is for code and data used to generate the results in the paper cited above. For the R-package with the richness estimator functions, please see the subfolder /Richness.

Main Matlab files:
Data_all.mat
RichnessEsts.m
RichnessEsts_ideal.m
biodiversitySamplingSim.m
bootRichnessEsts.m
bootRichnessEsts_all.m
bootRichnessEsts_all_indepIndivAssumed.m
runEmpiricalRichnessCorrection.m

Third Party Matlab files (by Kelly Kearney 2010, available at https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m):
boundedline.m
catuneven.m
inpaint_nans.m
outlinebounds.m
singlepatch.m

R files:
RichnessEsts.R
bootRichnessEsts.R

1. Run the "biodiversitySamplingSim.m" script to evaluate species richness estimators using simulated scenarios. The script calls "RichnessEsts.m", which is the function that takes in spatial abundance data and returns richness estimates, including Richness_raw, Chao1, GP, Chao2, ACE, S_aj2 (Jackknife_abundance), S_ij2 (Jackknife_incidence), Richness_apx (the proposed estimator 𝛺), Richness_taylor (the proposed Taylor expansion version 𝛺T), and Richness_0 (the proposed one-term approximate version 𝛺o). The script also calls "RichnessEsts_ideal.m", which is the analogous function to the previous that uses unobserved true quantities in the simulated data to return an idealized richness estimate (𝛺T, the theoretical best version of 𝛺o).

2. Run the "runEmpiricalRichnessCorrection.m" script to analyze a real multi-year, multi-site seaweed dataset. The script first loads the empirical data contained in "Data_all.mat". The script calls the function "bootRichnessEsts_all.m" to obtain bootstrapped samples of all estimators mentioned above assuming dependencies of individuals within species and sites, and otherwise calls the function "bootRichnessEsts.m" when only Richness_raw, Chao1, and Richness_apx are returned. Both functions calls "RichnessEsts.m" to estimate richness for every bootstrapped dataset. The bootstrap functions can optionally be replaced by "bootRichnessEsts_all_indepIndivAssumed.m" if bootstrapping should be done on individuals assuming they are independent. Plotting uses "boundedline.m" and associated files.

3. Analogous functions to "RichnessEsts.m" and "bootRichnessEsts.m" are provided in an R-package "Richness" at https://github.com/EWTekwa/RichnessEstimator/Richness.
