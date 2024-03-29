# RichnessEstimator
Code and data for:
An improved species richness estimator using spatial abundance data
Eden W. Tekwa1,2*, Matthew A. Whalen2,3,4, Patrick T. Martone3, Mary I. O’Connor1

1Department of Zoology, University of British Columbia, Vancouver, BC, Canada.
2Hakai Institute, Heriot Bay, BC, Canada
3Department of Botany, University of British Columbia, Vancouver, BC, Canada.
4Department of Biology, Virginia State University, Petersburg, VA, USA
*Corresponding author email: ewtekwa@gmail.com, ORCID: 0000-0003-2971-6128
Emails for MW: mawhal@gmail.com, PM: patrick.martone@botany.ubc.ca, MO: oconnor@zoology.ubc.ca.


Note: The following documentation is for code and data used to generate the results in the paper cited above. For the R-package with the richness estimator functions, please see the subfolder /Richness (https://github.com/EWTekwa/RichnessEstimator/Richness).

Main Matlab files:
Data_all.mat
RichnessEsts.m
RichnessEsts_ideal.m
biodiversitySamplingSim.m
bootRichnessEsts.m
bootRichnessEsts_all.m
bootRichnessEsts_all_indepIndivAssumed.m
runEmpiricalRichnessCorrection.m

Third Party Matlab files are need for plotting and copies are included in this repository (by Kelly Kearney 2010, available at https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m):
boundedline.m
catuneven.m
inpaint_nans.m
outlinebounds.m
singlepatch.m

R files:
RichnessEsts.R
bootRichnessEsts.R

1. Run the "biodiversitySamplingSim.m" script to evaluate species richness estimators using simulated scenarios. The script calls "RichnessEsts.m", which is the function that takes in spatial abundance data and returns richness estimates, including Richness_raw, Chao1, GP, Chao2, ACE, S_aj2 (Jackknife_abundance), S_ij2 (Jackknife_incidence), Richness_apx (the proposed estimator 𝛺), Richness_taylor (the proposed Taylor expansion version 𝛺T), and Richness_0 (the proposed one-term approximate version 𝛺o). The script also calls "RichnessEsts_ideal.m", which is the analogous function to the previous that uses unobserved true quantities in the simulated data to return idealized richness estimates (𝛺c, 𝛺Tc, and 𝛺oc corresponding to the operational 𝛺, 𝛺T, and 𝛺o). Simulation results used in manuscript are contained in "runEmpiricalRichnessCorrection_BCI.mat" for reference.

2. Run the "runEmpiricalRichnessCorrection_BCI.m" script to analyze the Barro Colorado Island tree census dataset. The script first loads the empirical data contained in "bci.tree_abundance.mat". Full data and spatial subsampling and local downsampling experiments are analyzed. The script calls the function "bootRichnessEsts.m" to obtain bootstrapped samples of Chao1, Chao2, and 𝛺 estimates assuming dependencies of individuals within species and sites. The function calls "RichnessEsts.m" to obtain point estimates. Analytical results used in manuscript are contained in "runEmpiricalRichnessCorrection_BCI.mat" for reference.

3. Run the "runEmpiricalRichnessCorrection_Seaweed.m" script to analyze the BC seaweed survey dataset. The script first loads the empirical data contained in "BC.seaweed_cover.mat". Full data and spatial subsampling and local downsampling experiments are analyzed. The script calls the function "bootRichnessEsts.m" to obtain bootstrapped samples of Chao1, Chao2, and 𝛺 estimates assuming dependencies of individuals within species and sites. The function calls "RichnessEsts.m" to obtain point estimates. Analytical results used in manuscript are contained in "runEmpiricalRichnessCorrection_Seaweed.mat" for reference.

4. Analogous functions to "RichnessEsts.m" and "bootRichnessEsts.m" are provided in an R-package "Richness" at https://github.com/EWTekwa/RichnessEstimator/Richness.
