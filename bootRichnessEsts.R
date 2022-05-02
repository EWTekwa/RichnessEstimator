bootRichnessEsts <- function(TransectAbundance,numBoot){

# bootRichnessEsts.m
# Ed Tekwa Feb 8, 2022
# bootRichnessEsts.R
# Matt Whalen rewrote original matlab script for R - Mar 27, 2022 
# function several estimates of richness:
  # - raw richness
  # - newly proposed method 
  # - Chao1 
  # - Chao2
  # - Abundance-based coverage estimator (ACE)
  # - Jackknife abundance estimator
  # - Jackknife incidence estimator
# Point estimates are calculated using script "RichnessEstsCov.R"
# This script calculates bootstrapped confidence intervals
  # the spatial TransectAbundance data: rows = transects, columns = species,
  #                                     values = individual counts

expectedRichness_raw = rep(0,numBoot)
expectedRichness_apx = rep(0,numBoot)
expectedChao1 = rep(0,numBoot)
expectedChao2 = rep(0,numBoot)
expectedACE   = rep(0,numBoot)
expectedS_aj2 = rep(0,numBoot)
expectedS_ij2 = rep(0,numBoot)

# RichnessEstsCov(TransectAbundance) # get point estimates from original dataset

maxSampleRichness = ncol(TransectAbundance)
for( resample in 1:numBoot) {
  sampleSetApx = TransectAbundance[,sample(1:maxSampleRichness,maxSampleRichness,replace = T)] # resample species with replacement
  for( j in 1:ncol(sampleSetApx) ){ # for each species
    sampleSetApx[,j] = sampleSetApx[ sample(1:nrow(sampleSetApx),nrow(sampleSetApx),replace = T),j ] # resample with replacement transects for each species
  }
  boots <- RichnessEstsCov(sampleSetApx)
  expectedRichness_raw[resample] = boots$Richness_raw
  expectedRichness_apx[resample] =  boots$Richness_apx
  expectedChao1[resample] =  boots$Chao1
  expectedChao2[resample] =  boots$Chao2
  expectedACE[resample] =  boots$ACE
  expectedS_aj2[resample] =  boots$S_aj2
  expectedS_ij2[resample] =  boots$S_ij2
}


return(list(as.list(RichnessEstsCov(TransectAbundance)),expectedRichness_raw=expectedRichness_raw,expectedRichness_apx=expectedRichness_apx,
            expectedChao1=expectedChao1,expectedChao2=expectedChao2,
            expectedACE=expectedACE,expectedS_aj2=expectedS_aj2,expectedS_ij2=expectedS_ij2))
}
