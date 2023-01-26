bootRichnessEsts <- function(TransectAbundance, numBoot = 100){

# bootRichnessEsts.m
# Ed Tekwa Feb 8, 2022
# bootRichnessEsts.R
# Matt Whalen rewrote original matlab script for R - Mar 27, 2022 
# function several estimates of richness:
#  - raw richness
#  - newly proposed method 
#  - Chao1 
#  - Chao2
#  - Abundance-based coverage estimator (ACE)
#  - Jackknife abundance estimator
#  - Jackknife incidence estimator
# Point estimates are calculated using script "RichnessEstsCov.R"
# This script calculates bootstrapped confidence intervals
# the spatial TransectAbundance data: rows = transects, columns = species,
#                                     values = individual counts
# Whalen update 16 January 2023
  

# Get point estimates from original dataset
TransectAbundance = TransectAbundance[ ,colSums(TransectAbundance) > 0 ] # remove empty species columns
# run function to calculate all point estimates - ignore correction terms states and approximated detection probabilities
pointests = RichnessEstsCov( TransectAbundance )[[1]]


# store bootstrapped estimates for raw richness and each estimator
expectedRichness_raw = rep(0,numBoot)    # raw
expectedChao1 = rep(0,numBoot)           # Chao1
expectedChao2 = rep(0,numBoot)           # Chao2
expectedACE   = rep(0,numBoot)           # ACE
expectedS_aj2 = rep(0,numBoot)           # Jackknife (abundance)
expectedS_ij2 = rep(0,numBoot)           # Jackknife (incidence)
expectedRichness_omega_0 = rep(0,numBoot)    # approximate Richness
expectedRichness_omega_T = rep(0,numBoot) # approximated Richness using Taylor expansion

# use point estimates from dataset as first bootstrap estimate
expectedRichness_raw[1] = pointests$Richness_raw
expectedChao1[1] = pointests$Chao1
expectedChao2[1] = pointests$Chao2
expectedACE[1] = pointests$ACE
expectedS_aj2[1] = pointests$S_aj2
expectedS_ij2[1] = pointests$S_ij2
expectedRichness_omega_0[1] = pointests$omega_0
expectedRichness_omega_T[1] = pointests$omega_T


# define variance function with normalizatin of n rather than n-1
varn <- function(x) mean((x-mean(x))^2)

# write observational process model
Dis <- expression( 1-(1-((1-exp(-C*nm))*P))^k )
# second-order derivatives 
d2Dis_dnm2 <- D(D( Dis, "nm" ), "nm")
d2Dis_dC2  <- D(D( Dis, "C" ), "C")
d2Dis_dP2  <- D(D( Dis, "P" ), "P")
d2Dis_dnmC  <- D(D( Dis, "nm" ), "C")
d2Dis_dnmP  <- D(D( Dis, "nm" ), "P")
d2Dis_dCP   <- D(D( Dis, "C" ), "P")


# get the number of sampling units (e.g., transects or quadrats)
numSamplingUnit = nrow(TransectAbundance)


# bootstrapping - scramble sampling units and resample

for( resample in 2:numBoot ){
  #  scramble sampling units (e.g., transects) first
    sampleSet_init = TransectAbundance[ sample(1:numSamplingUnit,numSamplingUnit,replace = T), ] 
    
    # using unit-scrambled data, next scramble the species
    sampleSetRaw = sampleSet_init[ ,sample(1:pointests$Richness_raw,pointests$Richness_raw,replace = T) ]
    sampleSetChao1 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Chao1),replace = T) ]
    sampleSetChao2 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$Chao2),replace = T) ]
    sampleSetACE = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$ACE),replace = T) ]
    sampleSetJK_a = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$S_aj2),replace = T) ]
    sampleSetJK_i = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$S_ij2),replace = T) ]
    sampleSetomega_0 = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$omega_0),replace = T) ]
    sampleSetomega_T = sampleSet_init[ ,sample(1:pointests$Richness_raw,round(pointests$omega_T),replace = T) ]
    
    # take out empty species columns
    sampleSetRaw = sampleSetRaw[ , colSums(sampleSetRaw) > 0 ] 
    sampleSetChao1 = sampleSetChao1[ , colSums(sampleSetChao1) > 0 ] 
    sampleSetChao2 = sampleSetChao2[ , colSums(sampleSetChao2) > 0 ]
    sampleSetACE = sampleSetACE[ , colSums(sampleSetACE) > 0 ] 
    sampleSetJK_a = sampleSetJK_a[ , colSums(sampleSetJK_a) > 0 ]
    sampleSetJK_i = sampleSetJK_i[ , colSums(sampleSetJK_i) > 0 ]
    sampleSetomega_0 = sampleSetomega_0[ , colSums(sampleSetomega_0) > 0 ]
    sampleSetomega_T = sampleSetomega_T[ , colSums(sampleSetomega_T) > 0 ]
    
    # compute raw estimate
    Richness_raw_boot = sum(colSums(sampleSetRaw)>0)
    
    # compute Chao1 estimate
    sampleSetChao1 = ceiling(sampleSetChao1) # convert data to discrete count
    raw = sum(colSums(sampleSetChao1) > 0)
    f1 = sum(colSums(sampleSetChao1) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetChao1) == 2) # number of doubleton species
    Chao1_boot = raw+f1*(f1-1)/(2*(f2+1)) # Chao1 richness estimator
    
    # compute Chao2 estimate
    raw = sum(colSums(sampleSetChao2) > 0)
    q1 = sum(colSums(sampleSetChao2 > 0) == 1) # number of species occurring in one transect only
    q2 = sum(colSums(sampleSetChao2 > 0) == 2) # number of species occurring in two transect only
    m = sum(colSums(sampleSetChao2 > 0))
    Chao2_boot = raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)) # Chao2 richness estimator
    
    # compute abundance-based coverage estimator (ACE)
    sampleSetACE = ceiling(sampleSetACE) # convert data to discrete count
    raw = sum(colSums(sampleSetACE) > 0)
    f1 = sum(colSums(sampleSetACE) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetACE) == 2) # number of doubleton species
    S_rare = sum(colSums(sampleSetACE) <= 10) # number of rare species (< = 10 individuals)
    S_abund = sum(colSums(sampleSetACE) > 10) # number of rare species (>10 individuals)
    n_rare = sum(sampleSetACE[,colSums(sampleSetACE) <= 10]) # total number of individuals in rare species
    C_ACE = 1-f1/n_rare # sample coverage estimate
    wsumfa = 0
    for( a in 1:10 ){
      wsumfa = wsumfa + a*(a-1)*sum(colSums(TransectAbundance) == a)
    }
    V2 = max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare - 1)) - 1),0) # coefficient of variation
    if( C_ACE > 0 ){
      ACE_boot = S_abund+S_rare/C_ACE+(f1/C_ACE)*V2
    } else {
      ACE_boot = Chao1
    }
    
    # compute jackknife abundance estimator
    sampleSetACE = ceiling(sampleSetACE) # convert data to discrete count
    raw = sum(colSums(sampleSetJK_a) > 0)
    f1 = sum(colSums(sampleSetJK_a) == 1) # number of singleton species
    f2 = sum(colSums(sampleSetJK_a) == 2) # number of doubleton species
    S_aj2_boot = raw+2*f1-f2
    
    # compute jackknife incidence estimator
    raw = sum(colSums(sampleSetJK_i) > 0)
    q1 = sum(colSums(sampleSetJK_i > 0) == 1) # number of species occurring in one transect only
    q2 = sum(colSums(sampleSetJK_i > 0) == 2) # number of species occurring in two transect only
    m = sum(colSums(sampleSetJK_i > 0))
    S_ij2_boot = raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)))
    
    # compute correction terms for proposed approximation method
    raw = sum(colSums(sampleSetomega_0) > 0)
    C_detected = rep(0,ncol(TransectAbundance)) # array of zeros to record clustering for each species
    P_detected = rep(0,ncol(TransectAbundance)) # array of zeros to record occupancy for each species
    for( species in 1:ncol(TransectAbundance) ) {
      if( numSamplingUnit == 1 ){
        C_detected[species] = 1+1/mean(TransectAbundance[,species]) # if only one sampled site, estimate spatial variance as poisson mean
      } else {
        C_detected[species] = 1+varn(TransectAbundance[,species])/(mean(TransectAbundance[,species])^2) # spatial variance normalized by N (not the normal N-1)
      }
      P_detected[species] = sum(TransectAbundance[,species] > 0)/numSamplingUnit # occupancy as number of transects occupied divided by number of transects
    }
    C_detected[C_detected==Inf] = NA
    P_detected[P_detected==0] = NA
    # compute on full dataset - means, variances, and covariances
    mean_C_detected = mean(C_detected, na.rm = TRUE )
    var_C_detected  = var(C_detected, na.rm = TRUE )
    mean_P_detected = mean(P_detected[P_detected>0], na.rm = TRUE )
    var_P_detected  = var(P_detected[P_detected>0], na.rm = TRUE )
    n_m_detected = colMeans(TransectAbundance[,colSums(TransectAbundance)>0])
    mean_n_m_detected = mean(n_m_detected)
    var_n_m_detected  = var(colMeans(TransectAbundance[,colSums(TransectAbundance)>0]))
    cov_nm_C_detected = cov( data.frame(x = colMeans(TransectAbundance[,colSums(TransectAbundance)>0]), y = na.omit(C_detected)), use = "complete.obs")
    cov_nm_P_detected = cov( data.frame(x = colMeans(TransectAbundance[,colSums(TransectAbundance)>0]), y = na.omit(P_detected)), use = "complete.obs")
    cov_C_P_detected = cov( data.frame(x = C_detected, y = P_detected), use = "complete.obs")
    if( length(cov_nm_P_detected[,2])>1 ){
      cov_nm_C_detected = cov_nm_C_detected[1,2]
      cov_nm_P_detected = cov_nm_P_detected[1,2]
      cov_C_P_detected = cov_C_P_detected[1,2]
    } else {
      cov_nm_C_detected = 0
      cov_nm_P_detected = 0
      cov_C_P_detected = 0
    }
    
    k = numSamplingUnit
    # k = 1
    C = mean_C_detected
    nm = mean_n_m_detected
    P = mean_P_detected
    # add "eval" before "(" below if using ke in place of k
    Apx_detectP_terms = (c( eval(Dis),
                            eval(d2Dis_dnm2)*var_n_m_detected/2,
                            eval(d2Dis_dC2)*var_C_detected/2,
                            eval(d2Dis_dP2)*var_P_detected/2,
                            eval(d2Dis_dnmC)*cov_nm_C_detected,
                            eval(d2Dis_dnmP)*cov_nm_P_detected,
                            eval(d2Dis_dCP)*cov_C_P_detected )) # Approximated detection probability in community
    
    # check for correction term relative to some threshold
    if( sum(Apx_detectP_terms, na.rm = T) > 0.1 ){ # if sum of correction terms is positive and greater than a threshold
      Ds_apx = sum(Apx_detectP_terms, na.rm = T) # use full correction
    } else {
      Ds_apx = Apx_detectP_terms[1] # else, use 0th order correction
    }
    if( sum(Apx_detectP_terms, na.rm = T) > 1 ){
      Ds_apx = 1
    }
    Richness_omega_T_boot = raw/Ds_apx
    
    # compute omega_o
    C = na.omit(C_detected)
    nm = n_m_detected
    P = na.omit(P_detected)
    Ds_means = eval(Dis)
    Ds_mean = mean(Ds_means, na.rm = T)
    Richness_omega_0_boot = raw/Ds_mean
    
    expectedRichness_raw[resample] = Richness_raw_boot
    expectedChao1[resample] = Chao1_boot
    expectedChao2[resample] = Chao2_boot
    expectedACE[resample] = ACE_boot
    expectedS_aj2[resample] = S_aj2_boot
    expectedS_ij2[resample] = S_ij2_boot
    expectedRichness_omega_0[resample] = Richness_omega_0_boot
    expectedRichness_omega_T[resample] = Richness_omega_T_boot
}




# return(list(as.list(RichnessEstsCov(TransectAbundance)),expectedRichness_raw=expectedRichness_raw,expectedRichness_apx=expectedRichness_apx,
#             expectedChao1=expectedChao1,expectedChao2=expectedChao2,
#             expectedACE=expectedACE,expectedS_aj2=expectedS_aj2,expectedS_ij2=expectedS_ij2))
return(data.frame(expectedRichness_raw = expectedRichness_raw,
            expectedRichness_omega_T = expectedRichness_omega_T,
            expectedRichness_omega_0 = expectedRichness_omega_0,
            expectedChao1 = expectedChao1, expectedChao2 = expectedChao2,
            expectedACE = expectedACE, expectedS_aj2 = expectedS_aj2, expectedS_ij2 = expectedS_ij2))
}
