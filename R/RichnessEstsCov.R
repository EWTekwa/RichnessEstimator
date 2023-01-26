RichnessEstsCov <- function(TransectAbundance){

# RichnessEsts.m
# Ed Tekwa Feb 8, 2022 
# RichnessEstsCov.R
# Matt Whalen rewrote original matlab script for R - Mar 27, 2022 
# function several estimates of richness:
  # - raw richness
  # - newly proposed method 
  # - Chao1 
  # - Chao2
  # - Abundance-based coverage estimator (ACE)
  # - Jackknife abundance estimator
  # - Jackknife incidence estimator
# The new method calculates richness based on the spatial TransectAbundance data: rows=transects, columns=species,
# values=individual counts
# Whalen update Oct 31, 2022
  # - return mean states
  # - fix a few errors
  # - return correction terms
  # - two versions of the estimate 
    # - Omega_not - uses average detection probabilities
    # - Omega_T - uses Taylor approximation



# define variance function with normalizatin of n rather than n-1
varn <- function(x) mean((x-mean(x))^2)

# # write observational process model
# Dix = (1-exp(-C*nm))*P # local detection probability of species i at sampled site x
# Dis = 1-(1-Dix)^k # detection probability of species i across all sampled sites x in community s
Dis <- expression( 1-(1-((1-exp(-C*nm))*P))^k )
# second-order derivatives 
d2Dis_dnm2 <- D(D( Dis, "nm" ), "nm")
d2Dis_dC2  <- D(D( Dis, "C" ), "C")
d2Dis_dP2  <- D(D( Dis, "P" ), "P")
d2Dis_dnmC  <- D(D( Dis, "nm" ), "C")
d2Dis_dnmP  <- D(D( Dis, "nm" ), "P")
d2Dis_dCP   <- D(D( Dis, "C" ), "P")


numTrans =  nrow(TransectAbundance) # get number of transects
Richness_raw = sum(colSums(TransectAbundance)>0) # get raw richness
C_detected = rep(0,ncol(TransectAbundance)) # array of zeros to record clustering for each species
P_detected = rep(0,ncol(TransectAbundance)) # array of zeros to record occupancy for each species


for( species in 1:ncol(TransectAbundance) ) {
  if( numTrans == 1 ){
  C_detected[species] = 1+1/mean(TransectAbundance[,species]) # if only one sampled site, estimate spatial variance as poisson mean
  } else {
  C_detected[species] = 1+varn(TransectAbundance[,species])/(mean(TransectAbundance[,species])^2) # spatial variance normalized by N (not the normal N-1)
  }
  P_detected[species] = sum(TransectAbundance[,species] > 0)/numTrans # occupancy as number of transects occupied divided by number of transects
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

# calculate singletons, doubletons, and the Chao1 richness estimator
f1 = sum(colSums(TransectAbundance) == 1) # number of singleton species
f2 = sum(colSums(TransectAbundance) == 2) # number of doubleton species
Chao1 = Richness_raw+f1*(f1-1)/(2*(f2+1)) # Chao1 richness estimator

# compute Chao2 estimate
q1 = sum(colSums(TransectAbundance >= 1) == 1) # number of species occurring in one sample only
q2 = sum(colSums(TransectAbundance >= 1) == 2) # number of species occurring in two samples only
m = sum(colSums(TransectAbundance >= 1)) # total number of samples
Chao2 = Richness_raw+((m-1)/m)*q1*(q1-1)/(2*(q2+1)) # Chao2 richness estimator

# compute abundance-based coverage estimator (ACE)
S_rare = sum(colSums(TransectAbundance) <= 10 & colSums(TransectAbundance) > 0 ) # number of rare species (<=10 individuals)
S_abund = sum(colSums(TransectAbundance) > 10) # number of abundant species (>10 individuals)
n_rare = sum(TransectAbundance[,colSums(TransectAbundance) <= 10 & colSums(TransectAbundance) > 0]) # total number of individuals in rare species
C_ACE = 1-f1/n_rare # sample coverage estimate
wsumfa = 0
for( a in 1:10 ){
  wsumfa = wsumfa + a*(a-1)*sum(colSums(TransectAbundance) == a)
}
V2 = max(((S_rare/C_ACE)*wsumfa/(n_rare*(n_rare - 1)) - 1),0) # coefficient of variation
if( C_ACE > 0 ){
  ACE = S_abund+S_rare/C_ACE+(f1/C_ACE)*V2
} else {
  ACE = Chao1
}

# compute jackknife abundance estimator
S_aj2 = Richness_raw+2*f1-f2

# compute jackknife incidence estimator
S_ij2 = Richness_raw+(q1*(2*m-3)/m-q2*((m-2)^2)/(m*(m-1)))



# compute correction terms for proposed approximation method
k <- numTrans
C <- mean_C_detected
nm <- mean_n_m_detected
P <- mean_P_detected
Apx_detectP_terms <- c( eval(Dis),
  eval(d2Dis_dnm2)*var_n_m_detected/2,
  eval(d2Dis_dC2)*var_C_detected/2,
  eval(d2Dis_dP2)*var_P_detected/2,
  eval(d2Dis_dnmC)*cov_nm_C_detected,
  eval(d2Dis_dnmP)*cov_nm_P_detected,
  eval(d2Dis_dCP)*cov_C_P_detected ) # Approximated detection probability in community

# check for correction term relative to some threshold
if( sum(Apx_detectP_terms, na.rm = T) > 0.1 ){ # if sum of correction terms is positive and greater than a threshold
  Ds_apx = sum(Apx_detectP_terms, na.rm = T) # use full correction
} else {
  Ds_apx = Apx_detectP_terms[1] # else, use 0th order correction
}
if( sum(Apx_detectP_terms, na.rm = T) > 1 ){
  Ds_apx = 1
}

omega_T <- Richness_raw/Ds_apx
 

# save the mean states of correction terms to return with the function
meanStates <- c(mean_n_m_detected,mean_C_detected,mean_P_detected,
                var_n_m_detected,var_C_detected,var_P_detected)


# compute Omega_o
C <- na.omit(C_detected)
nm <- n_m_detected
P <- na.omit(P_detected)
Ds_means <- eval(Dis)
Ds_mean <- mean(Ds_means, na.rm = T)
omega_0 <- Richness_raw/Ds_mean


return(list(data.frame(Richness_raw=Richness_raw,
                  Chao1=Chao1, Chao2=Chao2, 
                  ACE=ACE, S_aj2=S_aj2, S_ij2=S_ij2,
                  omega_0=omega_0,
                  omega_T=omega_T),meanStates=meanStates,
       Apx_detectP_terms=Apx_detectP_terms))
}
