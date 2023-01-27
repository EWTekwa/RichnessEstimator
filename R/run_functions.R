# code to run functions with Martone seaweed data

# pacakges
library(tidyverse)
library(R.matlab)
library(vegan)
data(BCI)

# read BCI data
bci <- readMat("bci.tree_abundance.mat")
bci.abundance <- bci[[2]]


# Martone Data
comm_2012 = read.csv("community_2012.csv")
comm_2013 = read.csv("community_2013.csv")
comm_2014 = read.csv("community_2014.csv")
comm_2015 = read.csv("community_2015.csv")
comm_2016 = read.csv("community_2016.csv")
comm_2017 = read.csv("community_2017.csv")
comm_2018 = read.csv("community_2018.csv")
comm_2019 = read.csv("community_2019.csv")
comm <- list( comm_2012, comm_2013, comm_2014, comm_2015, comm_2016, comm_2017, comm_2018, comm_2019 )

# TransectAbundance = TransectAbundance[,-1]
TransectAbundancePrepare <- function(comm) comm %>% 
  pivot_longer( Acrosiphonia:Unknown.crust, names_to = "taxon", values_to = "cover" ) %>% 
  mutate( cover = ifelse( cover < 0.5 & cover > 0, 0.5, cover) ) %>% 
  mutate( abundance = cover/0.5 ) %>% 
  separate( UID, c("site","beach","zone","year","quadrat"), sep = " " ) %>% 
  group_by( site, zone, taxon ) %>% 
  summarize( abundance = round(sum(abundance)) ) %>% # sum here because we need to track singleton and doubleton occurences
  ungroup %>% 
  pivot_wider( id_cols = 1:2, names_from = "taxon", values_from = "abundance") %>% 
  select( -site, -zone ) %>% 
  as.matrix()
comms <- lapply( comm, TransectAbundancePrepare )




# Point estimates and mean states
source("richnessEstsCov.R")
# using a data.frame or matrix
richnessEstsCov( TransectAbundance = comms[[1]] )
richnessEstsCov( TransectAbundance = bci.abundance[,,1] )
# using a list
estsall <- do.call( rbind,lapply( comms, RichnessEstsCov ) )
ests <- do.call( rbind, estsall[,1] )
estslong <- ests %>% 
  mutate( id = 1:length(comms)) %>% 
  pivot_longer( Richness_raw:omega_T, names_to = "estimate", values_to = "richness" )
# using an array
estsall <- do.call( rbind, apply( bci.abundance,3, RichnessEstsCov ) )
ests <- do.call( rbind, estsall[,1] )
estslong <- ests %>% 
  mutate( id = as.numeric(gl( n = dim(bci.abundance)[3], k = 1 )) ) %>% 
  pivot_longer( Richness_raw:omega_T, names_to = "estimate", values_to = "richness" )

# Bootstrapped estimates - note that the first bootstrap is the point estimate
source("bootRichnessEsts_all.R")
# using a data.frame or matrix
boot_test <- bootRichnessEsts( TransectAbundance = comms[[1]], numBoot = 100 )
boot_test <- bootRichnessEsts( TransectAbundance = bci.abundance[,,1], numBoot = 100 )
# plot(boot_test)
# using a list
ests <- do.call( rbind,lapply( comms, bootRichnessEsts, numBoot = 100 ) )
# define numBoot
numBoot = 100
estslong <- ests %>% 
  mutate( id = as.numeric(gl( n = length(comms), k = numBoot )) ) %>% 
  pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
# using an array
ests <- do.call( rbind,apply( bci.abundance, 3, bootRichnessEsts, numBoot = 100 ) )
estslong <- ests %>% 
  mutate( id = as.numeric(gl( n = dim(bci.abundance)[3], k = numBoot )) ) %>% 
  pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )


# run the wrapper function
source("estimateRichness.R")
estimateRichness( comms[[1]], boot = T )
estimateRichness( comms[[1]], boot = F )
estimateRichness( comms[[1]], boot = F, meanStates = T, Apx_detectP_terms = T  )
estimateRichness( BCI, boot = F  )
estimateRichness( comms, boot = F )
estimateRichness( comms, boot = T )
estimateRichness( bci.abundance, boot = F )
estimateRichness( bci.abundance, boot = T )
estimateRichness( comm = comms[[1]][,1], boot = F )




#####
## Plotting
##
# plot just richness_raw, chao1, and new richness_apx
ggplot( filter(estslong, estimate %in% c('Richness_raw','Chao2','Chao1','omega_T')), aes(x=id, y=richness, color=estimate) ) +
  geom_point( alpha = 0.5 ) +
  # geom_line() + 
  geom_smooth( se = F, alpha = 0.1 ) +
  # coord_cartesian( ylim = c(0,400)) +
  scale_color_manual( values = c("blue","red","grey","black") ) +
  theme_classic()
  
ggplot( estslong, aes(x=year, y=richness, color=estimate) ) +
  geom_line() +
  theme_classic()

# calcualting confidence intervals from bootstrapped estimates
bootRichnessEsts( comms[[8]], numBoot = 25 )
ests <- lapply( comms, bootRichnessEsts, numBoot = 25 ) 
point_ests <- do.call( rbind, lapply( ests, function(z) z[[1]]) )
boots_list <- lapply( ests, function(z) z[2:8] )

CI_calc <- function( boots ){
  do.call( cbind, lapply( boots, quantile, c(0.025,0.975) ) )
}

boots_CI_long <- do.call( rbind, lapply( boots_list, CI_calc )) 
boots_CI_long$year