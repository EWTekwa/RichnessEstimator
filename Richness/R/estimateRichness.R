# a wrapper function allowing mutliple data input types and choice of wheter to output point estimes or bootstrapped estimates
# created by Matt Whalen
# last update 27 Jan 2023
estimateRichness <- function( comm, boot = T, numBoot = 100, meanStates = F, Apx_detectP_terms = F ){
  if( any(class(comm) %in% c("data.frame","matrix")) ){ # note that mean states and approximated detection probability terms are only currently available when estimating single communities (i.e., no time series or comparison of communities across space)
    if( boot == T ){
      est = bootRichnessEsts( comm )
    } else if( meanStates == T & Apx_detectP_terms == T ){
        est = richnessEstsCov( comm )
      } else if( meanStates == T & Apx_detectP_terms == F ){
        est = richnessEstsCov( comm )[1:2]
      } else if( meanStates == F & Apx_detectP_terms == T ){
        est = richnessEstsCov( comm )[c(1:3)]
      } else if( meanStates == F & Apx_detectP_terms == F ){
        est = richnessEstsCov( comm )[[1]]
      }
  } else if( class(comm) == "list" ){ # each element of the list should be a community with taxa as columns and samples as rows
    if( boot == F ){
      estall = do.call( rbind,lapply( comm, richnessEstsCov ) )
      est = do.call( rbind, estall[,1] )
    } else {
      est = do.call( rbind,lapply( comm, bootRichnessEsts, numBoot = 100 ) )
      estlong = est %>% 
        mutate( id = as.numeric(gl( n = length(comm), k = numBoot )) ) %>% 
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else  if( any(class(comm) == "array" & length(dim(comm)) > 2) ){ # arrays should be arrange with taxa as columns, samples as rows, and time/space in the 3rd dimension
    if( boot == F ){
      estall = do.call( rbind, apply( comm,3, richnessEstsCov ) )
      est = do.call( rbind, estall[,1] )
      estlong <- est %>% 
        mutate( id = as.numeric(gl( n = dim(comm)[3], k = 1 )) ) %>% 
        pivot_longer( Richness_raw:omega_T, names_to = "estimate", values_to = "richness" )
      est = estlong
    } else {
      est = do.call( rbind,apply( comm, 3, bootRichnessEsts, numBoot = 100 ) )
      estlong = est %>% 
        mutate( id = as.numeric(gl( n = dim(comm)[3], k = numBoot )) ) %>% 
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else if( !(any(class(comm) %in% c("data.frame","matrix","list","array"))) ){
    est = c("sorry, this community data does not appear to be in a supported format",
            "try coercing the data to a matrix, data.frame, list, or array with species as columns and spatial sampling units as rows" )
  }
      return(est)
} 
