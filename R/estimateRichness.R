richest <- function( comm, boot = T, numBoot = 100, meanStates = F, Apx_detectP_terms = F ){
  if( any(class(comm) %in% c("data.frame","matrix")) ){ # note that mean states and approximated detection probability terms are only currently available when estimating single communities (i.e., no time series or comparison of communities across space)
    if( boot == T ){
      ests = bootRichnessEsts( comm )
    } else if( meanStates == T & Apx_detectP_terms == T ){
        ests = RichnessEstsCov( comm )
      } else if( meanStates == T & Apx_detectP_terms == F ){
        ests = RichnessEstsCov( comm )[1:2]
      } else if( meanStates == F & Apx_detectP_terms == T ){
        ests = RichnessEstsCov( comm )[c(1:3)]
      } else if( meanStates == F & Apx_detectP_terms == F ){
        ests = RichnessEstsCov( comm )[[1]]
      }
  } else if( class(comm) == "list" ){ # each element of the list should be a community with taxa as columns and samples as rows
    if( boot == F ){
      estsall = do.call( rbind,lapply( comm, RichnessEstsCov ) )
      ests = do.call( rbind, estsall[,1] )
    } else {
      ests = do.call( rbind,lapply( comm, bootRichnessEsts, numBoot = 100 ) )
      estslong = ests %>% 
        mutate( id = as.numeric(gl( n = length(comm), k = numBoot )) ) %>% 
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      ests = estslong
    }
  } else  if( any(class(comm) == "array" & length(dim(comm)) > 2) ){ # arrays should be arrange with taxa as columns, samples as rows, and time/space in the 3rd dimension
    if( boot == F ){
      estsall = do.call( rbind, apply( comm,3, RichnessEstsCov ) )
      ests = do.call( rbind, estsall[,1] )
      estslong <- ests %>% 
        mutate( id = as.numeric(gl( n = dim(comm)[3], k = 1 )) ) %>% 
        pivot_longer( Richness_raw:omega_T, names_to = "estimate", values_to = "richness" )
      ests = estslong
    } else {
      ests = do.call( rbind,apply( comm, 3, bootRichnessEsts, numBoot = 100 ) )
      estslong = ests %>% 
        mutate( id = as.numeric(gl( n = dim(comm)[3], k = numBoot )) ) %>% 
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      ests = estslong
    }
  } else if( !(any(class(comm) %in% c("data.frame","matrix","list","array"))) ){
    ests = "sorry, this community data does not appear to be in the correct format"
  }
      return(ests)
} 
