# a wrapper function allowing mutliple data input types and choice of wheter to output point estimes or bootstrapped estimates
# created by Matt Whalen
# last update 27 Jan 2023
estimateRichness <- function( Community, boot = T, numBoot = 100, meanStates = F, Apx_detectP_terms = F ){
  if( any(class(Community) %in% c("data.frame","matrix")) ){ # note that mean states and approximated detection probability terms are only currently available when estimating single Communities (i.e., no time series or comparison of independent Communities across space)
    if( boot == T ){
      est = bootRichnessEsts( Community )
    } else if( meanStates == T & Apx_detectP_terms == T ){
        est = richnessEstsCov( Community )
      } else if( meanStates == T & Apx_detectP_terms == F ){
        est = richnessEstsCov( Community )[1:2]
      } else if( meanStates == F & Apx_detectP_terms == T ){
        est = richnessEstsCov( Community )[c(1:3)]
      } else if( meanStates == F & Apx_detectP_terms == F ){
        est = richnessEstsCov( Community )[[1]]
      }
  } else if( class(Community) == "list" ){ # each element of the list should be a Community with taxa as columns and samples as rows
    if( boot == F ){
      estall = do.call( rbind,lapply( Community, richnessEstsCov ) )
      est = do.call( rbind, estall[,1] )
    } else {
      est = do.call( rbind,lapply( Community, bootRichnessEsts, numBoot = 100 ) )
      estlong = est %>%
        mutate( id = as.numeric(gl( n = length(Community), k = numBoot )) ) %>%
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else  if( any(class(Community) == "array" & length(dim(Community)) > 2) ){ # arrays should be arrange with taxa as columns, samples as rows, and time/space in the 3rd dimension
    if( boot == F ){
      estall = do.call( rbind, apply( Community,3, richnessEstsCov ) )
      est = do.call( rbind, estall[,1] )
      estlong <- est %>%
        mutate( id = as.numeric(gl( n = dim(Community)[3], k = 1 )) ) %>%
        pivot_longer( Richness_raw:Omega_T, names_to = "estimate", values_to = "richness" )
      est = estlong
    } else {
      est = do.call( rbind,apply( Community, 3, bootRichnessEsts, numBoot = 100 ) )
      estlong = est %>%
        mutate( id = as.numeric(gl( n = dim(Community)[3], k = numBoot )) ) %>%
        pivot_longer( Richness_raw:S_ij2, names_to = "estimate", values_to = "richness" )
      est = estlong
    }
  } else if( !(any(class(Community) %in% c("data.frame","matrix","list","array"))) ){
    est = c("sorry, this Communityunity data does not appear to be in a supported format",
            "try coercing the data to a matrix, data.frame, list, or array with species as columns and spatial sampling units as rows" )
  }
      return(est)
}
