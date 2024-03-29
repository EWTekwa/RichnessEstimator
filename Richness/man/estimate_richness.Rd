\name{estimate_richness}
\alias{estimate_richness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Richness in Community Data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper function that runs \code{RichnessEsts} or  \code{bootRichnessEsts} using a variety of data input types, including time series community data.
}
\usage{
estimate_richness(Community, boot = FALSE, numBoot = 100,
                             meanStates = FALSE, Apx_detectP_terms = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Community}{ Community data with species as columns and spatial sampling units as rows. Data can take the form of a data frame, matrix, list, or array of 3 dimensions.  List and array inputs are useful when calculating richness for community time series or for different communities who richness is to be compared. }
  \item{boot}{ logical. Should bootstrapped samples be computed? }
  \item{numBoot}{ The number of bootstrapped samples to return. }
  \item{meanStates}{ logical. Should mean states of correction terms bootstrapped samples be returned? Only used when \code{boot = FALSE} }
  \item{Apx_detectP_terms}{ logical. Should approximated detection probabilities be returned? Only used when \code{boot = FALSE} }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Tekwa, E.W., Whalen, M.W., Martone, P.T. and O'Connor, M.I. (2023) Theory and application of an improved species richness estimator. in revision.
%% \emph{Philosophical Transactions of the Royal Society B.} \doi{}
}
\author{
E.W. Tekwa and Matt Whalen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{RichnessEsts}}, \code{\link{bootRichnessEsts}}.
}

\examples{
## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function for a single community data frame without bootstrapping
# point estimates only
estimate_richness(BCI)
# with mean states for correction terms and approximated detection probabilities
estimate_richness(BCI, meanStates = TRUE, Apx_detectP_terms = TRUE)

## run the function for a single community data frame with bootstrapping
estimate_richness(BCI, boot = TRUE, numBoot = 10)

## run the function on a list
# Note: this may not be the most sensible comparison
require(vegan)
data("pyrifos")
week <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
communitylist <- split(pyrifos, f = week)
estimate_richness(communitylist)
# bootstrapping on a list returns a tidy data frame
require(tidyverse)
estimate_richness(communitylist, boot = TRUE, numBoot = 10)
}

