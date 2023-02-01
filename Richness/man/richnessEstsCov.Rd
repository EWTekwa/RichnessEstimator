\name{richnessEstsCov}
\alias{richnessEstsCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{richnessEstsCov
%%  ~~function to do ... ~~
}
\description{
This function takes spatial abundance data and returns several richness estimates, including two defined by Tekwa et al. 2023.  In addition to richness point estimates, the function returns mean states of correction terms and approximated detection probabilities used in Tekwa's estimation method.
}
\usage{
richnessEstsCov(Community)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Community}{ Data frame or matrix of community data with species as columns and spatial sampling units as rows.}
}
\details{
% can include plenty of details about the method introduced and the other richness estimators.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of 3 elements.
\item{comp1 }{point estimates for all richness estimators.}
\item{comp2 }{mean states of the correction terms for \eqn{\Omega{T}} and \eqn{\Omega{0}}.}
\item{comp2 }{approximated detection probabilities for \eqn{\Omega{T}} and \eqn{\Omega{0}}.}
}
\references{
Tekwa, E.W., Whalen, M.W., Martone, P.T. and O'Connor, M.I. (2023) Theory and application of an improved species richness estimator. \emph{Philosophical Transactions of the Royal Society B.} \doi{}
}
\author{
E.W. Tekwa and Matt Whalen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
 \code{\link{estimateRichness}}, \code{\link{bootRichnessEsts}}.
}
\examples{
## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function, which returns a list
richnessEstsCov(BCI)

## isolate the richness estimates
richnessEstsCov(BCI)[[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
