\name{bootstrapClusterClass}
\alias{bootstrapClusterClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap analysis of clustering algorithms
%%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Performs classification with \code{kmeans}, \code{hc} or \code{pam},
  then estimates uncertainty via  bootstrap resampling. 

  Bootstrap resampling is stratified such that the proportion of short
  and long latency in each bootstrap remains constant.
  
}
\usage{
bootstrapClusterClass(latency, its = 1000, restratify = FALSE,
algorithm=c("kmeans", "hc", "pam"), emModelNames = "E") 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{latency}{
    %%     ~~Describe \code{latency} here~~
    A vector of average latency values.
  }
  \item{its}{
    Number of bootstrap iterations.
    %%     ~~Describe \code{its} here~~
  }
  \item{restratify}{
    If TRUE, classify the data with \code{Mclust} and resample the
    original classification probabilistically. By default, the
    proportion of SL / LL in each bootstrap is fixed, based on the
    original classification (the prior). Optionally, we can resample the
    groups by using the prior probabilities from an EM clustering.  
    
    This allows us to vary the stratification according to uncertainty
    in the prior classification.
    
    For example, consider the subjects (A,B,C,D,E) with prior
    probabilities (0, 0.05, 0.45, 0.95, 1). This contains 3 SL and 2 LL.
    But if we resample the groups using these probabilities, we will
    classify A as LL with probability 0, B with probability 0.05, C with
    probability 0.45, and so on.
    
    %%     ~~Describe \code{restratify} here~~
  }
  \item{algorithm}{
    The clustering algorithm. By default, calls \code{kmeans} from the
    \code{stats} package. Alternatively, \code{pam} (\code{cluster}) or
    \code{hc} (\code{mclust}) can be used.
  }
  \item{emModelNames}{
    Only used if \code{restratify = TRUE}. Passed to Mclust to allow or
    disallow different models. The default "E" means that the Gaussian
    distributions inside a particular model will have equal standard
    deviation. If \code{modelNames = "V"}, then the SD of each Gaussian
    is variable. If \code{modelNames = c("E", "V")}, then \code{Mclust}
    will choose the appropriate model.
  }

}
\details{

  This function is used to classify average latency scores into two
  groups, the "short latency" (SL) with low stress resilience and the "long
  latency" (LL) with high stress resilience.
  
  The bootstrap resampling is done by sampling, with replacement, from
  the SL and LL groups defined by the initial call to the clustering
  algorithm on the original data. The proportion of SL and LL in each
  bootstrap remains fixed, unless restratify=TRUE, in which case an EM
  algorithm is used to define prior probabilities, and the
  stratification is recomputed by sampling these priors at each
  bootstrap.

  The classification boundary for each bootstrap is the mid-point
  between the two centroids / means. This boundary is used to
  the classify the original latency data. 

  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  A list with the following components

% stuff <- list(bootProbLL = bootProbLL, its = its, latency = latency, kBoundary = kBoundary, k_boundary_boot = k_boundary_boot, kCenters = kCenters, k_class_boot = k_class_boot, kClusters = kClusters)
 
  \item{bootProbLL }{The probability that a subject is classified as LL,
    defined as the number of times this subject was classified LL over all
    bootstraps.}
  \item{its }{The number of bootstrap iterations.}
  \item{latency }{The original latency data used for the original
    classification and resampled for the bootstrap.}
  \item{boundary }{Approximate cluster boundary between SL and LL means
    (from kmeans or hc) or medoids (from pam) computed on the original
    data.}
  \item{boundary_boot}{Approximate boundary point from all
    bootstraps.}
  \item{centers }{Cluster centers or medoids fit to the original data.}
  \item{class_boot }{A matrix containing the classification of sample
    data at each bootstrap.} 
  \item{clusters }{The classification of the original data: 1 (SL), 2
    (LL).}
 
 The ordering of the per-subject values (such as bootProbLL) is the same
 as in the latency vector.

 The cluster boundaries are defined as halfway between the two cluster
 means at each bootstrap.

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Philip A Cook <cookpa@pennmedicine.upenn.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{plotBootstrapClusterClass}}
  \code{\link{kmeans}}
  \code{\link{pam}}
  \code{\link{hc}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

set.seed(20140123)

sl <- rnorm(60, 200, 100)
sl[which(sl < 0)] <- 0

ll <- rnorm(15, 600, 100)
ll[which(ll > 900)] <- 900

boot = bootstrapClusterClass(c(sl, ll), its = 100, algorithm = "pam")


## Not run
\dontrun{plotBootstrapClusterClass(boot)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bootstrap }
\keyword{ ~clustering }% __ONLY ONE__ keyword per line
