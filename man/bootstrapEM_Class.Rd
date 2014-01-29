\name{bootstrapEM_Class}
\alias{bootstrapEM_Class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrap analysis of EM clustering
%%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Performs EM classification with \code{mclust}, then estimates uncertainty via
  bootstrap resampling. 

  The classification model is a Gaussian mixture model with two
  components and equal variance. Bootstrap resampling is stratified such
  that the proportion of short and long latency in each bootstrap
  remains constant.
}
\usage{
bootstrapEM_Class(latency, its = 1000, restratify = FALSE, modelNames = "E")
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
  If TRUE, resample the original classification probabilistically. By
  default, the proportion of SL / LL in each bootstrap is fixed, based
  on the original classification (the prior). Optionally, we can
  resample the groups by using the prior probabilities.
  
  This allows us to vary the stratification according to uncertainty in
  the prior classification.
  
  For example, consider the subjects (A,B,C,D,E) with prior
  probabilities (0, 0.05, 0.45, 0.95, 1). This contains 3 SL and 2 LL.
  But if we resample the groups using these probabilities, we will
  classify A as LL with probability 0, B with probability 0.05, C with
  probability 0.45, and so on.
  
%%     ~~Describe \code{restratify} here~~
}
\item{modelNames}{
  Passed to \code{Mclust} to allow or disallow different models. The
  default "E" means that the Gaussian distributions inside a particular
  model will have equal standard deviation. The SD will be different for
  each bootstrap, and the value for each bootstrap is recorded in the
  returned vector \code{em_std}. If \code{modelNames = "V"}, then the SD
  of each Gaussian is variable, and \code{em_std} contains two values
  for each bootstrap. If \code{modelNames = c("E", "V")}, then
  \code{Mclust} will choose the appropriate model for each bootstrap.

  Equal variance models imply a monotonically increasing probablity of
  LL classification with increasing latency. With unequal variance, the
  tails of the distributions might be quite different, lead to
  nonsensical results such as increasing probability of LL at very low 
  latencies. Be careful about interpreting the results if you allow
  unequal variance.
}
}
\details{

  This function is used to classify average latency scores into two
  groups, the "short latency" (SL) with low stress resilience and the "long
  latency" (LL) with high stress resilience.
  
  The bootstrap resampling is done by sampling, with replacement, from
  the SL and LL groups defined by the initial call to Mclust on the
  original data. The proportion of SL and LL in each bootstrap remains
  fixed, unless restratify=TRUE.

  The model for each bootstrap is fixed as a mixture of two Gaussian
  distributions with equal variance. Mclust is run on the resampled
  data, which uses hierarchical clustering for initialization followed
  by EM.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  A list with the following components
 %% bootProbLL = bootProbLL, clusters = clusters, em_mean = em_mean, em_mix
 %% = em_mix, em_std = em_std, its = its, latency = latency, priorProbLL =
 %% priorProbLL, r_boot = r_boot, r_curve_boot = r_curve_boot
  
  \item{bootProbLL }{The probability that a subject is classified as LL,
    defined as the number of times this subject was classified LL over all
    bootstraps. This is distinct from the EM probability derived from the
    Gaussian mixture model.}
  \item{bootThreshLL }{The smallest integer value of the average latency
    that would be classified as LL using the model from each bootstrap.}
  \item{clusters }{The object returned from the intial call to Mclust on
    the entire data set. This contains the baseline classification that
    used to initialize the bootstraps.}
  \item{em_mean }{A matrix containing the cluster means from all
    bootstraps.}
  \item{em_mix }{A matrix containing the Gaussian mixing parameters from
    all bootstraps.}
  \item{em_std }{A matrix containing the estimated variance of the
    Gaussian distributions from all bootstraps.}
  \item{its }{The number of bootstraps.}
  \item{latency }{The latency vector passed to the function.}
  \item{priorProbLL }{The probability of LL classification from Mclust
    on the original data. This is used to initialize the bootstraps and
    optionally to restratify the data before each bootstrap.}
  \item{r_boot }{A matrix containing the EM probability for each subject
    in the input latency vector, from all bootstraps.}
  \item{r_curve_boot }{A matrix containing the EM probability for all
    integers 1:900, from all bootstraps. Used for deriving percentiles
    for the classification across the spectrum of possible latency
    values (0-900 seconds).}
 
 The ordering of the per-subject values (such as bootProbLL) is the same
 as in the latency vector passed to bootstrapClassification. 

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Philip A Cook <cookpa@mail.med.upenn.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{plotBootstrapEM_Class}}
  \code{\link{Mclust}}
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

boot = bootstrapEM_Class(c(sl, ll), its = 100)

## Not run
\dontrun{plotBootstrapEM_Class(boot)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bootstrap }
\keyword{ ~clustering }% __ONLY ONE__ keyword per line
