\name{socialdefeat-package}
\alias{socialdefeat-package}
\alias{socialdefeat}
\docType{package}
\title{
  Package for analyzing social defeat data
}
\description{
  Contains code for analysis of social defeat data.

  Currently just a bunch of functions for doing a cluster analysis to
  classify animals as short latency (SL) or long latency (LL) based on
  their average latency.
}
\details{

  \tabular{ll}{
    Package: \tab socialdefeat\cr
    Type: \tab Package\cr
    Version: \tab 0.1\cr
    Date: \tab 2014-01-23\cr
    License: \tab GPL v3 or later\cr
  }
  
  The package contains functions for classifying animals into SL and LL
  groups based on average latency. The classification can be done with a
  model-based EM algorithm, or with a clustering algorithm such as
  k-means, pam, or hierarchical clustering.

  Uncertainty in classification can be estimated via bootstrap
  resampling. Functions to visualize the results of these experiments
  are also provided.

  See \code{bootstrapClusterClass} and \code{bootstrapEM_Class} for
  executable examples. These are just examples to play with using
  dummy data, see the function documentation for recommended
  parameters.
  

}
\author{
  Philip A Cook
  
  Maintainer: Philip A Cook <cookpa@pennmedicine.upenn.edu>
}
\references{
  
}

\keyword{ socialdefeat }
\seealso{
  \code{\link[mclust]{hc}}
  \code{\link[mclust]{Mclust}}
  \code{\link[stats]{kmeans}}
  \code{\link[cluster]{pam}}
}
\examples{

}
