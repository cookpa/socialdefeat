\name{emProbLL}
\alias{emProbLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Calculate the probability of LL classification from an EM Gaussian
  mixture model.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Classifies new data using a pre-existing Gaussian mixture from an
  object of type \code{Mclust}. 
}
\usage{
emProbLL(clusters, latency)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusters}{
    An object of class \code{Mclust}, which should contain a two-class
    Gaussian mixture model fit to average latency data.
}
  \item{latency}{
    A vector of average latencies. The probability of LL classification
    is computed for each entry.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A vector containing LL probabilities for each input latency.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  %% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  Philip A Cook <cookpa@mail.med.upenn.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

set.seed(20140123)

sl <- rnorm(60, 200, 100)
sl[which(sl < 0)] <- 0

ll <- rnorm(15, 600, 100)
ll[which(ll > 900)] <- 900

latency = c(sl, ll)

library(mclust)

clusters = Mclust(latency, G=2)

emProbLL(clusters, c(300, 400, 500))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
