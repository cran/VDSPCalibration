\name{sampletot}
\alias{sampletot}
\title{Samples Selection}
\description{Selects samples used in a calibration study}
\usage{
sampletot(x, index, n0, K)
}
\arguments{
  \item{x}{the old sample measurements needing calibration}
  \item{index}{ the ID list of the old sample measurements needing calibration}
  \item{n0}{the required sample size}
  \item{K}{the number of quantiles, it is 4 if we use quartiles (recommended)}
}
\details{The function selectes samples used in the calibration study}
\value{

\item{x }{the selected sample measurements to be used in the calibration study}
\item{index }{the id list of the selected samples to be used in the calibration study}
}
\references{
Tian L., Durazo-Arvizu R. A., Myers G., Brooks S., Sarafin K., and Sempos C. T. (2014), The estimation of calibration equations for variables with heteroscedastic measurement errors, Statist. Med., 33, pages 4420-4436
}
\author{Durazo-Arvizu, Ramon,   Sempos, Chris and Tian, Lu}

\examples{

VD.value= 60 + 25*rnorm(1000)
VD.index=1:1000

### x:     the VD value
### index: the index for VD value, it can be 1, 2, 3,....
### n0:    the number of samples we want to select
### K:     the number of quantiles, it is 4 if we use quartiles

sampletot(x=VD.value, index=VD.index, n0=100, K=4)
}

