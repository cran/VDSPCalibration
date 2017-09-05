\name{samplesize}
\alias{samplesize}
\title{ Uniform Sampling Within Quartiles }
\description{Estimates the sample size to achived the specified precision in the estimated calibration equation. }
\usage{
samplesize(x0, d0, cutpts = c(7.5, 42.5, 57.5, 72.5, 200), CVx, CVy)
}
\arguments{
  \item{x0}{
The value at which calibration will be esitmated (e.g., 30 nmol/L)
}
  \item{d0}{
Targeted width of the 95\% confidence interval of the calibrated value (e.g. 5nmol/L)
}
  \item{cutpts}{
Cut points used to define intervals, within which samples would be selected uniformly
}
  \item{CVx}{
CV\% of the old method (e.g. 12\%)
}
  \item{CVy}{
CV\% of the reference (new) method (e.g. 5\%)
}
}
\details{The function estimates the sample size to achived the specified precision in the estimated calibration equation. The precision is defined via x0 and d0}
\value{
Required sample size to achived the specified precision in the estimated calibration equation.
}
\references{
Tian L., Durazo-Arvizu R. A., Myers G., Brooks S., Sarafin K., and Sempos C. T. (2014), The estimation of calibration equations for variables with heteroscedastic measurement errors, Statist. Med., 33, pages 4420-4436
}
\author{Durazo-Arvizu, Ramon,   Sempos, Chris and Tian, Lu }



\examples{

samplesize(30, 5, cutpts=c(7,42,57,72,200),0.12, 0.05)

}

