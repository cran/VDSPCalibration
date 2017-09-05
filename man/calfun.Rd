\name{calfun}
\alias{calfun}
\title{Estimating the Calibration Equation}
\description{ Estimates the calibration equation based on CV information}
\usage{
calfun(x, y, CVx, CVy = CVx, lambda0 = 1)
}

\arguments{
  \item{x}{old VD measurements}
  \item{y}{reference (new) VD measurements}
  \item{CVx}{CV\% of the old VD measurements}
  \item{CVy}{CV\% of the new VD measurements}
  \item{lambda0}{the CV ratio of the new vs old measurements}
}
\details{ Estimation of the calibration equation. It covers 4 scenarios: Only CVx is known; only CVy is known;  both CVx and CVy are known; and  Only the ratio of CVy to CVx is known.}
\value{
 \item{coef }{estimated coefficients of the linear function}
 \item{se }{standard errors of the estimated coefficients}
 \item{lower CI}{ the lower end of the 95\% CI of the regression coefficients}
 \item{upper CI}{ the upper end of the 95\% CI of the regression coefficients}
}
\references{Tian L., Durazo-Arvizu R. A., Myers G., Brooks S., Sarafin K., and Sempos C. T. (2014), The estimation of calibration equations for variables with heteroscedastic measurement errors, Statist. Med., 33, pages 4420-4436}
\author{Durazo-Arvizu, Ramon; Sempos, Chris;  Tian, Lu}

\examples{
n=100
sigma0=10

beta0=5
beta1=1.2
CVx=0.15
CVy=0.07

lambda0=CVy^2/CVx^2

x0=runif(n, 20, 200)
y0=beta0+beta1*x0+rnorm(n)*sigma0
x=x0+x0*CVx*rnorm(n)
y=y0+y0*CVy*rnorm(n)

fit=calfun(x, y, CVx, CVy, lambda0)
fit

}
