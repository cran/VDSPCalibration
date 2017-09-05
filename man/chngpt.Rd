\name{chngpt}
\alias{chngpt}
\title{Piecewise Regression Estimation}
\description{
Estimate a piecewise linear regression equation}
\usage{
chngpt(x, y, start = quantile(x, probs = 0.1,

na.rm = "TRUE"),finish = quantile(x, probs = 0.9, na.rm = "TRUE"),

NbrSteps = 500)
}

\arguments{
  \item{x}{
old VD measurements
}
  \item{y}{
reference (new) VD measurements
}
  \item{start}{
lower bound of the changing point
}
  \item{finish}{
upper bound of the changing point
}
  \item{NbrSteps}{
number of points used in grid search
}
}
\details{
This function uses grid search method to fit a piecewise linear regression model with one changing point
}
\value{
\item{x}{old VD levels}
\item{y}{new VD levels}
\item{yfitted}{calibrated VD levels based on the fitted piecewise linear regression}
\item{chngpt}{the estimated chang point}
\item{coefficients}{the estimated regression coefficients for the piecewise linear regression}
}
\references{Tian L., Durazo-Arvizu R. A., Myers G., Brooks S., Sarafin K., and Sempos C. T. (2014), The estimation of calibration equations for variables with heteroscedastic measurement errors, Statist. Med., 33, pages 4420-4436
}
\author{Durazo-Arvizu, Ramon and   Sempos, Chris }

\examples{
### Generate equally spaced TEST VALUES in the interval [20,200]
set.seed(123456789)
x= 20 + 180*1:100/100
x2= (x - 95)*(x>=95)

# Generate REFERENCE VALUES
y = -8 + 1.5*x - 0.85*x2 + 15*rnorm(100)

#Plot test values versus reference values along with fitted piecewise model
plot(x,y)
fit.chngpt = chngpt(x,y)
plot(fit.chngpt$x[order(fit.chngpt$yfitted)],
     fit.chngpt$y[order(fit.chngpt$yfitted)],
     xlim=c(0,200), ylim=c(0,200), xlab="25-Hydroxyvitamin D (nmol/mL), IDS",
     ylab="25-Hydroxyvitamin D (nmol/mL), LC/MS", bty="n", las=1)
lines(fit.chngpt$x[order(fit.chngpt$yfitted)],
      fit.chngpt$yfitted[order(fit.chngpt$yfitted)], lty=2,col=2, lwd=2)
abline(v=fit.chngpt$chngpt, lty=2,col=3, lwd=2)
arrows(fit.chngpt$chngpt+20 ,15, fit.chngpt$chngpt,-8, length=0.1, lwd=2, col=4)
legend(fit.chngpt$chngpt + 5,30, legend=round(fit.chngpt$chngpt, digits=1),
       bty="n", col=4)


}
