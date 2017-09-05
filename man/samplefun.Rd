\name{samplefun}
\alias{samplefun}
\title{Uniformly Sampling}
\description{Draws samples uniformly (for internal use only)}
\usage{
samplefun(x, index, n0)
}
\arguments{
  \item{x}{ The VD values}
  \item{index}{the index for VD value, it can be 1, 2, 3,....}
  \item{n0}{Sample size}
}
\details{Uniform sampling (internal use only)}
\value{

\item{index }{selected ids }
\item{x }{selected VD levels}
}

\author{Durazo-Arvizu, Ramon,   Sempos, Chris and Tian, Lu}

\seealso{ \code{\link{sampletot}}}



\examples{
x=rnorm(100)
index=1:100
samplefun(x, index, 40)

}
