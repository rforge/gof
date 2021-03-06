\name{cumres.lvmfit}
\alias{cumres.lvmfit}
\title{Cumulative residual processes for structural equation models}
\usage{
  \method{cumres}{lvmfit} (model, y, x, full = FALSE,
    data = model.frame(model), p, R = 1000, b = 0,
    plots = min(R, 50), seed = round(runif(1, 1, 1e+09)),
    ...)
}
\arguments{
  \item{model}{\code{lvm} object}

  \item{y}{A formula specifying the association to be
  checked. Alternatively the outcome specified as a
  function or a string with the name of the outcome in the
  model.}

  \item{x}{Predictor. A function, vector or character}

  \item{full}{If FALSE the prediction, Pr, of the variable
  that are ordered after is only calculated based on the
  conditional distribution given covariates. If TRUE the
  conditional expectation is based on the largest set of
  covariates and endogenous variables such that the
  residual and Pr are uncorrelated.}

  \item{data}{data.frame (default is the model.frame of the
  model)}

  \item{p}{Optional parameter vector}

  \item{R}{Number of processes to simulate}

  \item{b}{Moving average parameter}

  \item{plots}{Number of processes to save for use with the
  plot method}

  \item{seed}{Random seed}

  \item{\dots}{Additional arguments parsed on to
  lower-level functions}
}
\value{
  Returns a \code{cumres} object with associated
  \code{plot},\code{print},\code{confint} methods
}
\description{
  Calculates GoF statistics based on cumulative residual
  processes for structural equation models fitted with the
  \code{lava} package.
}
\details{
  With \code{y} and \code{x} given as functions the user
  can decide which variables to use in the prediction of
  the outcome and predictor (use the \code{predict} method
  as below).
}
\examples{
\donttest{
library(lava)
m <- lvm(list(c(y1,y2,y3)~eta,eta~x)); latent(m) <- ~eta
## simulate some data with non-linear covariate effect
functional(m,eta~x) <- function(x) 0.3*x^2
d <- sim(m,100)

e <- estimate(m,d)
## Checking the functional form of eta on x
g <- cumres(e,eta~x,R=1000)
plot(g)

x <- function(p) predict(e,x=~y2+y3,p=p)[,"eta"]
## Checking the functional form of y1 on eta
cumres(e,y1~eta,R=1000)
g <- cumres(e,"y1",x=x,R=1000)
plot(g)

}
}
\author{
  Klaus K. Holst
}
\references{
  B.N. Sanchez and E. A. Houseman and L. M. Ryan (2009)
  \emph{Residual-Based Diagnostics for Structural Equation
  Models}. Biometrics Volume 65 (1), pp 104-115.
}
\keyword{models}
\keyword{regression}

