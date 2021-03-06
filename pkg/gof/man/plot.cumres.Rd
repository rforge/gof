\name{plot.cumres}
\alias{plot.cumres}
\title{Plot cumulative residuals from a 'cumres' object}
\usage{
  \method{plot}{cumres} (x, idx = 1:length(x$variable),
    col = c("grey"), ci = TRUE, col.ci = "darkblue",
    col.alpha = 0.3, lty.ci = 0, level = 0.95,
    legend = c("type1", "type2", "none"), xlab, ylab,
    vs = TRUE, ylim = NULL, title, ...)
}
\arguments{
  \item{x}{Object produced by the function \code{cumres}.}

  \item{idx}{vector of numbers (or variable names)
  indicating which processes from the \code{x} to plot.}

  \item{col}{Color of the sample processes. By setting this
  parameter to "none" or \code{NULL} no realizations will
  be drawn. The number of realizations is determined by the
  \code{cumres}-object.}

  \item{ci}{Type of prediction bands to plot. Defaults to
  none. Set to \code{TRUE} to obtain simultaneous
  prediction bands under the null (pointwise can be
  obtained by setting to "pointwise").}

  \item{col.ci}{Color of prediction band.}

  \item{col.alpha}{Degree of transparency (0-1) of the
  prediction bands.}

  \item{lty.ci}{Line type of prediction band.}

  \item{level}{The required prediction level.}

  \item{legend}{Type of legend where "type1" gives p-values
  of GoF-tests and "type2" gives usual type of legends.}

  \item{xlab}{Optional label of x-axis}

  \item{ylab}{Optional label of y-axis}

  \item{vs}{Label of predictor}

  \item{ylim}{Range of y axis}

  \item{title}{Main title}

  \item{...}{Additional arguments passed to the
  plot-routine.}
}
\description{
  \code{plot} displays the observed cumulative residual
  process with realizations under the null. 95\% prediction
  bands
}
\examples{
n <- 500; x <- abs(rnorm(n,sd=0.2))+0.01; y <- sqrt(x) + rnorm(n,sd=0.2)
l <- lm(y ~ x)
g <- cumres(l, R=500)
plot(g, idx=1, ci="sim", col=NULL, col.ci="purple", legend="type2")
}
\author{
  Klaus K. Holst <kkho@biostat.ku.dk>
}
\seealso{
  \code{\link[gof]{cumres}}
}
\keyword{hplot}
\keyword{regression}

