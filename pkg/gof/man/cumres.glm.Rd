\name{cumres.glm}
\alias{cumres}
\alias{cumres.glm}
\alias{cumres.lm}
\title{Calculates GoF statistics based on cumulative residual processes}
\usage{
  \method{cumres}{glm} (model,
    variable = c("predicted", colnames(model.matrix(model))),
    data = data.frame(model.matrix(model)), R = 1000,
    b = 0, plots = min(R, 50), breakties = 1e-12,
    seed = round(runif(1, 1, 1e+09)), ...)
}
\arguments{
  \item{model}{Model object (\code{lm} or \code{glm})}

  \item{variable}{List of variable to order the residuals
  after}

  \item{data}{data.frame used to fit model (complete
  cases)}

  \item{R}{Number of samples used in simulation}

  \item{b}{Moving average bandwidth (0 corresponds to
  infinity = standard cumulated residuals)}

  \item{plots}{Number of realizations to save for use in
  the plot-routine}

  \item{breakties}{Add unif[0,breakties] to observations}

  \item{seed}{Random seed}

  \item{...}{additional arguments}
}
\value{
  Returns an object of class 'cumres'.
}
\description{
  Given the generalized linear models model
  \deqn{g(E(Y_i|X_{i1},...,X_{ik})) = \sum_{i=1}^k
  \beta_jX_{ij}} the \code{cumres}-function calculates the
  the observed cumulative sum of residual process,
  cumulating the residuals, \eqn{e_i}, by the jth
  covariate: \deqn{W_j(t) = n^{-1/2}\sum_{i=1}^n
  1_{\{X_{ij}<t\}}e_i} and Kolmogorov-Smirnov and
  Cramer-von-Mises test statistics are calculated via
  simulation from the asymptotic distribution of the
  cumulative residual process under the null (Lin et al.,
  2002).
}
\note{
  Currently linear (normal), logistic and poisson
  regression models with canonical links are supported.
}
\examples{
sim1 <- function(n=100, f=function(x1,x2) {10+x1+x2^2}, sd=1, seed=1) {
  if (!is.null(seed))
    set.seed(seed)
  x1 <- rnorm(n);
  x2 <- rnorm(n)
  X <- cbind(1,x1,x2)
  y <- f(x1,x2) + rnorm(n,sd=sd)
  d <- data.frame(y,x1,x2)
  return(d)
}
d <- sim1(100); l <- lm(y ~ x1 + x2,d)
system.time(g <- cumres(l, R=100, plots=50))
g
\donttest{plot(g)}
g1 <- cumres(l, c("y"), R=100, plots=50)
g1
g2 <- cumres(l, c("y"), R=100, plots=50, b=0.5)
g2
}
\author{
  Klaus K. Holst
}
\references{
  D.Y. Lin and L.J. Wei and Z. Ying (2002)
  \emph{Model-Checking Techniques Based on Cumulative
  Residuals}. Biometrics, Volume 58, pp 1-12.

  John Q. Su and L.J. Wei (1991) \emph{A lack-of-fit test
  for the mean function in a generalized linear model}.
  Journal. Amer. Statist. Assoc., Volume 86, Number 414, pp
  420-426.
}
\seealso{
  \code{\link[timereg]{cox.aalen}} in the
  \code{timereg}-package for similar GoF-methods for
  survival-data.
}
\keyword{models}
\keyword{regression}

