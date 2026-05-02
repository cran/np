## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(np.messages = FALSE)

## -----------------------------------------------------------------------------
library(np)
data(cps71, package = "np")

bw <- npregbw(logwage ~ age, data = cps71)
summary(bw)

fit <- npreg(bws = bw)
summary(fit)

## ----fig.width = 6, fig.height = 4--------------------------------------------
plot(cps71$age, cps71$logwage, cex = 0.25, col = "grey")
lines(cps71$age, fitted(fit), col = 2, lwd = 2)

## -----------------------------------------------------------------------------
set.seed(42)
mydat <- data.frame(
  y = rnorm(200),
  x_cont = runif(200),
  x_unordered = factor(sample(c("a", "b", "c"), 200, replace = TRUE)),
  x_ordered = ordered(sample(1:4, 200, replace = TRUE))
)

bw_mixed <- npregbw(y ~ x_cont + x_unordered + x_ordered, data = mydat)
fit_mixed <- npreg(bws = bw_mixed)
summary(fit_mixed)

## -----------------------------------------------------------------------------
if (requireNamespace("crs", quietly = TRUE) &&
    utils::packageVersion("crs") >= package_version("0.15-41")) {
  set.seed(7)
  n <- 120
  x <- runif(n, -1, 1)
  y <- x + 0.4 * x^2 + rnorm(n, sd = 0.18)

  fit_nomad <- npreg(y ~ x, nomad = TRUE, degree.max = 1L, nmulti = 1L)
  fit_nomad$bws$nomad.shortcut

  # Tune one component explicitly while leaving the rest of the preset in place.
  fit_nomad_direct <- npreg(
    y ~ x,
    nomad = TRUE,
    search.engine = "nomad",
    degree.max = 1L,
    nmulti = 1L
  )
}

