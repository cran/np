## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(np.messages = FALSE)

## -----------------------------------------------------------------------------
library(np)
set.seed(42)

n <- 250
x <- rnorm(n)
y <- rnorm(n)

npunitest(x, y, bootstrap = TRUE)

