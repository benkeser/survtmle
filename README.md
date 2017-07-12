
<!-- README.md is generated from README.Rmd. Please edit that file -->
R/`survtmle`
============

[![Travis-CI Build Status](https://travis-ci.org/benkeser/survtmle.svg?branch=master)](https://travis-ci.org/benkeser/survtmle) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/benkeser/survtmle?branch=master&svg=true)](https://ci.appveyor.com/project/benkeser/survtmle) [![Coverage Status](https://img.shields.io/codecov/c/github/benkeser/survtmle/master.svg)](https://codecov.io/github/benkeser/survtmle?branch=master) [![CRAN](http://www.r-pkg.org/badges/version/survtmle)](http://www.r-pkg.org/pkg/survtmle) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted Minimum Loss-Based Estimation (TMLE) for Survival Analysis with Competing Risks

**Authors:** [David Benkeser](https://www.benkeserstatistics.com/) and [Nima Hejazi](http://nimahejazi.org)

------------------------------------------------------------------------

Description
-----------

`survtmle` is an R package designed to use targeted minimum loss-based estimation (TMLE) to compute covariate-adjusted marginal cumulative incidence estimates in right-censored survival settings with and without competing risks. The estimates can leverage ensemble machine learning via the `SuperLearner` package.

------------------------------------------------------------------------

Installation
------------

You can install a stable release of `survtmle` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) with:

``` r
devtools::install_github("benkeser/survtmle")
```

For standard use, we recommend installing the package from [CRAN](https://cran.r-project.org/) (an initial release will be coming soon) via

``` r
install.packages("survtmle")
```

------------------------------------------------------------------------

Issues
------

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/benkeser/survtmle/issues).

------------------------------------------------------------------------

Example
-------

This minimal example shows how to use `survtmle` to obtain cumulative incidence estimates with a very simple, simulated data set.

``` r
library(survtmle)
#> survtmle: Targeted Learning for Survival Analysis
#> Version: 0.1.2
set.seed(341796)

# simulate data
n <- 100
t_0 <- 10
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, 0.5)
T <- rgeom(n,plogis(-4 + W$W1 * W$W2 - A)) + 1
C <- rgeom(n, plogis(-6 + W$W1)) + 1
ftime <- pmin(T, C)
ftype <- as.numeric(ftime == T)

# apply survtmle for estimation
fit <- survtmle(ftime = ftime, ftype = ftype,
                trt = A, adjustVars = W,
                glm.trt = "1",
                glm.ftime = "I(W1*W2) + trt + t",
                glm.ctime = "W1 + t", method = "hazard", 
                t0 = t_0)

# extract cumulative incidence at each time point
tpfit <- timepoints(fit, times = seq_len(t_0))

# examine output object produced by the timepoints function
print(tpfit)
#> $est
#>              t1         t2         t3         t4         t5         t6
#> 0 1 0.032997470 0.06492788 0.09582530 0.12572293 0.15465313 0.18264737
#> 1 1 0.008014555 0.01603567 0.02406256 0.03209448 0.04013064 0.04817027
#>             t7         t8         t9        t10
#> 0 1 0.20973629 0.23594966 0.26131640 0.28586459
#> 1 1 0.05621257 0.06425675 0.07230203 0.08034761
#> 
#> $var
#>               t1           t2           t3           t4           t5
#> 0 1 4.565496e-04 0.0005099637 0.0004951975 0.0004852183 0.0005862762
#> 1 1 2.111604e-06 0.0003345475 0.0003208485 0.0003079889 0.0006317931
#>               t6           t7           t8           t9         t10
#> 0 1 0.0012883360 0.0013743743 0.0012847990 0.0020298371 0.003137739
#> 1 1 0.0006229943 0.0009761833 0.0009544762 0.0009534062 0.000956650

# examine plot of cumulative incidences
plot(tpfit)
```

<img src="README-example-1.png" style="display: block; margin: auto;" />

------------------------------------------------------------------------

License
-------

© 2016-2017 [David C. Benkeser](http://www.benkeserstatistics.com)

The contents of this repository are distributed under the MIT license. See below for details:

    The MIT License (MIT)

    Copyright (c) 2016-2017 David C. Benkeser

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
