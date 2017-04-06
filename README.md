# R/`survtmle`

[![Travis-CI Build Status](https://travis-ci.org/benkeser/survtmle.svg?branch=master)](https://travis-ci.org/benkeser/survtmle)
[![AppVeyor Build  Status](https://ci.appveyor.com/api/projects/status/github/benkeser/survtmle?branch=master&svg=true)](https://ci.appveyor.com/project/benkeser/survtmle)
[![Coverage Status](https://img.shields.io/codecov/c/github/benkeser/survtmle/master.svg)](https://codecov.io/github/benkeser/survtmle?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/survtmle)](http://www.r-pkg.org/pkg/survtmle)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Targeted Minimum Loss-Based Estimation (TMLE) for Survival Analysis with
> Competing Risks

---

## Description

`survtmle` is an R package designed to use Targeted Minimum Loss-Based
Estimation (TMLE) to compute marginal cumulative incidence estimates in
right-censored survival settings with and without competing risks, including
vaccine sieve analysis. This tool also provides facilities for computing and
obtaining inference for the use of data adaptive target parameters in such
settings.

---

## Installation

- Install the most recent _stable release_:
  `devtools::install_github("benkeser/survtmle")`

- To contribute, install the _development version_:
  `devtools::install_github("benkeser/survtmle", ref = "develop")`

---

## License

&copy; 2016-2017 [David C. Benkeser](http://www.benkeserstatistics.com)

The contents of this repository are distributed under the MIT license. See
below for details:
```
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
```
