# perform coverage test locally
library(covr)
cov <- package_coverage(type = "all", combine_types = FALSE,
                        line_exclusions = list("R/plots.R")
                       )
cov
