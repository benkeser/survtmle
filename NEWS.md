# survtmle 1.1.1 
* Minor bug fixes and documentation updates. 

# survtmle 1.1.0
* Adds support for the use of `speedglm` to fit the numerous regressions fit in
    the estimation procedure. Users may see warnings when `speedglm` fails, in
    which case the code defaults back to standard `glm`.
* Fixes problems with the `plot.tp.survtmle` method induced, by changes in the
    inner working of `tidyr` as of `tidyr` v0.8.0.
* Adds a method `confint.tp.survtmle` that computes and provides output tables
    for statistical inference directly from objects of class `tp.survtmle`. This
    provides information equivalent to that output by `confint.survtmle`.

# survtmle 1.0.0
* The first public release made available on CRAN.
