# survtmle 1.1.3
* Extensive improvement to documentation throughout the package, including the
  use of Roxygen macros like `\link{}`.
* Addition of a new argument `returnCall` to `survtmle` that allows disabling
  of the return of `call` objects. This is useful for resolving the efficiency
  issues described in https://github.com/benkeser/survtmle/issues/40.

# survtmle 1.1.2
* Minor functionality updates, bug fixes and documentation edits:
  * Removal of `geom_point` when `plot.tp.survtmle` with more than 1 time point.
  * Add `CITATION` functionality for using of `citation("survtmle")`.
  * Add argument `cv` to allow altering cross-validation in `SuperLearner`.

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
