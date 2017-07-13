#' Mock RTSS/AS01 data set
#'
#' A dataset containing data that is similar in structure to the RTSS/AS01
#' malaria vaccine trial. Privacy agreements prevent the sharing of the real
#' data, so please note THAT THIS IS NOT THE REAL RTSS,S DATA. The data set is a
#' \code{list} of 10 simulated multiple outputation draws. The covariate data,
#' \code{ftime}, and \code{vaccine} stay the same across the data sets; however,
#' the \code{ftype} variable changes, simulating output data sets of multiply
#' infected trial participants.
#'
#' @format A \code{list} with 10 entries, each a \code{data.frame} with 6,890
#'         rows and 18 columns.
#' \describe{
#'   \item{ftime}{number of months until first recorded malaria disease}
#'   \item{ftype}{the genotype of sampled malaria parasite
#'   (0 = censored, 1 = CSP matched, 2 = CSP mismatched)}
#'   \item{vaccine}{vaccine assignment (0 = control vaccine, 1 = vaccine)}
#'   \item{ageWeeks}{participant's age in weeks at trial enrollment}
#'   \item{weightForAgeZscore}{WHO weight-for-age Z-score}
#'   \item{sex}{participant's sex (0 = male, 1 = female)}
#'   \item{site1-5}{Indicator of study site}
#'   \item{heightForAgeZscore}{WHO height-for-age Z-score}
#'   \item{weightForHeightZscore}{WHO weight-for-height Z-score}
#'   \item{armCircumZscore}{WHO arm circumference Z-score}
#'   \item{hemog}{hemoglobin}
#'   \item{distInpatient}{distance from nearest inpatient clinic}
#'   \item{distOutpatient}{distance from nearest outpatient clinic}
#'   \item{startMonthCat}{study site-specific indicator of rainy (=1) versus dry
#'                        (=0) season}
#'   ...
#' }

"rtss"
