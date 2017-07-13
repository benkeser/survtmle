#' Mock RV144 data set
#'
#' A dataset containing data that is similar in structure to the RV144 "Thai
#' trial" of the ALVAC/AIDSVAX vaccine. Privacy agreements prevent the sharing
#' of the real data, so please note THAT THIS IS NOT THE REAL RV144 DATA.
#'
#' @format A data frame with 15,955 rows and 10 columns:
#' \describe{
#'   \item{ftime}{number of six month visit windows until first recorded incidence of HIV}
#'   \item{ftype}{the genotype of HIV (0 = censored, 1 = amino acid site 169
#'                matched, 2 = amino acid site 169 mismatched)}
#'   \item{vax}{vaccine assignment (0 = placebo, 1 = vaccine)}
#'   \item{male}{male gender (0 = no, 1 = yes)}
#'   \item{year04}{trial enrollment year 2004 (0 = no, 1 = yes)}
#'   \item{year05}{trial enrollment year 2005 (0 = no, 1 = yes)}
#'   \item{medRisk}{medium category of risk behaviors (0 = no, 1 = yes)}
#'   \item{highRisk}{high category of risk behaviors (0 = no, 1 = yes)}
#'   \item{medAge}{medium category for age (0 = no, 1 = yes)}
#'   \item{highAge}{high category for age (0 = no, 1 = yes)}
#'   ...
#' }

"rv144"
