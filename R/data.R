#' Demographic data for Thailand 1960--1995
#'
#' A list containing population counts in 1960 and fertility rates,
#' survival proportions, migration proportions, and sex ratios at
#' birth at five year intervals between 1960 and 1995
#' (inclusive). Population counts, fertility rates, and survival and
#' migration proportions are by five-year age group.
#'
#' @format A list with five elements of various types:
#' \describe{
#'   \item{thai_base_pop_counts}{A list with components named
#'   \dQuote{male} and \dQuote{female} containing 17*1 matrices
#'   of age-sex-specific population counts by five-year age group in 1960.}
#'
#'   \item{thai_fert_rates}{A 17*8 matrix containing annual age-specific annual fertility rates. These can be interpreted as the average annual age-specific fertility rates which prevailed in the five-year intervals \eqn{[1960, 1965)}, \eqn{[1965, 1970)}, \ldots{}, \eqn{[1995, 2000)}.}
#'
#'   \item{thai_mig_props}{list with components named \dQuote{male} and \dQuote{female} containing 17*8 matrices of age-sex-specific annual migration proportions. These indicate the average annual net number of international migrants as a proportion of the population of Thailand in the five-year intervals \eqn{[1960, 1965)}, \eqn{[1965, 1970)}, \ldots{}, \eqn{[1995, 2000)}.}
#'
#'   \item{thai_surv_props}{list with components named \dQuote{male} and \dQuote{female} containing 18*8 matrices of age-sex-specific 5-year survival proportions. Rows 1--17 indicate the proportion of those aged \eqn{[a, a+5)} at year \eqn{t} that survive to be aged \eqn{[a+5, a+10)} at year \eqn{t+5}. The 18th row of these matrices is sex-specific survival in the open age interval, that is proportion of those aged \eqn{[\omega, \infty)}{w+} at year \eqn{t} that survive five more years.}
#'
#'   \item{thai_srb}{A 1*8 matrix containing sex ratios at birth prevailaing over the five-year intervals \eqn{[1960, 1965)}, \eqn{[1965, 1970)}, \ldots{}, \eqn{[1995, 2000)}.}
#'
#' }
#' @source Wheldon, M. C., Raftery, A. E., Clark, S. J., and Gerland, P. (2015), \emph{Bayesian Reconstruction of Two-Sex Populations by Age: Estimating Sex Ratios at Birth and Sex Ratios of Mortality}, Journal of the Royal Statistical Society: Series A (Statistics in Society), 178, 977--1007. \url{https://doi.org/10.1111/rssa.12104}
"Thailand_demog"
