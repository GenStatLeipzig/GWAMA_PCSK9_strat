#' @param beta 	Vector of effect size estimates.
#' @param se 		Vector of associated standard errors.
#' @param priorsd 	Scalar specifying the standard deviation of the prior on true effect sizes.
#' @param log 		Whether to return results on a natural log scale.
#'
#' @return
#'    A vector of approximate Bayes factors, on a log scale if
#'    \code{log=TRUE}.  Higher values indicate stronger support for
#'    association (which is inverted relative to the original definition).
#'
#'
#' @examples
#' data(agtstats)
#' agtstats$pval <- with(agtstats, pchisq((beta/se.GC)^2, df = 1, lower.tail = FALSE))
#' max1 <- function(bf) return(bf/max(bf, na.rm = TRUE))
#' agtstats$BF.normal <- with(agtstats, max1(abf.Wakefield(beta, se.GC, 0.05)))
#' agtstats$BF.t <- with(agtstats, max1(abf.t(beta, se.GC, 0.0208)))
#' with(agtstats, plot(-log10(pval), log(BF.normal)))
#' with(agtstats, plot(-log10(pval), log(BF.t)))
#'
#' @author   Toby Johnson \email{Toby.x.Johnson@gsk.com}
#'
#' @export
## wakefield ABF inverted to be BF for alternate model relative to null
## assumes prior mean is zero
## taken from https://rdrr.io/github/tobyjohnson/gtx/src/R/abf.R
abf.Wakefield <- function(beta, se, priorsd, log = FALSE) {
  if (log) {
    return(log(sqrt(se^2/(se^2 + priorsd^2))) + 
             (beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2))
  } else {
    return(sqrt(se^2/(se^2 + priorsd^2)) * 
             exp((beta/se)^2/2 * priorsd^2/(se^2 + priorsd^2)))
  }
}