#' Summarize differential methylation analysis
#'
#' Create a brief numerical based on the output of a \code{limma} analysis of differentially methylated CpG sites
#'
#' @param tt Data.frame top-table object output from \code{\link{limma.tt}}. Must have column names 'adj.P.Val' and "betaFC"
#'
#' @return named vector of length 2.
#' "Signif":  number of CpG sites with FDR q-value < 0.05
#' "Sig.Eff" number of CpG sites with FDR q-value < 0.05 and mean effect size of \eqn{|\beta| > 0.2}
#'
#'  @export
summarize.tt <- function(tt){
  stopifnot(all(c('adj.P.Val', 'betaFC') %in% names(tt)))
  sig <- tt$adj.P.Val < 0.05
  sig.eff <- sum(abs(tt$betaFC[sig]) > 0.2 )
  c(Signif = sum(sig), Sig.Eff = sig.eff)
}

