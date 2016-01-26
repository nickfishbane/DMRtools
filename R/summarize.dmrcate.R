#' Summarize differential methylation analysis
#'
#' Create a brief numerical based on the output of a \code{DMRcate} analysis of differentially methylated regions (DMRs)
#'
#' @param res.obj \code{dmrcate.output} object output from \code{\link{DMRcate::dmrcate}}.
#'
#' @return named vector of length 3.
#' "min3":  number of DMRs with a minimum of 3 CpG sites
#' "min3_beta02": number of DMRs with at least 3 CpG sites and a max effect size of \eqn{|\beta| > 0.2}
#' "min3_anno": number of DMRs with at least 3 CpG sites and at least one gene annotated to the region.
#'
#'  @export

summarize.dmrcate <- function(res.obj){
  stopifnot(is(res.obj, 'dmrcate.output'))
  tab <- subset(res.obj$results, no.probes >= 3)
  sig.eff <- sum(abs(tab$maxbetafc) > 0.2)
  ann <- sum(tab$gene_assoc != '')
  c(min3 = nrow(tab), min3_beta02 = sig.eff, min3_anno = ann)
}
