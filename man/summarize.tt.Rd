% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.tt.R
\name{summarize.tt}
\alias{summarize.tt}
\title{Summarize differential methylation analysis}
\usage{
summarize.tt(tt)
}
\arguments{
\item{tt}{Data.frame top-table object output from \code{\link{limma.tt}}. Must have column names 'adj.P.Val' and "betaFC"}
}
\value{
named vector of length 2.
"Signif":  number of CpG sites with FDR q-value < 0.05
"Sig.Eff" number of CpG sites with FDR q-value < 0.05 and mean effect size of \eqn{|\beta| > 0.2}

 @export
}
\description{
Create a brief numerical based on the output of a \code{limma} analysis of differentially methylated CpG sites
}

