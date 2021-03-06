% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.genes.R
\name{write.genes}
\alias{write.genes}
\title{write gene list from results}
\usage{
write.genes(res.obj, file, thresh = 0.05)
}
\arguments{
\item{res.obj}{Either a data.frame (limma top-table) with column names "adj.P.Val" and "Gene" or an object of class \code{dmrcate.output} created by the \code{dmrcate} function in the \code{DMRcate} package}

\item{file}{name (with suggested '.txt' extension) for output of gene list}

\item{thresh}{Numerical value for p-value cut-off when writing genes from top-table}
}
\description{
writes all the associated genes found in a \code{dmrcate.output} object or a
top-table generated by \code{limma.tt} that satisfy passing a certain threshold
}
\details{
if \code{res.obj} is a \code{dmrcate.output} object, it will print all genes with minimum p-value (across probes in the region) less than the supplied threshold.
If a threshold is desired, this should be specified in the creation of the
object, using \code{pcutoff} in the \code{DMRcate::dmrcate()} function. If
\code{res.obj} is a top-table, it will print all gene symbols where the FDR
q-value is below the supplied threshold
}

