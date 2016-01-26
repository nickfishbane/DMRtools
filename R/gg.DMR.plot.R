### This is a wrapper for gene.plot (above) which takes in either a region e.g.
### 'chr5:122433740-122435550' like in DMRcate results OR a gene symbol
#' Plot methylation pattern in genomic region
#'
#' @seealso \code{\link{gene.plot}}
#'
#' @param region A character value formatted e.g. "chr5:122433740-122435550"
#' @param gene A gene symbol to plot
#' @param f A data.frame of CpG site annotations for the 450K annotation, which is subsetted for the \code{gene.plot} function
#' @param ... Arguments to pass to \code{gene.plot()} function
#'
#' @return ggplot object, same as output of \code{gene.plot}
gg.DMR.plot <- function(region = NULL, gene = NULL, f = fDat, ...){
  stopifnot(!is.null(region) | !is.null(gene))
  if(!is.null(region)){
    chrom <- gsub(':.*', '', region)
    sta.pos <- as.numeric(gsub('.*:', '', gsub('-.*', '', region)))
    sto.pos <- as.numeric(gsub('.*-', '', region))
    f <- f %>% filter(chr == chrom, pos <= sto.pos, pos >= sta.pos)
  } else{
    ix <- laply(strsplit(f$GeneSymbol, ';'), function(x) gene %in% x)
    f <- f[ix, ]
  }
  gene.plot(f, ...)

}
