#' Plot for a DMRcate output
#'
#' Volcano-style plot of effect-size vs. significance with size of DMR represented
#'
#' Each DMR in the input object will be colored according to whether they are considered hypermethylated in the reference group, hypomethylated, or neither. This will be determined by user-specified p-value and effect-size thresholds
#'
#' @param res.obj A \code{dmrcate.output} object
#' @param pval.fld Which p-value from DMRcate analysis to plot, either the mean p-value in the DMR (\code{'mean'}) or the minimum p-value (\code{'min'}).
#' @param pthresh Numeric value signifying the maximum p-value at which DMRs should be colored
#' @param contrast.thresh Numeric value signifying the threshold for maximum absolute effect size at which DMRs should be colored.
#' @param min.p Numeric minimum value for plotted p-values. Gives an upper limit to the y-axis
#' @param xlabel,ylabel,title Text labels for axis and main titles in plot
#' @param label.thresh Value for maximum p-value threshold with which to label points with gene names
#' @param colors for neutral, hypomethylated, and hypermethylated regions, respectively
#'
#' @return A \code{ggplot} object.
#' @export
#'
#'
## plot for a DMRcate result object:
## volcano plot but have variable size per point, based on number of probes in gene
## volcano plot function
### Label points with label.fld (if not NULL)
### use label.thresh and label.thresh.fld to determine what gets labeled.
### automatically replaces p-values of 0 with 1e-300
volcano_plot.DMRcate <- function(res.obj, pval.fld = c('mean', 'min'),
                                 pthresh = 0.01, contrast.thresh = 0.25, min.p = 1e-100,
                                 xlabel = 'max beta difference',
                                 ylabel = '-log10 p-value',
                                 title = '',
                                 label.thresh = 0.001,
                                 colors = c(neut = 'darkgrey', low = 'red', high = 'blue')){
  require(ggplot2)

  stopifnot(is(res.obj, 'dmrcate.output'))
  p.fld <- match.arg(pval.fld)
  df <- res.obj$results[c('gene_assoc', paste0(p.fld, 'pval'), 'maxbetafc', 'no.probes')]
  names(df) <- c('Gene', 'pval', 'contrast', 'no.probes')
  df$colorClass <- 'neut'
  df$colorClass[df$pval < pthresh & df$contrast < -contrast.thresh] <- 'hypomethylated'
  df$colorClass[df$pval < pthresh & df$contrast > contrast.thresh] <- 'hypermethylated'
  names(colors)[2:3] <- c('hypomethylated', 'hypermethylated')

  ### labeling:
  df$labl <- ''
  df$hjust <- as.numeric(df$contrast < 0)

  ix <- which(df$pval < label.thresh)
  df$labl[ix] <- gsub(';.*', '', df$Gene[ix])

  ## add p-value to label for points below min.p
  ix <- which(df$pval < min.p)
  df$labl[ix] <- sprintf('%s(%.0e)', df$labl[ix], df$pval[ix])
  # set low p-values to min.p + jitter
  ix <- ix[order(df$pval[ix])]
  df$pval[ix] <- 10^(log10(min.p) + sort(rnorm(length(ix), 0, 3)))

  ### add label to p-value threshold line:
  pval.line.label <- paste('p =', pthresh)
  pval.x <- min(df$contrast)

  ## for x-limit:
  xmax <- max(abs(df$contrast))*1.2

  ggplot(df, aes(x = contrast, y = -log10(pval), color = colorClass, size = no.probes)) +
    geom_point(alpha = 0.33) +
    xlab(xlabel) + ylab(ylabel) +
    ylim(c(0, max(-log10(df$pval)) + 1)) + xlim(c(-xmax, xmax)) +
    scale_color_manual(name = 'effect direction', values = colors,
                       limits = c('hypomethylated', 'hypermethylated')) +
    scale_size_area(name = '#probes', max_size = 15) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text(aes(x = contrast, y = -log10(pval), label = labl, hjust = hjust),
              data = df, vjust = 0, color = 'black', size = 3) +
    theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'),
          plot.title = element_text(size = 15, face = 'bold'),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = 'right',
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
             hjust = 0, vjust = 0, size = rel(4), fontface = 3)
  #guides(fill = guides_legend())
}
