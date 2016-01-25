#' Plot methylation pattern in genomic region
#'
#' \code{gene.plot} returns a plot of \eqn{\beta}-values for a particular genomic region
#' and optionally according to a grouping variable for the samples
#'
#' To add: gene names and genomic locations (e.g. 5'UTR, Body, etc.)
#' This function should be called by \code{gg.DMR.plot()}
#' @seealso \code{\link{gg.DMR.plot}}
#'
#' @param fSub A data.frame of CpG site 450K annotation.
#' @param betas A matrix of beta-values containing all CpG sites contained in fSub.
#' @param pheno A data.frame of phenotypes with rows in same order as columns of betas.
#' @param groupVar (optional) Character value corresponding to a column of pheno, indicating grouping variable with which to color points in the plot.
#' @param gene.fld,cpg.fld,pos.fld,chr.fld Character values for names of fSub denoting, respectively: gene symbol/gene name, CpG identifier, genomic position (b.p.) and chromosome.
#'
#' @examples
#' f <- data.frame(GeneSymbol = c(rep('', 2), rep('BLAH', 3), rep('', 1)),
#'                 cpgSite = paste0('cg0000', 1:6),
#'                 pos = 17500000 + c(1, 147, 1458, 1500, 2157, 2699),
#'                 chr = 'chr19')
#'
#' p <- data.frame(ID = paste0('S', 1:10), Group = rep(LETTERS[1:2], each = 5))
#' set.seed(123)
#' means <- c(2.5, 2, 3, 3, 0, -1, 0, -1.2, 0.2, -1, 0.5, -0.5)
#' SD <- c(0.5, 0.5, 0.7, 0.8, 1, 0.5)
#' mvals <- matrix(rnorm(60, rep(means, each = 5), rep(SD, each = 10)), ncol = 6,
#'             dimnames = list(p$ID, f$cpgSite))
#' b <- t(2^mvals/(1+2^mvals))
#' gene.plot(fSub = f, betas = b, pheno = p, groupVar = 'Group')
#'
#' @return A \code{ggplot} object for the plot of interest
#'
#' @export
#'
gene.plot <- function(fSub, betas, pheno, groupVar = NULL,
                      gene.fld = 'GeneSymbol', cpg.fld = 'cpgSite',
                      pos.fld = 'pos', chr.fld = 'chr'){
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  stopifnot(ncol(betas) == nrow(pheno),
            groupVar %in% names(pheno),
            all(c(gene.fld, cpg.fld, pos.fld, chr.fld) %in% names(fSub)))
  pca.theme <- function(i = 1){
    require(ggplot2)
    theme(text = element_text(colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = NA),
          panel.border = element_rect(fill = NA, colour = 'black'),
          panel.grid = element_blank(),
          axis.text = element_text(size = rel(i), colour = 'black'),
          axis.title = element_text(size = rel(i*1.5)),
          plot.title = element_text(size = rel(i*1.5)))
  }

  ### get CpG sites relating to gene, using features
  feat <- fSub[c(gene.fld, cpg.fld, pos.fld, chr.fld)]
  names(feat) <- c('Gene', 'CpG', 'POS', 'CHR')
  feat <- filter(feat, CpG %in% rownames(betas))
  stopifnot(n_distinct(feat$CHR) == 1)
  if(!is.null(groupVar)){
    plotDat <- cbind(as.data.frame(t(betas[feat$CpG, ])),
                     varb = pheno[[groupVar]])

    plotDat <- melt(plotDat, id.var = 'varb', variable.name = 'cgSite',
                    value.name = 'Beta')
    plotDat <- merge(plotDat, feat, all.x = T, by.x = 'cgSite', by.y = 'CpG')
    plotDat.med <- plotDat %>% group_by(cgSite, varb) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))

    d <- min(dist(unique(plotDat$POS)))
    ggplot(plotDat, aes(x = POS, y = Beta, color = varb, fill = varb)) +
      geom_point(alpha = 0.5,
                 position = position_jitterdodge(jitter.width = 0,
                                                 dodge.width = d*0.8)) +
      geom_line(aes(x = POS, y = medianBeta, color = varb), plotDat.med,
                position = position_jitterdodge(jitter.width = 0,
                                                dodge.width = d*0.8)) +
      scale_color_discrete(name = groupVar) +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme() +
      scale_fill_discrete(guide = 'none')
  } else{
    plotDat <- t(betas[feat$CpG, ])
    plotDat <- plotDat %>%
      melt(varnames = c('subj', 'cgSite'), value.name = 'Beta') %>%
      merge(feat, all.x = T, by.x = 'cgSite', by.y = 'CpG') %>%
      mutate(POSeven = factor(match(POS, sort(unique(POS))) %% 2))
    plotDat.med <- plotDat %>% group_by(cgSite) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))

    ggplot(plotDat, aes(x = POS, y = Beta, color = POSeven)) +
      geom_point(alpha = 0.5) +
      geom_line(aes(x = POS, y = medianBeta), plotDat.med, color = 'black') +
      scale_color_manual(values = c('red', 'blue')) +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme() +
      theme(legend.position = 'none') +
      geom_hline(yintercept = c(-1,0,1), color = 'grey20') +
      ylim(range(plotDat$Beta))

  }
  ###############
  ####
  ###!!!!!!!!!! TO BE COMPLETED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}
