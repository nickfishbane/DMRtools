#' Summarize top-table for various cutoffs
#'
#' Get counts from a limma analysis of how many features pass statistical significance and effect-size relevance
#'
#' @param tt A top-table output from \code{limma}
#' @param pval.thresh Numerical vector of p-values in (0,1)
#' @param beta.thresh Numerical vector of contrast (effect) sizes
#'
#' @return Matrix of integers for each contrast and p-value pair, how many features have a lower p-value and contrast of greater magnitude
#'
#'

summarize.tt.cutoffs <- function(tt, pval.thresh,
                                 beta.thresh = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)){
  library(dplyr)
  library(reshape2)
  stopifnot(is.data.frame(tt),
            all(c('betaFC', 'P.Value') %in% names(tt)),
            all(pval.thresh < 1), all(pval.thresh > 0))
  df <- data.frame(P = rep(sort(pval.thresh), each = length(beta.thresh)),
                   B = rep(sort(beta.thresh, decreasing = T), length(pval.thresh)))
  df <- df %>% group_by(P, B) %>%
    summarise(N = nrow(filter(tt, P.Value < P, abs(betaFC) > B))) %>% ungroup
  acast(df, P ~ B)
}
